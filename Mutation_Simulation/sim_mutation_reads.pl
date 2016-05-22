#!/usr/bin/perl -w
#
#   sim_mutation_reads.pl -- Simulate reads with synthetic mutations use the method describled in Keightley et al. 2014.
#
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2016-05-22
#   Version: 2.0.0
#
#   Change logs:
#   Version 1.0.0 15/01/21: The initial version.
#   Version 1.0.1 15/07/30: Skip scaffold while sampling.
#   Version 2.0.0 16/05/22: Update: add options "--exclude", "--max-shared" and "--min-depth";
#                                   add support for simulate shared mutations;
#                                   add check for whether mutation is successfully generated;
#                                   now the number of simulated mutations were weighted according
#                                   to the chromosome length.




use strict;

use Data::Dumper;
use Data::Random qw(:all);
use Getopt::Long;
use List::Util qw(shuffle);
use List::Util::WeightedChoice qw(choose_weighted);


use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my %options             = ();
   $options{rand_size}  = 100;
   $options{out_format} = 'fastq';
   $options{max_shared} = 1;
   $options{min_depth}  = 1;
my ($no_rc, $use_rg_id,
    $samtools_opts, $min_seq_len, $max_clipping, $min_insert_size, $max_insert_size);
GetOptions(
            "fasta=s"          => \$options{fasta_file},
            "depth=s"          => \$options{depth_file},
            "random-size=i"    => \$options{rand_size},
            
            "exclude=s{,}"     => \@{$options{exclude_ids}},
            
            "bams=s{,}"        => \@{$options{bam_files}},
            
            "output=s"         => \$options{output},
            
            "samtools=s"       => \$samtools_opts,
            "min-insert=i"     => \$min_insert_size,
            "max-insert=i"     => \$max_insert_size,
            "min-len=i"        => \$min_seq_len,
            "max-clipping=i"   => \$max_clipping,
            "use-rg"           => \$use_rg_id,
            
            "max-shared=i"     => \$options{max_shared},
            "min-depth=i"      => \$options{min_depth},
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless($options{fasta_file} && $options{depth_file} && (@{$options{bam_files}} > 0) && $show_help ) {
    print <<EOF;

$0  -- Extract all reads with the 

Version: $VERSION

Usage:   perl $0 [options]

Input Options:

    -f, --fasta       <filename>
        reference sequences in fasta format, required
    
    -d, --depth       <filename>
        input file contain empirical distribution of read depths, required

    -b, --bams        <filename>
        bam file(s), at least one bam file should be specified


Output Options:

    -o, --output <filename>
        output filename, default to STDOUT
        
    -n, --no-rc
        do not reverse complement sequence with negtive strand
    -u, --use-rg
        add read group id to extracted records
        
Simulation Options:
    --random-size <int>
        number of mutations to be simulated, default: 100

    -s, --samtools <string>
        directly pass samtools view options to this script, e.g.
        "-f 4 -F 8"

    --min-len      <int>
        minium sequence length

    --max-clipping <int>
        maximum allowed clipping length, include both soft and hard
        clipping bases
    
    --min-insert   <int>
    --max-insert   <int>
        screen out records with insert size wihtin this range

    --max-shared   <int>
        choose how many samples could shared the mutation loci, the samples
        were also choosed by random [default: 1]
    --min-depth    <int>
        minimum required depth of "mutated" reads, loci with a depth smaller
        than this threshold would be set to this threshold [default: 1]
    
    -e, --exclude  <strings>
        exclude unwanted chromosomes or scaffolds while simulating, all
        chromosomes with ids match strings specified here would be ignored 

EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


if ($options{output}) {
    open (STDOUT, "> $options{output}") || die $!;
}


##
## sampling the numbers of nonreference reads from realistic distributions
##
print STDERR ">> Start sampling read depths ...";
my $rh_emp_depths = get_empirical_dist($options{depth_file}, $options{rand_size});
print STDERR "done!\n";



##
## parse reference genome infos
##
print STDERR ">> Start parsing reference sequences ... ";
my %CHROM_SEQs = ();
parse_fasta_SEQs(\%CHROM_SEQs, $options{fasta_file});


my %CHROM_LENGTHs = ();
my $GENOME_SIZE   = 0;
for my $chrom (sort keys %CHROM_SEQs)
{
    if (@{$options{exclude_ids}} > 0) {
        my $exclude_str = join '|', @{$options{exclude_ids}};
        next if ($chrom =~ /($exclude_str)/);
    }
    
    $CHROM_LENGTHs{$chrom} = length $CHROM_SEQs{$chrom};
    
    $GENOME_SIZE += $CHROM_LENGTHs{$chrom};
}

##
## weight each chromosome by its length
##
my @CHROM_IDs = sort keys %CHROM_LENGTHs;
my @chrom_weights = ();
for my $chrom (@CHROM_IDs)
{
    my $weight = $CHROM_LENGTHs{$chrom} / $GENOME_SIZE;
    push @chrom_weights, $weight;
}
print STDERR "done!\n";



##
## randomly generate mutation loci
##
print STDERR ">> Start generating simulated reads ... ";
print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "##Tag(SAM): replaced reads\n";
print STDOUT "##Tag(MUT): simulated loci\n";
print STDOUT "#Tag\tChrom\tStart\tEnd\tSamples\tRef\tMut\tSimulated_depths(Total_depth:Sampled_depth:Actually_replaced_depth)\n";
my $rand_index = 1;
while ($rand_index <= $options{rand_size})
{
    print STDERR "\r>> Start generating random mutations .. $rand_index \/ " . $options{rand_size};

    my $chrom = choose_weighted([@CHROM_IDs], [@chrom_weights]);
    my $pos   = int(rand($CHROM_LENGTHs{$chrom}+1));

    my $ref_base   = uc(substr($CHROM_SEQs{$chrom}, $pos-1, 1));
    my @alt_bases  = grep { $_ ne $ref_base } qw(A T G C);
    my $mut_base   = $alt_bases[int(rand(3))];

    my @samples_selected = rand_set(set => $options{bam_files}, min => 1, max => $options{max_shared});
       @samples_selected = sort @samples_selected;
       
    my @sim_samples = ();
    my @sim_depths  = ();
    for my $bam_file (@samples_selected)
    {
        my $sim_info = gen_mutated_reads($bam_file, $chrom, $pos, $ref_base, $mut_base, $rh_emp_depths);
        
        if ($sim_info ne -2) {
            my ($sample, $depth) = (split /\t/, $sim_info);
            
            push @sim_samples, $sample;
            push @sim_depths,  $depth;
        }
    }
    
    if (@sim_samples >= 1) {
        $rand_index ++;
        
        my $mut_samples = join ";", @sim_samples;
        my $mut_depths  = join ";", @sim_depths;
        
        print STDOUT "MUT\t$chrom\t$pos\t$pos\t$mut_samples\t$ref_base\t$mut_base\t$mut_depths\n";
    }
}
print STDERR "\tdone!\n";



print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################


=head2 get_empirical_dist

    About   : Sample numbers of nonreference reads from realistic distributions
    Usage   : get_empirical_dist($file);
    Args    : File contains numbers of nonreference reads for each depth.
    Returns : Hash reference of empirical distribution.

=cut
sub get_empirical_dist
{
    my ($in) = @_;
    
    my %emp_depths = ();
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/^\#/ || /^\s+$/);
        
        my ($total_depth, $alt_depth, $count) = (split /\s+/);
        
        push @{$emp_depths{$total_depth}->{depth}}, $alt_depth;
        push @{$emp_depths{$total_depth}->{count}}, $count;
    }
    
    return \%emp_depths;
}



=head2 gen_mutated_reads

    About   : Generate reads with mutations.
    Usage   : gen_mutated_reads($bam_file, $chrom, $pos, $ref_base, $mut_base, $rh_emp_depths);
    Args    : Bam file to extract reads;
              Chromosome of mutation to be synthesized;
              Position of mutation to be synthesized;
              Reference base for consistent check;
              Alternative base used as synthesized mutation allele;
              Hash reference to empirical distributions of read depths.
    Returns : Null.

=cut
sub gen_mutated_reads
{
    my ($in, $mut_chrom, $mut_pos, $ref_base, $mut_base, $rh_emp_depths) = @_;


    ##
    ## get the sample id from bam header
    ##
    open (my $header_fh, "samtools view -H $in |") || die $!;
    my $mut_sample = '';
    while (<$header_fh>)
    {
        next unless(/^\@RG.*SM:(.*?)(\s|$)/);
        $mut_sample = $1;
        last;
    }
    close $header_fh;
    
    if ($mut_sample eq '') {
        print STDERR "Error: No SM tag found in bam header!\n"; exit(2);
    }
    
    ##
    ## extract reads from bam file
    ##
    my $mut_region = "$mut_chrom:$mut_pos-$mut_pos";
    my $pipe_str = "samtools view $in $mut_region |";
    
    if ($samtools_opts) {
        $pipe_str = "samtools view $samtools_opts $in $mut_region |";
    }
    
    open (my $fh, $pipe_str) || die $!;
    my @reads = shuffle <$fh>;
    close $fh;
    
    
    ##
    ## randomly generate mutation depth based on empirical distribution,
    ## note some regions would have total depth not present in previous sampled
    ## distributions, a random number is used despite the empirical distribution
    ##
    my $total_depth = scalar @reads;
    
    return -2 if $total_depth < $options{min_depth};
    
    my $mut_depth   = $rh_emp_depths->{$total_depth} ?
      choose_weighted($rh_emp_depths->{$total_depth}->{depth}, $rh_emp_depths->{$total_depth}->{count}) : int(rand($total_depth));
    
    $mut_depth = $options{min_depth} if $mut_depth < $options{min_depth};
    

    ##
    ## get the relative position in each read, if this position contain a reference
    ## base, replace it with the mutation base, replicate this process until the
    ## desired detph of reads contain mutations is generated
    ##
    my $replaced_count = 0;
    my @replaced_reads = ();
    for (my $i=0; $i < @reads; $i++)
    {
        next if ($reads[$i] =~ /^@/ || $reads[$i] =~ /^\s+$/); ## skip header
        
        chomp($reads[$i]);
        
        my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
            $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/, $reads[$i]);
        
        ##
        ## reads contain indels may cause the replace position exceeds the read length,
        ## such case is simply skipped
        ##
        my $read_mut_pos = $mut_pos - $POS + 1;
        if ($read_mut_pos > 100) {
            next;
        }
        
        my $read_base = uc(substr($SEQ, $read_mut_pos-1, 1));
        my $mut_read  = $reads[$i];
        if ($read_base eq $ref_base && $replaced_count < $mut_depth) {
            substr($SEQ, $read_mut_pos-1, 1, $mut_base);
            
            $mut_read = join "\t", ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
                                    $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT);
            
            $replaced_count ++;
            
            print STDOUT "SAM:$mut_sample:$mut_region\t$mut_read\n";
        }
    }
    
    return -2 if $replaced_count < $options{min_depth};
    
    return "$mut_sample\t$total_depth:$mut_depth:$replaced_count";
}
