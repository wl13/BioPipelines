#!/usr/bin/perl -w
#
#   sim_mutation_reads.pl -- Simulate reads with synthetic mutations use the method describled in Keightley et al. 2014.
#
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2015-07-30
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 15/01/21: The initial version.
#   Version 1.0.1 15/07/30: Skip scaffold while sampling.



use strict;

use Data::Dumper;
use Getopt::Long;
use List::Util qw(shuffle);
use List::Util::WeightedChoice qw(choose_weighted);


use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my %options             = ();
   $options{rand_size}  = 100;
   $options{out_format} = 'fastq';
my ($no_rc, $extend_size, $use_rg_id,
    $samtools_opts, $min_seq_len, $max_clipping, $min_insert_size, $max_insert_size);
GetOptions(
            "fasta=s"          => \$options{fasta_file},
            "depth=s"          => \$options{depth_file},
            "random-size=i"    => \$options{rand_size},
            
            "extend=i"         => \$extend_size,
            
            "bams=s{,}"        => \@{$options{bam_files}},
            
            "output=s"         => \$options{output},
            
            "samtools=s"       => \$samtools_opts,
            "min-insert=i"     => \$min_insert_size,
            "max-insert=i"     => \$max_insert_size,
            "min-len=i"        => \$min_seq_len,
            "max-clipping=i"   => \$max_clipping,
            "use-rg"           => \$use_rg_id,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $options{fasta_file} && $options{depth_file} && (@{$options{bam_files}} > 0) && $show_help ) {
    print <<EOF;

$0  -- Extract all reads with the 

Version: $VERSION

Usage:   perl $0 [options]

Input Options:

    --fasta  <filename>
        reference sequence file in fasta format, required
    
    --depth  <filename>
        input file contain empirical distribution of read depths, required

    -b, --bams   <filename>
        bam file(s), at least one bam file should be specified


Output Options:

    -o, --output <filename>
        output filename, default to STDOUT
        
    -n, --no-rc
        do not reverse complement sequence with negtive strand
    -u, --use-rg
        add read group id to extracted records
        
Filtering Options:
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

EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


if ($options{output}) {
    open (STDOUT, "> $options{output}") || die $!;
}


##
## step1: generate random mutation loci in the genome
##
print STDERR ">> Start generating random mutations ...";
my $rh_mut_loci = gen_mut_loci($options{fasta_file}, $options{rand_size});
print STDERR "done!\n";


##
## step2: sampling the numbers of nonreference reads from realistic distributions
##
print STDERR ">> Start sampling read depths ...";
my $rh_emp_depths = get_empirical_dist($options{depth_file}, $options{rand_size});
print STDERR "done!\n";


##
## step3: change reads in bam files
##
print STDERR ">> Start generating simulated reads ... ";
print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "##Tag(SAM): replaced reads\n";
print STDOUT "##Tag(MUT): simulated loci\n";
print STDOUT "#Tag\tChrom\tStart\tEnd\tSample\tRef\tMut\tTotal_depth\tSim_depth\tReplaced_reads\n";
my $index = 0;
my $samples_count = scalar @{$options{bam_files}};
for my $chrom (sort keys %{$rh_mut_loci})
{
    for my $pos (sort {$a <=> $b} keys %{$rh_mut_loci->{$chrom}})
    {
        ++$index;
        
        my ($ref_base, $mut_base) = (split /\t/, $rh_mut_loci->{$chrom}->{$pos});
        
        my $sample_selected = int(rand($samples_count));
        
        gen_mutated_reads($options{bam_files}->[$sample_selected], "$chrom:$pos-$pos",
                          $ref_base, $mut_base, $rh_emp_depths);
    }
}
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################


=head2 gen_mut_loci

    About   : Generate random positions.
    Usage   : gen_mut_loci($file, $rand_size);
    Args    : Reference sequences in fasta format;
              Number of positions needed to pick out.
    Returns : Randomly generated mutations loci.

=cut
sub gen_mut_loci
{
    my ($in, $rand_size) = @_;
    
    ##
    ## read chromosome sequences
    ##
    my %chrom_seqs = ();
    parse_fasta_SEQs(\%chrom_seqs, $in);
    
    my %chrom_lengths = ();
    
    for my $chrom (sort keys %chrom_seqs)
    {
        next if ($chrom =~ /(mitochondria|chloroplast|scaffold)/);
        
        $chrom_lengths{$chrom} = length $chrom_seqs{$chrom};
    }
    
    
    my @chrom_ids = sort keys %chrom_lengths;
    my $chrom_num = scalar @chrom_ids;
    
    ##
    ## randomly generate mutation loci
    ##
    my %rand_mutations = ();
    for (my $i=0; $i<$rand_size; $i++)
    {
        my $chrom = @chrom_ids[int(rand($chrom_num))];
        my $pos   = int(rand($chrom_lengths{$chrom}+1));
        
        my $ref_base  = uc(substr($chrom_seqs{$chrom}, $pos-1, 1));
        my @alt_bases = grep { $_ ne $ref_base } qw(A T G C);
        
        ###print STDERR "$ref_base\t#@alt_bases#\n";exit;
        
        my $mut_base  = $alt_bases[int(rand(3))];
        
        ###print "$chrom\t$pos\t$ref_base\t$mut_base\n";
        
        $rand_mutations{$chrom}->{$pos} = "$ref_base\t$mut_base";
    }

    return \%rand_mutations;
}


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
    Usage   : gen_mutated_reads($bam_file, $region, $ref_seq, $mut_seq, $rh_emp_depths);
    Args    : Bam file to extract reads;
              Region to be replaced;
              Reference sequence in this region;
              Replace sequence in this region;
              Hash reference to empirical distributions of read depths.
    Returns : Null.

=cut
sub gen_mutated_reads
{
    my ($in, $region, $ref_seq, $mut_seq, $rh_emp_depths) = @_;
    
    my ($mut_chrom, $mut_start, $mut_end) = ($region =~ /(.*?)\:(\d+)\-(\d+)/);
    
    my $mut_len = $mut_end - $mut_start + 1;
    
    
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
    my $pipe_str = "samtools view $in $region |";
    
    if ($samtools_opts) {
        $pipe_str = "samtools view $samtools_opts $in $region |";
    }
    
    open (my $fh, $pipe_str) || die $!;
    my @reads = shuffle <$fh>;
    close $fh;
    
    ###print "$in\t$region\t@reads\n";exit;
    
    ##
    ## randomly generate mutation depth based on empirical distribution,
    ## note some regions would have total depth not present in previous sampled
    ## distributions, a random number is used despite the empirical distribution
    ##
    my $total_depth = scalar @reads;
    my $mut_depth   = $rh_emp_depths->{$total_depth} ?
      choose_weighted($rh_emp_depths->{$total_depth}->{depth}, $rh_emp_depths->{$total_depth}->{count}) : int(rand($total_depth));
    
    ###print join '', @reads;
    
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
        
        ###$reads[$i] =~ /RG:Z:(.*?)\s+/;
        ###$QNAME =~ s/\/\d$//;
        
        
        ##
        ## reads contain indels may cause the replace position exceeds the read length,
        ## such case is simply skipped
        ##
        my $read_region_start = $mut_start - $POS + 1;
        if ($read_region_start > 100) {
            #print STDERR "ERROR\t$mut_chrom\t$mut_start\t$mut_sample\t$ref_seq\t$mut_seq\t$total_depth\t$mut_depth\n";
            #print STDERR "$read_region_start\t$reads[$i]\n";exit;
            next;
        }
        
        my $read_region_seq   = substr($SEQ, $read_region_start-1, $mut_len);
        
        ###print "$reads[$i]\n";

        my $mut_read = $reads[$i];
        if ($read_region_seq eq $ref_seq && $replaced_count < $mut_depth) {
            substr($SEQ, $read_region_start-1, $mut_len, $mut_seq);
            
            $mut_read = join "\t", ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
                                    $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT);
            
            $replaced_count ++;
            
            ###print "$mut_read\n";
            ###print "$read_region_start\n";
            
            ###push @replaced_reads, "$QNAME\t$FLAG\t$SEQ";
            
            print STDOUT "SAM:$mut_sample:$region\t$mut_read\n";
        }
    }
    
    print STDOUT "MUT\t$mut_chrom\t$mut_start\t$mut_end\t$mut_sample\t$ref_seq\t$mut_seq\t$total_depth\t$mut_depth\t$replaced_count\n";
    
    ###return [@replaced_reads];
}




__END__

=head2 gen_mut_pos

    About   : Generate random positions.
    Usage   : gen_mut_pos($file, $rand_size);
    Args    : Bam file contain SQ tags in header;
              Number of positions needed to pick out.
    Returns : Randomly picked read depths.

=cut
sub gen_mut_pos
{
    my ($in, $rand_size) = @_;
    
    
    my %chrom_lengths = ();
    my $pipe_str = "samtools view -H $in |";
    open (my $fh, $pipe_str) || die $!;
    while (<$fh>)
    {
        if (/^\@SQ/) {
            my ($SQ, $SN, $LN) = (split /\s+/);
            
            next if ($SN =~ /(mitochondria|chloroplast)/);
            
            my $chrom  = (split /\:/, $SN)[1];
            my $length = (split /\:/, $LN)[1];
            
            ###print "$chrom\t$length\n";
            
            $chrom_lengths{$chrom} = $length;
        }
    }
    
    my @chrom_ids = sort keys %chrom_lengths;
    my $chrom_num = scalar @chrom_ids;
    
    my @rand_pos = ();
    
    for (my $i=0; $i<$rand_size; $i++)
    {
        my $chrom = @chrom_ids[int(rand($chrom_num))];
        my $pos   = int(rand($chrom_lengths{$chrom}+1));
        
        ###print "$chrom\t$pos\n";
        
        push @rand_pos, "$chrom:$pos";
    }

    return [@rand_pos];
}


=head2 get_ref_base

    About   : Get the reference base in certain positions
    Usage   : get_ref_base($file, $ra_pos);
    Args    : Reference sequences in fasta format;
              Array reference to query regions in the format "chr:pos"
    Returns : Array reference of extracted reference bases.

=cut
sub get_ref_base
{
    my ($in, $ra_pos) = @_;
    
    my %ref_seqs = ();
    parse_fasta_SEQs(\%ref_seqs, $in);
    
    my %ref_bases = ();
    for my $region ($ra_pos)
    {
        my ($chrom, $pos) = (split /\:/, $region);
        
        my $ref_base = substr($ref_seqs{$chrom}, $pos - 1, 1);
        
        $ref_bases{$region} = $ref_base;
        
        print "$chrom\t$pos\t$ref_base\n";
    }
    
    return \%ref_bases;
}
