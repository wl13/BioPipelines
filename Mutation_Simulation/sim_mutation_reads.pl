#!/opt/nfs/bin/perl -w
#
#   sim_mutation_reads.pl -- Simulate reads with synthetic mutations use the method describled in Keightley et al. 2014, 2015.
#
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2019-07-12
#   Version: 2.4.0
#
#   Change logs:
#   Version 1.0.0 15/01/21: The initial version.
#   Version 1.0.1 15/07/30: Skip scaffold while sampling.
#   Version 2.0.0 16/05/22: Update: add options "--exclude", "--max-shared" and "--min-depth";
#                                   add support for simulate shared mutations;
#                                   add check for whether mutation is successfully generated;
#                                   now the number of simulated mutations were weighted according
#                                   to the chromosome length.
#   Version 2.0.1 16/05/23: Update: stop generate mutations on reads with indels or soft-clips;
#                                   rearrange output infos.
#   Version 2.1.0 16/07/16: Update: set the default value of "--min-depth" to 0 to suppress success
#                                   test by default; add option "--no-ref-n" to skip reference N sites;
#                                   fix various issues due to indels; remove test for nucleotide 
#                                   consistent compared to reference allele.
#   Version 2.1.1 16/07/18: Update: add option "--no-deletion"; now the deletions were not ignored by default.
#   Version 2.2.0 16/09/30: Update: add option "--group-file" to simulate group-specific mutations.
#   Version 2.2.1 17/01/11: Update: add option "--valid-depth" to output numbers of valid synthetic mutations.
#   Version 2.2.2 17/01/13: Bug fixed: use no less than instead of greater in valid depth test.
#   Version 2.3.0 19/07/09: Update: add options "--min-group-size" and "--max-group-size" to specify number of samples
#                                   could share the simulated mutations within each group.
#   Version 2.4.0 19/07/12: Update: add option "--use-filename".



use strict;

use Data::Dumper;
use Data::Random qw(:all);
use Getopt::Long;
use List::Util qw(shuffle);
use List::Util::WeightedChoice qw(choose_weighted);


use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.4.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my %options                  = ();
   $options{rand_size}       = 100;
   $options{out_format}      = 'fastq';
   $options{max_shared}      = 1;
   $options{min_depth}       = 0;
   $options{max_frac_in_del} = 0.8;
   $options{skip_deletions}  = 0;
   $options{valid_depth}     = 5;
   $options{min_group_size}  = 1;
   $options{max_group_size}  = -1;
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
            "no-ref-N"         => \$options{skip_ref_n_sites},
            "max-frac-del=f"   => \$options{max_frac_in_del},
            "no-deletion"      => \$options{skip_deletions},
            
            "group-file=s"     => \$options{group_file},
            "min-group-size=i" => \$options{min_group_size},
            "max-group-size=i" => \$options{max_group_size},
            
            "valid-depth=i"    => \$options{valid_depth},
            
            "use-filename"     => \$options{use_file_name},
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless($options{fasta_file} && $options{depth_file} && (@{$options{bam_files}} > 0 || $options{group_file}) && $show_help ) {
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
        bam file(s), at least one bam file should be specified unless the
        "--group-file" option is specified

    -g, --group-file  <filename>
        simulate group-specific mutations, e.g. mutations shared among samples
        within a certain group, this file should contain the absolute path to
        each bam file and the relevant group infos in the format:
            
            bam_file1 group1
            bam_file2 group1
            bam_file3 group2
            ...

    
Output Options:

    -o, --output <filename>
        output filename, default to STDOUT
        
    -n, --no-rc
        do not reverse complement sequence with negtive strand
    --use-rg
        add read group id to extracted records
    --use-filename
        use filename as sample id instead of SM tag in SAM file, only the
        string between "/" and the first "." of a filename would be used here

Simulation Options:
    --random-size <int>
        number of mutations to be simulated, default: 100

    --max-shared   <int>
        set this option over 1 to simulate shared mutations different from
        the "--group-file" option, where all samples would be chosen by random.
        If "--group-file" option is specified, this option would be ignored.
        [default: 1]
        
    --min-group-size <int>
    --max-group-size <int>
        Restrict the number of "mutated" samples within this range for each
        group, default: [1, all members]
    
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

        
    -e, --exclude  <strings>
        exclude unwanted chromosomes or scaffolds while simulating, all
        chromosomes with ids match strings specified here would be ignored 

    --min-depth    <int>
        minimum required depth of "mutated" reads, loci with a depth smaller
        than this threshold would be set to this threshold. This option would
        also cause those loci where no mutated reads were successful generated
        to be skipped [default: 0]
    
    --valid-depth  <int>
        add a LowDepth tag to sites where no sample meet this depth threshold,
        [default: 5]
    
    --no-ref-N
        do not generate mutations in reference N sites
        
    --no-deletion
        skip deletions while simulating
    --max-frac-del  <float>
        skip deletion loci in the simulated samples, as generating mutations in
        deletions is meaningless, this option determines whether a locus would
        be considered as deletion, the fraction was calculated by
        
            replaced reads with mutation in deletion / total replaced reads
        
        default: 0.8

EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


if ($options{output}) {
    open (STDOUT, "> $options{output}") || die $!;
}

if ($options{max_shared} > scalar @{$options{bam_files}}) {
    $options{max_shared} = scalar @{$options{bam_files}};
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
## read group infos
##
my %group_infos = ();
if ($options{group_file}) {
    get_group_info($options{group_file}, \%group_infos);
}

my @group_ids = sort keys %group_infos;


##
## randomly generate mutation loci
##
print STDERR ">> Start generating simulated reads ... ";
print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "##Tag(SAM): replaced reads\n";
print STDOUT "##Tag(MUT): simulated loci\n";
print STDOUT "#Tag\tChrom\tStart\tEnd\tSamples\tRef\tMut\tSimulated_depths(Non_Replaced,Replaced)\tSamples_Count(All,Depth>=$options{valid_depth})\n";
my $rand_index = 1;
while ($rand_index <= $options{rand_size})
{
    print STDERR "\r>> Start generating random mutations .. $rand_index \/ " . $options{rand_size};

    my $chrom = choose_weighted([@CHROM_IDs], [@chrom_weights]);
    my $pos   = int(rand($CHROM_LENGTHs{$chrom}+1));

    my $ref_base   = uc(substr($CHROM_SEQs{$chrom}, $pos-1, 1));
    
    next if ($options{skip_ref_n_sites} && $ref_base eq 'N');
    
    my @alt_bases  = grep { $_ ne $ref_base } qw(A T G C);
    my $mut_base   = $alt_bases[int(rand(3))];

    my @samples_selected = ();
    
    ##
    ## generate group-specific mutations
    ##
    if ($options{group_file}) {
        my @group_selected   = rand_set(set => \@group_ids, size => 1);
        
        my $max_group_size = scalar @{$group_infos{$group_selected[0]}};
        
        if ($options{max_group_size} >= $options{min_group_size}) {
            $max_group_size = $options{max_group_size};
        }
        
        if ($max_group_size > scalar @{$group_infos{$group_selected[0]}}) {
            $max_group_size = scalar @{$group_infos{$group_selected[0]}};
        }
        
        if ($max_group_size > $options{min_group_size}) {
            @samples_selected = rand_set(set => $group_infos{$group_selected[0]},
                                        min => $options{min_group_size},
                                        max => $max_group_size);
        }
        else {
            @samples_selected = rand_set(set => $group_infos{$group_selected[0]}, size => $max_group_size);
        }

    }
    else {
        @samples_selected = rand_set(set => $options{bam_files}, min => 1, max => $options{max_shared});
    }
    
    my @sim_samples     = ();
    my @sim_depths      = ();
    my @replaced_depths = ();
    for my $bam_file (@samples_selected)
    {
        my $sim_info = gen_mutated_reads($bam_file, $chrom, $pos, $ref_base, $mut_base, $rh_emp_depths);
        
        if ($sim_info ne -2) {
            my ($sample, $depth) = (split /\t/, $sim_info);
            
            push @sim_samples, $sample;
            push @sim_depths,  $depth;
            
            my $replaced_depth = (split /\,/, $depth)[1];
            push @replaced_depths, $replaced_depth;
        }
    }
    
    my $total_samples = scalar @sim_depths;
    my $valid_replace = grep {$_ >= $options{valid_depth}} @replaced_depths;
    
    if (@sim_samples >= 1) {
        $rand_index ++;
        
        
        ## sort output sampless
        my @sorted_index   = sort {$sim_samples[$a] cmp $sim_samples[$b]} (0..$#sim_samples);
        my @sorted_samples = @sim_samples[@sorted_index];
        my @sorted_depths  = @sim_depths[@sorted_index];
        
        my $mut_samples = join ";", @sorted_samples;
        my $mut_depths  = join ";", @sorted_depths;
        
        print STDOUT "MUT\t$chrom\t$pos\t$pos\t$mut_samples\t$ref_base\t$mut_base\t$mut_depths\t$total_samples,$valid_replace\n";
    }
}
print STDERR "\tdone!\n";



print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################

=head2 get_group_info

    About   : Get group infos of each sample
    Usage   : get_group_info($group_file);
    Args    : File contain group infos
    Returns : Null

=cut
sub get_group_info
{
    my ($in, $rh_group_infos) = @_;
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/\#/ || /^\s+$/);
        
        my ($sample_file, $group_id) = (split /\s+/);
        
        push @{$rh_group_infos->{$group_id}}, $sample_file;
    }
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
    
    if ($options{use_file_name}) {
        $mut_sample = $in;
        $mut_sample =~ s/.*\///;
        $mut_sample =~ s/\..*//;
    }
    
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
    
    if ($options{min_depth} > 0 && $total_depth < $options{min_depth}) {
        return -2;
    }
    
    my $mut_depth   = $rh_emp_depths->{$total_depth} ?
      choose_weighted($rh_emp_depths->{$total_depth}->{depth}, $rh_emp_depths->{$total_depth}->{count}) : int(rand($total_depth));
    
    if ($options{min_depth} > 0 && $mut_depth < $options{min_depth}) {
        $mut_depth = $options{min_depth};
    }
    
    ##
    ## get the relative position in each read, if this position contain a reference
    ## base, replace it with the mutation base, replicate this process until the
    ## desired detph of reads contain mutations is generated
    ##
    my $replaced_count    = 0;
    my @replaced_reads    = ();
    my $reads_in_deletion = 0;
    for (my $i=0; $i < @reads; $i++)
    {
        next if ($reads[$i] =~ /^@/ || $reads[$i] =~ /^\s+$/); ## skip header
        
        chomp($reads[$i]);
        
        my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
            $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/, $reads[$i]);
        
        my $read_mut_pos = $mut_pos - $POS + 1;
        
        ## reads contain indels or soft-clips may cause a failure in mapping original replace positions
        if ($CIGAR =~ /(I|D)/) {
            my $full_read_str  = '';
            my @ins_counts_all = ();
            
            while ($CIGAR =~ m/(\d+)(\w)/g)
            {
                $full_read_str .= $2 x $1;
                
                if ($2 eq 'I') {
                    push @ins_counts_all, $1;
                }
            }
            
            my $base_str = substr($full_read_str, $read_mut_pos-1, 1);
            my $bef_str  = substr($full_read_str, 0, $read_mut_pos-1);
            my $aft_str  = substr($full_read_str, $read_mut_pos);
            
            my $bef_del_count = ($bef_str =~ tr/D/D/);
            my $bef_ins_count = ($bef_str =~ tr/I/I/);
            
            ## recalculate the mutation position in the reads
            if ($base_str eq 'M' && ($bef_str =~ /(I|D)/)) {
                ###print "##corrected1a:$CIGAR\t$mut_pos\t$POS\t$read_mut_pos#\n";
                
                $read_mut_pos = $read_mut_pos - $bef_del_count + $bef_ins_count;
                
                ###print "##corrected1b:$read_mut_pos\t$full_read_str\t$bef_str\t$base_str\t$aft_str#\n";
            }
            elsif ($base_str eq 'I') {
                ###print "##corrected2a:$CIGAR\t$mut_pos\t$POS\t$read_mut_pos#\n";
                
                $aft_str =~ /(^I+)/;
                my $aft_ins_count = $1 ? length($1) : 0;
                
                $read_mut_pos = $read_mut_pos - $bef_del_count + $bef_ins_count + 1 + $aft_ins_count;
                
                ###print "##corrected2b:$read_mut_pos\t$full_read_str\t$bef_str\t$base_str\t$aft_str#\n";
            }
            elsif ($base_str eq 'D') {
                ###print "##skiped:$CIGAR\t$mut_pos\t$POS\t$read_mut_pos#\n";
                
                $reads_in_deletion ++; next;
            }
        }
        
        my $read_base = uc(substr($SEQ, $read_mut_pos-1, 1));
        my $mut_read  = $reads[$i];
        ###print "#1#$mut_read\n";
        if ($replaced_count < $mut_depth) {
            substr($SEQ, $read_mut_pos-1, 1, $mut_base);
            
            $mut_read = join "\t", ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
                                    $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT);
            
            $replaced_count ++;
            
            ###print STDOUT "#2#$mut_read\n";
            
            print STDOUT "SAM:$mut_sample:$mut_region\t$mut_read\n";
        }
    }
    
    if ($options{skip_deletions} && $reads_in_deletion >= 1 && $reads_in_deletion >= $options{max_frac_in_del} * $replaced_count) {
        return -2; ## skip deletions
    }
    
    if ($options{min_depth} > 0 && $replaced_count < $options{min_depth}) {
        return -2; ## skip low depth
    }
    
    my $non_replaced_depth = $total_depth - $replaced_count;
    
    return "$mut_sample\t$non_replaced_depth,$replaced_count";
}
