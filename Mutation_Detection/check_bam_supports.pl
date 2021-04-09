#!/usr/bin/perl -w
#
#   check_bam_supports.pl -- Annotate mutation allele with the mapping details of its supporting reads.
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2021-04-09
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 21/04/01: The initial version.
#   Version 1.0.1 21/04/09: A few updates on input and output.



use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my %threshold                 = ();
   $threshold{max_clip_ratio} = 1;
my $min_supp_depth            = 5;
my ($input, $bam_list, $output, $samtools_opts);
GetOptions(
            "input=s"          => \$input,
            "bam-list=s"       => \$bam_list,

            "output=s"         => \$output,
            
            "samtools=s"       => \$samtools_opts,
            
            "no-both-clips"    => \$threshold{no_both_clips},
            "max-clip-ratio=f" => \$threshold{max_clip_ratio},
            
            "min-supp-depth=i" => \$min_supp_depth,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $input && $bam_list && $show_help ) {
    print <<EOF;

$0  -- Annotate mutation allele with the mapping details of its supporting reads.

This script checks whether a mutation allele was supported by reliable reads
by excluding reads failed various thresholds here, multiple thresholds will
be applied together if specified

*Note: current version does not handle reads with indels well, so if reads
supporting the mutation alleles also contain indels, these reads could be
missed

Version: $VERSION

Usage:   perl $0 [options]

Options:

    -i, --input    <filename>
        vcf-like file, perferred format by "detect_mutations.pl", expect 10
        rows with only one sample row, the 3rd row is the sample ids, while
        the INFO field contains the "MA" tag for mutation allele, required

    -b, --bam-list <filename>
        list of bam files correspond to each sample id, in the format:
        "Sample_id BAM_file_location", required
        
    -o, --output   <filename>
        output filename, default to STDOUT

    -s, --samtools <string>
        options passed to samtools view, e.g. "-f 4 -F 8", use this option to
        remove reads that do not want to confirm the filters below, e.g.
        unmapped/supplementary/secondary reads

    --max-clip-ratio <float>
        maximum ratio of soft/hard clipped length, calculated as
        (soft-clipped-bases + hard-clipped-bases) / read_length
        default: 1 (i.e., no filtering)
    
    -n, --no-both-clips  <string>
        only filter reads with both ends clipped when exceeding
        the --max-clip-ratio, default: no filtering

        
    --min-supp-depth <int>
        minimum requirement of realible supporting read-depth for each sample,
        default: 5

EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}



print STDERR ">> Start parsing $bam_list ... ";
my %bam_files   = ();
my $bam_list_fh = getInputFilehandle($bam_list);
while (<$bam_list_fh>)
{
    next if (/^\#/ || /^\s+$/); ## skip header and blank lines
    
    my ($sample_id, $file_path) = (split /\s+/)[0,1];
    
    $bam_files{$sample_id} = $file_path;
}
print STDERR "done!\n";


print STDOUT  "#CHROM\tPOS\tSAMPLES\tREF\tMUT\tSupporting_Reads\n";

print STDERR ">> Start checking mutation allele in $input ... ";
my $mut_fh       = getInputFilehandle($input);
while (<$mut_fh>)
{
    next if (/^\#/ || /^\s+$/); ## skip header and blank lines
    
    my ($chrom, $pos, $sample_ids, $ref_allele, $alt_alleles, $qual, $filter, $info, $format, $mut_sample, $append) = (split /\s+/);
    
    
    ## if the mutation allele from the INFO field
    next unless ($info =~ /MA=(\w+)/);
    my $mut_allele = $1;  
    
    
    my @samples = split /\;/, $sample_ids;
    
    my @s_reads_cnts = ();
    my $reliable_cnt = 0;
    for my $sample (@samples)
    {
        my $s_reads_cnt = count_support_reads($sample, $chrom, $pos, $mut_allele, \%bam_files, \%threshold);
        push @s_reads_cnts, $s_reads_cnt;
        
        $reliable_cnt++ if $s_reads_cnt >= $min_supp_depth;
    }
    
    my $s_reads_cnts = join "\;", @s_reads_cnts;
    
    my $state = ($reliable_cnt < 1) ? "POOR_SUPPORT" : "QUALIFIED_SUPPORT";
    
    print STDOUT "$chrom\t$pos\t$sample_ids\t$mut_allele\t$state\:$s_reads_cnts\n";
}
print STDERR "done!\n";



print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 count_support_reads

    About   : Count post-filtered reads which support the mutation allele.
    Usage   : count_support_reads($sample, $chrom, $var_pos, $var_allele, $rh_bam_files, $rh_threshold, $ra_read_filters);
    Args    : File in SAM format;
              Prefix of output files.
    Returns : Null

=cut
sub count_support_reads
{
    my ($sample, $chrom, $var_pos, $var_allele, $rh_bam_files, $rh_threshold) = @_;

    my $pipe_str = "samtools view $rh_bam_files->{$sample} $chrom:$var_pos-$var_pos |";
    
    if ($samtools_opts) {
        $pipe_str = "samtools view $samtools_opts $rh_bam_files->{$sample} $chrom:$var_pos-$var_pos |";
    }
    
    my $supp_reads = 0;
    
    open (my $fh, $pipe_str) || die $!;
    while (<$fh>)
    {
        next if (/^@/ || /^\s+$/); ## skip header
        
        chomp(my $record = $_);
        
        my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
            $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/, $record);
        
        
        ## parse soft/hard clips in CIGAR string
        ## original codes from: https://github.com/tseemann/samclip/blob/master/samclip
        if (($rh_threshold->{max_clip_ratio} < 1) && ($CIGAR =~ /\d[SH]/)) {
            my($HL, $SL, undef, $SR, $HR) = ($CIGAR =~ m/^(?:(\d+)H)?(?:(\d+)S)?(.*?)(?:(\d+)S)?(?:(\d+)H)?$/x);
            
            $HL ||= 0;
            $SL ||= 0;
            $SR ||= 0;
            $HR ||= 0;
            
            my $clip_ratio = ($HL + $SL + $SR + $HR) / length($SEQ);
            
            if ($clip_ratio > $rh_threshold->{max_clip_ratio}) {           ## check if the read is overclipped
                if ($threshold{no_both_clips}) {
                    next if ((($HL + $SL) > 0) && (($SR + $HR) > 0));  ## check if the clips happen on both ends
                }
                else {
                    next;
                }
            }
        }
        
        ###print STDERR "$var_pos\t$RNAME\t$POS\n" if ($var_pos - $POS > length($SEQ));
        if ($var_pos - $POS > length($SEQ)) { ## skip reads with preceding INDELs
            next;
        }
        
        ## test whether this read supports the mutation allele
        my $bases = substr($SEQ, $var_pos - $POS, length($var_allele));
        
        ###my $flank = substr($SEQ, $var_pos - $POS-2, length($var_allele)+4);
        ###print STDOUT "$sample\t$QNAME\t$RNAME\t$POS\t$flank\t$bases\t$var_allele\n";
        
        next unless($bases eq $var_allele);
        
        $supp_reads ++;
        
        ###next if ($rh_threshold->{min_insert_size} && abs($TLEN) < $rh_threshold->{min_insert_size});
        ###next if ($rh_threshold->{max_insert_size} && abs($TLEN) > $rh_threshold->{max_insert_size});
    }
    
    return $supp_reads;
}


