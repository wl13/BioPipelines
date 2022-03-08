#!/usr/bin/perl -w
#
#   annotate_mut_vcf.pl -- Annotate mutation sites with various metrics to predict their confidences.
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2021-11-06
#   Version: 2.2.0
#
#   Change logs:
#   Version 1.0.0 21/04/01: The initial version.
#   Version 2.0.0 21/04/06: Change this script to a more flexible mutation confidence annotating script.
#   Version 2.0.1 21/04/08: Fixed a bug when missing allele detail in some pre-exist variant sites.
#   Version 2.0.2 21/04/09: Add support for threading when processig bam files.
#   Version 2.1.0 21/04/09: Moved the bam processing to another script as it is too slow; change the inaccurate
#                           bam evidence to warn level.
#   Version 2.1.1 21/04/09: Revised the approach to estimate the sample distance.
#   Version 2.2.0 21/11/06: Add options to test co-exist variants.


use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.2.0';
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my %annotates                    = ();
my %threshold                    = ();
   $threshold{max_clip_ratio}    = 1;
   $threshold{min_var_score}     = 30;
   $threshold{min_called_frac}   = 0;
   $threshold{apply_filters}     = '';
   $threshold{max_pre_exist}     = 5;
   $threshold{max_cnts}          = "-1+-1";
   $threshold{remove_rules}      = "fatal-first";
   $threshold{max_sample_dist}   = -1;
my @appendix_outputs             = ();
my $min_supp_depth               = 5;
my ($input, $output, $out_stats, $use_indel_tag);
GetOptions(
            "input=s"           => \$input,
            
            "trf=s"             => \$annotates{tandem_regions},
            "surround-indel=s"  => \$annotates{surround_indel},
            "preexist-var=s"    => \$annotates{pre_exist_vars},
            "bam-supp=s"        => \$annotates{bam_support},
            "kinship=s"         => \$annotates{kinship},
            "coexist-vars=s"    => \$annotates{co_exist_vars},
            
            "output=s"          => \$output,
            "out-stats=s"       => \$out_stats,
            
            "use-indel-tag"     => \$use_indel_tag,
            
            "max-pre-exist=i"   => \$threshold{max_pre_exist},
            "min-co-exist=i"    => \$threshold{min_co_exist},
            "min-var-score=f"   => \$threshold{min_var_score},
            "min-call-frac=f"   => \$threshold{min_called_frac},
            "max-sample-dist=i" => \$threshold{max_sample_dist},
            "incl-grps"         => \$threshold{incl_grps},
            
            "out-appends=s{,}"  => \@appendix_outputs,
            
            "remove=s"          => \$threshold{max_cnts},
            "apply-rule=s"      => \$threshold{remove_rules},
            "apply-set=s"       => \$threshold{applied_sets},
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;


unless( $input && $show_help ) {
    print <<EOF;

$0 -- Annotate mutation sites with various metrics to predict their confidences.

Version: $VERSION

Usage:   perl $0 [options]

Options:

    --input    <filename>
        vcf-like file, perferred format by "detect_mutations.pl", expect 10
        rows with only one sample row, required
        

    --output       <filename>
        output filename, default to STDOUT
    
    --out-stats    <filename>
        output statistics of used annotates to this file
    
    --out-appends  <strings>
        write additional details of different annotations, supported tags are
        
        "trf", "bam_support", "sets_combined", "pre-exist", "sample_distance"


Options for INFO-based annotations
    
    By default, this script output below annotations:
        (1) allele frequencies in "called mutated group (Shared)" and "uncalled
            mutated group (GRPS)";
        (2) calling methods, e.g., UG+HC, or HC-only;
        (3) calling strategy, e.g., Topology-based (TP) or Frequency-based (FQ);
        (4) variant quality score;
        (5) combined sets of sources;
        (6) a confidence set will be marked if all annotates pass
        (7) user-defined combination of sources, i.e., "SET_*" tags
        
        
    --min-var-score <float>
        minimum require variant score, sites with score lower than this will
        be marked as LowQual [default: 30]
    
    --min-call-frac <float>
        minimum fractions of called samples among all focal-group samples with
        "mutated reads", calculated as "nCalled / (nCalled + nGRPS)", lower
        fraction suggests higher calling bias [default: 0]
    
      
Options for other annotations
        
    -b, --bam-supp <filename>
        file contains info of supporting reads for mutation allele in the
        format
        "CHROM	POS	Sample_ids	Mut_allele	Status:Readcounts"
        


    --preexist-var  <filename>
        give a file of pre-existing variants, and check whether the mutation
        allele is already present there, the basic format for this file is
        
        #CHROM  POS     Overall_Allele_Depth
        Pp01    139282  [GT]:5320,[G]:478
        Pp01    179607  [CCT]:3749,[C]:218
        Pp01    185182  [CAT]:1813,[C]:74
        ...
        
        The numbers in third low, e.g., "[GT]:5320,[G]:478" stand for the
        read-depths of each allele, the reference allele (or pre-mutated
        allele) should also be present here
    
    --max-pre-exist <int>
        check whether the mutation allele is identical to the prior existed
        variants provided in the appendix row, identical allele in appendix row
        with read-depth over this value will be considered as "pre-existing"
        [default: 5]


    --coexist-vars <filename>
        give a file of co-existing variants, and check whether the mutation
        allele is truly present there, the basic format for this file is
        
        #CHROM  POS     Overall_Allele_Depth
        Pp01    139282  [GT]:5320,[G]:478
        Pp01    179607  [CCT]:3749,[C]:218
        Pp01    185182  [CAT]:1813,[C]:74
        ...
        
        The numbers in third low, e.g., "[GT]:5320,[G]:478" stand for the
        read-depths of each allele, the reference allele (or pre-mutated
        allele) should also be present here
        
    --min-co-exist <int>
        check whether the mutation allele is identical to the co-existed
        variants provided in the appendix row, identical allele in appendix row
        with read-depth below this value will NOT be considered as "co-exist"
        [default: 5]
        
        
    --trf            <filename>
        a file contains the slippage information of each mutation predicted
        by trf, in the format:
        
            Pp01	453800	Pp01,453787,453809,GTGGGGTGTGTGGAGGTGTGTG([GTGGGGTGT]x2.3)
            Pp01	39762	N/A (for sites with no predicted tandem repeat)
            ...

    --surround-indel <filename>
        a file contains the surrounding information for each SNV mutation in
        the format:
        
            Pp01	7092	NO  (means no indel within a defined range)
            Pp01	20190	YES (means indel exists within a defined range)
            ...
            
    
    --kinship        <filename>
        a file contains belonging groups of all sequenced sample ids to infer
        the relationship of each sample, require formated group ids such as
        "B1-1-1- ..." which contains the branching information, e.g.,
        
        #Sample_ID      Group_ID
        DHQ1_B1-10-L1	B1-10
        DHQ1_B1-10-L2	B1-10
        DHQ1_B1-1-10-L1	B1-1-10
        DHQ1_B1-1-11-L1	B1-1-11
        ...
        
        so DHQ1_B1-10-L1 and DHQ1_B1-10-L2 both belong to B1-10, while
        DHQ1_B1-1-10-L1 and DHQ1_B1-1-11-L1 both belong to B1-1, the four
        samples belong to B1
        
    
    --max-sample-dist <int>
        maximum allowed distance among samples carry the same mutation, the
        distance represents for the relatedness of shared samples, which is
        estimated by sum up the numbers of non-mutated close-related samples
        untill reach the furthest-related sample, for example:
        
        for four samples B1-1-1, B1-1-2, B1-2-3, B2 sequenced, one mutation is
        called to be shared by B1-1-1 and B1-1-2, the distance is 0; for
        another mutation shared by B1-1-1 and B1-2-3, the distance is 1 (i.e.,
        since this mutation is earlier than B1-1 and B1-2, it is more likely
        B1-1-2 also carry this mutation)
        
        the larger the distance, the higher chance that far-related samples are
        more likely to share the mutation than the close-related samples, which
        indicates non-proper co-occurence
        [default: -1, i.e., no filtering]
    
    --incl-grps
        also considering low-reads supported samples from the same focal group
        (GRPS) when estimating the sample distance, and choose the smaller one
         
    
Options for applying filters based on the confidence annotations

    Annotations are assigned to different categories based on prior
    experiences as shown below:
        
    [FATAL level]
        (1) ref_allele: a reference allele is most likely a common
            allele;
        (2) proper_reads_only: mutation is only seen when calling without
            anomalous mapped reads;
        (3) high_mapq_only: mutation is only seen when calling with high
            mapping quality;
        (4) UG_misalign: variant only called by UnifideGenotyper in
            slippage/tandem region with nearby INDEL;
        (5) pre_exist: mutation allele is found in pre-existing variants;
        (6) unexpected_cooccur: mutation co-occurred preferentially in
            non-closely related samples rather than closely related samples
        (7) lack_co_exist: mutation allele is absent in external samples
            where it is expected to be present
        
    [WARN level]
        (1) anomalous_reads_only: mutation is only seen when calling with
            anomalous mapped reads;
        (2) low_mapq_only: mutation is only seen when calling with low mapping
            quality;
        (3) missing_calls: presence of one or more missing calls;
        (4) caller_bias: variant only called by part of callers;
        (5) low_qual: low variant quality score;
        (6) biased_scall: many samples in focal group have the "mutated
            reads" but only very few have sufficient read-depth;
        (7) indel_nearby: SNVs called with nearby INDELs are prone to
            mis-alignemnt;
        (8) unexpected_presence: one or more "mutated-reads" seen in
            non-mutated group;
        (9) clustered: mutations clustered within a small range is suspect
            to mapping/calling artefacts especially with "high_mapq_only",
            require "Cluster" tag in FILTER field
        (10) lack_reliable_reads: mutation allele is only supported by
            non-realible (e.g. clipped) reads;
            
    --remove     <string>
        remove sites with number of fatal- and warn- level annotates
        exceed these numbers, given in the format
        N_FATAL+N_WARN, e.g., 1+3 means at most 1 fatal annotate and 3
        warn annotates, set to -1 if no filtering is needed
        [default: -1+-1, i.e., no removing]
    
    --apply-rule <string>
        "fatal-first": when removing, first test the number of fatals,
            only when fatals equal the threshold then test warns,
            so 1+3 will not remove sites with 4 warns but 0 fatals
        "independent": independently compare two metrics, 1+3 will remove
            sites with either over 1 fatal or over 3 warns
            
        [default: fatal-first]
    
    --apply-set  <string>
        apply each annotate only to certain set(s), require "SET_*" tags,
        for example, "SET_1:unexpected_cooccur;SET_2:missing_calls" will
        only annotate variants from SET_1 with "unexpected_cooccur",
        variants from SET_2 with missing_calls, all other annotates are
        applied on all sets
        [default: apply all annotates to all variants]
        
EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

## parse BAM supportings
my %bam_supports   = ();
if ($annotates{bam_support}) {
    print STDERR ">> Start parsing $annotates{bam_support} ... ";
    
    my $bam_list_fh = getInputFilehandle($annotates{bam_support});
    while (<$bam_list_fh>)
    {
        next if (/^\#/ || /^\s+$/); ## skip header and blank lines
        
        my ($chrom, $pos, $sample_ids, $mut_allele, $supp_infos) = (split /\s+/);
        
        $bam_supports{$chrom}->{$pos} = $supp_infos;
    }
    
    print STDERR "done!\n";
}

## parse pre-existing variants
my %pre_exist_vars = ();
if ($annotates{pre_exist_vars}) {
    print STDERR ">> Start parsing $annotates{pre_exist_vars} ... ";
    
    my $pre_exist_fh = getInputFilehandle($annotates{pre_exist_vars});
    while (<$pre_exist_fh>)
    {
        next if (/^\#/ || /^\s+$/); ## skip header and blank lines
        
        my ($chrom, $pos, $allele_details) = (split /\s+/);
        
        next unless ($allele_details =~ /\:/);
        
        $pre_exist_vars{$chrom}->{$pos}->{info} = $allele_details;
        
        my @cmp_alleles = (split /\,/, $allele_details); ## parsing pre-existing alleles such as [G]:4123,[A]:1741
        
        for (@cmp_alleles)
        {
            $_ =~ /\[(\w+)\]\:(\d+)/;
            $pre_exist_vars{$chrom}->{$pos}->{depth}->{$1} = $2;
        }
    }
    
    print STDERR "done!\n";
}



## parse co-existing variants
my %co_exist_vars = ();
if ($annotates{co_exist_vars}) {
    print STDERR ">> Start parsing $annotates{co_exist_vars} ... ";
    
    my $co_exist_fh = getInputFilehandle($annotates{co_exist_vars});
    while (<$co_exist_fh>)
    {
        next if (/^\#/ || /^\s+$/); ## skip header and blank lines
        
        my ($chrom, $pos, $allele_details) = (split /\s+/);
        
        next unless ($allele_details =~ /\:/);
        
        $co_exist_vars{$chrom}->{$pos}->{info} = $allele_details;
        
        my @cmp_alleles = (split /\,/, $allele_details); ## parsing pre-existing alleles such as [G]:4123,[A]:1741
        
        for (@cmp_alleles)
        {
            $_ =~ /\[(\w+)\]\:(\d+)/;
            $co_exist_vars{$chrom}->{$pos}->{depth}->{$1} = $2;
        }
    }
    
    print STDERR "done!\n";
}



## parse tandem repeats
my %trf_infos = ();
if ($annotates{tandem_regions}) {
    print STDERR ">> Start parsing $annotates{tandem_regions} ... ";
    
    my $trf_fh = getInputFilehandle($annotates{tandem_regions});
    while (<$trf_fh>)
    {
        next if (/^\#/ || /^\s+$/); ## skip header and blank lines
        
        my ($chrom, $pos, $trf_predict) = (split /\s+/);
        
        $trf_predict =~ s/.*\,//;
        $trf_predict =~ s/\;.*//g;
        $trf_predict =~ s/\w+\(//;
        $trf_predict =~ s/\)//;
        
        $trf_infos{$chrom}->{$pos} = $trf_predict;
    }
    
    print STDERR "done!\n";
}

## parse surrounding indels
my %surround_indel = ();
if ($annotates{surround_indel}) {
    print STDERR ">> Start parsing $annotates{surround_indel} ... ";
    
    my $surr_indel_fh = getInputFilehandle($annotates{surround_indel});
    while (<$surr_indel_fh>)
    {
        next if (/^\#/ || /^\s+$/); ## skip header and blank lines
        
        my ($chrom, $pos, $is_indel_around) = (split /\s+/);
        
        $surround_indel{$chrom}->{$pos} = $is_indel_around;
    }
    
    print STDERR "done!\n";
}

## parse sample relationships
my %sample_kinships = ();
if ($annotates{kinship}) {
    get_group_info($annotates{kinship}, \%sample_kinships);
}

my $out_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION\tVar_type\tCallers\tMut_Type\tCalled_fq\tGRPS_sum\tConf_level";


## check how many additional infos should we add to the final results
my %appendix_outputs = ();
   $appendix_outputs{$_} = 1 for @appendix_outputs;



if ($appendix_outputs{surround_indel}) {
    $out_header .= "\tIs_indel_nearby";
}
if ($appendix_outputs{trf}) {
    $out_header .= "\tTRF_predicts";
}
if ($annotates{bam_support} && $appendix_outputs{bam_support}) {
    $out_header .= "\tBAM_supports";
}
if ($appendix_outputs{"pre-exist"}) {
    $out_header .= "\tPre-exist_alleles";
}
if ($appendix_outputs{"lack-co-exist"}) {
    $out_header .= "\tCo-exist_alleles";
}
if ($appendix_outputs{sample_distance}) {
    $out_header .= "\tSample_relatedness";
}
if ($appendix_outputs{sets_combined}) {
    $out_header .= "\tCombined_sets";
}

my ($max_fatal_cnt, $max_warn_cnt) = (split /\+/, $threshold{max_cnts});


## parse set rules
my %set_rules = ();
if ($threshold{applied_sets}) {
    my @applied_sets = (split /\;/, $threshold{applied_sets});
    
    for my $set_ann (@applied_sets)
    {
        my ($set, $ann) = (split /\:/, $set_ann);
        $set_rules{$ann}->{$set} = 1;
    }
}


print STDOUT "##source=$SOURCE $CMDLINE\n";
print STDOUT "$out_header\n";

my %ann_stats = ();

print STDERR ">> Start checking mutation allele in $input ... ";
my $mut_fh       = getInputFilehandle($input);
while (<$mut_fh>)
{
    next if (/^\#/ || /^\s+$/); ## skip header and blank lines
    
    chomp(my $inline = $_);
    
    my ($chrom, $pos, $sample_ids, $ref_allele, $alt_allele, $qual, $filter, $info, $format, $mut_sample, $append) = (split /\s+/, $inline);
    
    my $user_set = ($info =~ m/(SET\_.*?)(\(|\;|$)/) ? $1 : 'all'; ###print STDERR "user:$user_set\n";
    
    my @samples = split /\;/, $sample_ids;
    
    my ($mut_allele) = ($info =~ /MA=(\w+)/);  ## obtain the mutation allele from the INFO field
    
    my $var_type = (length($ref_allele) eq length($mut_allele)) ? "SNV" : "INDEL";
    
    if ($use_indel_tag && $append) {  ## if the appendix defines the variant type then use it
        $var_type = ($append =~ /INDEL/) ? "INDEL" : "SNV";
        $append =~ s/\:INDEL//;
    }
    
    my @warn_annotates  = ();
    my @fatal_annotates = ();
    my @append_infos    = ();
    
    ## check for clustered mutations
    if ($filter =~ /Cluster/) {
        if (!$threshold{applied_sets} || !$set_rules{clustered} || $set_rules{clustered}->{$user_set}) {
            push @warn_annotates, "clustered";
        }
    }

    ###print STDERR Dumper(%set_rules); exit;
    
    my $called_fq = ($info =~ /Shared\=(\d+)/)      ? $1 : 1;   ## number of samples called with the mutation allele
    my $GRPS_fqs  = ($info =~ /GRPFQ\=(.*?)(\;|$)/) ? $1 : 0;   ## number of samples have too few reads supporting the mutation allele
    my @GRPS_fqs  = (split /\,/, $GRPS_fqs);
    my $GRPS_sum  = 0; if (@GRPS_fqs > 0) {$GRPS_sum += $_ for @GRPS_fqs;}
    
    my $call_frac = $called_fq / ($called_fq + $GRPS_sum);
    
    if ($call_frac < $threshold{min_called_frac}) {
        if (!$threshold{applied_sets} || !$set_rules{biased_scall} || $set_rules{biased_scall}->{$user_set}) {
            push @warn_annotates, "biased_scall";
        }
    }
    
    ## check whether its an "backward mutation to reference allele"
    my $mut_type  = ($mut_allele eq $ref_allele) ? "REF" : "ALT";  

    if ($mut_type eq "REF") {
        ## a reference allele is often treated as a common allele
        if (!$threshold{applied_sets} || !$set_rules{ref_allele} || $set_rules{ref_allele}->{$user_set}) {
            push @fatal_annotates, "ref_allele";
        }
    }
    
    ## check source of callers
    my $caller_src = "UG_Single";
    if ($info =~ /UG_Single\+HC_GVCF/){
        $caller_src = "UG_Single+HC_GVCF";
    } elsif ($info =~ /HC_GVCF/){
        $caller_src = "HC_GVCF";
    }
    
    if ($caller_src ne "UG_Single+HC_GVCF") {
        ## variant only called by a single caller
        if (!$threshold{applied_sets} || !$set_rules{caller_bias} || $set_rules{caller_bias}->{$user_set}) {
            push @warn_annotates, "caller_bias";
        }
    }
    

    
    ## check for SNV mutations with surrounding indels
    if ($annotates{surround_indel}) {
        if ($surround_indel{$chrom}->{$pos} && ($surround_indel{$chrom}->{$pos} eq "YES")) {
            if (!$threshold{applied_sets} || !$set_rules{indel_nearby} || $set_rules{indel_nearby}->{$user_set}) {
                push @warn_annotates, "indel_nearby";
            }
        }
    }
    
    ## check for variant called in slippage/tandem region regions
    if ($annotates{tandem_regions}) {
        ###if ($trf_infos{$chrom}->{$pos}) {
        ###    if (!$threshold{applied_sets} || !$set_rules{slippage} || $set_rules{slippage}->{$user_set}) {
        ###        push @warn_annotates, "slippage";
        ###    }
        ###}
        
        if ($appendix_outputs{trf}) {
            if ($trf_infos{$chrom}->{$pos}) {
                push @append_infos, $trf_infos{$chrom}->{$pos};
            }
            else {
                push @append_infos, "N/A";
            }
        }
    }
    
    if (($caller_src eq "UG_Single") && ($trf_infos{$chrom}->{$pos}) && $surround_indel{$chrom}->{$pos}) {
        ## variant only called by UnifideGenotyper has higher false positive rate around INDELs/slippage regions
        if (!$threshold{applied_sets} || !$set_rules{UG_misalign} || $set_rules{UG_misalign}->{$user_set}) {
            push @fatal_annotates, "UG_misalign";
        }
    }

    
    ## check varaint score
    if ($qual < $threshold{min_var_score}) {
        if (!$threshold{applied_sets} || !$set_rules{low_qual} || $set_rules{low_qual}->{$user_set}) {
            push @warn_annotates, "low_qual";
        }
    }
    
    ## check for anomalous reads
    if ($info =~ /Combine=AR(\;|$)/) {
        ## mutation is only seen when calling with the anomalous mapped reads, indicating the mutation allele is supporting by anomalous reads only
        if (!$threshold{applied_sets} || !$set_rules{anomalous_reads_only} || $set_rules{anomalous_reads_only}->{$user_set}) {
            push @warn_annotates, "anomalous_reads_only";
        }
    }
    if ($info =~ /Combine=NAR(\;|$)/) {
        ## mutation is only seen when calling without anomalous mapped reads, indicating a putative false-positive call
        if (!$threshold{applied_sets} || !$set_rules{proper_reads_only} || $set_rules{proper_reads_only}->{$user_set}) {
            push @fatal_annotates, "proper_reads_only";
        }
    }
    
    ## check for low-quality mapping
    if ($info =~ /Combine=MQ0(\;|$)/) {
        ## mutation is only seen when calling with the low MAPQ reads, indicating the mutation allele is supporting by reads with low MAPQ only
        if (!$threshold{applied_sets} || !$set_rules{low_mapq_only} || $set_rules{low_mapq_only}->{$user_set}) {
            push @warn_annotates, "low_mapq_only";
        }
    }
    if ($info =~ /Combine=MQ20(\;|$)/) {
        ## mutation is only seen when calling with high MAPQ (>=20) reads, indicating a putative false-positive call
        if (!$threshold{applied_sets} || !$set_rules{high_mapq_only} || $set_rules{high_mapq_only}->{$user_set}) {
            push @fatal_annotates, "high_mapq_only";
        }
    }

    ## check for missing rate
    if ($info !~ /NMISS=0/) {
        ## if there is any sample missing, then mark it
        if (!$threshold{applied_sets} || !$set_rules{missing_calls} || $set_rules{missing_calls}->{$user_set}) {
            push @warn_annotates, "missing_calls";
        }
    }
    
    ## check for unexpected cooccurence
    if (($info !~ /FPD=0/) && ($info !~ /FPD=1;FPFQ=1(\;|$)/)) {
        ## if over one "mutated-reads" seen in unexpected group, then mark it
        if (!$threshold{applied_sets} || !$set_rules{unexpected_presence} || $set_rules{unexpected_presence}->{$user_set}) {
            push @warn_annotates, "unexpected_presence";
        }
    }

    
    ## check for mapping status in BAM file
    if ($annotates{bam_support}) {
        my $bam_support = $bam_supports{$chrom}->{$pos} ? $bam_supports{$chrom}->{$pos} : "N/A";

        if ($bam_support =~ "POOR_SUPPORT") {
            if (!$threshold{applied_sets} || !$set_rules{lack_reliable_reads} || $set_rules{lack_reliable_reads}->{$user_set}) {
                push @warn_annotates, "lack_reliable_reads";
            }
        }
        
        if ($appendix_outputs{bam_support}) {
            push @append_infos, $bam_support;
        }
    }
    
    
    ## check for pre-existing alleles
    if ($annotates{pre_exist_vars} && $pre_exist_vars{$chrom}->{$pos}->{depth}) {
        ## also confirm the reference allele in pre-existing alleles to make sure the "change" is the same
        if ($pre_exist_vars{$chrom}->{$pos}->{depth}->{$ref_allele} &&
            $pre_exist_vars{$chrom}->{$pos}->{depth}->{$mut_allele} &&
            $pre_exist_vars{$chrom}->{$pos}->{depth}->{$ref_allele} > 0 &&
            $pre_exist_vars{$chrom}->{$pos}->{depth}->{$mut_allele} >= $threshold{max_pre_exist}) {
            if (!$threshold{applied_sets} || !$set_rules{pre_exist} || $set_rules{pre_exist}->{$user_set}) {
                push @fatal_annotates, "pre_exist";
            }
        }
    }
    
    if ($appendix_outputs{"pre-exist"}) {
        if ($pre_exist_vars{$chrom}->{$pos}->{info}) {
            push @append_infos, $pre_exist_vars{$chrom}->{$pos}->{info};
        }
        else {
            push @append_infos, "N/A";
        }
    }
    
    
    ## check for co-existing alleles
    if ($annotates{co_exist_vars}) {
        if ($co_exist_vars{$chrom}->{$pos}->{depth} &&
            $co_exist_vars{$chrom}->{$pos}->{depth}->{$ref_allele} &&
            $co_exist_vars{$chrom}->{$pos}->{depth}->{$mut_allele} &&
            $co_exist_vars{$chrom}->{$pos}->{depth}->{$ref_allele} > 0 &&
            $co_exist_vars{$chrom}->{$pos}->{depth}->{$mut_allele} >= $threshold{min_co_exist}) {
        }
        else {
            if (!$threshold{applied_sets} || !$set_rules{lack_co_exist} || $set_rules{lack_co_exist}->{$user_set}) {
                push @fatal_annotates, "lack_co_exist";
            }
        }
    }
    
    if ($appendix_outputs{"lack-co-exist"}) {
        if ($co_exist_vars{$chrom}->{$pos}->{info}) {
            push @append_infos, $co_exist_vars{$chrom}->{$pos}->{info};
        }
        else {
            push @append_infos, "N/A";
        }
    }
    
    
    
    ## check for non-proper co-occurence of mutations
    if ($annotates{kinship}) {  ## hierarchic group counts
        my ($sample_dist, $out_grp_counts) = estimate_sample_distance(\@samples, \%sample_kinships);

        if ($threshold{incl_grps} && ($info =~ /GRPS\=(.*?)(\;|$)/)) {
            my $grps_samples = $1;
               $grps_samples =~ s/\(//g;
               $grps_samples =~ s/\)//g;
            my @add_samples = (split /\,/, $grps_samples);
            
            my ($sample_dist2, $out_grp_counts2) = estimate_sample_distance([@samples, @add_samples], \%sample_kinships);
            
            if ($sample_dist2 < $sample_dist) {  ## choose the smaller one as the true distance
                $sample_dist    = $sample_dist2;
                $out_grp_counts = $out_grp_counts2;
            }
        }
        
        if ($threshold{max_sample_dist} >= 0 && $sample_dist > $threshold{max_sample_dist}) {
            if (!$threshold{applied_sets} || !$set_rules{unexpected_cooccur} || $set_rules{unexpected_cooccur}->{$user_set}) {
                push @fatal_annotates, "unexpected_cooccur";
            }
        }
        
        if ($appendix_outputs{sample_distance}) {
            push @append_infos, "$sample_dist|$out_grp_counts";
        }
    }
    
    ## check calling strategy and combined set
    my $conf_level = "TP";
    if ($info =~ /Combine=Grouped\+NonGrouped/) {
        $conf_level = "TP+FQ";
    } elsif ($info =~ /Combine=NonGrouped/) {
        $conf_level = "FQ";
    }
    
    my $warn_cnt  = scalar @warn_annotates;
    my $fatal_cnt = scalar @fatal_annotates;
    
    if ($warn_cnt + $fatal_cnt == 0) {
        $conf_level .= "(Confidence)";
    } else {
        my $fatal_tags = join ",", @fatal_annotates;
        my $warn_tags  = join ",", @warn_annotates;
        
        if ($fatal_cnt > 0 && $warn_cnt > 0) {
            $conf_level .= "(FATAL:$fatal_cnt|$fatal_tags;WARN:$warn_cnt|$warn_tags)";
        }
        elsif ($fatal_cnt > 0) {
            $conf_level .= "(FATAL:$fatal_cnt|$fatal_tags)";
        }
        elsif ($warn_cnt > 0) {
            $conf_level .= "(WARN:$warn_cnt|$warn_tags)";
        }
    }
    
    $ann_stats{"$user_set\t$var_type\t$conf_level"}->{"non-filtered"} ++;
    
    if ($threshold{remove_rules} eq "fatal-first") {
        next if ($max_fatal_cnt >= 0 && $fatal_cnt > $max_fatal_cnt);
        next if ($fatal_cnt == $max_fatal_cnt && $max_warn_cnt  >= 0 && $warn_cnt  > $max_warn_cnt);
    }
    elsif ($threshold{remove_rules} eq "independent") {
        next if ($max_fatal_cnt >= 0 && $fatal_cnt > $max_fatal_cnt);
        next if ($max_warn_cnt  >= 0 && $warn_cnt  > $max_warn_cnt);
    }
    
    $ann_stats{"$user_set\t$var_type\t$conf_level"}->{"post-filtered"} ++;
    
    if ($appendix_outputs{sets_combined}) {
        my @combines = ($info =~ m/(SET\_.*?)(\(|\;|$)/);
           @combines = grep {/SET/} @combines;
           
        my $combine = join ",", @combines;
        
        ###print STDERR "#@combines#$info#$combine#\n"; exit;
        
        push @append_infos, $combine;
    }
    
    my $append_infos = join "\t", @append_infos;
    
    my $outline = "$inline\t$var_type\t$caller_src\t$mut_type\t$called_fq\t$GRPS_sum\t$conf_level\t$append_infos";
    
    print STDOUT "$outline\n";
}
print STDERR "done!\n";


if ($out_stats) {
    open (STATS, "> $out_stats") || die $!;
    print STATS "##source=$SOURCE $CMDLINE\n";
    print STATS "#user_set\tvar_type\tconf_level\tnon-filtered\tpost-filtered\n";
    
    for my $state (sort keys %ann_stats)
    {
        my $non_flt_cnt  = $ann_stats{$state}->{"non-filtered"}  ?
                           $ann_stats{$state}->{"non-filtered"}  : 0;
        my $post_flt_cnt = $ann_stats{$state}->{"post-filtered"} ?
                           $ann_stats{$state}->{"post-filtered"} : 0;
        
        print STATS "$state\t$non_flt_cnt\t$post_flt_cnt\n";
    }
    close STATS;
}



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
        
        my ($sample_id, $group_id) = (split /\s+/);
        
        $rh_group_infos->{ID}->{$sample_id} = $group_id;
        
        if ($annotates{kinship}) {
            my @group_levels = (split /\-/, $group_id);
            
            for (my $i=0; $i<@group_levels; $i++)
            {
                my $new_gid = join "\-", @group_levels[0..$i];
                
                $rh_group_infos->{NUM}->{$new_gid} ++;
                $rh_group_infos->{LVL}->{$i}->{$new_gid} ++;
            }
        }
        else {
            $rh_group_infos->{NUM}->{$group_id} ++;
        }
    }
}


=head2 estimate_sample_distance

    About   : Estimate distances between shared samples
    Usage   : estimate_sample_distance($ra_samples, $rh_sample_kinships);
    Args    : Reference to array of shared samples;
              Reference to hash of sample kinship.
    Returns : Estiamted sample distance;
              Details of traversed branch levels.

=cut
sub estimate_sample_distance
{
    my ($ra_samples, $rh_sample_kinships) = @_;
    
    my %group_counts = ();
    my $sample_num   = scalar @{$ra_samples};
    
    for my $id (@{$ra_samples})
    {
        ## parse the belonging groups of the sample
        my $group_id = $rh_sample_kinships->{ID}->{$id};
        
        ###unless($group_id) {print STDERR "$id\n"; exit;}
        
        my @group_levels = (split /\-/, $group_id);  
        
        for (my $i=0; $i<@group_levels; $i++)
        {
            my $new_gid = join "\-", @group_levels[0..$i];
            
            $group_counts{$i}->{$new_gid} ++;
        }
    }
    
    my @grp_counts  = ();
    my @grp_lvls    = sort {$b <=> $a} keys %group_counts;
    for my $lvl (@grp_lvls)   ## count from higher to lower levels
    {
        my @lvl_counts = ();
        my $lvl_sum    = 0;
        my @gids       = sort keys %{$group_counts{$lvl}};
        for my $gid (@gids)                                       ## count different branches at the same level
        {
            my $total_cnt = $rh_sample_kinships->{NUM}->{$gid};   ## obtain the total samples in this branch level
            
            $lvl_sum += $group_counts{$lvl}->{$gid};
            
            push @lvl_counts, "$gid:$group_counts{$lvl}->{$gid}/$total_cnt";
        }
        
        next unless(@lvl_counts > 0);
        
        my $lvl_counts = join "\,", @lvl_counts;
        
        push @grp_counts, "Lv" . ($lvl+1) . "($lvl_counts)";
        
        last if (($lvl_sum == $sample_num) && (@lvl_counts == 1));  ## stop climbing down if all mutated samples are collapsed to the most recent common level
    }

    
    my $sample_dist     = 0;
    my @out_grp_counts  = ();
    if (scalar @grp_counts > 1) {
        ## calculate the sample distance as the samples missed before collapse to the most recent common level
        my @mrc_groups = (split /\,/, $grp_counts[-2]);
        
        for my $cnt (@mrc_groups)
        {
            $cnt =~ /.*\:(\d+)\/(\d+)/;
            $sample_dist += ($2 - $1);
        }
        
        ## eliminate saturated level
        for (my $i=0; $i+1<@grp_counts; $i++)
        {
            my $lost_cnt = 0;
            my @groups   = (split /\,/, $grp_counts[$i]);
            
            for my $cnt (@groups)
            {
                $cnt =~ /.*\:(\d+)\/(\d+)/;
                $lost_cnt += ($2 - $1);
            }
            
            if ($lost_cnt > 0) {
                push @out_grp_counts, $grp_counts[$i];
            }
        }
    }
    
    push @out_grp_counts, $grp_counts[-1];
    
    my $out_grp_counts = join "\;", @out_grp_counts;
    
    return ($sample_dist, $out_grp_counts);
}