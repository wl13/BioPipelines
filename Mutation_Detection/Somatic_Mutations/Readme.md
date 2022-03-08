## Pipelines and scripts to identify and analysis somatic mutations in plants.

<br />

### somatic-mutation-detection.md

An adaption of previous pipeline to deal with somatic mutations

<br />


<br />

## Miscellaneous scripts

<br />

#### fillVcfDepth-split.pl

> An alternative version of fillVcfDepth.pl which proceeds chromosome by chromosome to save memory (can be slow if the reference draft contain many scaffolds/contigs)

* Options:

        -v, --vcf    <filename>
            input vcf file, required
        -o, --output <filename>
            output filename, default to STDOUT

        -l, --list   <filename>
            file contains list of depth files and related sample names in the
            format:

            sample_id library_id filename

            delimited by tab(s) or space(s), required

            *Note:
            All files specified in this list should be compressed with bgzip and
            indexed using tabix:
            bgzip example.readcounts
            tabix -S1 -s1 -b2 -e2 example.readcounts.gz

        -f, --filter <strings>
            skip filter loci, can have multiple values, separate by blanks, e.g.
            "LowQual SNPFilter" ...    [default: no filtering]
            -p, --phased
            skip unphased sites

        -t, --type   <string>
            set "snp" to process snp sites only, or set "indel" to process indels
            only

        -u, --update-AD
            update AD field according to the results of readcounts

        -m, --minimum-vcf
            remove original INFO and FORMAT fields



<br />

#### test_uniformity.pl
> D (Variance/Mean) test by comparing observed D and simulated D\
> Code for random simulation was borrowed from the post in https://stackoverflow.com/a/22381248 by ikegami (the second approach) and modified to allow for zero-values

* **Options:**   

            -q, --query-file  <string>
                file contain list of numbers, required

            -o, --output <filename>
                output file, default to STDOUT

            -t, --times
                random times [default: 1]
<br />




#### check_bam_supports.pl
> This script checks whether a mutation allele was supported by reliable reads by excluding reads failed various thresholds, multiple thresholds will be applied together if specified
> Note: current version does not handle reads with indels well, so if reads supporting the mutation alleles also contain indels, these reads could be missed

* **Options:**   

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
<br />




#### annotate_mut_vcf.pl
> Annotate mutation sites with various metrics to predict their confidences

* **Options:**   

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

<br />

## Publications

**More detailed descriptions of this pipeline could also be found in the below publications:**

[1] L. Wang, Y. Ji, Y. Hu, H. Hu, X. Jia, M. Jiang, X. Zhang, L. Zhao, Y. Zhang, Y. Jia, C. Qin, L. Yu, J. Huang, S. Yang, L.D. Hurst, D. Tian, The architecture of intra-organism mutation rate variation in plants, <i>PLOS Biology</i>. 17 (2019) e3000191. doi:10.1371/journal.pbio.3000191.

[2] Ren, Y., He, Z., Liu, P., Traw, B., Sun, S., Tian, D., Yang, S., Jia, Y., and Wang, L. (2021). Somatic Mutation Analysis in Salix suchowensis Reveals Early-Segregated Cell Lineages. <i>Molecular Biology and Evolution</i> 38: 5292–5308.

