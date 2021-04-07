### Pipelines to detect somatic mutations among collected samples

<br />

Some details which have been descibed in previous section are not mentioned here. Please refer to the README.md in parental Mutation_Detection folder for those details.

<br />

## Workflow

### Step1: Generate initial candidate targets

* Use vcf_process.pl to screen out **all candidate sites**. Theoretically, every genomic position with a non-fixed variant could be a somatic mutation (leave out any putative somatic recombinations). Assume we collected ***n*** samples from a single individual, a site of allele frequency = ***n-1*** indicates non-fixed and we would choose this as a start point. To take an example, we use ***n*** = 20 here.

* Two variant callers **UnifiedGenotyper** and **HaplotypeCaller** (http://www.broadinstitute.org/gatk/) were used here, one can add his own favorite caller

#### Step1-1: UnifiedGenotyper
      vcf_process.pl --vcf samples.ug.vcf.gz \
          --quality 30 --rare-only 19 | \
          bgzip -c > ug_single/samples.ug.fq19.vcf.gz && \
          tabix -p vcf ug_single/samples.ug.fq19.vcf.gz

      bgzip -dc ug_single/samples.ug.fq19.vcf.gz | \
          vcf-annotate --fill-type | \
          perl -ne 'next if(/\#/); next unless(/snp/); 
              my ($chrom, $pos)=(split /\s+/)[0,1]; 
              print "$chrom\t",($pos-1),"\t$pos\n";' \
          > ug_single/samples.ug.fq19.snp.bed

      bgzip -dc ug_single/samples.ug.fq19.vcf.gz | \
          vcf-annotate --fill-type | \
          perl -ne 'next if(/\#/); next unless(/ins/ || /del/ || /\*/); 
              my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
              my $start = $pos-1; my $end = $pos+length($ref); 
              print "$chrom\t$start\t$end\n";' | uniq \
          > ug_single/samples.ug.fq19.indel.bed



#### Step1-2: HaplotypeCaller
      vcf_process.pl --vcf samples.hc.vcf.gz \
          --quality 30 --rare-only 19 | \
          bgzip -c > hc_gvcf/samples.hc.fq19.vcf.gz && \
          tabix -p vcf hc_gvcf/samples.hc.fq19.vcf.gz

      bgzip -dc hc_gvcf/samples.hc.fq19.vcf.gz | \
          vcf-annotate --fill-type | \
          perl -ne 'next if(/\#/); next unless(/snp/); 
            my ($chrom, $pos)=(split /\s+/)[0,1]; 
            print "$chrom\t",($pos-1),"\t$pos\n";' \
          > hc_gvcf/samples.hc.fq19.snp.bed

      bgzip -dc hc_gvcf/samples.hc.fq19.vcf.gz | \
          vcf-annotate --fill-type | \
          perl -ne 'next if(/\#/); next unless(/ins/ || /del/ || /\*/); 
              my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
              my $start = $pos-1; my $end = $pos+length($ref); 
              print "$chrom\t$start\t$end\n";' | uniq \
          > hc_gvcf/samples.hc.fq19.indel.bed



#### Step1-3: Merge candidate SNV target regions
      cat ug_single/samples.ug.fq19.snp.bed \
          hc_gvcf/samples.hc.fq19.snp.bed | \
          sort -k1,1 -k2,3n | bedtools merge -i - \
          > combined/samples.fq19.snp.bed

<br />

### Step2: Generate more accurate allele depths for further filtering

* For SNVs, we obtain more accurate allele depths for each site using samtools mpileup (version 0.1.19, https://github.com/samtools/samtools) and Varscan (version 2.3.6, http://dkoboldt.github.io/varscan/)

* 3 sets of results will be generated in this step by combining slightly different criteria - MQ0 (mapping quality >= 0), MQ20 (mapping quality >= 20), AR (include anomalous reads) and NAR (exclude anomalous reads). Explanation of these criteria could be found in previous section.

* "xargs" is used here to do it in parallel

      find . -name "*.bam" -print | sed 's/.bam$//' | xargs -n 1 -P 12 -I PREFIX \
          sh -c '
              sample=`basename PREFIX | cut -d"." -f1`

              echo "[`date`]: Start processing ${sample} ... "

              ## count in anomalous reads, mapping quality >= 20
              samtools mpileup -Ad100000 -q20 -f reference_genome.fasta \
                  -l combined/samples.fq19.snp.bed \
                  PREFIX.bam | grep -vP "\t0\t" \
                  > readcounts/samples.fq19.snp.${sample}.MQ20.AR.mpileup

              java -jar VarScan/VarScan.v2.3.6.jar readcounts \
                  readcounts/samples.fq19.snp.${sample}.MQ20.AR.mpileup \
                  --min-base-qual 20 --min-coverage 1 \
                  --output-file readcounts/samples.fq19.snp.${sample}.MQ20.AR.readcounts


              ## count in anomalous reads, mapping quality >= 0
              samtools mpileup -Ad100000 -q0 -f reference_genome.fasta \
                  -l combined/samples.fq19.snp.bed \
                  PREFIX.bam | grep -vP "\t0\t" \
                  > readcounts/samples.fq19.snp.${sample}.MQ0.AR.mpileup

              java -jar VarScan/VarScan.v2.3.6.jar readcounts \
                  readcounts/samples.fq19.snp.${sample}.MQ0.AR.mpileup \
                  --min-base-qual 20 --min-coverage 1 \
                  --output-file readcounts/samples.fq19.snp.${sample}.MQ0.AR.readcounts

              ## only use proper pairs, mapping quality >= 20
              samtools mpileup -d100000 -q20 -f reference_genome.fasta \
                  -l combined/samples.fq19.snp.bed \
                  PREFIX.bam | grep -vP "\t0\t" \
                  > readcounts/samples.fq19.snp.${sample}.MQ20.NAR.mpileup

              java -jar VarScan/VarScan.v2.3.6.jar readcounts \
                  readcounts/samples.fq19.snp.${sample}.MQ20.NAR.mpileup \
                  --min-base-qual 20 --min-coverage 1 \
                  --output-file readcounts/samples.fq19.snp.${sample}.MQ20.NAR.readcounts

              echo "[`date`]: Finished processing ${sample}"
          '

<br />
    
* The **fillVcfDepth.pl** needs to read the whole file into memory and could be rather memory-consuming. An alternative script **fillVcfDepth-split.pl** only read each chromosome per time and is thus more preferable for large files. **fillVcfDepth-split.pl** requires a bgzip-compressed file and a tabix index (usually bundled with samtools). However, **fillVcfDepth-split.pl** could be extremely slow if the genome contains many sequences, e.g., a draft genome with tens of thousands of scaffolds.
    
        find readcounts/ -name "samples.fq19.snp.*.MQ20.NAR.readcounts" | \
            xargs -n 1 -P 6 -I {} bgzip {}
            
        find readcounts/ -name "samples.fq19.snp.*.MQ20.NAR.readcounts.gz" | \
            xargs -n 1 -P 6 -I {} tabix -S1 -s1 -b2 -e2 {}

        ## count in anomalous reads, MQ>=20
        for f in `find readcounts/ \
            -name "samples.fq19.snp.*.MQ20.AR.readcounts.gz" | sort`;
        do
            library=`basename $f | cut -d"." -f9`
            sample=${library}

            echo "${sample} ${library} ${f}"
        done > combined/samples.fq19.snp.MQ20.AR.readcounts.list

        ## count in anomalous reads, MQ>=0
        for f in `find readcounts/ \
            -name "samples.fq19.snp.*.MQ0.AR.readcounts.gz" | sort`;
        do
            library=`basename $f | cut -d"." -f9`
            sample=${library}

            echo "${sample} ${library} ${f}"
        done > combined/samples.fq19.snp.MQ0.AR.readcounts.list

        ## count proper mapped reads only, MQ>=20
        for f in `find readcounts/ \
            -name "samples.fq19.snp.*.MQ20.NAR.readcounts.gz" | sort`;
        do
            library=`basename $f | cut -d"." -f9`
            sample=${library}

            echo "${sample} ${library} ${f}"
        done > combined/samples.fq19.snp.MQ20.NAR.readcounts.list


        find combined/ -name "samples.fq19.snp.*.readcounts.list" | \
            xargs -n 1 -P 4 -I INPUT \
            sh -c '
                out_set=`echo INPUT | sed "s/.*.fq19.snp.//" | sed "s/.readcounts.list//"`

                echo "${out_set}"

                ## UnifiedGenotyper
                fillVcfDepth-split.pl --vcf ug_single/samples.ug.fq19.vcf.gz \
                    --list INPUT --minimum-vcf --update-AD | bgzip -c \
                    > ug_single/samples.ug.fq19.snp.${out_set}.vcf.gz && \
                    tabix -p vcf ug_single/samples.ug.fq19.snp.${out_set}.vcf.gz

                ## HaplotypeCaller
                fillVcfDepth-split.pl --vcf hc_gvcf/samples.hc.fq19.vcf.gz \
                    --list INPUT --minimum-vcf --update-AD | bgzip -c \
                    > hc_gvcf/samples.hc.fq19.snp.${out_set}.vcf.gz && \
                    tabix -p vcf hc_gvcf/samples.hc.fq19.snp.${out_set}.vcf.gz
            '

<br />


* For small insertion/deleletions (indels), we run HaplotypeCaller in **multiple-sample calling** mode around previous candidate indel sites

      ## UnifiedGenotyper
      Samples=`find . -name "*.bam" -print | \
          xargs -I BAM_FILE echo -n "-I BAM_FILE "`
      java -jar GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
          -R reference_genome.fasta \
          -T HaplotypeCaller -nct 1 -stand_call_conf 30.0 -stand_emit_conf 30.0 -ip 20 \
          -L ug_single/samples.ug.fq19.indel.bed \
          -o ug_single/samples.ug.fq19.indel.hc_multi.vcf \
          ${Samples} 2>& 1 | \
          tee ug_single/samples.ug.fq19.indel.hc_multi.log

      bgzip ug_single/samples.ug.fq19.indel.hc_multi.vcf && \
          tabix -p vcf ug_single/samples.ug.fq19.indel.hc_multi.vcf.gz


      ## HaplotypeCaller
      Samples=`find . -name "*.bam" -print | \
          xargs -I BAM_FILE echo -n "-I BAM_FILE "`
      java -jar GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
          -R reference_genome.fasta \
          -T HaplotypeCaller -nct 1 -stand_call_conf 30.0 -stand_emit_conf 30.0 -ip 20 \
          -L hc_gvcf/samples.hc.fq19.indel.bed \
          -o hc_gvcf/samples.hc.fq19.indel.hc_multi.vcf \
          ${Samples} 2>& 1 | \
          tee hc_gvcf/samples.hc.fq19.indel.hc_multi.log

      bgzip hc_gvcf/samples.hc.fq19.indel.hc_multi.vcf && \
          tabix -p vcf hc_gvcf/samples.hc.fq19.indel.hc_multi.vcf.gz
    

<br />

### Step3: Screen out candidate mutations

* Parallel comparison is extreme-efficient in removing false positive calls but rely heavily on whether a proper "control" is available. For germline analysis, we usually only focused on sample-specific mutations, and all other non-focal samples could be used as controls. For somatic analysis, such a control is not so clear as many mutations present in multiple samples.

* For somatic analysis, one expectation is that samples with closer physical or genealogical relations are more likely to share true mutations (mutations raised before they seperating), thus we could first set up groups for each sample if we know these details. Take a tree as an example, leaves from the same branch are more likely to share true somatic mutations, while "mutations" found in two leaves from different branches are more likely to be false positives. Therefore, we could group samples by branches, and then compare each branch. This method is termed as "topology-based" here. This approach remain the notion of a parallel compare and is robust to various sequencing/mapping/calling artifacts.

* However, it's not clear whether the pre-assumption of topology-based approach is always true. A more direct approach is to just look at the allele frequency of each variant as mentioned at the beginning, and filtering out those non-fixed variant sites, termed as the "frequency-based" approach. In practice, this approach is subjected to various artifacts (as the focal samples and control samples are not well distinguishable) and is of high false positive rate as well as high false negative rate. Even so, it serve as a good complementary to test whether the topology-based approach really lose many true mutations. 

* When using the frequency-based approach, to deal with high false positive discover rate, a good practice is to test several criteria (e.g., --max-shared-freq, --max-cmp-miss, and etc.), from most strict to most loose, and evaluate the results by inspecting those differences between different criteria.

<br />

#### Step3-1: filtering substitution mutations

    find . -name "samples.*.fq19.snp.*.vcf.gz" | \
        xargs -n 1 -P 6 -I VCFIN \
        sh -c '
            out_prefix=`echo VCFIN | sed "s/vcf.gz/mut/"`

            echo "${out_prefix}"

            ## screen according to topology
            detect_mutations.pl -v VCFIN --max-cmp-depth 2 --max-cmp-total 3 \
                --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
                -g groups.txt | \
                vcf-annotate -f c=3,150 --fill-type \
                --output ${out_prefix}.c2t3d5m5.topology.vcf

            ## screen according to frequency
            detect_mutations.pl -v VCFIN \
                --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 \
                --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 15 | \
                vcf-annotate -f c=3,150 --fill-type \
                --output ${out_prefix}.s15c2t3d5m0.frequency.vcf
        '

<br />

#### Step3-2: filtering indel mutations

* For indels, the filtering is based on AD field. Very few sites from the VCF file generated by HaplotypeCaller do not have AD field (which usually with low confidence) need to be filtered use the script "flt_vcf_AD"

      find  -name "samples.*.fq19.indel.*.vcf.gz" | \
          xargs -n 1 -P 0 -I VCFIN \
          sh -c '
              out_prefix=`echo VCFIN | sed "s/hc_multi.vcf.gz/mut/"`

              ## screen according to topology
              vcf_process.pl --vcf VCFIN --quality 30 --var-type indel | flt_vcf_AD | \
                  detect_mutations.pl -v - --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 \
                  --max-cmp-miss 5 -g groups.txt | \
                  vcf-annotate -f c=3,150 --fill-type \
                  > ${out_prefix}.c2t3d5m5.topology.vcf

              ## screen according to frequency
              vcf_process.pl --vcf VCFIN --quality 30 --var-type indel | flt_vcf_AD | \
                  detect_mutations.pl -v - --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 \
                  --max-cmp-miss 0 --max-shared-freq 15 | \
                  vcf-annotate -f c=3,150 --fill-type \
                  > ${out_prefix}.s15c2t3d5m0.frequency.vcf
          '

<br />

### Step4: Collecting candidate mutations from different set

* A "softmask" strategy is highly recommended to collect all candidate mutations from different criteria/approaches/methods, and rank each site by tagging it with those details. Through inspection of those details, one can easily distinguish false positive calls while minimize false negative rate.

<br />

#### Step4-1: Combine results from topology-based approach

#### 1) Combine results from UnifiedGenotyper

    vcf_process.pl --vcf ug_single/samples.ug.fq19.snp.MQ20.NAR.mut.c2t3d5m5.topology.vcf \
        --secondary-vcf ug_single/samples.ug.fq19.snp.MQ20.AR.mut.c2t3d5m5.topology.vcf \
        --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
        vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
        --secondary-vcf ug_single/samples.ug.fq19.snp.MQ0.AR.mut.c2t3d5m5.topology.vcf \
        > ug_single/samples.ug.fq19.snp.mut.c2t3d5m5.topology.combined.vcf

#### 2) Combine results from HaplotypeCaller

    vcf_process.pl --vcf hc_gvcf/samples.hc.fq19.snp.MQ20.NAR.mut.c2t3d5m5.topology.vcf \
        --secondary-vcf hc_gvcf/samples.hc.fq19.snp.MQ20.AR.mut.c2t3d5m5.topology.vcf \
        --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
        vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
        --secondary-vcf hc_gvcf/samples.hc.fq19.snp.MQ0.AR.mut.c2t3d5m5.topology.vcf \
        > hc_gvcf/samples.hc.fq19.snp.mut.c2t3d5m5.topology.combined.vcf

<br />

#### Step4-2: Combine results from frequency-based approach

#### 1) Combine results from UnifiedGenotyper

    vcf_process.pl --vcf ug_single/samples.ug.fq19.snp.MQ20.NAR.mut.s15c2t3d5m0.frequency.vcf \
        --secondary-vcf ug_single/samples.ug.fq19.snp.MQ20.AR.mut.s15c2t3d5m0.frequency.vcf \
        --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
        vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
        --secondary-vcf ug_single/samples.ug.fq19.snp.MQ0.AR.mut.s15c2t3d5m0.frequency.vcf \
        > ug_single/samples.ug.fq19.snp.mut.s15c2t3d5m0.frequency.combined.vcf

#### 2) Combine results from HaplotypeCaller

    vcf_process.pl --vcf hc_gvcf/samples.hc.fq19.snp.MQ20.NAR.mut.s15c2t3d5m0.frequency.vcf \
        --secondary-vcf hc_gvcf/samples.hc.fq19.snp.MQ20.AR.mut.s15c2t3d5m0.frequency.vcf \
        --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
        vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
        --secondary-vcf hc_gvcf/samples.hc.fq19.snp.MQ0.AR.mut.s15c2t3d5m0.frequency.vcf \
        > hc_gvcf/samples.hc.fq19.snp.mut.s15c2t3d5m0.frequency.combined.vcf

<br />

#### Step4-3: Combine different methods

    vcf_process.pl --vcf hc_gvcf/samples.hc.fq19.snp.mut.c2t3d5m5.topology.combined.vcf \
        --secondary-vcf ug_single/samples.ug.fq19.snp.mut.c2t3d5m5.topology.combined.vcf \
        --combine-rows 0 1 --compare-rows 2 3 4 \
        --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
        > combined/samples.fq19.snp.mut.c2t3d5m5.topology.combined.vcf

    vcf_process.pl --vcf hc_gvcf/samples.hc.fq19.snp.mut.s15c2t3d5m0.frequency.combined.vcf \
        --secondary-vcf ug_single/samples.ug.fq19.snp.mut.s15c2t3d5m0.frequency.combined.vcf \
        --combine-rows 0 1 --compare-rows 2 3 4 \
        --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
        > combined/samples.fq19.snp.mut.s15c2t3d5m0.frequency.combined.vcf

<br />

#### Step4-4: Combine topology and frequency

    vcf_process.pl --vcf combined/samples.fq19.snp.mut.c2t3d5m5.topology.combined.vcf \
        --secondary-vcf combined/samples.fq19.snp.mut.s15c2t3d5m0.frequency.combined.vcf \
        --combine-rows 0 1 --compare-rows 2 3 4 \
        --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
        > combined/samples.fq19.snp.mut.combined.vcf


<br />


#### Similar for indels, simply combine different all results

    vcf_process.pl --vcf hc_gvcf/samples.hc.fq19.indel.mut.c2t3d5m5.topology.vcf \
        --secondary-vcf ug_single/samples.ug.fq19.indel.mut.c2t3d5m5.topology.vcf \
        --combine-rows 0 1 --compare-rows 2 3 4 \
        --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
        > combined/samples.fq19.indel.mut.c2t3d5m5.topology.combined.vcf

    vcf_process.pl --vcf hc_gvcf/samples.hc.fq19.indel.mut.s15c2t3d5m0.frequency.vcf \
        --secondary-vcf ug_single/samples.ug.fq19.indel.mut.s15c2t3d5m0.frequency.vcf \
        --combine-rows 0 1 --compare-rows 2 3 4 \
        --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
        > combined/samples.fq19.indel.mut.s15c2t3d5m0.frequency.combined.vcf

    vcf_process.pl --vcf combined/samples.fq19.indel.mut.c2t3d5m5.topology.combined.vcf \
        --secondary-vcf combined/samples.fq19.indel.mut.s15c2t3d5m0.frequency.combined.vcf \
        --combine-rows 0 1 --compare-rows 2 3 4 \
        --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
        > combined/samples.fq19.indel.mut.combined.vcf


<br />

### Step5: Evidence-collecting and filtering

* The original lengthy "inline" PERL codes is now replaced by a new script "annotate_mut_vcf.pl";
  
* Step1~4 generally collecting all possible candidates and only remove (1) obvious common varaints; (2) apparent artefacts repeatedly show up in less related samples. Step5 does the remaining filtering works to obtain a confidence set of mutations. It's highly recommend to further inspec (manually or experimentally) the results from this step, including the ones of higher confidence (assess false positive discover rate) as well as with lower confidence (assess false negative discover rate), especially for a pilot survey. 

* A range of evidences (annotations) which likely suggest false positives are collected here based on prior experiences. The current avaliable metrics includes:

* **[FATAL level: most likely false positive calls]**
        1. ref_allele: a reference allele is most likely a common allele;
        2. proper_reads_only: mutation is only seen when calling without anomalous mapped reads (e.g., NAR);
        3. high_mapq_only: mutation is only seen when calling with high mapping quality (e.g., MQ20);
        4. UG_misalign: variant only called by UnifideGenotyper in slippage/tandem region with nearby INDEL (UG-only + near_indel + trf_predict);
        5. clipped_only: mutation allele is only supported by clipped reads (no read-support after excluding the clipped ones);
        6. pre_exist: mutation allele is found in pre-existing variants (require a list of pre-existing alleles that should not be taken as mutations);
        7. unexpected_cooccur: mutation co-occurred preferentially in non-closely related samples rather than closely related samples (based on the experimental design on the kinship of samples)
	        
* **[WARN level: need cautious]**
        1. anomalous_reads_only: mutation is only seen when calling with anomalous mapped reads (e.g. AR);
        2. low_mapq_only: mutation is only seen when calling with low mapping quality (e.g., MQ0);
        3. missing_calls: presence of one or more missing calls;
        4. caller_bias: variant only called by part of callers;
        5. low_qual: low variant quality score;
        6. biased_scall: many samples in focal group have the "mutated reads" but only very few have sufficient read-depth;
        7. indel_nearby: SNVs called with nearby INDELs are prone to mis-alignemnt (need to provide the information of whether there is INDEL nearby);
        8. unexpected_presence: one or more "mutated-reads" seen in non-mutated group (more strict requirement of maximum FPD);
        9. clustered: mutations clustered within a small range is suspect to mapping/calling artefacts especially with "high_mapq_only" (require "Cluster" tag in FILTER field, can be generated by "vcf-annotate -f" as shown above)
	    

		annotate_mut_vcf.pl --input combined/samples.fq19.snp.mut.combined.vcf \
		    --surround-indel combined/samples.fq19.snp.mut.combined.indels_nearby.csv \
		    --trf combined/samples.fq19.snp.mut.combined.trf.csv \
		    --kinship samples.kinship.txt \
		    --min-call-frac 0.25 --apply-set SET_1:unexpected_cooccur \
			--max-sample-dist 3 --incl-grps --remove "1+3" \
		    --out-appends sets_combined trf sample_distance \
		    --out-stats combined/samples.fq19.snp.mut.combined.ann.stats.csv \
		    --output combined/samples.fq19.snp.mut.combined.ann.csv


<br />


### Accompany scripts and files used in this step are given below:

#### annotate_mut_vcf.pl

> Annotate mutation sites with various metrics to predict their confidences.

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
	
	
* **Options for INFO-based annotation:**
	    
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
	    
	      
* **Options for BAM-based annotation:**
	        
	        This script checks whether a mutation allele was supported by reliable
	        reads by excluding reads failed various thresholds here, multiple
	        thresholds will be applied together if specified
	        
	        *Note: current version does not handle reads with indels well, so if
	        reads supporting the mutation alleles also contain indels, these reads
	        could be missed
	        
	    -b, --bam-list <filename>
	        check supporting reads of a mutation allele by looking back into bam
	        files, a list of bam files correspond to each sample id should be
	        given in the format:
	        "Sample_id BAM_file_location"
	        
	    -s, --samtools <string>
	        options passed to samtools view, e.g. "-f 4 -F 8", use this option to
	        remove reads that do not want to confirm the filters below, e.g.
	        unmapped/supplementary/secondary reads
	
	    --max-clip-ratio <float>
	        maximum ratio of soft/hard clipped length, calculated as
	        (soft-clipped-bases + hard-clipped-bases) / read_length
	        [default: 1 (i.e., no filtering)]
	    
	    -n, --no-both-clips  <string>
	        only filter reads with both ends clipped when exceeding
	        the --max-clip-ratio [default: no filtering]
	        
	    --min-supp-depth <int>
	        minimum requirement of realible supporting read-depth for each sample
	        [default: 5]
	
	        
* **Options for other annotation:**
	
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
	         
	    
* **Options for applying filters based on the confidence annotations:**
	
	    See above explanations on evidences/annotations used here.
	        
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


#### Annotate nearby INDELs

* vcf-annotate in vcftools (http://vcftools.sourceforge.net/) and bcftools (http://samtools.github.io/bcftools/bcftools.html) are used here, "samples.fq19.snp.mut.combined.pos.csv" contains the positions (two rows: "CHROM" and "POS") of the combined/samples.fq19.snp.mut.combined.vcf

		bgzip -dc samples.hc.vcf.gz | \
		    vcf-annotate --fill-type | \
		    perl -ne 'next if(/\#/); next unless(/ins/ || /del/); 
				my ($chrom, $pos, $ref, $alt)=(split /\s+/)[0,1,3,4];
		        my $len = abs(length($ref) - length($alt));
		        my $start = $pos-1; my $end = $pos+length($ref); 
				print "$chrom\t$start\t$end\n";' | uniq \
		    > samples.hc.indel.bed
		
		perl -e 'print "##INFO=<ID=GapDistance,Number=1,Type=String,Description=\"Gaps Nearby\">\n";' \
		    > indel.ex100.hdr
		
		perl -ne 'my ($chrom, $start, $end) = (split /\t/)[0..2]; 
			$start -= 100; $end += 100; $start = 1 if $start < 1;
		    print "$chrom\t$start\t$end\tex100bp\n";' \
		   samples.hc.indel.bed \
		    > samples.hc.indel.ex100.tab
		
		bgzip samples.hc.indel.ex100.tab && \
		    tabix -s1 -b2 -e3 samples.indel.ex100.tab.gz
		
		
		awk 'BEGIN{OFS="\t"} {if(!/\#/){$8 = ".";} print;}' \
		    combined/samples.fq19.snp.mut.combined.vcf | \
		    bcftools annotate -a samples.hc.indel.ex100.tab.gz \
		    -h indel.ex100.hdr \
		    -c CHROM,FROM,TO,GapDistance - | \
		    awk 'BEGIN{OFS="\t"} !/\#/ {if(/GapDistance=ex100bp/){print $1,$2,"YES";}else{print $1,$2,"NO";}}' | \
		    map_records.pl --query combined/samples.fq19.snp.mut.combined.pos.csv \
		    --subject - --rows1 0 1 --rows2 0 1 --merge-dups \
		    > combined/samples.fq19.snp.mut.combined.indels_nearby.csv


#### Annotate slippage/tandem regions

* require trf (https://tandem.bu.edu/trf/trf.html) and bedtools (https://bedtools.readthedocs.io/en/latest/), "convert_trf.pl" and "map_pos2intervals.pl" are in the BioScripts repository


		trf reference.fasta 2 7 7 80 10 20 150 -m -f -d -h
		
		convert_trf.pl -i reference.fasta.2.7.7.80.10.20.150.dat \
		    -f bed | sort -k1,1 -k2,3n | bedtools merge -i - -nms \
		    > reference.trf.w150.bed
			
		map_pos2intervals.pl --query combined/samples.fq19.snp.mut.combined.pos.csv \
		    --subject reference.trf.w150.bed \
		    > combined/samples.fq19.snp.mut.combined.trf.csv


#### 


### Step6: Manual inspection and criteria adjusting

* Please refer to the parental folder for codes to generate the figures and alignments in helping manual inspection.

<br />
