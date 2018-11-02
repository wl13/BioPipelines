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

### Step5: Formating results

* As a rule of thumb, a site pass following criteria is of highest confidence and will be marked as "Confidence" here.

  1. Pass "AR", "NAR" filters - meaning no anamlous reads;
  2. Pass "MQ0", "MQ20" filters - meaning high mapping confidence;
  3. Pass "Missing" filter - the higher the missing rate the lower informative of a site;
  4. Could be called by HaplotypeCaller (HC_GVCF);
  5. "FPD<=1" - One or more reads carrying identical "mutation allele" in control samples usually suggest a false positive call;
  6. Variant quality >= 50, this criterion also depends on the sequencing depth, and could be lowered, but generally has less influence
  
* All other sites could be viewed as "Evaluation set". Further inspecition (manually or experimentally) of the confidence set (assess false positive discover rate) and evalution set (assess false negative discover rate) is highly recommended, and is usually necessary for a pilot survey. Some criteria could be fine-tuned through this trial-and-error process.


      perl -ne 'next if (/^\#\#/); if (/\#/) {
          print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION\tLength\tMethods\tMut_Type\tFrequency\tMut_Set\tCombine\n";
          next;} my @line = (split /\s+/); my $info = $line[7]; $info =~ /MA=(\w+)/; my $mut_allele = $1;
          my $var_type = 0; if ($info =~ /MATYPE\=INDEL/) {$var_type = "INDEL";}
          my $fq_sum = 1; if($info =~ /Shared\=(\d+)/){$fq_sum = $1;}
          my $type = ($mut_allele eq $line[3]) ? "REF" : "ALT"; my $out_line = join "\t", @line;
          my $method = "UG_Single"; 
          if ($info =~ /UG_Single\+HC_GVCF/){$method = "UG_Single+HC_GVCF";}elsif($info =~ /HC_GVCF/){$method = "HC_GVCF";}
          my @filters = (); if ($info =~ /Combine=AR;/) {push @filters, "AR";} if ($info =~ /Combine=NAR;/) {push @filters, "NAR";}
          if ($info =~ /Combine=MQ20;/) {push @filters, "MQ20";} if ($info =~ /Combine=MQ0;/) {push @filters, "MQ0";}
          if ($info !~ /NMISS=0/) {push @filters, "MISSING";} 
          if (($info !~ /FPD=0/) && ($info !~ /FPD=1;FPFQ=1;/)) {push @filters, "FPD";}
          if ($method eq "UG_Single") {push @filters, "UG";} if ($line[5] < 50) {push @filters, "LowQual";}
          my $set = "TP"; 
          if ($info =~ /Combine=Grouped\+NonGrouped/){$set = "TP+FQ";}elsif($info =~ /Combine=NonGrouped/){$set = "FQ";}
          if (scalar @filters == 0) {$set .= "(Confidence)";} else {my $filters = join ",", @filters; $set .= "($filters)";}
          print "$out_line\t$var_type\t$method\t$type\t$fq_sum\t$set\t$combine\n";' \
          combined/samples.mut.combined.vcf \
          > combined/samples.mut.combined.csv
    
    
<br />

* Tips: For manually inspection, the previous section contains codes for generating figures and alignments which ease the whole process.

<br />
