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
    
    
### Step3: Screen out candidate mutations

* Parallel comparison is extreme-efficient in removing false positive calls but rely heavily on whether a proper "control" is available. For germline analysis, we usually only focused on sample-specific mutations, and all other non-focal samples could be used as controls. For somatic analysis, such a control is not so clear as many mutations present in multiple samples.

* For somatic analysis, one expectation is that samples with closer physical or genealogical relations are more likely to share true mutations (mutations raised before they seperating), thus we could first set up groups for each sample if we know these details. Take a tree as an example, leaves from the same branch are more likely to share true somatic mutations, while "mutations" found in two leaves from different branches are more likely to be false positives. Therefore, we could group samples by branches, and then compare each branch. This method is termed as "topology-based" here. This approach remain the notion of a parallel compare and is robust to various sequencing/mapping/calling artifacts.

* However, it's not clear whether the pre-assumption of topology-based approach is always true. A simple solution is to just look into the allele frequency, and filtering out those truely non-fixed mutations, termed as a "frequency-based" approach. Ideally, this approach should work 

#### Step3-1: substitution mutations
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
                --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 50 | \
                vcf-annotate -f c=3,150 --fill-type \
                --output ${out_prefix}.s50c2t3d5m0.frequency.vcf
        '

#### Step3-2: indel mutations
find  -name "samples.*.fq19.indel.*.vcf.gz" | \
    xargs -n 1 -P 0 -I VCFIN \
    sh -c '
        out_prefix=`echo VCFIN | sed "s/hc_multi.vcf.gz/mut/"`
        
        ## screen according to topology
        vcf_process.pl --vcf VCFIN --quality 30 --var-type indel | flt_vcf_AD | \
            detect_mutations.pl -v - --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 \
            -g samples.3_groups.txt | \
            vcf-annotate -f c=3,150 --fill-type \
            > ${out_prefix}.c2t3d5m5.topology.vcf
        
        ## screen according to frequency
        vcf_process.pl --vcf VCFIN --quality 30 --var-type indel | flt_vcf_AD | \
            detect_mutations.pl -v - --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 0 --max-shared-freq 50 | \
            vcf-annotate -f c=3,150 --fill-type \
            > ${out_prefix}.s50c2t3d5m0.frequency.vcf
    '


