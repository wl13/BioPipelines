###Pipelines to detect "mutations" (rare variants) among a cohort of samples


Assume the sequencing data were already mapped and pre-processed (e.g. mark PCR duplicates, etc.)


###Pipelines

####Step1: Generate initial candidate targets

* Use vcf_process.pl to screen out **rare variants**, set a frequency higher than the expected frequency of the mutations to tolerant genotyping errors

		vcf_process.pl --vcf samples.vcf.gz --quality 30 --rare-only 3 | \
		    bgzip -c > samples.fq3.vcf.gz && tabix -p vcf samples.fq3.vcf.gz
    
    
* Extract candidate positions for **base substitutions(single nucleotide variants)**

		bgzip -dc samples.fq3.vcf.gz | vcf-annotate --fill-type | \
		    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1];
		        print "$chrom\t",($pos-1),"\t$pos\n";' > samples.fq3.snp.bed


* Extract candidate regions for **indels**

		bgzip -dc samples.fq3.vcf.gz | vcf-annotate --fill-type | \
		    perl -ne 'next if(/\#/); next unless(/ins/ || /del/); my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
		        my $start = $pos-1; my $end = $pos+length($ref); print "$chrom\t$start\t$end\n";' | uniq | \
		    bedtools slop -i - -b 10 -g reference.fasta.fai > samples.fq3.indel.bed

	**Note for extracting indel target regions:**  
 1. the region need to **cover almost the whole deletion** for HaplotypeCaller to re-generate the same call;  
 2. extend calling regions will possibly change the emitted quality values or lost several calls;  
 3. the extend positions **should not exceed the chromosme length**, otherwise will cause the HaplotypeCaller to fail. 
	


####Step2: Count accurate allele depths for each locus and each sample


* Use VarScan (version 2.3.6, http://dkoboldt.github.io/varscan/) and samtools (version 0.1.19, https://github.com/samtools/samtools) to generate read counts for each SNP

		find bam_files/ -name "*.bam" -print | sed 's/.bam$//' | xargs -n 1 -P 6 -I PREFIX \
		    sh -c '
		        sample=`basename PREFIX | cut -d"." -f1`
		        
		        echo "[`date`]: Start processing ${sample} ... "
		        
		        samtools mpileup -Ad100000 -q20 -f reference.fasta -l samples.fq3.snp.bed \
		            PREFIX.bam | grep -vP "\t0\t" \
		            > readcounts/samples.fq3.snp.${sample}.mpileup
		        
		        java -jar VarScan.v2.3.6.jar readcounts readcounts/samples.fq3.snp.${sample}.mpileup \
		            --min-base-qual 20 --min-coverage 1 \
		            --output-file readcounts/samples.fq3.snp.${sample}.readcounts
		        
		        echo "[`date`]: Finished processing ${sample}"
		    '

	**Note for reads counted:**  
	1. Requires samtools version >= 0.1.9, as older versions did not deploy Base Alignment Quality (BAQ), which aims to provide an efficient and effective way to rule out false SNPs caused by nearby INDELs.  
	2. The default settings for mpileup is to filter reads with bitwise flag 0X704. So for mpileup generation the following reads will not been considered:  
	  1. 0x0400 (aka 1024 or "d"), duplicate;   
	  2. 0x0200 (aka 512 or "f"), failed QC;  
	  3. 0x0100 (aka 256 or "s"), non primary alignment;  
	  4. 0x0004 (aka 4 or "u"), unmapped.   
	3. Apply -A to use anomalous read pairs in mpileup, which are not used by default (requring r874+).


* Get the file list of all counting results

		for f in `readcounts/ -name "samples.fq3.snp.*.readcounts" | sort`;
		do
		    library=`basename $f | cut -d"." -f4`
		    sample=${library}
		    
		    echo "${sample} ${library} ${f}"
		done > samples.fq3.snp.readcounts.list


* Fill depth fields (RC) use fillVcfDepth.pl

		fillVcfDepth.pl --vcf samples.fq3.vcf.gz --list samples.fq3.snp.readcounts.list \
		    --minimum-vcf --update-AD | bgzip -c \
		    > samples.fq3.snp.vcf.gz && tabix -p vcf samples.fq3.snp.vcf.gz


* Run HaplotypeCaller in multiple-sample calling mode around indels to generate more accurate AD values for each sample

		Samples=`find bam_files/ -name "*.bam" -print | xargs -I BAM_FILE echo -n "-I BAM_FILE "`
		java -jar GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
		    -R reference.fasta \
		    -T HaplotypeCaller -nct 1 -stand_call_conf 30.0 -stand_emit_conf 30.0 \
		    -L samples.fq3.indel.bed -o samples.fq3.indel.hc_multi.vcf \
		    ${Samples} 2>& 1 | tee samples.fq3.indel.hc_multi.log
		
		bgzip samples.fq3.indel.hc_multi.vcf && tabix -p vcf samples.fq3.indel.hc_multi.vcf.gz



#####Step3: Screen out candidate mutations


## filtering point mutations
detect_mutations.pl -v samples.fq3.snp.vcf.gz --max-cmp-depth 2 --max-cmp-total 3 \
    --max-cmp-miss 5 --min-supp-depth 5 --min-supp-plus 1 --min-supp-minus 1 \
    --mask-only LowDepth -g samples.group.txt | \
    vcf-annotate -f c=3,150 --fill-type > samples.fq3.snp.mut.c2t3d5m5.vcf


## filtering indel mutations
vcf_process.pl --vcf samples.fq3.indel.hc_multi.vcf.gz --quality 30 --var-type indel | \
    awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 \
    --max-cmp-miss 5 --mask-only LowDepth -g samples.group.txt | \
    vcf-annotate -f c=2,150 --fill-type > samples.fq3.indel.mut.c2t3d5m5.vcf


##
## Step4: rerun all Step1~3 with one or more different callers, preferentially those
## implemented with a different algorithm
##



##
## Step5: combine results with different callers
##

vcf_process.pl --vcf caller1.mut.vcf --secondary-vcf caller2.mut.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag CALLER1 --secondary-tag CALLER2 \
    --intersect-tag CALLED_BOTH > mut.combined.vcf



##
## Step6: generate alignments and figures for manually inspections
##

## local re-align reads to reference with another aligner like ClustalW2 (http://www.clustal.org/clustal2/)
for record in `cat mut.vcf | perl -ne 'next if (/^\#/); my ($chrom, $pos, $sample) = (split /\s+/)[0,1,2,7];
    if($sample =~ /\;/) {$sample = "Shared"} print "$sample;$chrom:$pos#$1\n";'`;
do
    sample=${record/;*}
    mutation=${record/*;}
    mut_base=${mutation/*\#}
    mutation=${mutation/\#*}
    chrom=${mutation/:*}
    mut_pos=${mutation/*:}

    start_pos=`echo "${mut_pos}-200" | bc`
    end_pos=`echo "${mut_pos}+200" | bc`

    echo "${sample} ${chrom} ${mut_pos} ${mut_base}"
    
    if [[ -n alignments/${sample} ]]; then
        mkdir -pv alignments/${sample}
    fi
    
    echo -e "${chrom} ${start_pos} ${end_pos}\n${chrom} ${start_pos} ${mut_pos}\n${chrom} ${mut_pos} ${end_pos}\n" | \
        fasta_process.pl --rows 0 1 2 --subset 1 2 --query - \
        --fasta reference.fasta > alignments/${sample}/${chrom}_${mut_pos}.fa
    
    for bam_file in `ls bam_files/*.bam`;
    do
        sample=`basename ${bam_file} | cut -d"." -f1`
        samtools view -X -F 3844 ${bam_file} ${chrom}:${mut_pos}-${mut_pos} | \
            awk -v name=${sample} 'BEGIN {OFS = FS = "\t"}; {print ">"name"|"$2"|"$1"\n"$10;}' \
            >> alignments/${sample}/${chrom}_${mut_pos}.fa
    done

    reference_align.pl -i alignments/${sample}/${chrom}_${mut_pos}.fa \
        > alignments/${sample}/${chrom}_${mut_pos}.aln.fas && \
        rm -v alignments/${sample}/${chrom}_${mut_pos}.fa
done




## generate figures of bam alignments in IGV (https://www.broadinstitute.org/igv/)
echo "snapshotDirectory alignments" > samples.run_igv.txt

for bam_file in `ls bam_files/*.bam`;
do
    echo "load ${bam_file}" >> samples.run_igv.txt
done

cat mut.vcf | perl -ne 'next if (/^\#/); my ($chrom, $pos, $sample) = (split /\s+/)[0,1,2];
    if($sample =~ /\;/) {$sample = "Shared"}
    my $ex_start=$pos-55; my $ex_end=$pos+55;
    print "goto $chrom:$ex_start-$ex_end\nsnapshot $sample\/$chrom\_", "$pos", ".ex55.png\n";
    $ex_start=$pos-250; $ex_end=$pos+250;
    print "goto $chrom:$ex_start-$ex_end\nsnapshot $sample\/$chrom\_", "$pos", ".ex250.png\n";
    $ex_start=$pos-3000; $ex_end=$pos+3000;
    print "goto $chrom:$ex_start-$ex_end\nsnapshot $sample\/$chrom\_", "$pos", ".ex3000.png\n";' \
    >> samples.run_igv.txt







