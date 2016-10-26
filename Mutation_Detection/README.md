###Pipelines to detect "mutations" (rare variants) among a cohort of samples use mapping results


Assume the sequencing data were already mapped and pre-processed (e.g. mark PCR duplicates, etc.)




##Pipelines   

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
	
<br />
   
####Step2: Count accurate allele depths for each locus and each sample


* Use VarScan (version 2.3.6, http://dkoboldt.github.io/varscan/) and samtools (version 0.1.19, https://github.com/samtools/samtools) to generate read counts for each **SNV site**

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
	1. Requires **samtools version >= 0.1.9**, as older versions did not deploy Base Alignment Quality (BAQ), which aims to provide an efficient and effective way to rule out false SNPs caused by nearby INDELs.  
	2. The default settings for mpileup is to filter reads with bitwise flag 0X704. So for mpileup generation the following reads will not been considered:  
	  1. 0x0400 (aka 1024 or "d"), duplicate;   
	  2. 0x0200 (aka 512 or "f"), failed QC;  
	  3. 0x0100 (aka 256 or "s"), non primary alignment;  
	  4. 0x0004 (aka 4 or "u"), unmapped.   
	3. Apply "-A" option would include anomalous read pairs, those are reads with mate unmapped or not proper paired (for the later case, I'm not sure the exactly definition of those "non-proper pairs" ...). For samples merely contaminated, apply "-A" usually works better than turn it off (as those "non-proper" reads in one sample could occasionally become "proper" in another sample, which would cause failure in removing those "non-specific" calls); however, for samples not "clean" enough, apply this option would result in more "contaminated" calls. Several workarounds could be used:
	  1. Run mpileup twice, **with or without "-A"**, and **only trust the intersection** of final mutation calls from both set;
	  2. Only remove reads with mate unmapped before feed the bam file to mpileup and apply "-A" in mpileup, as those contaminations usually diverged from reference genome, it's less likely for both pairs to be mapped, this could be done by  
		* 

			    samtools view -b -F 3852 -L samples.fq3.snp.bed PREFIX.bam \
			    	> samples.fq3.snp.PF.bam
			    
			    samtools mpileup -Ad100000 -q20 -f reference.fasta -l samples.fq3.snp.bed \
			    	samples.fq3.snp.PF.bam | grep -vP "\t0\t" \
			    	> readcounts/samples.fq3.snp.${sample}.mpileup 
     
	4. It's also suggested to generate a third set with **"-q0"** instead of **"-q20"** if the sample size is limited, an **intersection between different mapping quality thresholds** could also help to reduce the false positive rate while only marginally increase the false negative rate. 

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


* Run HaplotypeCaller in **multiple-sample calling** mode around **indels** to generate more accurate AD values for each sample

		Samples=`find bam_files/ -name "*.bam" -print | xargs -I BAM_FILE echo -n "-I BAM_FILE "`
		java -jar GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
		    -R reference.fasta \
		    -T HaplotypeCaller -nct 1 -stand_call_conf 30.0 -stand_emit_conf 30.0 \
		    -L samples.fq3.indel.bed -o samples.fq3.indel.hc_multi.vcf \
		    ${Samples} 2>& 1 | tee samples.fq3.indel.hc_multi.log
		
		bgzip samples.fq3.indel.hc_multi.vcf && tabix -p vcf samples.fq3.indel.hc_multi.vcf.gz

<br />

####Step3: Screen out candidate mutations

* Use detect_mutations.pl to screen out candidate **point mutations**

		detect_mutations.pl -v samples.fq3.snp.vcf.gz --max-cmp-depth 2 --max-cmp-total 3 \
		    --max-cmp-miss 5 --min-supp-depth 5 --min-supp-plus 1 --min-supp-minus 1 \
		    --mask-only LowDepth -g samples.group.txt | \
		    vcf-annotate -f c=3,150 --fill-type > samples.fq3.snp.mut.c2t3d5m5.vcf


* Use detect_mutations.pl to screen out candidate **indel mutations**

		vcf_process.pl --vcf samples.fq3.indel.hc_multi.vcf.gz --quality 30 --var-type indel | \
		    awk '/\#/ || ($9 ~ /AD/)' | \
		    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 \
		    --max-cmp-miss 5 --mask-only LowDepth -g samples.group.txt | \
		    vcf-annotate -f c=2,150 --fill-type > samples.fq3.indel.mut.c2t3d5m5.vcf

<br />

####Step4: Rerun Step1~3 with another variant sets from different callers or even mappers, preferentially those implemented with a different algorithm

<br />

####Step5: Collect all cadidate mutations from different mappers, callers or different parameters

* Combine various results into a single vcf file use vcf_process.pl  

		vcf_process.pl --vcf caller1.mut.vcf --secondary-vcf caller2.mut.vcf \
		    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag CALLER1 --secondary-tag CALLER2 \
		    --intersect-tag CALLED_BOTH > mut.combined.vcf

<br />

####Step6: generate alignments and figures for manually inspections


* **Local re-align reads to reference** with another aligner like ClustalW2 (http://www.clustal.org/clustal2/)   

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




* Generate figures of **bam alignments** in IGV (https://www.broadinstitute.org/igv/)

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

<br />

##Scripts


####fillVcfDepth.pl   
> Add read numbers counted using VarScan readcounts to vcf file 

* **Options:**   

		-v, --vcf    <filename>  
			input vcf file, required
		-o, --output <filename>   
			output filename, default to STDOUT
		
		-l, --list   <filename>
			file contains list of depth files and related sample names in the format:
				"sample_id library_id filename"
			delimited by tab(s) or space(s), required
		
		-f, --filter <strings>
			skip filter loci, can have multiple values, separate by blanks, 
			e.g. "LowQual SNPFilter" ...    [default: no filtering]
		-p, --phased
			skip unphased sites
		
		-t, --type   <string>
			set "snp" to process snp sites only, or set "indel" to process indels only
		
		-u, --update-AD
			update AD field according to the results of readcounts
		
		-m, --minimum-vcf
			remove original INFO and FORMAT fields

<br />

####detect_mutations.pl 
> Screen out candidate sample/group-specific mutations

* **Options:**   

		    -v, --vcf     <filename>
		        input vcf file, required
		    -o, --output  <filename>
		        output filename, default to STDOUT
		
		    -f, --filter  <strings>
		        skip filter loci, can have multiple values, separate by space, e.g.
		        "LowQual SNPFilter ..."
		    -M, --match   <strings>
		        only retain loci matches, can have multiple values, separate by space,
		        e.g. "PASS ..."
		    -q, --quality     <float>
		        loci with quality smaller than this value will filtered
		    
		    -g, --group-file  <file>
		        file contain group infos of each sample, each sample per line, e.g.
		        sample1 group1
		        sample2 group1
		        sample3 group2
		        ...
		        set this option to screen group-specific mutation, only samples belong
		        to different groups would be used as compare samples
		        
		    --max-shared-freq <int>
		        locus with allele frequency below this value will be considering as a
		        shared mutation locus
		    
		    --min-supp-depth  <int>
		        minimum number of supporting reads [default: 1]
		    --min-supp-plus   <int>
		        minimum number of supporting reads in plus strand
		    --min-supp-minus  <int>
		        minimum number of supporting reads in minus strand
		    --min-lib-depth   <int>
		        minimum number of supporting reads in each library
		    --min-lib-cnt     <int>
		        minimum number of supporting libraries
		    --max-cmp-miss    <int>
		        maximum allowed missing alleles in compare samples
		    --no-ref-mut
		        remove mutations with reference allele
		    
		    --max-cmp-depth   <int>
		        maximum allowed number of reads contain same mutation base (termed as
		        mutation-like base here, possible from sequencing or mapping errors)
		        in each compared sample, set --max-cmp-depth to 2 means none of compare
		        samples could contain more than 2 mutation-like reads (e.g., FPD<=2 for
		        each compared sample), a value of 3 indicates no more than 3, etc.
		    --max-cmp-perc    <float>
		        maximum allowed percentage of reads contain mutation-like base in
		        each compared sample, this is an alternative option of "--max-cmp-depth"
		    
		    --max-cmp-total   <int>
		        maximum allowed number of total reads containing mutation-like base
		        across all compared samples
		        
		    --max-cmp-freq    <string>
		        detailed settings of maximum allowed number of reads with the mutation
		        -like base in compare samples at each depth, this option should be
		        setted according to "--max-cmp-depth" option, e.g.
		        if "--max-cmp-depth" is set to 2, this option should have 2 values
		        seperated by comma, like "3,1", which means at most 3 samples can have
		        1 read, while only 1 sample can have 2 reads
		    
		    --controls  <strings>
		        specify samples served as controls where no missing calls is allowed,
		        and shared mutations contain those samples will be filtered
		    
		    --mask-only <strings>
		        set proper FILTER field for those records failed given criteria rather
		        than remove them, can have multiple values, support filtering types:
		        LowDepth (--min-supp-depth)
		        StrandBias (--min-supp-plus or --min-supp-minus)
		        HighMissing (--max-cmp-miss)
		        NonSpecific (--max-cmp-total)
		        NoControl (--controls)
		        Shared

	**Note:**  
	1. This script is designed to process vcf files with AD (Allele Depth) field for each sample, and mainly used for processing output file from another script fillVcfDepth.pl, which could give all required informations used in this script;
	2. "--max-cmp-depth" or "--max-cmp-perc" are cutoff values which actually determines whether an sample will be treated as contain "candidate mutation" (contain mutation-like reads more than these thresholds, a true mutation allele) or belong to compare samples (contain mutation-like reads below this thresholds, possible sequencing errors). A lower value will considering more samples exceed this thresholds as "candidate mutation" samples, e.g. --max-cmp-depth 0 indicates samples contain any mutation-like reads are possible candidates, in other words compared samples could not have any mutation-like reads. A lower thresholds would give a higher false negative due to less tolerant of sequencing or mapping errors, while a higher thresholds could lead to slightly more false positives. In practice,"--max-cmp-detph 2" will mostly works fine.


<br />

##Publications

**More detailed descriptions of this pipeline could be also be found in the below publications:**

1. Yang, S., Wang, L., Huang, J., Zhang, X., Yuan, Y., Chen, J.-Q., Hurst, L.D., and Tian, D. (2015). Parent-progeny sequencing indicates higher mutation rates in heterozygotes. *Nature* 523: 463–467.

2. Xie, Z., Wang, L., Wang, L., Wang, Z., Lu, Z., Tian, D., Yang, S., and Hurst, L.D. (2016). Mutation rate analysis via parent–progeny sequencing of the perennial peach. I. A low rate in woody perennials and a higher mutagenicity in hybrids. *Proc. R. Soc. B* 283: 20161016.


