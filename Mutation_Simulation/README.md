####Generate synthetic (in silico) mutations from real sequencing reads. 
Those synthetic mutations could be used to test the false-negative rate of a certain mutation-detection pipeline.

####Original methods described in:

* 1) Keightley, P. D., Ness, R. W., Halligan, D. L. & Haddrill, P. R. Estimation of the Spontaneous Mutation Rate per Nucleotide Site in a Drosophila melanogaster Full-Sib Family. Genetics 196, 313–320 (2014).

* 2) Keightley, P. D. et al. Estimation of the Spontaneous Mutation Rate in Heliconius melpomene. Mol Biol Evol 32, 239–243 (2015).


####Pipelines

* **Step1:** get the empirical distributions of the depth of non-reference alleles for all heterozygous genotypes in 20 sequenced samples

		bcftools filter -i 'QUAL>=50 && TYPE="snp" && MAF>=0.05' \       ## try to reduce false positives
		    20_samples.hc.vcf.gz | \                                     ## more filtering could be added
		    perl -e 'my %alt_depths = ();                                ## sampling all heterozygous loci
		      while(<>){
		        next if (/^\#/);
		        my ($chrom, $pos, $id, $ref, $alt, $qual, 
		        	  $filter, $info, $format, @samples) = (split /\t/);
		        for my $sample (@samples) {
		            next unless($sample =~ /(\d+\/\d+):(\d+),(\d+)/);    ## make use of the AD field
		            my ($GT, $REF_DP, $ALT_DP) = ($1, $2, $3);
		            next unless($GT eq "0/1"); 
		            my $total_DP = $REF_DP + $ALT_DP; 
		            $alt_depths{$total_DP}->{$ALT_DP}++;
		          }
		      }
		      for my $depth (sort {$a <=> $b} keys %alt_depths) { 
		        for my $alt_depth (sort {$a <=> $b} keys %{$alt_depths{$depth}}) {
		          print "$depth\t$alt_depth\t$alt_depths{$depth}->{$alt_depth}\n";
		        }
		      }' > 20_samples.hc.depths.csv


* **Step2:** generate synthetic reads (reads with synthetic mutations) for all 20 sequenced samples from mapping results

		## Get all bam files
		BAM_FILES=`find examples/ -name "*.bam" -print | \
				xargs -I BAM echo -n "BAM "`
		
		## Random select a genome position with a read depth x, then replace
		## y reads with a random picked nucleotide (different from the original
		## one), the y was determined according to the empirical distribution.
		## Reads with following situations were not counted in x (FLAG 3844)
		##    read unmapped
		##    not primary alignment
		##    read fails platform/vendor quality checks
		##    read is PCR or optical duplicate
		##    supplementary alignment
		sim_mutation_reads.pl --fasta reference.fasta \
		    --depth 20_samples.hc.depths.csv \
		    --bams ${BAM_FILES} --random-size 1000 --samtools "-F 3844" \
		    > 20_samples.simulated.dat


* **Step3:** extract not only the synthetic reads, but also its paired reads

		## Extract synthetic nucleotide changes
		awk 'BEGIN{OFS="\t"} /^\#Tag/ || $1 == "MUT" {if(/\#Tag/){$2 = "#Chrom";} print;}' \
		    20_samples.simulated.dat | cut -f 2- | sort -k1,1 -k2,2n > 20_samples.simulated.vars.csv
		
		## Extract reads carry the synthetic mutations
		awk '$1 ~ /SAM:/' 20_samples.simulated.dat | cut -f 2- > 20_samples.simulated.reads.sam
		
		## Generate synthetic reads
		find examples/ -name "*.bam" -print | xargs -n 1 -P 4 -I BAM_FILE sh -c '
		    sample=`basename BAM_FILE | cut -d"." -f1`
		    
		    echo "${sample}"
		    
		    ## Extract the read pairs cover all synthetic mutations
		    extract_bam_pairs.pl --samtools "-F 3844" --extend 10000 --bam BAM_FILE \
		        --input 20_samples.simulated.vars.csv --patches 20_samples.simulated.reads.sam \
		        > reads/${sample}_ex10k.sam 2> reads/${sample}_ex10k.log
		    
		    ## Convert sam file to fastq reads
		    sam2fastq.pl -i reads/${sample}_ex10k.sam \
		        -o reads/${sample}_ex10k >> reads/${sample}_ex10k.log 2>&1
		'


####Known Issues:
1) A certain portion of synthetic mutations could appear upon a indel region in the sequenced sample, those "mutations" could appear in a different position nearby.

