### Generate synthetic (in silico) mutations from real sequencing reads. 
Those synthetic mutations could be used to test the false-negative rate of a certain mutation-detection pipeline.

<br />

#### Original methods described in:

* 1) Keightley, P. D., Ness, R. W., Halligan, D. L. & Haddrill, P. R. Estimation of the Spontaneous Mutation Rate per Nucleotide Site in a Drosophila melanogaster Full-Sib Family. Genetics 196, 313–320 (2014).

* 2) Keightley, P. D. et al. Estimation of the Spontaneous Mutation Rate in Heliconius melpomene. Mol Biol Evol 32, 239–243 (2015).

<br />


## Pipelines

#### Step1: get the empirical distributions of the depth of non-reference alleles

* sampling non-reference read-depth of heterozygous loci (this requires the AD field)

		bcftools filter -i 'QUAL>=50 && TYPE="snp" && MAF>=0.05' \
		    20_samples.hc.vcf.gz | \
		    perl -e 'my %alt_depths = ();
		      while(<>){
		        next if (/^\#/);
		        my ($chrom, $pos, $id, $ref, $alt, $qual, 
		        	  $filter, $info, $format, @samples) = (split /\t/);
		        for my $sample (@samples) {
		            next unless($sample =~ /(\d+\/\d+):(\d+),(\d+)/);
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

<br />

#### Step2: generate synthetic reads (reads with synthetic mutations) for all 20 sequenced samples from mapping results

* Get all bam files

		BAM_FILES=`find examples/ -name "*.bam" -print | \
				xargs -I BAM echo -n "BAM "`

* Random select a genome position with a read depth x, then replace y reads with a random picked nucleotide (different from the original one), the y was determined according to the empirical distribution. Generate reads with synthesized mutations, excluding non-informative reads (-F 3844):    
	1. unmapped;   
	2. not primary alignment;   
	3. read fails platform/vendor quality checks;   
	4. read is PCR or optical duplicate;   
	5. supplementary alignment;      

* 
   
		sim_mutation_reads.pl --fasta reference.fasta \
		    --depth 20_samples.hc.depths.csv \
		    --bams ${BAM_FILES} --random-size 1000 --samtools "-F 3844" \
		    > 20_samples.simulated.dat

<br />

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

<br />


## Scripts

#### sim_mutation_reads.pl   
> Script used to generate synthetic mutated sites from original bam files 

* **Options:**   

		Input Options:
		
		    -f, --fasta       <filename>
		        reference sequences in fasta format, required
		    
		    -d, --depth       <filename>
		        input file contain empirical distribution of read depths, required
		
		    -b, --bams        <filename>
		        bam file(s), at least one bam file should be specified
		
		
		Output Options:
		
		    -o, --output <filename>
		        output filename, default to STDOUT
		        
		    -n, --no-rc
		        do not reverse complement sequence with negtive strand
		    -u, --use-rg
		        add read group id to extracted records
		        
		Simulation Options:
		    --random-size <int>
		        number of mutations to be simulated, default: 100
		
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
		
		    --max-shared   <int>
		        choose how many samples could shared the mutation loci, the samples
		        were also choosed by random [default: 1]
		
		    -e, --exclude  <strings>
		        exclude unwanted chromosomes or scaffolds while simulating, all
		        chromosomes with ids match strings specified here would be ignored 
		
		    --min-depth    <int>
		        minimum required depth of "mutated" reads, loci with a depth smaller
		        than this threshold would be set to this threshold. This option would
		        also cause those loci where no mutated reads were successful generated
		        to be skipped [default: 0]
		        
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

<br />

