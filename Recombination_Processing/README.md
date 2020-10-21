### Pipelines to detect recombination events (Crossover or Gene conversions) from HTS data


Here mainly focused on F2 samples, should also be easily adapted to F3, F4 ... or even inbred lines.


## Pipelines   

#### Step1: Screen out confident bi-allelic markers distinguishable between parental ("crossed") samples

* Use some hard filtering (e.g. variant quality, depth, missing rate, etc.) to remove any putative false variant calls

      vcf_process.pl --vcf F1_and_F2s.vcf.gz --quality 50 --min-alleles 2 --max-alleles 2 \
           --min-sample-depth 10 --max-sample-depth 80 --max-missing 18 > F1_and_F2s.flt.vcf
           
* A "marker" should be able to distinguish the parental source in each F2 sample, it's prefered to also have the parental samples or F1 sample sequenced. The example assumes the F1 sample was also sequenced, and appears at the first sample in the vcf file. "A" and "B" stands for different parental source. Fill with "./." if no information was available for certain samples, those alleles will ignored in those samples during downstream process.

      cat F1_and_F2s.flt.vcf | \
          perl -ne 'if(/\#/){print;next;} my @line = (split /\s+/);
          my ($f1_gt) = ($line[9] =~ /(\d\/\d):/); next unless(($f1_gt eq "0/1"));
          my @new_fields = (); $line[7] = ".";
          for (my $i=9;$i<@line;$i++) { my $sc = "./."; my $gt = "./.";
              if ($line[$i] =~ /((\d)\/(\d)):/) { $gt = $1;
                  if ($gt eq "0/1") { $sc = "A/B"; } 
                  elsif ($gt eq "0/0") { $sc = "A/A"; } 
                  else { $sc = "B/B"; }
              } push @new_fields, "$gt:$sc";
          } $line[8] = "GT:SC"; my $out_line = join "\t", @line[0..8];
          my $out_fields = join "\t", @new_fields; print "$out_line\t$out_fields\n";
      ' > F1_and_F2s.markers.vcf

* The output vcf file would have a new field "SC" represents different parental source of each sampleExample of the output vcf file. Here is an example:

		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
		chr01   161     .       C       A       54.12   .       .       GT:SC   0/0:A/A
		chr01   431     .       C       T       44.2    .       .       GT:SC   0/0:A/A
		chr01   1641    .       G       A       64.15   .       .       GT:SC   1/1:B/B
		chr01   4165    .       C       A       34.31   .       .       GT:SC   1/1:B/B
		...
  
	
<br />
   
#### Step2: Clustering markers to form recombinant blocks

This was a subfunction of vcf_process.pl, also see "Clustering variants" section in vcf_process.pl

Clustering was done by a "seeding-and-extension" approach (originally described in "Wijnker, E. et al. The genomic landscape of meiotic crossovers and gene conversions in Arabidopsis thaliana. eLife 2, e01426 (2013)").


* It's recommended to first generate a "raw" blocks by simply link consecutive markers with same source type, this would give an overview of all putative recombinants 

		vcf_process.pl --vcf markers.vcf.gz --out-blocks --source-tag "SC" > markers.blocks.csv

* Clustering would reduce noises due to false variant calls or putative gene conversion / mutation events, and generate purified recombinant backgrounds. Here use two criteria to determine a reliable "seed". However, the results largely depend on the quality of markers, the criteria matter less. In practice, use a minimum length threshold "--min-frag-length 10000" should mostly be fine, "--min-frag-markers" depends on the marker density, around 15 to 25 often give good results for a divergent over 0.3% (e.g. more than 3 markers per kilo bp in average).

		vcf_process.pl --vcf markers.vcf.gz --out-blocks --source-tag "SC" --fill-gaps \
		    --min-frag-length 10000 --min-frag-markers 25 > markers.blocks.l10km25.csv

* It's also recommended to test different criteria or compared to the "raw" blocks, a plot would be useful to visualize the differences

		awk 'BEGIN{OFS="\t";} !/\#/ {$1 = $1"-original"; print;}' markers.blocks.csv | \
		    cat markers.blocks.l10km25.csv - | \
		    paintGenomeBlocks.pl --input - \
		    --width 1600 --height 3000 --thickness 10 --chrom-distance 20 --block-distance 2 \
		    --output markers.blocks.cmp --length reference_genome.fasta.fai \
		    --colors "type1:strong_red2;B:strong_blue2" --sort-blocks sample-original sample --format png

* A much nicer plot could be generated using [scripts](https://github.com/jiaxianqing/Ploting_scripts) written by [Xianqing Jia](https://github.com/jiaxianqing/), below is an example taken from there:

```
<img src="https://github.com/jiaxianqing/Ploting_scripts/blob/master/examples/plot.block.png"  div align = "center" width="75%" height="75%" />

<br />

#### Step3: Refine the boundaries

* One may need some manual works to get the most accuracy boundaries needed, especially where markers were scarce. Take a look at the informations give in the block results, then went back to those filtered variants or even dig into the original mapping results. 


<br />

#### Step4: Detect crossover events

* Once the recombinat blocks were formed, then we knwo where the crossovers have occured

      detect_recomb_events.pl --rows 0 1 2 6 7 --min-co-len 500000 \
          --blocks markers.blocks.l10km25.csv \
          > markers.blocks.l10km25.co_500k.csv

<br />

#### Step5: Detect gene conversions

* Find the converted genotypes differs to its background, most times one need to further purge the markers (i.e., markers.flt.vcf.gz) to reduce any possible false positives

        detect_recomb_events.pl --blocks markers.blocks.l10km25.csv \
            --vcf markers.flt.vcf.gz \
            --source-tag SC --min-markers 1 --min-gc-len 1 --max-gc-len 10000 --min-bg-len 10000 \
            --output gene_convs.csv


<br />

#### Additional Step: Phasing blocks

* If the parental samples were not sequenced, one may need a phasing step to assign the source of blocks for each sample. This could be done by a first-round clustering to only distinguish "homozygous" and "heterozygous" regions, and then phasing each chromosome based on "least occurrence of COs". 


<br />

## Publications

**More details about the whole pipeline and the phasing stage could be found in the below publications:**

1. Wang, L., Zhang, Y., Qin, C., Tian, D., Yang, S., and Hurst, L.D. (2016). Mutation rate analysis via parent–progeny sequencing of the perennial peach. II. No evidence for recombination-associated mutation. Proc. R. Soc. B 283: 20161785.



