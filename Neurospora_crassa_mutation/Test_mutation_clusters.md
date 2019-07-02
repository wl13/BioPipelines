### Simulation 1: Random simulate mutations within duplicates

#### Step 1) Prepare a genome file with non-duplicates been masked 

        bedtools subtract -a neurospora_crassa.chrom.bed \
                -b neurospora_crassa.duplicates.bed \
                > neurospora_crassa.non-duplicates.bed


        bedtools maskfasta -fi neurospora_crassa.fasta \
            -bed neurospora_crassa.non-duplicates.bed \
            -fo neurospora_crassa.non-duplicates.mask_non-dup.fasta


#### Step 2) Generate same amount of mutations within duplicates to see how many of them were clustered

        cat C_site_mutations_in_duplicates.vcf | \
            test_var_clusters.pl --query-file - -d 1000 -s 1 -t 10000 \
            -f neurospora_crassa.non-duplicates.mask_non-dup.fasta \
            > C_site_mutations_in_duplicates.test.csv




### Simulation 2:  Shuffling clustered regions within duplicates

#### Step 1) Check mutation clusters in observed results

        cat C_site_mutations_in_duplicates.vcf | \
            test_var_clusters.pl --query-file - -d 1000 -s 2 \
            > C_site_mutations_in_duplicates.clusters.csv

        awk 'BEGIN{OFS="\t";} !/\#/ && $7 > 1 {print $3,$4-1,$5,$7,"Observed";}' \
            C_site_mutations_in_duplicates.clusters.csv | \
            sort -k1,1 -k2,3n \
            > C_site_mutations_in_duplicates.clusters.bed

#### Step 2) Check overlapped clusters in observed results

        cat C_site_mutations_in_duplicates.clusters.bed | \
            bedtools merge -c 4,5 -o sum,collapse -i - | \
            perl -e 'my $olp_cnt = 0; $olp_lengths = 0; my $olp_mut = 0;
                print "#source\tolp_cnt\tolp_lengths\tolp_mut\n";
                while(<>){
                    my @line = (split /\s+/); my $freq = ($line[4] =~ tr/\,/\,/) + 1; next unless($freq > 1);
                    my $length = $line[2] - $line[1]; $olp_lengths += $length;
                    $olp_cnt++; $olp_mut += $line[3];
                }
                print "Observed\t$olp_cnt\t$olp_lengths\t$olp_mut\n";' \
            > C_site_mutations_in_duplicates.clusters.olp.csv


#### Step 3) Simulate clusters by shuffling clustered regions within duplicates

        for i in {1..10000};
        do
            echo -ne "${i}\r"
            bedtools shuffle -chrom \
                -excl neurospora_crassa.non-duplicates.bed \
                -i C_site_mutations_in_duplicates.clusters.bed \
                -g neurospora_crassa.chrom.fasta.fai | \
                awk -v time=$i 'BEGIN{OFS="\t";} {print $1,$2,$3,$4,"Rand"time}' | \
                sort -k1,1 -k2,3n | bedtools merge -c 4,5 -o sum,collapse -i - | \
                perl -e 'my $olp_cnt = 0; $olp_lengths = 0; my $olp_mut = 0;
                    while(<>){
                        my @line = (split /\s+/); my $freq = ($line[4] =~ tr/\,/\,/) + 1; next unless($freq > 1);
                        my $length = $line[2] - $line[1]; $olp_lengths += $length;
                        $olp_cnt++; $olp_mut += $line[3];
                    }
                    print "Rand\t$olp_cnt\t$olp_lengths\t$olp_mut\n";' \
                >> C_site_mutations_in_duplicates.clusters.olp.csv
        done


