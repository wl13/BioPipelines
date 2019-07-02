## random simulation to test mutation clusters

awk 'BEGIN{OFS="\t";} {print $1,0,$2;}' \
    /opt/nfs/share/data/ref/neurospora/Broadinstitute/NC12_rename/neurospora_crassa.fasta.fai \
    > /opt/nfs/share/data/ref/neurospora/Broadinstitute/NC12_rename/neurospora_crassa.bed

bedtools subtract -a /opt/nfs/share/data/ref/neurospora/Broadinstitute/NC12_rename/neurospora_crassa.bed \
    -b /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.bed \
    > /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.out.bed

##
## site simulation
##
bedtools maskfasta -fi /opt/nfs/share/data/ref/neurospora/Broadinstitute/NC12_rename/neurospora_crassa.fasta \
    -bed /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.out.bed \
    -fo /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.masked.fasta

parseSEQs.pl --fasta /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.masked.fasta \
    --rate --sum-all \
    > /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.masked.nt_rate.csv


time awk '/\#/ || (/RepeatState=Repeat;/ && ($4 == "G" || $4 == "C"))' \
    /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.vcf | \
    test_var_clusters.pl --query-file - -d 1000 -s 1 -t 100 \
    -f /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.masked.fasta \
    > /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/clusters/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.c1k.sim2.csv



##
## region simulation
##

awk 'BEGIN{OFS="\t";} !/\#/ && $7 > 1 {print $3,$4-1,$5,$7,"Observed";}' \
    /home/wl/Data/data/neurospora/03.analysis/02.mutations/NC_all/clusters/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.c1k.csv | \
    sort -k1,1 -k2,3n \
    > /home/wl/Data/data/neurospora/03.analysis/02.mutations/NC_all/clusters/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.c1k.o2.bed

cat /home/wl/Data/data/neurospora/03.analysis/02.mutations/NC_all/clusters/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.c1k.o2.bed | \
    bedtools merge -c 4,5 -o sum,collapse -i - | \
    perl -e 'my $olp_cnt = 0; $olp_lengths = 0; my $olp_mut = 0;
        print "#source\tolp_cnt\tolp_lengths\tolp_mut\n";
        while(<>){
            my @line = (split /\s+/); my $freq = ($line[4] =~ tr/\,/\,/) + 1; next unless($freq > 1);
            my $length = $line[2] - $line[1]; $olp_lengths += $length;
            $olp_cnt++; $olp_mut += $line[3];
        }
        print "Observed\t$olp_cnt\t$olp_lengths\t$olp_mut\n";' \
    > /home/wl/Data/data/neurospora/03.analysis/02.mutations/NC_all/clusters/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.c1k.olp.csv


for i in {1..10000};
do
    echo -ne "${i}\r"
    bedtools shuffle -chrom \
        -excl /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/duplications/blast/neurospora_crassa.chrom.w200.self_blastn.o65.out.bed \
        -i /home/wl/Data/data/neurospora/03.analysis/02.mutations/sexual/clusters/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.c1k.o2.bed \
        -g /opt/nfs/share/data/ref/neurospora/Broadinstitute/NC12_rename/neurospora_crassa.fasta.fai | \
        awk -v time=$i 'BEGIN{OFS="\t";} {print $1,$2,$3,$4,"Rand"time}' | \
        sort -k1,1 -k2,3n | bedtools merge -c 4,5 -o sum,collapse -i - | \
        perl -e 'my $olp_cnt = 0; $olp_lengths = 0; my $olp_mut = 0;
            while(<>){
                my @line = (split /\s+/); my $freq = ($line[4] =~ tr/\,/\,/) + 1; next unless($freq > 1);
                my $length = $line[2] - $line[1]; $olp_lengths += $length;
                $olp_cnt++; $olp_mut += $line[3];
            }
            print "Rand\t$olp_cnt\t$olp_lengths\t$olp_mut\n";' \
        >> /home/wl/Data/data/neurospora/03.analysis/02.mutations/NC_all/clusters/NC_s268.NC12.bwa.sort.dedup.realn.mut_hyw_wr.c1k.olp.csv
done
