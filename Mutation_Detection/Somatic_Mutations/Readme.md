## Pipelines and scripts to identify and analysis somatic mutations in plant tissues e.g. leaves, root and etc.

<br />

### somatic-mutation-detection.md

An adaption of previous pipeline to deal with somatic mutations

<br />

### somatic-inheritance-analysis.md

Test inheritance of somatic mutations

<br />

## Miscellaneous scripts

<br />

#### fillVcfDepth-split.pl

> An alternative of fillVcfDepth.pl

* Options:

        -v, --vcf    <filename>
            input vcf file, required
        -o, --output <filename>
            output filename, default to STDOUT

        -l, --list   <filename>
            file contains list of depth files and related sample names in the
            format:

            sample_id library_id filename

            delimited by tab(s) or space(s), required

            *Note:
            All files specified in this list should be compressed with bgzip and
            indexed using tabix:
            bgzip example.readcounts
            tabix -S1 -s1 -b2 -e2 example.readcounts.gz

        -f, --filter <strings>
            skip filter loci, can have multiple values, separate by blanks, e.g.
            "LowQual SNPFilter" ...    [default: no filtering]
            -p, --phased
            skip unphased sites

        -t, --type   <string>
            set "snp" to process snp sites only, or set "indel" to process indels
            only

        -u, --update-AD
            update AD field according to the results of readcounts

        -m, --minimum-vcf
            remove original INFO and FORMAT fields



<br />

#### test_uniformity.pl
> D (Variance/Mean) test by comparing observed D and simulated D\
> Code for random simulation was borrowed from the post in https://stackoverflow.com/a/22381248 by ikegami (the second approach) and modified to allow for zero-values

* **Options:**   

            -q, --query-file  <string>
                file contain list of numbers, required

            -o, --output <filename>
                output file, default to STDOUT

            -t, --times
                random times [default: 1]
<br />


<br />

## Publications

**More detailed descriptions of this pipeline could also be found in the below publications:**

[1] L. Wang, Y. Ji, Y. Hu, H. Hu, X. Jia, M. Jiang, X. Zhang, L. Zhao, Y. Zhang, Y. Jia, C. Qin, L. Yu, J. Huang, S. Yang, L.D. Hurst, D. Tian, The architecture of intra-organism mutation rate variation in plants, <i>PLOS Biology</i>. 17 (2019) e3000191. doi:10.1371/journal.pbio.3000191.


