02-methylation-summaries
================
Kathleen Durkin
2025-10-27

- [1 Create a Bash variables file](#1-create-a-bash-variables-file)
- [2 Download BAM files](#2-download-bam-files)
  - [2.0.1 Group 1, Library 4, MinION](#201-group-1-library-4-minion)
  - [2.0.2 Group 4, Library 1, MinION](#202-group-4-library-1-minion)
  - [2.0.3 Group 4, Library 2, MinION](#203-group-4-library-2-minion)
- [3 Methylation summaries from original BAM
  files](#3-methylation-summaries-from-original-bam-files)
  - [3.1 modkit](#31-modkit)
    - [3.1.1 Group 1, Library 4, MinION](#311-group-1-library-4-minion)
    - [3.1.2 Group 4, Library 1, MinION](#312-group-4-library-1-minion)
    - [3.1.3 Group 4, Library 2, MinION](#313-group-4-library-2-minion)
- [4 Summary](#4-summary)

# 1 Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export nanopore_dir=/home/shared/8TB_HDD_02/shedurkin/SIFP-nanopore'
echo 'export output_dir_top=${nanopore_dir}/M-multi-group/output/02-methylation-summaries'
echo 'export data_dir_top=${nanopore_dir}/M-multi-group/data/02-methylation-summaries'

echo "# Paths to programs"
echo 'export modkit=/home/shared/dist_modkit_v0.5.1_8fa79e3/modkit'
echo 'export samtools=/home/shared/samtools-1.12/samtools'
echo 'export conda=/home/shared/8TB_HDD_02/shedurkin/.local/share/r-miniconda/bin/conda'
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=20'
echo ""

echo "# Input/output files"
echo 'export raw_checksums=checksums.md5'
echo 'export trimmed_checksums=trimmed_fastq_checksums.md5'
echo ""
} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export nanopore_dir=/home/shared/8TB_HDD_02/shedurkin/SIFP-nanopore
    export output_dir_top=${nanopore_dir}/M-multi-group/output/02-methylation-summaries
    export data_dir_top=${nanopore_dir}/M-multi-group/data/02-methylation-summaries
    # Paths to programs
    export modkit=/home/shared/dist_modkit_v0.5.1_8fa79e3/modkit
    export samtools=/home/shared/samtools-1.12/samtools
    export conda=/home/shared/8TB_HDD_02/shedurkin/.local/share/r-miniconda/bin/conda

    # Set number of CPUs to use
    export threads=20

    # Input/output files
    export raw_checksums=checksums.md5
    export trimmed_checksums=trimmed_fastq_checksums.md5

# 2 Download BAM files

### 2.0.1 Group 1, Library 4, MinION

Barcodes 12, 13, and 14

``` bash
source .bashvars

wget \
--directory-prefix ${data_dir_top}/bam_pass/G1L4_MinION \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 6 \
--no-host-directories \
--no-parent \
--quiet \
--accept *.bam \
https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group1_MinION/Library4/20250819_1337_MD-101223_FBD09922_51407a5d/bam_pass
```

### 2.0.2 Group 4, Library 1, MinION

Barcodes 27, 28, 29

Note that G4L1 MinION ended up split into 2 nanopore runs, due to
erroring out between washes

``` bash
source .bashvars

wget \
--directory-prefix ${data_dir_top}/bam_pass/G4L1_MinION \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 6 \
--no-host-directories \
--no-parent \
--quiet \
--accept *.bam \
https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group4_MinION/Library1/20250905_1158_MD-101223_FBD08455_b87a1f92/bam_pass

wget \
--directory-prefix ${data_dir_top}/bam_pass/G4L1_MinION_wash3 \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 6 \
--no-host-directories \
--no-parent \
--quiet \
--accept *.bam \
https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group4_MinION/Library1_wash3/20250907_1505_MD-101223_FBD08455_50a8f15e/bam_pass
```

### 2.0.3 Group 4, Library 2, MinION

Barcodes 30, 31, 32

``` bash
source .bashvars

wget \
--directory-prefix ${data_dir_top}/bam_pass/G4L2_MinION \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 6 \
--no-host-directories \
--no-parent \
--quiet \
--accept *.bam \
https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group4_MinION/Library2/20250909_1146_MN47571_FBD36396_e03e7ede/bam_pass
```

# 3 Methylation summaries from original BAM files

Using BAM files output from sequencing runs, *not* from recalling the
pod5 files

## 3.1 modkit

### 3.1.1 Group 1, Library 4, MinION

Try merging all bam files for each barcode to obtain unified summaries

``` bash
source .bashvars

### Barcode 12 ###
mergedbam=${output_dir_top}/G1L4_MinION_barcode12_merged.bam
# Merge all bam files from the given barcode into one
$samtools merge -u $mergedbam ${data_dir_top}/bam_pass/G1L4_MinION/barcode12/*.bam
# Total number of reads included in the merged file
$samtools view -c $mergedbam
# Sample reads and summarize methylation
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G1L4_MinION_barcode12_modkit_summary.txt
# Delete merged bam file (space-saving)
rm $mergedbam

### Barcode 13 ###
mergedbam=${output_dir_top}/G1L4_MinION_barcode13_merged.bam
$samtools merge -u $mergedbam ${data_dir_top}/bam_pass/G1L4_MinION/barcode13/*.bam
$samtools view -c $mergedbam
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G1L4_MinION_barcode13_modkit_summary.txt
rm $mergedbam

### Barcode 14 ###
mergedbam=${output_dir_top}/G1L4_MinION_barcode14_merged.bam
$samtools merge -u $mergedbam ${data_dir_top}/bam_pass/G1L4_MinION/barcode14/*.bam
$samtools view -c $mergedbam
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G1L4_MinION_barcode14_modkit_summary.txt
rm $mergedbam
```

    2065657
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: C: 0.8046875
    2489352
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: C: 0.8125
    4090390
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: C: 0.82421875

``` bash
source .bashvars

cat ${output_dir_top}/G1L4_MinION_barcode12_modkit_summary.txt
cat ${output_dir_top}/G1L4_MinION_barcode13_modkit_summary.txt
cat ${output_dir_top}/G1L4_MinION_barcode14_modkit_summary.txt
```

    # bases             C 
    # total_reads_used  100000 
    # count_reads_C     100000 
    # pass_threshold_C  0.8046875 
     base  code  pass_count  pass_frac    all_count  all_frac 
     C     -     1556909     0.92845637   1664842    0.8941322 
     C     h     18106       0.01079744   58699      0.031525314 
     C     m     101864      0.060746185  138423     0.07434247 
    # bases             C 
    # total_reads_used  100000 
    # count_reads_C     100000 
    # pass_threshold_C  0.8125 
     base  code  pass_count  pass_frac     all_count  all_frac 
     C     -     1471527     0.92578214    1573397    0.8910261 
     C     h     16036       0.0100887325  54004      0.030582855 
     C     m     101933      0.06412913    138425     0.078391075 
    # bases             C 
    # total_reads_used  100000 
    # count_reads_C     100000 
    # pass_threshold_C  0.82421875 
     base  code  pass_count  pass_frac    all_count  all_frac 
     C     -     1621732     0.93584156   1733758    0.9014899 
     C     h     14749       0.008511103  53815      0.027981805 
     C     m     96432       0.05564734   135641     0.07052829 

### 3.1.2 Group 4, Library 1, MinION

Try merging all bam files for each barcode to obtain unified summaries

``` bash
source .bashvars

### Barcode 27 ###
mergedbam1=${output_dir_top}/prewash3_G4L1_MinION_barcode27_merged.bam
mergedbam2=${output_dir_top}/wash3_G4L1_MinION_barcode27_merged.bam
mergedbam=${output_dir_top}/G4L1_MinION_barcode27_merged.bam
# Merge all bam files from the given barcode into one
$samtools merge -u $mergedbam1 ${data_dir_top}/bam_pass/G4L1_MinION/barcode27/*.bam
$samtools merge -u $mergedbam2 ${data_dir_top}/bam_pass/G4L1_MinION_wash3/barcode27/*.bam
$samtools merge -u $mergedbam ${output_dir_top}/*G4L1_MinION_barcode27_merged.bam
# Total number of reads included in the merged file
$samtools view -c $mergedbam
# Sample reads and summarize methylation
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G4L1_MinION_barcode27_modkit_summary.txt
# Delete merged bam file (space-saving)
rm $mergedbam1
rm $mergedbam2
rm $mergedbam

### Barcode 28 ###
mergedbam1=${output_dir_top}/prewash3_G4L1_MinION_barcode28_merged.bam
mergedbam2=${output_dir_top}/wash3_G4L1_MinION_barcode28_merged.bam
mergedbam=${output_dir_top}/G4L1_MinION_barcode28_merged.bam
$samtools merge -u $mergedbam1 ${data_dir_top}/bam_pass/G4L1_MinION/barcode28/*.bam
$samtools merge -u $mergedbam2 ${data_dir_top}/bam_pass/G4L1_MinION_wash3/barcode28/*.bam
$samtools merge -u $mergedbam ${output_dir_top}/*G4L1_MinION_barcode28_merged.bam
$samtools view -c $mergedbam
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G4L1_MinION_barcode28_modkit_summary.txt
rm $mergedbam1
rm $mergedbam2
rm $mergedbam

### Barcode 29 ###
mergedbam1=${output_dir_top}/prewash3_G4L1_MinION_barcode29_merged.bam
mergedbam2=${output_dir_top}/wash3_G4L1_MinION_barcode29_merged.bam
mergedbam=${output_dir_top}/G4L1_MinION_barcode29_merged.bam
$samtools merge -u $mergedbam1 ${data_dir_top}/bam_pass/G4L1_MinION/barcode29/*.bam
$samtools merge -u $mergedbam2 ${data_dir_top}/bam_pass/G4L1_MinION_wash3/barcode29/*.bam
$samtools merge -u $mergedbam ${output_dir_top}/*G4L1_MinION_barcode29_merged.bam
$samtools view -c $mergedbam
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G4L1_MinION_barcode29_modkit_summary.txt
rm $mergedbam1
rm $mergedbam2
rm $mergedbam
```

    484845
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: C: 0.6171875 A: 0.7324219
    385251
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: A: 0.7324219 C: 0.6113281
    943264
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: A: 0.7363281 C: 0.6015625

``` bash
source .bashvars

cat ${output_dir_top}/G4L1_MinION_barcode27_modkit_summary.txt
cat ${output_dir_top}/G4L1_MinION_barcode28_modkit_summary.txt
cat ${output_dir_top}/G4L1_MinION_barcode29_modkit_summary.txt
```

    # bases             A,C 
    # total_reads_used  100000 
    # count_reads_A     100000 
    # count_reads_C     100000 
    # pass_threshold_C  0.6171875 
    # pass_threshold_A  0.7324219 
     base  code   pass_count  pass_frac    all_count  all_frac 
     C     -      4709514     0.9370899    5037398    0.90344435 
     C     m      183430      0.036498543  288506     0.051742807 
     C     21839  132736      0.02641155   249866     0.044812825 
     A     -      6998666     0.9448404    7497537    0.91100186 
     A     a      408581      0.05515963   732454     0.088998154 
    # bases             C,A 
    # total_reads_used  100000 
    # count_reads_C     100000 
    # count_reads_A     100000 
    # pass_threshold_C  0.6113281 
    # pass_threshold_A  0.7324219 
     base  code   pass_count  pass_frac    all_count  all_frac 
     C     -      4599511     0.9394284    4927496    0.90604055 
     C     m      149560      0.030546924  244914     0.04503342 
     C     21839  147003      0.030024668  266084     0.048926044 
     A     -      7187994     0.9458192    7701594    0.912195 
     A     a      411761      0.054180827  741331     0.087804995 
    # bases             A,C 
    # total_reads_used  100000 
    # count_reads_A     100000 
    # count_reads_C     100000 
    # pass_threshold_A  0.7363281 
    # pass_threshold_C  0.6015625 
     base  code   pass_count  pass_frac    all_count  all_frac 
     A     -      7643002     0.94833636   8188903    0.91499233 
     A     a      416377      0.05166366   760793     0.08500769 
     C     -      5183410     0.9336648    5546190    0.90081483 
     C     m      207997      0.037465584  321663     0.05224466 
     C     21839  160275      0.028869629  289006     0.046940494 

### 3.1.3 Group 4, Library 2, MinION

Try merging all bam files for each barcode to obtain unified summaries

``` bash
source .bashvars

### Barcode 30 ###
mergedbam=${output_dir_top}/G4L2_MinION_barcode30_merged.bam
# Merge all bam files from the given barcode into one
$samtools merge -u $mergedbam ${data_dir_top}/bam_pass/G4L2_MinION/barcode30/*.bam
# Total number of reads included in the merged file
$samtools view -c $mergedbam
# Sample reads and summarize methylation
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G4L2_MinION_barcode30_modkit_summary.txt
# Delete merged bam file (space-saving)
rm $mergedbam

### Barcode 31 ###
mergedbam=${output_dir_top}/G4L2_MinION_barcode31_merged.bam
$samtools merge -u $mergedbam ${data_dir_top}/bam_pass/G4L2_MinION/barcode31/*.bam
$samtools view -c $mergedbam
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G4L2_MinION_barcode31_modkit_summary.txt
rm $mergedbam

### Barcode 32 ###
mergedbam=${output_dir_top}/G4L2_MinION_barcode32_merged.bam
$samtools merge -u $mergedbam ${data_dir_top}/bam_pass/G4L2_MinION/barcode32/*.bam
$samtools view -c $mergedbam
$modkit summary --num-reads 100000 $mergedbam > ${output_dir_top}/G4L2_MinION_barcode32_modkit_summary.txt
rm $mergedbam
```

    274217
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: A: 0.7402344 C: 0.59375
    [E::bgzf_read_block] Invalid BGZF header at offset 7768473
    [E::bgzf_read] Read block operation failed with error 2 after 0 of 4 bytes
    samtools merge: "/home/shared/8TB_HDD_02/shedurkin/SIFP-nanopore/M-multi-group/data/02-methylation-summaries/bam_pass/G4L2_MinION/barcode31/FBD36396_pass_barcode31_e03e7ede_3381b7e8_0.bam" is truncated
    [W::bam_hdr_read] EOF marker is absent. The input is probably truncated
    8117
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: A: 0.7480469 C: 0.6308594
    [E::bgzf_read_block] Invalid BGZF header at offset 17718219
    [E::bgzf_read] Read block operation failed with error 2 after 0 of 4 bytes
    samtools merge: "/home/shared/8TB_HDD_02/shedurkin/SIFP-nanopore/M-multi-group/data/02-methylation-summaries/bam_pass/G4L2_MinION/barcode32/FBD36396_pass_barcode32_e03e7ede_3381b7e8_0.bam" is truncated
    [W::bam_hdr_read] EOF marker is absent. The input is probably truncated
    17282
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: C: 0.625 A: 0.7519531

``` bash
source .bashvars

cat ${output_dir_top}/G4L2_MinION_barcode30_modkit_summary.txt
cat ${output_dir_top}/G4L2_MinION_barcode31_modkit_summary.txt
cat ${output_dir_top}/G4L2_MinION_barcode32_modkit_summary.txt
```

    # bases             A,C 
    # total_reads_used  100000 
    # count_reads_A     100000 
    # count_reads_C     100000 
    # pass_threshold_C  0.59375 
    # pass_threshold_A  0.7402344 
     base  code   pass_count  pass_frac    all_count  all_frac 
     A     -      7663789     0.9469944    8207808    0.9136961 
     A     a      428961      0.05300559   775275     0.08630389 
     C     -      5260505     0.9280016    5633731    0.8956567 
     C     m      233652      0.041218366  354310     0.05632859 
     C     21839  174481      0.030780056  302015     0.048014674 
    # bases             A,C 
    # total_reads_used  8117 
    # count_reads_A     8117 
    # count_reads_C     8117 
    # pass_threshold_A  0.7480469 
    # pass_threshold_C  0.6308594 
     base  code   pass_count  pass_frac    all_count  all_frac 
     A     -      604043      0.95120835   646373     0.91786575 
     A     a      30984       0.048791625  57840      0.08213424 
     C     -      387650      0.9449186    415593     0.9117501 
     C     m      11550       0.028153772  19071      0.041838974 
     C     21839  11047       0.02692768   21155      0.046410967 
    # bases             C,A 
    # total_reads_used  17282 
    # count_reads_A     17282 
    # count_reads_C     17282 
    # pass_threshold_C  0.625 
    # pass_threshold_A  0.7519531 
     base  code   pass_count  pass_frac    all_count  all_frac 
     C     -      941964      0.93773484   1009393    0.9045421 
     C     m      35732       0.03557157   54957      0.04924833 
     C     21839  26814       0.026693612  51566      0.04620957 
     A     -      1434663     0.95336837   1538782    0.92035663 
     A     a      70173       0.04663166   133159     0.07964336 

# 4 Summary

These are super interesting! In the Group 4 (oldest) sequences, there
are notable levels of 6mA methylation (indicated by base A with code a),
in fact the 6mA methylation level (`pass_frac`) is higher than the 5mC
cytosine methylation (base C, code m). I’m very skeptical of this, since
6mA is very uncommon in eukaryotes, but it hasn’t been evaluated much in
Cnidarians.

It seems most likely that the high 6mA levels are the result of either
contamination (eg. present in contaminant fungi) or an artifact of
age-related degradation. I can’t actually tell right now whether the
modern samples (Group 1) have any 6mA methylation, because that
sequencing was done before the MinKnow software update that allowed for
real-time 6mA modification calling during sequencing. The Group 1 stuff
was only sequenced with 5mC modification calling.

I’m also curious why Group 1 shows notable 5hmC levels, but 5hmC is
*not* present in either of the Group 4 runs. I believe all of the
sequencing runs were performed with softare that is 5hmC-sensitive, so
if its present in one run I’d expect to find it in others…

Finally, in the Group 4 runs, which used the updated basecaller
software, there are also codes for Cytosine with the `21839`
modification code. According to the Nanopore `Dorado` basecaller
[documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/mods/),
this modification code indicates 4mC (N(4)-methylcytosine). 4mC is
believed to be essentially absent outside of prokaryotes (with a few
very rare exceptions, associated with horizontal gene transfer from a
prokaryote). As such, it seems *extremely* likely that the occurance of
4mC in the older samples is spurious.

I’ll have to re-basecall everything with the updated models, including
adapter/barcode trimming and and genome alignment, before proceeding
with methylation summaries.
