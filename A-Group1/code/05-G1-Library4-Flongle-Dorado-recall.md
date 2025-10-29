05-G1-Library4-Flongle-Dorado-recall
================
Kathleen Durkin
2025-10-27

- [1 Create a Bash variables file](#1-create-a-bash-variables-file)
- [2 Raw reads](#2-raw-reads)
  - [2.1 Download raw pod5 files](#21-download-raw-pod5-files)
  - [2.2 Verify raw read checksums](#22-verify-raw-read-checksums)
- [3 Trimming and Demultiplexing](#3-trimming-and-demultiplexing)
  - [3.1 modkit](#31-modkit)

# 1 Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export nanopore_dir=/home/shared/8TB_HDD_02/shedurkin/SIFP-nanopore'
echo 'export output_dir_top=${nanopore_dir}/A-Group1/output/05-G1-Library4-Flongle-Dorado-recall'
echo 'export data_dir_top=${nanopore_dir}/A-Group1/data/05-G1-Library4-Flongle-Dorado-recall'
echo 'export raw_pod5_dir=${data_dir_top}/pod5'
echo 'export raw_pod5_url="https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group1_Flongle/Group1/Library4/20250818_1549_MD-101223_AYL248_54e75371/pod5/"'
echo 'export dorado_models_dir=${nanopore_dir}/dorado-models'
echo 'export genome_dir=${nanopore_dir}/data/GCA_965233905.1_jaEunKnig1.1/'

echo "# Paths to programs"
echo 'export modkit=/home/shared/dist_modkit_v0.5.1_8fa79e3/modkit'
echo 'export samtools=/home/shared/samtools-1.12/samtools'
echo 'export dorado=/home/shared/dorado-1.2.0-linux-x64/bin/dorado'
echo ""

echo "# Set pod5 filename patterns"
echo "export pod5_pattern='*.pod5'"
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
    export output_dir_top=${nanopore_dir}/A-Group1/output/05-G1-Library4-Flongle-Dorado-recall
    export data_dir_top=${nanopore_dir}/A-Group1/data/05-G1-Library4-Flongle-Dorado-recall
    export raw_pod5_dir=${data_dir_top}/pod5
    export raw_pod5_url="https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group1_Flongle/Group1/Library4/20250818_1549_MD-101223_AYL248_54e75371/pod5/"
    export dorado_models_dir=${nanopore_dir}/dorado-models
    export genome_dir=${nanopore_dir}/data/GCA_965233905.1_jaEunKnig1.1/
    # Paths to programs
    export modkit=/home/shared/dist_modkit_v0.5.1_8fa79e3/modkit
    export samtools=/home/shared/samtools-1.12/samtools
    export dorado=/home/shared/dorado-1.2.0-linux-x64/bin/dorado

    # Set pod5 filename patterns
    export pod5_pattern='*.pod5'

    # Set number of CPUs to use
    export threads=20

    # Input/output files
    export raw_checksums=checksums.md5
    export trimmed_checksums=trimmed_fastq_checksums.md5

# 2 Raw reads

## 2.1 Download raw pod5 files

Reads are downloaded from:
<https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group1_Flongle/Group1/Library4/20250818_1549_MD-101223_AYL248_54e75371/pod5/>

Note that this directory contains multiple subdirectories, each
representing one barcode (specimen) and containing the fastq.gz files
associated with that barcode

The `--cut-dirs 7` command cuts the preceding directory structure
(i.e. `nightingales/P_evermanni/30-789513166/`) so that we just end up
with the reads.

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${raw_pod5_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 6 \
--no-host-directories \
--no-parent \
--quiet \
--accept ${pod5_pattern} ${raw_pod5_url}
```

``` bash
# Load bash variables into memory
source .bashvars

ls -lh "${raw_pod5_dir}"
```

    total 4.0K
    drwxr-xr-x 2 shedurkin labmembers 4.0K Oct 27 18:03 pod5

## 2.2 Verify raw read checksums

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${raw_pod5_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 6 \
--no-host-directories \
--no-parent \
--quiet \
--accept 'checksums.md5' ${raw_pod5_url}

cd "${raw_pod5_dir}"

# Recursively verify checksums in all subdirectories
find $(pwd) -type d | while read -r DIR; do
    # Check if checksums.md5 exists in this directory
    if [[ -f "$DIR/checksums.md5" ]]; then
        echo "Verifying checksums in $DIR"
        (
            cd "$DIR" || exit 1
            md5sum -c checksums.md5
        )
        echo ""
    fi
done
```

# 3 Trimming and Demultiplexing

We want to demultiplex and trim adapter and barcode sequences

Notes:

- When selecting multiple modification models, only one modification per
  canonical base may be active at once. For example, sup,4mC,5mC is
  invalid as both modification models operate on the C canonical base
  context.

- Since basecalling is GPU-intensive, but alignment is CPU-intensive,
  including alignment during basecalling shouldn’t have a significant
  impact on throughput

- Dorado uses `minimap2` for alignment by default, which seems to be
  more of a k-mer-based pseudoalignment

- When using alignment during basecalling, no reads are excluded
  following alignment. All reads are basecalled. Then, reads that do not
  align to the provided reference are included in the output SAM/BAM as
  *unmapped* reads, with the appropriate tags. If you want to exclude
  unmapped reads, you’ll have to do that *after* basecalling

``` bash
source .bashvars

# --trim 'all': will trim adapter, primer (if present) and barcode sequences
# --kit-name <kit-name>: allows simultaneous demultiplexing, with which the trimming does not interfere
# --reference: allows basecalling with alignment
# --modified-bases: allows for modified basecalling, space-separated. Note that, for cytosine modifications, there are options for calling modifications on *all* cytosines, and for calling modifications that occur *specifically* in a CG context. Cannot specify modifications that occur on same base (e.g., cannot simultaneously call 4mC and 5mCG)

# --trim: 'all', trim adapter and barcodes
# --kit-name: specifying the library prep i used, Native Barcoding V14 96
# --reference: aligning to Eunicea knighti genome (closest published relative of Eunicea tourneforti)
# -- modified-bases: 5mCG_5hmCG calls 5mC and 5hmC occuring in CG contexts; 6mA calles 6mA modifications.

$dorado basecaller \
hac \
-r ${raw_pod5_dir}/ \
--kit-name SQK-NBD114-96 \
--trim 'all' \
--reference ${genome_dir}/GCA_965233905.1_jaEunKnig1.1_genomic.fna \
--modified-bases 5mCG_5hmCG 6mA \
> ${output_dir_top}/AYL248_recalled.bam
```

Separate out mapped reads, and separate by barcode

``` bash
source .bashvars

cd ${output_dir_top}

# Alignment summary
$samtools flagstat AYL248_recalled.bam > alignment_summary.txt
head alignment_summary.txt

echo ""

# Summarize by barcode
$samtools view AYL248_recalled.bam | \
awk '{for(i=12;i<=NF;i++) if($i ~ /^BC:Z:/) {print substr($i,6)}}' | \
sort | uniq -c | sort -nr
```

    134052 + 0 in total (QC-passed reads + QC-failed reads)
    77531 + 0 secondary
    7462 + 0 supplementary
    0 + 0 duplicates
    120219 + 0 mapped (89.68% : N/A)
    0 + 0 paired in sequencing
    0 + 0 read1
    0 + 0 read2
    0 + 0 properly paired (N/A : N/A)
    0 + 0 with itself and mate mapped

      54749 SQK-NBD114-96_barcode14
      37131 SQK-NBD114-96_barcode13
      34373 SQK-NBD114-96_barcode12
          1 SQK-NBD114-96_barcode46

``` bash
source .bashvars

cd ${output_dir_top}

# Separate out mapped reads
$samtools view -b -F 4 AYL248_recalled.bam > AYL248_recalled_mapped.bam

# Separate by barcode
samtools view -h AYL248_recalled_mapped.bam | awk '/^@/ || /BC:Z:SQK-NBD114-96_barcode12/' | samtools view -b > AYL248_recalled_mapped_barcode12.bam
samtools view -h AYL248_recalled_mapped.bam | awk '/^@/ || /BC:Z:SQK-NBD114-96_barcode13/' | samtools view -b > AYL248_recalled_mapped_barcode13.bam
samtools view -h AYL248_recalled_mapped.bam | awk '/^@/ || /BC:Z:SQK-NBD114-96_barcode14/' | samtools view -b > AYL248_recalled_mapped_barcode14.bam
```

## 3.1 modkit

``` bash
source .bashvars
$samtools view -c ${output_dir_top}/AYL248_recalled.bam
$samtools view -c ${output_dir_top}/AYL248_recalled_mapped.bam

$modkit summary -n 100000 ${output_dir_top}/AYL248_recalled.bam
$modkit summary -n 100000 ${output_dir_top}/AYL248_recalled_mapped.bam
```

    134052
    120219
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: A: 0.8417969 C: 0.7675781
    # bases             C,A 
    # total_reads_used  49044 
    # count_reads_A     49039 
    # count_reads_C     46384 
    # pass_threshold_C  0.7675781 
    # pass_threshold_A  0.8417969 
     base  code  pass_count  pass_frac    all_count  all_frac 
     A     -     6647558     0.99483615   7249998    0.9773442 
     A     a     34505       0.005163824  168062     0.022655789 
     C     -     548284      0.9152803    587667     0.88322896 
     C     h     7433        0.012408311  22748      0.034188908 
     C     m     43317       0.07231142   54947      0.082582116 
    [0;32m>[0m sampling 100000 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: C: 0.78515625 A: 0.8613281
    # bases             A,C 
    # total_reads_used  35226 
    # count_reads_C     34304 
    # count_reads_A     35226 
    # pass_threshold_C  0.78515625 
    # pass_threshold_A  0.8613281 
     base  code  pass_count  pass_frac    all_count  all_frac 
     A     -     5414591     0.99571174   5906424    0.97891504 
     A     a     23319       0.004288228  127219     0.02108494 
     C     -     470935      0.9489109    506862     0.920187 
     C     h     4587        0.00924258   15921      0.028903916 
     C     m     20768       0.041846503  28042      0.05090909 

In both unmapped and mapped reads, we observe negligible 6mA (~0.5%) and
low 5hmCG (~1%). 5mCG is low in both, but differs somewhat from unmapped
(~7%) to mapped (~4%)
