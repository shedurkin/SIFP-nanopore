02-G4-Library1-Flongle-Dorado-recall
================
Kathleen Durkin
2025-10-24

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
echo 'export output_dir_top=${nanopore_dir}/D-Group4/output/02-G4-Library1-Flongle-Dorado-recall'
echo 'export data_dir_top=${nanopore_dir}/D-Group4/data/02-G4-Library1-Flongle-Dorado-recall'
echo 'export raw_pod5_dir=${data_dir_top}/pod5_pass'
echo 'export raw_pod5_url="https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group4_Flongle/Library1/20250904_2022_MN47571_AXT657_1336ce16/pod5_pass/"'
echo 'export dorado_models_dir=${nanopore_dir}/dorado-models'
echo 'export genome_dir=${nanopore_dir}/data/GCA_965233905.1_jaEunKnig1.1/'

echo 'export trimmed_fastqc_dir=${output_dir_top}/trimmed-fastqc'
echo 'export trimmed_reads_dir=${nanopore_dir}/D-Group4/output/02-G4-Library1-Flongle-Dorado-recall/trimmed-reads'
echo 'export trimmed_reads_url=""'
echo ""

echo "# Paths to programs"
echo 'export modkit=/home/shared/dist_modkit_v0.5.1_8fa79e3/modkit'
echo 'export samtools=/home/shared/samtools-1.12/samtools'
echo 'export dorado=/home/shared/dorado-1.2.0-linux-x64/bin/dorado'
echo ""

echo "# Set FastQ filename patterns"
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
    export output_dir_top=${nanopore_dir}/D-Group4/output/02-G4-Library1-Flongle-Dorado-recall
    export data_dir_top=${nanopore_dir}/D-Group4/data/02-G4-Library1-Flongle-Dorado-recall
    export raw_pod5_dir=${data_dir_top}/pod5_pass
    export raw_pod5_url="https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group4_Flongle/Library1/20250904_2022_MN47571_AXT657_1336ce16/pod5_pass/"
    export dorado_models_dir=${nanopore_dir}/dorado-models
    export genome_dir=${nanopore_dir}/data/GCA_965233905.1_jaEunKnig1.1/
    export trimmed_fastqc_dir=${output_dir_top}/trimmed-fastqc
    export trimmed_reads_dir=${nanopore_dir}/D-Group4/output/02-G4-Library1-Flongle-Dorado-recall/trimmed-reads
    export trimmed_reads_url=""

    # Paths to programs
    export modkit=/home/shared/dist_modkit_v0.5.1_8fa79e3/modkit
    export samtools=/home/shared/samtools-1.12/samtools
    export dorado=/home/shared/dorado-1.2.0-linux-x64/bin/dorado

    # Set FastQ filename patterns
    export pod5_pattern='*.pod5'

    # Set number of CPUs to use
    export threads=20

    # Input/output files
    export raw_checksums=checksums.md5
    export trimmed_checksums=trimmed_fastq_checksums.md5

# 2 Raw reads

## 2.1 Download raw pod5 files

Reads are downloaded from:
<https://gannet.fish.washington.edu/kdurkin1/SIFP_2025/Group4_Flongle/Library1/20250904_2022_MN47571_AXT657_1336ce16/pod5_pass/>

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

    total 12K
    drwxr-xr-x 2 shedurkin labmembers 4.0K Oct 24 17:37 barcode27
    drwxr-xr-x 2 shedurkin labmembers 4.0K Oct 24 17:37 barcode28
    drwxr-xr-x 2 shedurkin labmembers 4.0K Oct 24 17:37 barcode29

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
> ${output_dir_top}/AXT657_pass_recalled.bam
```

Separate out mapped reads, and separate by barcode

``` bash
source .bashvars

cd ${output_dir_top}

# Alignment summary
$samtools flagstat AXT657_pass_recalled.bam > alignment_summary.txt
head alignment_summary.txt

echo ""

# Summarize by barcode
$samtools view AXT657_pass_recalled.bam | \
awk '{for(i=12;i<=NF;i++) if($i ~ /^BC:Z:/) {print substr($i,6)}}' | \
sort | uniq -c | sort -nr
```

    [E::bgzf_read_block] Invalid BGZF header at offset 12177310
    [E::bgzf_read] Read block operation failed with error 2 after 0 of 4 bytes
    samtools flagstat: error reading from "AXT657_pass_recalled.bam"

    [E::bgzf_read_block] Invalid BGZF header at offset 12177310
    [E::bgzf_read] Read block operation failed with error 2 after 0 of 4 bytes
    samtools view: error reading file "AXT657_pass_recalled.bam"
    samtools view: error closing "AXT657_pass_recalled.bam": -1
      27713 SQK-NBD114-96_barcode29
      13028 SQK-NBD114-96_barcode27
      11196 SQK-NBD114-96_barcode28

``` bash
source .bashvars

cd ${output_dir_top}

# Separate out mapped reads
$samtools view -b -F 4 AXT657_pass_recalled.bam > AXT657_pass_recalled_mapped.bam

# Separate by barcode
samtools view -h AXT657_pass_recalled_mapped.bam | awk '/^@/ || /BC:Z:SQK-NBD114-96_barcode27/' | samtools view -b > AXT657_recalled_mapped_barcode27.bam
samtools view -h AXT657_pass_recalled_mapped.bam | awk '/^@/ || /BC:Z:SQK-NBD114-96_barcode28/' | samtools view -b > AXT657_recalled_mapped_barcode28.bam
samtools view -h AXT657_pass_recalled_mapped.bam | awk '/^@/ || /BC:Z:SQK-NBD114-96_barcode29/' | samtools view -b > AXT657_recalled_mapped_barcode29.bam
```

## 3.1 modkit

``` bash
source .bashvars
$samtools view -c ${output_dir_top}/AXT657_pass_recalled.bam
$samtools view -c ${output_dir_top}/AXT657_pass_recalled_mapped.bam

$modkit summary ${output_dir_top}/AXT657_pass_recalled.bam
$modkit summary ${output_dir_top}/AXT657_pass_recalled_mapped.bam
```

    [E::bgzf_read_block] Invalid BGZF header at offset 12177310
    [E::bgzf_read] Read block operation failed with error 2 after 0 of 4 bytes
    samtools view: error reading file "/home/shared/8TB_HDD_02/shedurkin/SIFP-nanopore/D-Group4/output/02-G4-Library1-Flongle-Dorado-recall/AXT657_pass_recalled.bam"
    samtools view: error closing "/home/shared/8TB_HDD_02/shedurkin/SIFP-nanopore/D-Group4/output/02-G4-Library1-Flongle-Dorado-recall/AXT657_pass_recalled.bam": -1
    37239
    [0;32m>[0m sampling 10042 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: A: 0.6738281 C: 0.6816406
    # bases             C,A 
    # total_reads_used  10042 
    # count_reads_C     8707 
    # count_reads_A     10035 
    # pass_threshold_A  0.6738281 
    # pass_threshold_C  0.6816406 
     base  code  pass_count  pass_frac    all_count  all_frac 
     C     -     37492       0.87994933   39715      0.83909065 
     C     h     1297        0.030441007  2715       0.057361983 
     C     m     3818        0.08960969   4901       0.10354736 
     A     -     543577      0.97744554   587318     0.95050347 
     A     a     12543       0.022554485  30584      0.04949652 
    [0;32m>[0m sampling 10042 reads from BAM
    [0;32m>[0m calculating threshold at 10(th) percentile
    [0;32m>[0m calculated thresholds: A: 0.7402344 C: 0.59765625
    # bases             C,A 
    # total_reads_used  10042 
    # count_reads_A     10042 
    # count_reads_C     9384 
    # pass_threshold_A  0.7402344 
    # pass_threshold_C  0.59765625 
     base  code  pass_count  pass_frac    all_count  all_frac 
     A     -     731233      0.984624     791874     0.961347 
     A     a     11419       0.015375976  31839      0.038653027 
     C     -     78316       0.9205849    83828      0.8884037 
     C     h     3180        0.037380103  5993       0.06351343 
     C     m     3576        0.042034984  4537       0.048082832 

Interesting…

So in both unmapped and mapped reads, there are bases called with 5mC,
5mhC, and 6mA modifications. The only modification which seems to differ
from unmapped to mapped reads is 5mC, which is present at ~9% in the
unmapped reads, but at ~4% in the mapped reads.

I find it *very* interesting that there is notable 6mA presence in the
mapped reads…

I can also check to see if 4mC modifications are present, but I’ll have
to re-basecall. Since 4mC modifications occur on the same base as 5mC
and 5hmC, it cannot be called simultaneously with them.
