# DNAseq_Analysis
# NGS DNAseq Pipeline for Cancer Mutation Detection 

This pipeline identifies specific, known cancer-driving mutations in a targeted set of 15 genes (including TP53) from medulloblastoma tumor samples.It uses a targeted sequencing approach with the Illumina TruSight Tumor 15 assay, enabling fast, cost-effective, and sensitive detection of variants, especially in the TP53 variant subgroup 4 of medulloblastoma.

# Overview
Objective: Detect known mutations in exons of 15 cancer-related genes, focusing on TP53 mutations in medulloblastoma.

Assay: Illumina TruSight Tumor 15 kit—targets exons only, improving detection sensitivity for low-frequency variants.

Cancer Type: Medulloblastoma; TP53 variant subgroup 4 samples highlight poor-outcome driver mutations.

Reference Guide: GenCore BioNYU Tutorial

Guide Used: https://learn.gencore.bio.nyu.edu/

Demo Data: SRA ERX2405312

# Prerequisites
Recommended: Mac OS (Homebrew used for installations)
Tools: curl, wget, gunzip, fastqc, fastp, bwa, samtools, picard, gatk, snpeff, veep, igv
Data: SRA access for downloading FASTQ files

# Pipeline Workflow
1. Setup Project Directory
```
mkdir -p Project_data/Rawdata
cd Project_data/Rawdata
```
2. Download Raw Data
```
curl -o ERR2356709_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/009/ERR2356709/ERR2356709_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/009/ERR2356709/ERR2356709_2.fastq.gz
```
3. Unzip the Fastq File
```
gunzip *.fastq.gz
```

4. Quality Check & Trimming
```
#Install tools
brew install fastqc fastp
fastqc ERR2356709_1.fastq ERR2356709_2.fastq
```
```
#Adapter trimming (with custom file)
cat > adapter_seq.fasta << 'EOF'
>Illumina_Universal_Adapter
AGATCGGAAGAG
>Illumina_Primer
GATCGGAAGAGC
>Illumina_Small_RNA_Adapter
TGGAATTCTCGG
>Nextera_Adapter
CTGTCTCTTATA
EOF
fastp -i ERR2356709_1.fastq -o Trim_ERR2356709_1.fastq -I ERR2356709_2.fastq -O Trim_ERR2356709_2.fastq --adapter_fasta adapter_seq.fasta --thread 2 --html fastp_report.html -q 30 -l 50
```
4. Reference Genome Preparation
```
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
gunzip chr17.fa.gz
```
5. Alignment (BWA)
```
brew install bwa
```
6. Load Chromosome fasta file
```
bwa index chr17.fa
bwa mem -t 2 chr17.fa Trim_ERR2356709_1.fastq Trim_ERR2356709_2.fastq > ERR2356709_aligned.sam
```
7. Samtools Processing
```
brew install samtools
```
```
samtools view -bS ERR2356709_aligned.sam > ERR2356709_aligned.bam
samtools sort ERR2356709_aligned.bam -o ERR2356709_sorted.bam
samtools index ERR2356709_sorted.bam
samtools rmdup -sS ERR2356709_sorted.bam rmdup_ERR2356709_sorted.bam
```
8. Picard Read Group & Dictionary
```
brew install picard-tools
```
```
#load Chromosome dictionary
picard CreateSequenceDictionary -R chr17.fa -O chr17.dict
```
```
picard AddOrReplaceReadGroups \
    I=rmdup_ERR2356709_sorted.bam \
    O=picard_output.bam \
    RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=ERR2356709 \
    SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
```
9. GATK Preprocessing & Variant Calling
```
wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
unzip gatk-4.3.0.0.zip
```
```
samtools faidx chr17.fa
gatk CreateSequenceDictionary -R chr17.fa -O chr17.dict
```
```
gatk HaplotypeCaller \
  -R chr17.fa \
  -I picard_output.bam \
  -O GATK_output.vcf
```
10. Variant Annotation & Filtering
SnpEff / SnpSift
```
curl -OL https://github.com/pcingola/SnpEff/releases/latest/download/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
java -jar snpEff.jar -version
```
```
# Or
brew install snpeff
snpEff -version
SnpSift -version
```
```
SnpSift filter "((QUAL >= 30) & (DP >= 10) & (MQ >= 30))" GATK_output.vcf > filtered_variants.vcf
```
11. VEP (Variant Effect Predictor)
```
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```
```
# Or method 2
wget https://github.com/Ensembl/ensembl-vep/archive/refs/heads/main.zip
unzip main.zip
cd ensembl-vep-main
perl INSTALL.pl
```
10. Visualization (IGV)
```
brew install igv
```
# Use IGV tools to index chr17 and visualize variants
Output
Key output files: Filtered VCF file containing detected cancer-driving mutations in TP53 and other targeted genes.

Reports: FastQC and Fastp HTML reports (quality/trimming), annotated variant files.

Final step: Integrative visualization and variant impact interpretation for medulloblastoma mutation research.

References & Guides
Data/Workflow example: NCBI SRA ERX2405312

Protocol Guide: GenCore NYU NGS Tutorials
