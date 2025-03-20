# lrsv_oysters

This file describes the subprocess of long-read structural variation (SV) analysis in oysters, with additional SV analysis included. For the complete SV analysis workflow, please visit the aforementioned website: [Github Repository](https://github.com/OZTaekOppa/lrsv_oysters?tab=readme-ov-file)

## Citation
Hyungtaek Jung et al. 2024: Unraveling the biotechnological potential of Crassostrea gigas: comparative genomics and structural variations, BioRxiv.

## Contents:
- [LICENSE](#license)
- [GETTING STARTED](#getting-started)
- [FAQ](#faq)
- [WIKI PAGE](#wiki-page)
- [AUTHORS](#authors)
- [COPYRIGHT](#copyright)

## GETTING STARTED
The files used in this study is also available in main brance: [Github Repository](https://github.com/OZTaekOppa/lrsv_oysters?tab=readme-ov-file)

### Tested Files
- C. gigas1: Genome assembly [GCA_011032805.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_011032805.1/), Project number [PRJNA598006](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA598006) and Raw data [SRA_598006](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=598006) in NCBI
- C. gigas2: Genome assembly [GCF_902806645.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/), Project number [PRJEB35351](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB35351) and Raw data [SRA_596972](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=596972) in NCBI

### Additional SV caller (Sniffles2)
We added another SV caller, [Sniffles2](https://www.nature.com/articles/s41587-023-02024-y), and compared its results with those of CuteSV. This process was applied to all five alignment tools used in this study: LRA, Minimap2, NGMLR, Vulcan, and WinnowMap.

An example code for LRA is provided below:
```
module load lra
module load samtools
module load cutesv
module load sniffles2

# LRA Indexing and Mapping for PacBio Data
lra index -CLR /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
lra align -CLR /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta /Cra_gig2/Cra_gig2_PBallRN.fastq -t 8 -p s > LRA_Cg1RCg2PB.sam

# Samtools: Conversion of SAM to BAM, Including Mapping Statistics
samtools faidx /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta
samtools view -bT /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta LRA_Cg1RCg2PB.sam > LRA_Cg1RCg2PB.bam
samtools sort -@8 -T $TMPDIR/LRA_Cg1RCg2PB.sorted -O bam -o LRA_Cg1RCg2PB.sorted.bam LRA_Cg1RCg2PB.bam
samtools index LRA_Cg1RCg2PB.sorted.bam
samtools flagstat LRA_Cg1RCg2PB.sorted.bam &> LRA_Cg1RCg2PB_stats.txt

# CuteSV: Aligning Cg2 PacBio Data Against the Cg1 Reference Genome
# BAM File Generation via LRA and Samtools
# Utilise --genotype with CuteSV to Generate a Genotype Score, a Prerequisite for SURVIVOR Analysis
cuteSV --genotype LRA_Cg1RCg2PB.sorted.bam /Cgig1/GCA_011032805.1_ASM1103280v1_genomic.fasta Cg1RCg2PB_LRACuteSV.vcf /Cgig1/LR_MappingSV/LRA --threads 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5

# Sniffles2: Aligning Cg2 PacBio Data Against the Cg1 Reference Genome
# BAM File Generation via LRA and Samtools
sniffles -i LRA_Cg1RCg2PB.sorted.bam --threads 8 -v Cg1RCg2PB_LRAsniffles.vcf
```
- For specific parameter details, refer to the official websites of [Sniffles2](https://github.com/fritzsedlazeck/Sniffles).

### Cross-validation 
To exclude low-quality SVs, we applied filtering criteria which remove variants with insufficient read depth or low confidence, and retaining only those flagged as ‘PASS’. Finally, cross-validation between the SV callers was performed using SURVIVOR, allowing us to select SVs consistently detected by both callers.

```
module load vcftools
module load SURVIVOR

# Obtain quality-controlled VCF files through vcftools
vcftools --vcf Cg1RCg2PB_LRACuteSV.vcf --recode --recode-INFO-all --out PASS_Cg1RCg2PB_LRACuteSV --remove-filtered-all
vcftools --vcf Cg1RCg2PB_LRAsniffles.vcf --recode --recode-INFO-all --out PASS_Cg1RCg2PB_LRASnF2 --remove-filtered-all

# Merge SVs Called by Both SV Calling Tools
SURVIVOR merge LRA_vcfs 400 2 2 1 0 30 LRASV_Cg1Rmerged_Def.vcf
SURVIVOR merge MM2_vcfs 400 2 2 1 0 30 MM2SV_Cg1Rmerged_Def.vcf
SURVIVOR merge NgM_vcfs 400 2 2 1 0 30 NgMSV_Cg1Rmerged_Def.vcf
SURVIVOR merge VLC_vcfs 400 2 2 1 0 30 VLCSV_Cg1Rmerged_Def.vcf
SURVIVOR merge WmM_vcfs 400 2 2 1 0 30 WmMSV_Cg1Rmerged_Def.vcf

# Checking for Overlapping VCF Files with a PERL Script
# Ensure to run the overlap.txt check each time to prevent overlapp.txt issues
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' sample_merged.vcf | sed -e 's/\(.\)/\1 /g' > sample_merged_overlap.txt 

# Plotting
SURVIVOR genComp CuteSV_Cg1Rmerged_Def.vcf CuteSV_Cg1Rmerged_Def.mat.txt
```

## FAQ
We encourage users to use the [Issues](../../issues).

## AUTHORS
Hyungtaek Jung and et al. 2024 Unraveling the biotechnological potential of Crassostrea gigas: comparative genomics & structural variations.

## COPYRIGHT
The full scripts are distributed under the MIT license.
