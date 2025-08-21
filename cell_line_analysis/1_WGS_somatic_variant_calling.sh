#!/bin/bash

source /home/m7huang/miniconda3/etc/profile.d/conda.sh

#################### UPSTREAM ####################
# wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
# wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

#path=/expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/
#conda activate gatk
#gatk SelectVariants \
#        -V $path/tools/af-only-gnomad.hg38.vcf.gz \
#        -select-type SNP -restrict-alleles-to BIALLELIC \
#        -select "AF > 0.05" \
#        -O $path/tools/af-only-gnomad-common-biallelic.grch38.main.vcf.gz \
#        --lenient

# fasta_formatter -i /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.fa -w 5000 > /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.5000chr.fa

#################### DNA-Seq Alignment Command Line Parameters ####################
### from: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/

## ensure that the average read length > 70 bps
# samtools stats $input/${sample}_R1.fastq.gz.bam  | grep "average length"

input="/expanse/lustre/projects/csd728/avwilliams/2024_25_wgs_concat/250611_celllines"
output="/expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/samples/"
sample=$1

echo $sample

cd /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/samples/
mkdir $sample
mkdir $sample/tmp

# run with 8 threads!
conda activate gatk
bwa mem \
    -t 8 \
    -T 0 \
    /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.fa \
    $input/${sample}_R1.fastq.gz $input/${sample}_R2.fastq.gz > $output/$sample/$sample.sam

conda activate samtools
samtools view $output/$sample/$sample.sam \
    -Shb \
    -o $output/$sample/$sample.bam

 rm $output/$sample/$sample.sam

##################### DNA-Seq Co-Cleaning Command Line Parameters ####################

conda activate gatk
gatk SortSam \
    CREATE_INDEX=true \
    INPUT=$output/$sample/$sample.bam \
    OUTPUT=$output/$sample/$sample.sort.bam \
    SORT_ORDER=coordinate \
    TMP_DIR=$output/$sample/tmp \
    VALIDATION_STRINGENCY=STRICT

gatk MarkDuplicates \
    CREATE_INDEX=true \
    INPUT=$output/$sample/$sample.sort.bam \
    O=$output/$sample/$sample.sort.marked_duplicates.bam \
    M=$output/$sample/$sample.marked_dup_metrics.txt \
    TMP_DIR=$output/$sample/tmp \
    VALIDATION_STRINGENCY=STRICT

#### add random read group ####
gatk --java-options "-Xmx8G" AddOrReplaceReadGroups \
   I=$output/$sample/$sample.sort.marked_duplicates.bam \
   O=$output/$sample/$sample.sort.marked_duplicates.rg.bam \
   TMP_DIR=$output/$sample/tmp \
   RGID=1 \
   RGLB=1 \
   RGPL=1 \
   RGPU=1 \
   RGSM=1

gatk --java-options "-Xmx8G" BaseRecalibrator \
   --input $output/$sample/$sample.sort.marked_duplicates.rg.bam \
   -R /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.fa \
   --tmp-dir $output/$sample/tmp \
   --known-sites /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/1000G_phase3_v4_20130502.sites.hg38.vcf \
   -O $output/$sample/$sample.recal_data.table

gatk --java-options "-Xmx8G" ApplyBQSR \
   -R /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.fa \
   -I $output/$sample/$sample.sort.marked_duplicates.rg.bam \
   --tmp-dir $output/$sample/tmp \
   --bqsr-recal-file $output/$sample/$sample.recal_data.table \
   -O $output/$sample/$sample.sort.marked_duplicates.bqsr.bam

################### Tumor-Only Variant Call Command-Line Parameters ####################
## updated workflow from: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
conda activate gatk

## run with 4 threads!
gatk --java-options "-Xmx8G" Mutect2 \
    -R /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.fa \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 \
    -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM \
    -I $output/$sample/$sample.sort.marked_duplicates.bqsr.bam \
    -O $output/$sample/$sample.mutect2.vcf \
    --tmp-dir $output/$sample/tmp \
    --f1r2-tar-gz $output/$sample/$sample.f1r2.tar.gz
    --germline-resource /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/af-only-gnomad.hg38.vcf.gz \
    -pon /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/1000g_pon.hg38.vcf.gz

gatk --java-options "-Xmx8G" LearnReadOrientationModel \
    -I $output/$sample/$sample.f1r2.tar.gz \
    -O $output/$sample/$sample.read.orientation.model.tar.gz \
    --tmp-dir $output/$sample/tmp

gatk --java-options "-Xmx8G" GetPileupSummaries \
    -I $output/$sample/$sample.sort.marked_duplicates.bqsr.bam \
    -O $output/$sample/$sample.targeted_sequencing.table \
    -R /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.fa \
    --tmp-dir $output/$sample/ \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 \
    -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM \
    -V /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/af-only-gnomad-common-biallelic.grch38.main.vcf.gz

gatk --java-options "-Xmx8G" CalculateContamination \
    --tmp-dir $output/$sample/ \
    -I $output/$sample/$sample.targeted_sequencing.table \
    -O $output/$sample/$sample.targeted_sequencing.contamination.table

gatk --java-options "-Xmx8G" FilterMutectCalls \
    -O $output/$sample/$sample.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz \
    -R /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/tools/GRCh38.d1.vd1.fa \
    -V $output/$sample/$sample.mutect2.vcf \
    --contamination-table $output/$sample/$sample.targeted_sequencing.contamination.table \
    --tmp-dir $output/$sample/tmp \
    --ob-priors $output/$sample/$sample.read.orientation.model.tar.gz \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 \
    -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

echo $sample > $output/$sample/$sample.sample_name.txt

conda activate bcftools
bcftools reheader -s $output/$sample/$sample.sample_name.txt $output/$sample/$sample.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz > $output/$sample/$sample.vcf.gz

#################### VEP Annotation ####################
conda activate vep

/home/m7huang/miniconda3/envs/vep/bin/perl /home/m7huang/software/ensembl-vep/vep \
    --input_file $output/$sample/$sample.vcf \
    --output_file $output/$sample/$sample.vep.2.vcf \
    --offline --dir /expanse/lustre/projects/csd728/kfisch/placental_variants/ \
    --fasta /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
    --force_overwrite --species homo_sapiens --assembly GRCh38 --no_progress --no_stats --buffer_size 50 --symbol \
    --numbers --show_ref_allele --polyphen b --variant_class --allele_number --failed 1 --vcf --format vcf --fork 4 --cache_version 107

perl /home/m7huang/software/vcf2maf-1.6.21/vcf2maf.pl \
    --input-vcf $output/$sample/$sample.vep.2.vcf \
    --output-maf $output/maf_min_vep/$sample.maf \
    --ref-fasta /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
     --tumor-id $sample --inhibit-vep --ncbi-build GRCh38 --retain-info DP,ECNT,PON,POPAF,TLOD --retain-fmt AF,AD,DP,GQ,PL --retain-ann QUAL,FILTER
     
exit 0

find $input -type f -name '*_R1.fastq.gz' -exec basename {} _R1.fastq.gz \;

while IFS= read -r sample;
do
    echo "$sample"
    sbatch 20250602_somatic_WGS_cell_line.sh $sample

done < /expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/vcf_names2.txt
