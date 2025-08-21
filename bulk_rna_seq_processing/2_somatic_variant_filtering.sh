#!/bin/bash

source /home/m7huang/miniconda3/etc/profile.d/conda.sh

#### UPSTREAM from https://github.com/claupomm/RNA-seq_snv_tumour_only?tab=readme-ov-file ####
### RNA editing:
## cd /expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/tools
## wget http://srv00.recas.ba.infn.it/webshare/ATLAS/download/TABLE1_hg38_v3.txt.gz

### Low complexity regions:
#wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
#unzip snpEff_latest_core.zip
#cd snpEff/
#curl -o LCR-hs38_with_chr.bed.gz -L https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true
#gunzip -c LCR-hs38_with_chr.bed.gz > LCR-hs38_with_chr.bed
#sed 's/^chr//g' LCR-hs38_with_chr.bed > LCR-hs38.bed

### 1000 Genomes
## wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase3_v4_20130502.sites.hg38.vcf
##bcftools filter -i 'INFO/AF>0.1' -o /expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/tools/1000G_phase3_v4_20130502.sites.hg38.af01.vcf /expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/tools/1000G_phase3_v4_20130502.sites.hg38.vcf

### dbSNP ###
## wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
## zgrep "^#" GCF_000001405.40.gz > GCF_000001405.40.common.vcf
## zgrep ";COMMON" GCF_000001405.40.gz >> GCF_000001405.40.common.vcf
## refseq=$(cut -f1 gcf_refseq2chr.tsv)
## IFS=' ' read -ra refseq <<< $(echo $refseq)
## chr=$(cut -f2 gcf_refseq2chr.tsv)
## IFS=$' ' read -ra chr <<< $(echo $chr)
## touch GCF_000001405.40.common.chr.tsv
##
##for i in ${!refseq[@]}; do
##date
##echo ${chr[$i]}
##zgrep ${refseq[$i]} GCF_000001405.40.common.vcf.gz | sed "s/${refseq[$i]}/${chr[$i]}/g" | cut -f1,2 >> GCF_000001405.40.common.chr.tsv
##done

#### End upstream ####

#### 1: Filtering ####

#path=/expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/
#
#cd $path
#
#mkdir filtered_vcfs
#mkdir samples
#mkdir annotated_samples
#mkdir maf_samples
#
## directory with samples
#ls /expanse/lustre/projects/csd728/kfisch/placental_variants/new_annotations/filtered_vcfs > vcf_names.txt
#
## path to tool folder
#tools/
#rna_edit=$path/tools/TABLE1_hg38_v3.txt.gz
#SNPEFF=$path/tools/snpEff
#filter=$path/tools/1000G_phase3_v4_20130502.sites.hg38.af01.vcf
#filter2=$path/tools/GCF_000001405.40.common.chr.tsv

#while IFS= read -r sample;
#do
#
#    echo "$sample"
#    conda activate vcftools
#    gzip -c samples/$sample.notGVCF.vcf > samples/$sample.notGVCF.vcf.gz
#
#    vcftools --gzvcf samples/$sample.notGVCF.vcf.gz \
#        --recode --exclude-positions $rna_edit --stdout | gzip -c > samples/$sample.pass2.vcf.gz
#
#    java -jar $SNPEFF/SnpSift.jar intervals \
#        -noLog -x -i samples/$sample.pass2.vcf.gz \
#        $SNPEFF/LCR-hs38.bed | gzip > samples/$sample.pass2.lcr.vcf.gz
#
#    conda activate vcftools
#    vcftools --gzvcf samples/$sample.pass2.lcr.vcf.gz \
#     --exclude-positions $filter --recode --recode-INFO-all --stdout | pigz -p4 > samples/$sample.hc.pass2.lcr.1kG.vcf.gz
#
#    vcftools --gzvcf samples/$sample.hc.pass2.lcr.1kG.vcf.gz \
#        --exclude-positions $filter2 --recode --recode-INFO-all --stdout | pigz -p4 > filtered_samples/$sample.hc.pass2.lcr.1kG.dbsnp.vcf.gz
#
#done < vcf_names.txt


#### 2: Annotation ####

orig_path=/expanse/lustre/projects/csd728/kfisch/somatic_RNAseq_variants/circulating_batch2_redo
path=/expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/all_samples/circulating_batch2
sample=$1

cd $path

## change by file name
#basename -a $orig_path/*.notGVCF.vcf.hc.pass2.lcr.1kG.dbsnp.vcf.gz > file_names.txt
#sed 's/\.notGVCF.vcf.hc.pass2.lcr.1kG.dbsnp.vcf.gz//g' file_names.txt > sample_names.txt

/home/m7huang/miniconda3/envs/vep/bin/perl /home/m7huang/software/ensembl-vep/vep \
    --input_file $orig_path/$sample.hc.pass2.lcr.1kG.dbsnp.vcf.gz \
    --output_file $path/annotated_samples/$sample.vep.vcf \
    --offline --dir /expanse/lustre/projects/csd728/kfisch/placental_variants/ \
    --fasta /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
    --plugin LoF,loftee_path:/home/m7huang/software/loftee,human_ancestor_fa:/home/m7huang/software/human_ancestor.fa.gz --dir_plugins /home/m7huang/software/loftee --force_overwrite --species homo_sapiens --assembly GRCh38 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --format vcf  --pubmed --fork 4 --cache_version 107 --polyphen b --everything --af_gnomadg --af_gnomade
    
perl /home/m7huang/software/vcf2maf-1.6.21/vcf2maf.pl \
    --input-vcf $path/annotated_samples/$sample.vep.vcf \
    --output-maf $path/maf_samples/$sample.maf \
    --ref-fasta /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
     --tumor-id $sample --inhibit-vep --ncbi-build GRCh38 --retain-info AC,AF,AN,BaseQRankSum,ClippingRankSum,DP,ExcessHet,FS,MLEAC,MLEAF,MQ,MQ0,MQRankSum,QD,ReadPosRankSum,SOR --retain-fmt GT,AD,DP,GQ,PL --retain-ann QUAL,FILTER,MAX_AF,MAX_AF_POPS,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_OTH_AF,gnomADg_SAS_AF,LoF,LoF_filter,LoF_flags,LoF_info,MANE_SELECT,MANE_PLUS_CLINICAL,APPRIS,UNIPROT_ISOFORM,miRNA

exit 0



while IFS= read -r sample;
do

    echo "$sample"
    bash 20250531_somatic_filtering.sh $sample

done < /expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/all_samples/circulating_batch1/sample_names.txt
