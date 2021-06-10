cd /hpcshare/genomics/yabili_pipeline_vda/yussline

reference="/hpcshare/genomics/references/gatk_bundle/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

merged_leftrim="/hpcshare/genomics/yabili_pipeline_vda/yussline/HKNPC_090_merged_leftalignandtrim.vcf"

output="/hpcshare/genomics/yabili_pipeline_vda/yussline/output.vcf"

gatk SelectVariants \
-R $reference  \
--exclude-filtered \
-V $merged_lefttrim \
--output $output


#VEP Annotation
cd /workspace/somatic

vep -i ~/workspace/somatic/exome.merged.norm.pass_only.vcf --cache --dir /opt/vep_cache/ --format vcf --vcf --plugin Downstream --plugin Wildtype --symbol --terms SO --flag_pick --transcript_version -o ~/workspace/somatic/exome.merged.norm.annotated.vcf

######################################################
#Adding Bam-readcounts to VCF file

cd /workspace/somatic

mkdir -p /workspace/somatic/bam_readcounts

# First activate the bam-readcount conda enironment

source activate bam-readcount

# Running bam-readcount on annotated vcf file, once using tumor bam and another time for normal bam file

python -u /usr/local/bin/bam_readcount_helper.py ~/workspace/somatic/exome.merged.norm.annotated.vcf 'TUMOR' ~/workspace/inputs/references/genome/ref_genome.fa /workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam ~/workspace/somatic/bam_readcounts/

python -u /usr/local/bin/bam_readcount_helper.py ~/workspace/somatic/exome.merged.norm.annotated.vcf 'NORMAL' ~/workspace/inputs/references/genome/ref_genome.fa /workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam ~/workspace/somatic/bam_readcounts/

#pt.2
git clone https://github.com/griffithlab/VAtools

cd ~/workspace/somatic/bam_readcounts/
cat TUMOR_bam_readcount_snv.tsv TUMOR_bam_readcount_indel.tsv > TUMOR_bam_readcount_combined.tsv
cat NORMAL_bam_readcount_snv.tsv NORMAL_bam_readcount_indel.tsv > NORMAL_bam_readcount_combined.tsv

mkdir -p ~/workspace/somatic/final/
cd ~/workspace/somatic/final/
vcf-readcount-annotator ~/workspace/somatic/exome.merged.norm.annotated.vcf ~/workspace/somatic/bam_readcounts/TUMOR_bam_readcount_combined.tsv DNA -s TUMOR -o ~/workspace/somatic/final/exome.tumordna_annotated.vcf.gz
vcf-readcount-annotator ~/workspace/somatic/final/exome.tumordna_annotated.vcf.gz ~/workspace/somatic/bam_readcounts/NORMAL_bam_readcount_combined.tsv DNA -s NORMAL -o ~/workspace/somatic/final/exome.annotated.vcf.gz
# Remove the intermediate file
rm exome.tumordna_annotated.vcf.gz
tabix -p vcf exome.annotated.vcf.gz
source deactivate

#Generating Table from VCF file
wget https://github.com/genome/docker-cle
wget https://github.com/genome/docker-cle/blob/master/add_annotations_to_table_helper.py

# Adjust the output fields accordingly

cd ~/workspace/somatic/final/

gunzip exome.annotated.vcf.gz

#java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar -T VariantsToTable -R ~/workspace/inputs/references/genome/ref_genome.fa --variant ~/workspace/somatic/final/exome.annotated.vcf -F CHROM -F POS -F ID -F REF -F ALT -F set -F AC -F AF -F DP -F AF -F CSQ -o variants.tsv

gatk VariantsToTable -V ~/workspace/somatic/final/exome.annotated.vcf -F CHROM -F POS -F ID -F REF -F ALT -F set -F AC -F AF -GF GT -GF AD -O variants.tsv

wget -c https://raw.githubusercontent.com/genome/docker-cle/master/add_annotations_to_table_helper.py

source activate bam-readcount

python -u add_annotations_to_table_helper.py variants.tsv exome.annotated.vcf Consequence,Gene ./

source deactivate


#Machine learning filtration
https://github.com/griffithlab/DeepSVR/
