#ILLUMINA_PREPROSS

echo "****************************STARTING*******************************"
#example files 
#strelka
#/hpcshare/genomics/yabili_pipeline_vda/yussline/HKNPC-090_PAIRED_ILLUMINA/results/VariantCall    ing/T_vs_N/Strelka/StrelkaBP_T_vs_N_somatic_snvs.vcf.gz
#/hpcshare/genomics/yabili_pipeline_vda/yussline/HKNPC-090_PAIRED_ILLUMINA/results/VariantCall    ing/T_vs_N/Strelka/StrelkaBP_T_vs_N_somatic_indels.vcf.gz
strelka_snv="/home/yabili/yussab/mergeVCF_MutectStrelka/test/StrelkaBP_T_vs_N_somatic_snvs.vcf.gz"
strelka_indels="/home/yabili/yussab/mergeVCF_MutectStrelka/test/StrelkaBP_T_vs_N_somatic_indels.vcf.gz"
#strelka_snv=""
#strelka_indels=""

#mutect
#/hpcshare/genomics/yabili_pipeline_vda/yussline/HKNPC-090_PAIRED_GATK/results/VariantCalling/Mutect2/filt_vcf/T_vs_N_filtered.vcf.gz
mutect="/home/yabili/yussab/mergeVCF_MutectStrelka/test/T_vs_N_filtered.vcf.gz"
#mutect=""

results="/home/yabili/yussab/mergeVCF_MutectStrelka/results"
##############################ILLUMINA_PREPROCESSING################################################
#In this step we add the GT field in ILLUMINA files

zcat ${strelka_snv} | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > ${results}/somatic.snvs.gt.vcf
zcat ${strelka_indels} | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - >  ${results}/somatic.indels.gt.vcf

bgzip  ${results}/somatic.snvs.gt.vcf
tabix  ${results}/somatic.snvs.gt.vcf.gz
bgzip  ${results}/somatic.indels.gt.vcf
tabix  ${results}/somatic.indels.gt.vcf.gz

#In this step we merge Illumina snvs and indels in one file
bcftools concat -a -o  ${results}/StrelkaBP_T_vs_N_somatic_snvs_indels.vcf.gz  -O z  ${results}/somatic.snvs.gt.vcf.gz  ${results}/somatic.indels.gt.vcf.gz

tabix  ${results}/StrelkaBP_T_vs_N_somatic_snvs_indels.vcf.gz


################################GATK_PREPROCESSING###############################################

#cd /hpcshare/genomics/yabili_pipeline_vda/yussline/HKNPC-090_PAIRED_GATK/results/VariantCalling/Mutect2/filt_vcf
echo ${mutect}  >  ${results}/wgs_vcf.fof

#WARNING - you need the .tbi inde of ${mutect}
bcftools concat --allow-overlaps --remove-duplicates --file-list  ${results}/wgs_vcf.fof --output-type z --output  ${results}/mutect_wgs.vcf.gz

#mv mutect_wgs.vcf.gz wgs.vcf.gz

tabix  ${results}/mutect_wgs.vcf.gz
gunzip  ${results}/mutect_wgs.vcf.gz

#WARNING - this step need to be coded to retrive sample name from ${mutect} and convert it to NORMAL/TUMOR
#sed -i 's/N/NORMAL/' mutect_wgs.vcf #create a variable for Normal sample
#sed -i 's/T/TUMOR/' mutect_wgs.vcf #create a variable for Tumor sample

#############################################################################################################################################
#COMBINE VCF
strelka="${results}/StrelkaBP_T_vs_N_somatic_snvs_indels.vcf.gz"
gatk="${results}/mutect_wgs.vcf"
reference="/home/yabili/storage/references/gatk_bundle/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
out="${results}/HKNPC_090_merged.vcf.gz"

#eseguire questo comando usando singularity shell -B /hpcshare  gatk3*  + GenomeAnalysisTK.jar 

java -Xmx24g -jar /home/yabili/miniconda3/envs/bioinfo/opt/gatk-3.8/GenomeAnalysisTK.jar -T CombineVariants \
	-R ${reference} \
	-genotypeMergeOptions PRIORITIZE \
	--rod_priority_list mutect,strelka \
	--variant:strelka ${strelka} \
	--variant:mutect ${gatk}  \
	-o ${out}

#Left Align and Trim
#Reference for explaining left align and trim: https://genome.sph.umich.edu/wiki/Variant_Normalization#Left_alignment

gatk --java-options "-Xmx24G" LeftAlignAndTrimVariants \
	-V  ${results}/HKNPC_090_merged.vcf.gz \
	-O  ${results}/HKNPC_090_merged_leftalignandtrim.vcf \
	-R $reference

#Splitting Multi-allelic 
#bcftools norm --threads 8 -m - HKNPC_090_merged_leftalignandtrim.vcf   > HKNPC_090_merged_leftalignandtrim.decomposed.vcf
vt decompose -s ${results}/HKNPC_090_merged_leftalignandtrim.vcf -o ${results}/HKNPC_090_merged_leftalignandtrim.decomposed.vcf
#vt decompose -s /workspace/somatic/exome.merged.leftalignandtrim.vcf -o /workspace/somatic/exome.merged.leftalignandtrim.decomposed.vcf

echo "****************************FINISHED*******************************"
