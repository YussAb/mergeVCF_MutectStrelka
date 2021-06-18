#!/usr/bin/env nextflow
/*
 ================================================================================
                          yabili/mergeVCF_MutectStrelka
 ================================================================================
 Started June 2021.
 --------------------------------------------------------------------------------
 yabili/mergeVCF_MutectStrelka:
   This pipeline merge the vcf files genereated by GATK-Mutect2 and Illumina-Strelka2
 --------------------------------------------------------------------------------
  @Homepage
  https://github.com/YussAb/mergeVCF_MutectStrelka
 --------------------------------------------------------------------------------
*/
 
def helpMessage() {
     log.info"""
     Usage:
     The typical command for running the pipeline is as follows:
     nextflow mergeVCF_MutectStrelka.nf --strelka_snv path/to/strelka_snv --strelka_indels path/to/strelka_indels --mutect path/to/mutect --tumor_id mutectTumorName  --normal_id mutectNormalName   -with-report report.html
     Mandatory arguments:
     --strelka_snv      [path/to/strelka_snv]
     --strelka_indels   [path/to/strelka_indels]
     --mutect           [path/to/mutect]
     --tumor_id         [mutectTumorName]
     --normal_id        [mutectNormalName]
     """.stripIndent()    
 }
 
if (params.help) exit 0, helpMessage()


params.reference="/home/yabili/storage/references/gatk_bundle/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
params.input_tsv=""

Channel
    .fromPath(params.input_tsv)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(file(row.strelka_snvs), file(row.strelka_indels), file(row.mutect), row.tumor_id , row.normal_id , row.set_name) }
    .set { samples_ch }

samples_ch.into {samples_ch1; samples_ch2; samples_ch3; samples_ch4 }

//ricordarsi di aggiungere un conda env tra le configuration


/*
params.strelka_snv="/home/yabili/yussab/mergeVCF_MutectStrelka/test/HKNPC-076T_vs_HKNPC-076N/Strelka/StrelkaBP_HKNPC-076T_vs_HKNPC-076N_somatic_snvs.vcf.gz"
params.strelka_indels="/home/yabili/yussab/mergeVCF_MutectStrelka/test/HKNPC-076T_vs_HKNPC-076N/Strelka/StrelkaBP_HKNPC-076T_vs_HKNPC-076N_somatic_indels.vcf.gz"
params.mutect="/home/yabili/yussab/mergeVCF_MutectStrelka/test/HKNPC-076T_vs_HKNPC-076N/Mutect2/Mutect2_filtered_HKNPC-076T_vs_HKNPC-076N.vcf.gz"
params.outdir="results"
params.tumor_id="HKNPC-076T"
params.normal_id="HKNPC-076N"
params.set_name="HKNPC-076"

Channel.fromPath( params.strelka_snv, checkIfExists:true )
            .set{strelka_snv_ch}

Channel.fromPath( params.strelka_indels, checkIfExists:true )
            .set{strelka_indels_ch}


Channel.fromPath(params.mutect, checkIfExists:true )
            .set{mutect_ch}
*/



process illuminaMerge_snvIndels {
    conda '/home/yabili/yussab/mergeVCF_MutectStrelka/nextflow.yml' 
    publishDir "${params.outdir}/illuminaMerge_snvIndels", mode:'copy'
 
    input:
    // file (strelka_snv) from  strelka_snv_ch
    // file (strelka_indels) from strelka_indels_ch
    set file(strelka_snv), file(strelka_indels), file(mutect), tumor_id , normal_id , set_name from samples_ch1

    output:
    set file("${strelka_indels.simpleName}.snvs_indels.vcf.gz"), file("${strelka_indels.simpleName}.snvs_indels.vcf.gz.tbi" ) into illuminamerge_ch
    
    shell:
    '''
    zcat !{strelka_snv} | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description= Genotype >\\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > !{strelka_snv.simpleName}.snvs.gt.vcf
    
    zcat !{strelka_indels} | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description= Genotype >\\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > !{strelka_indels.simpleName}.indels.gt.vcf

    bgzip  !{strelka_snv.simpleName}.snvs.gt.vcf 
    tabix  !{strelka_snv.simpleName}.snvs.gt.vcf.gz
    bgzip  !{strelka_indels.simpleName}.indels.gt.vcf
    tabix  !{strelka_indels.simpleName}.indels.gt.vcf.gz
    bcftools concat -a -o  !{strelka_indels.simpleName}.snvs_indels.vcf.gz  -O z  !{strelka_snv.simpleName}.snvs.gt.vcf.gz  !{strelka_indels.simpleName}.indels.gt.vcf.gz
    tabix  !{strelka_indels.simpleName}.snvs_indels.vcf.gz
    '''
}


process mutect_preprocess {
    conda '/home/yabili/yussab/mergeVCF_MutectStrelka/nextflow.yml' 
    publishDir "${params.outdir}/mutect_preprocess", mode:'copy'
 

    input:
    // file (mutect) from mutect_ch
    set file(strelka_snv), file(strelka_indels), file(mutect), tumor_id , normal_id , set_name from samples_ch2

    output:
    file("${mutect.simpleName}_wgs.vcf") into mutectprocessed_ch

    shell:
    '''
    echo !{mutect}  >  wgs_vcf.fof

    #WARNING - you need the .tbi inde of ${mutect}
    
    tabix !{mutect}
    bcftools concat --allow-overlaps --remove-duplicates --file-list  wgs_vcf.fof --output-type z --output  !{mutect.simpleName}_wgs.vcf.gz

    tabix   !{mutect.simpleName}_wgs.vcf.gz
    gunzip  !{mutect.simpleName}_wgs.vcf.gz

    sed -i 's/!{normal_id}/NORMAL/' !{mutect.simpleName}_wgs.vcf   #create a variable for Normal sample
    sed -i 's/!{tumor_id}/TUMOR/' !{mutect.simpleName}_wgs.vcf     #create a variable for Tumor sample
    '''
} 


process combineVCF {
    conda '/home/yabili/miniconda3/envs/bioinfo' 
    publishDir "${params.outdir}/combineVCF", mode:'copy'
 
    input:
    set file(strelka_snv), file(strelka_indels), file(mutect), tumor_id , normal_id , set_name from samples_ch3
    file (mutect) from mutectprocessed_ch
    set file (strelka), file(tbi) from illuminamerge_ch

    output:
    set file("${set_name}_merged.vcf.gz"), file("${set_name}_merged.vcf.gz.tbi") into combinevcf_ch

    script:
    """
    java -Xmx24g -jar /home/yabili/miniconda3/envs/bioinfo/opt/gatk-3.8/GenomeAnalysisTK.jar -T CombineVariants \
	-R ${params.reference} \
	-genotypeMergeOptions PRIORITIZE \
	--rod_priority_list mutect,strelka \
	--variant:strelka ${strelka[0]} \
	--variant:mutect ${mutect}  \
	-o ${set_name}_merged.vcf.gz

    """
}


process dataNormalization {
    conda '/home/yabili/miniconda3/envs/bioinfo' 
    publishDir "${params.outdir}/dataNormalizedfinal", mode:'copy'
 
    input:
    set file(merged), file(tbi)  from combinevcf_ch
    set file(strelka_snv), file(strelka_indels), file(mutect), tumor_id , normal_id , set_name from samples_ch4

    
    output:
    file("${set_name}_merged.norm.pass_only.vcf") into normalized_ch

    script:
    """
    gatk --java-options "-Xmx24G" LeftAlignAndTrimVariants \
	-V  ${merged[0]} \
	-O  merged_leftalignandtrim.vcf \
	-R  ${params.reference}

    vt decompose -s merged_leftalignandtrim.vcf -o merged_leftalignandtrim.decomposed.vcf
 
    java -Xmx24g -jar /home/yabili/miniconda3/envs/bioinfo/opt/gatk-3.8/GenomeAnalysisTK.jar -T SelectVariants\
        -R ${params.reference}\
        --excludeFiltered --variant merged_leftalignandtrim.decomposed.vcf\
        -o ${set_name}_merged.norm.pass_only.vcf
    """
}

