# ----------------------- #
# 设置全局变量（路径等）
# ----------------------- #

# 软件路径
gatk = "/home/CJM/seq_tools/gatk/GenomeAnalysisTK.jar"  ##3.8版本的gatk
java="/home/CJM/seq_tools/java/jre1.8.0_241/bin/java"   ##1.8版本的java
beagle = "/home/CJM/seq_tools/beagle/beagle.22Jul22.46e.jar"  
bcftools = "/home/CJM/seq_tools/bcftools/bcftools-1.8/bcftools" 
# 参考基因组
REF = "/home/CJM/sequence/reference/ref.fasta"  

# ----------------------- #
# 定义样本列表和染色体区域
# ----------------------- #

SAMPLES = ["1","2"]  # 需要分析的样本列表

# ----------------------- #
# # 定义目标规则
# # ----------------------- #

rule all:
     input:
         "./vcf/combined.SNP.clean.vcf.gz",
         "./vcf/combined.INDEL.clean.vcf.gz"

# 1. generate gvcf file for indi in chr region

rule generate_gvcf:
    input:
        fa=REF,  
        bam="{sample}.rmdup.bam",  # 这里更新为 clean.bam
        bai="{sample}.rmdup.bam.bai"  # 这里更新为 clean.bam.bai 
    output:
        gvcf="./gvcf/{sample}.g.vcf.gz"  # 每个样本每条染色体的 GVCF
    threads: 50
    shell:
        '''
        {java} -Xmx200g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar {gatk} -T HaplotypeCaller -ERC GVCF \
            -R {input.fa} \
            -I {input.bam} \
            -o {output.gvcf} \
            -nct {threads}
        '''
# 2. CombineGVCFs for all samples
rule combine_gvcfs:
    input:
        ref=REF,
        gvcfs=expand("./gvcf/{sample}.g.vcf.gz", sample=SAMPLES)
    output:
        combined_gvcf="./vcf/combined.g.vcf.gz"
    params:
        gvcf_inputs=lambda wildcards, input: " ".join(f"-V {g}" for g in input.gvcfs)
    shell:
        '''
        {java} -Xmx8g -jar {gatk} -T CombineGVCFs \
            -R {input.ref} \
            {params.gvcf_inputs} \
            -o {output.combined_gvcf} 
        '''

# 3. GenotypeGVCFs on the combined GVCF
rule genotype_gvcfs:
    input:
        ref=REF,
        combined_gvcf="./vcf/combined.g.vcf.gz"
    output:
        combined_vcf="./vcf/combined.vcf.gz"
    threads: 16
    shell:
        '''
        {java} -Xmx8g -jar {gatk} -T GenotypeGVCFs \
            -R {input.ref} \
            -V {input.combined_gvcf} \
            -o {output.combined_vcf} \
	    -nt {threads}
        '''

# 4.combine genotype in single chr and SelectVariants

rule select_snp_variants:
    input:
        ref=REF,
        combined_vcf="./vcf/combined.vcf.gz"
    output:
        snp_vcf="./vcf/combined.SNP.vcf.gz",
        indel_vcf="./vcf/combined.INDEL.vcf.gz"
    threads: 50
    shell:
        '''
        {java} -jar -Xmx64g {gatk} -T SelectVariants \
            -R {input.ref} \
            -V {input.combined_vcf} \
            -selectType SNP \
            -o {output.snp_vcf} \
            -nt {threads}
			
	{java} -jar -Xmx64g {gatk} -T SelectVariants \
            -R {input.ref} \
            -V {input.combined_vcf} \
            -selectType INDEL \
            -o {output.indel_vcf} \
            -nt {threads}	
        '''

# 5. Mark Hard Filter for SNPs
rule mark_hard_filter_snp:
    input:
        ref=REF,
        snp_vcf="./vcf/combined.SNP.vcf.gz",
        indel_vcf="./vcf/combined.INDEL.vcf.gz"
    output:
        snp_filtered_vcf=temp("./vcf/combined.SNP.filtered.vcf.gz"),
        indel_filtered_vcf=temp("./vcf/combined.INDEL.filtered.vcf.gz"),
        snp_clean_vcf="./vcf/combined.SNP.clean.vcf.gz",
        indel_clean_vcf="./vcf/combined.INDEL.clean.vcf.gz"
    shell:
        '''
        {java} -Xmx64g -jar {gatk} -T VariantFiltration \
            -R {input.ref} \
            -V {input.snp_vcf} \
	    -filter "QD<2.0 || MQ<40.0 || FS>60.0 || SOR> 3.0 || MQRankSum< -12.5 || ReadPosRankSum< -8.0" \
            -filterName  "Filter" \
            -o {output.snp_filtered_vcf}
			
	{java} -Xmx8g -jar {gatk} -T VariantFiltration \
            -R {input.ref} \
            -V {input.indel_vcf} \
	    -filter "QD<2.0 || FS>200.0 || SOR> 10.0 || MQRankSum< -12.5 || ReadPosRankSum< -8.0" \
            -filterName "Filter" \
            -o {output.indel_filtered_vcf}

        zgrep -v -P "\\tFilter" {output.snp_filtered_vcf} | bgzip > {output.snp_clean_vcf}
        zgrep -v -P "\\tFilter" {output.indel_filtered_vcf} | bgzip > {output.indel_clean_vcf}
       '''
