# ----------------------- #
# 设置全局变量（路径等）
# ----------------------- #

# 软件路径
picard = "/home/CJM/miniconda3/bin/picard"
samtools = "/home/CJM/seq_tools/samtools/samtools-1.9/samtools" 
# ----------------------- #
# 定义样本列表
# ----------------------- #

SAMPLES = ["SRR9651232","SRR9651233"] # 需要分析的样本列表

# ----------------------- #
# # 定义目标规则
# ----------------------- #

rule all:
    input:
        expand("./{sample}.rmdup.bam", sample=SAMPLES)

# 1. Picard CleanSam

rule picard_clean_bam:
    input:
        rmdup_bam="{sample}.rmdup.bam"
    output:
        rmdup_bam="{sample}.rmdup.bam",
        bam_index="{sample}.rmdup.bam.bai"
    threads: 4
    shell:
        '''
        {picard} CleanSam I={input.rmdup_bam} O=temp.{wildcards.sample}.clean.bam && \
        mv temp.{wildcards.sample}.clean.bam {output.rmdup_bam} && \
        samtools index {output.rmdup_bam}
        '''





		
