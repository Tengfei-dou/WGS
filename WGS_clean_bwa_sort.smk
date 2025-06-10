# ----------------------- #
# 设置全局变量（路径等）
# ----------------------- #

# 软件路径
fastp="/home/CJM/seq_tools/fastp/fastp"
bwa = "/home/CJM/seq_tools/bwa/bwa-0.7.17/bwa"  
samtools = "/home/CJM/seq_tools/samtools/samtools-1.9/samtools"  
sambamba = "/home/CJM/miniconda3/bin/sambamba"

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
             expand("{sample}.rmdup.bam.bai", sample=SAMPLES)

# 1. clean data

rule clean_data:
    input:
        r1="{sample}_1.fastq.gz",
        r2="{sample}_2.fastq.gz"
    output:
        r1_fq="{sample}_1.fq.fastq.gz",
        r2_fq="{sample}_2.fq.fastq.gz"
    threads: 16
    shell:
        "{fastp} -i {input.r1} -I {input.r2} -o {output.r1_fq} -O {output.r2_fq} -w {threads} "

		
# 2. bwa align and sorted (bwa v0.7.5a-r405, samtools v1.9)

rule bwa_align_and_sorted:
    input:
        r1_fq="{sample}_1.fq.fastq.gz",
        r2_fq="{sample}_2.fq.fastq.gz",
       ref=REF
    output:
       unsorted_bam=temp("{sample}.unsort.bam"),
       sorted_bam=temp("{sample}.sort.bam"),
       rmdup_bam="{sample}.rmdup.bam",
       rmdup_index = "{sample}.rmdup.bam.bai"
    threads: 50
    shell:
        '''
        {bwa} mem -t {threads} -M -R  "@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:Illumina\\tPU:{wildcards.sample}\\tSM:{wildcards.sample}" \
		{input.ref} {input.r1_fq} {input.r2_fq} | \
		{samtools} view -b -o  {output.unsorted_bam} - && \
        {samtools} sort -@ {threads} -O bam -o {output.sorted_bam} {output.unsorted_bam} && \
	    {sambamba} markdup  --remove-duplicates -t {threads} --overflow-list-size 100000000 --tmpdir=/tmp {output.sorted_bam} {output.rmdup_bam} &&\
        {samtools} index {output.rmdup_bam}
        '''
