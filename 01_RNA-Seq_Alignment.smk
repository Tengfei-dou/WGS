# ----------------------- #
# 设置全局变量（路径等）
# ----------------------- #

# 软件路径
trimmomatic = "/home/CJM/miniconda3/bin/trimmomatic" 
STAR="/home/CJM/miniconda3/bin/STAR"
stringtie="/home/CJM/miniconda3/bin/stringtie"
featureCounts="/home/CJM/miniconda3/bin/featureCounts"
fastp="/home/CJM/seq_tools/fastp/fastp"
# 参考基因组
REF = "/home/CJM/sequence/reference/ref.fasta"  
GTF = "/home/CJM/sequence/reference/Sus_scrofa.Sscrofa11.1.113.gtf"
# ----------------------- #
# 定义样本列表和染色体区域
# ----------------------- #

SAMPLES = ["A1","A2","A3","B4","B5","B6"]  # 需要分析的样本列表

# ----------------------- #
# # 定义目标规则
# # ----------------------- #

rule all:
     input:
           expand("{sample}/STARgenome/{sample}.gtf", sample=SAMPLES),
           expand("{sample}/STARgenome/{sample}.tsv", sample=SAMPLES),
           expand("{sample}/STARgenome/{sample}.featureCounts.txt", sample=SAMPLES),
           "gene_count_matrix.csv",
           "transcript_count_matrix.csv"
# 1.构建索引
rule build_reference:
    input:
        fa=REF,  # 猪基因组序列（FASTA）
        gtf=GTF        # 基因注释文件（GTF）
    output:
        index_dir=directory("ref_genome")                   
    threads: 50
    shell:
        '''
        mkdir -p ref_genome

        {STAR} --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output.index_dir} \
             --genomeFastaFiles {input.fa} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 100
        '''
# 2.质控比对
rule qc_and_map:
    input:
        r1="{sample}_1.fq.gz",
        r2="{sample}_2.fq.gz",
        index_dir="ref_genome"
    output:
        r1_clean="{sample}/{sample}_1.clean.fq.gz",
        r2_clean="{sample}/{sample}_2.clean.fq.gz",
        bam="{sample}/Aligned.sortedByCoord.out.bam"  
    threads: 10
    shell:
        '''
        mkdir -p {wildcards.sample}

            # 双端数据：去除低质量序列和接头
             {fastp} -i {input.r1} -I {input.r2} -o {output.r1_clean} -O  {output.r2_clean} -w {threads}

            # 双端数据比对
            {STAR} --runThreadN {threads} \
                --genomeDir ref_genome \
                --sjdbGTFfile {GTF} \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesCommand zcat \
                --outFilterMismatchNmax 3 \
                --readFilesIn {output.r1_clean} {output.r2_clean} \
                --outFileNamePrefix {wildcards.sample}/
        '''
# 3.基因表达量的定量分析
rule expr_gene:
    input:
        bam="{sample}/Aligned.sortedByCoord.out.bam"  # STAR 比对输出的 BAM 文件
    output:
        "{sample}/STARgenome/{sample}.gtf",  # stringtie 输出的 GTF 文件
        "{sample}/STARgenome/{sample}.tsv",  # stringtie 输出的 FPKM 表
        "{sample}/STARgenome/{sample}.featureCounts.txt"  # featureCounts 输出的计数文件
    threads: 50
    shell:
        '''       
        mkdir -p {wildcards.sample}/STARgenome

        # 使用 stringtie 进行 FPKM 定量
        {stringtie} -p {threads} -e -B -G {GTF} -o {output[0]} -A {output[1]} {input.bam}
        
        # 使用 featureCounts 进行基因计数
        {featureCounts} -T {threads} -p -t exon -g gene_id -a {GTF} -o {output[2]} {input.bam}
        '''
# ----------------------- #
# 4. 自动生成 sample_list.txt 文件
# ----------------------- #

rule generate_sample_list:
    input:
        gtf_files = expand("{sample}/STARgenome/{sample}.gtf", sample=SAMPLES)
    output:
        "sample_list.txt"
    run:
        with open(output[0], "w") as f:
            for sample in SAMPLES:
                gtf_path = f"{sample}/STARgenome/{sample}.gtf"
                f.write(f"{sample} {gtf_path}\n")

# ----------------------- #
# 5. 使用 prepDE.py 生成表达矩阵
# ----------------------- #

rule prepDE:
    input:
        sample_list = "sample_list.txt"
    output:
        gene_matrix = "gene_count_matrix.csv",
        transcript_matrix = "transcript_count_matrix.csv"
    params:
        script = "/home/CJM/DTF/RNA_Temp/prepDE.py"  # 改成你实际的 prepDE.py 路径
    shell:
        """
        python2 {params.script} -i {input.sample_list} -g {output.gene_matrix} -t {output.transcript_matrix}
        """
