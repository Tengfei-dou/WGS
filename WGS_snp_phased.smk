# ----------------------- #
# 设置全局变量（路径等）
# ----------------------- #

# 软件路径
plink = "/home/CJM/seq_tools/plink1.9/plink" ##1.9版本的plink
beagle = "/home/CJM/seq_tools/beagle/beagle.22Jul22.46e.jar"   

# ----------------------- #
# # 定义目标规则
# # ----------------------- #

rule all:
     input:
        fam = "./vcf/snp.fam",
        bim = "./vcf/snp.bim",
        bed = "./vcf/snp.bed"

# 1. extract_chr（1-18, X, Y）
rule extract_chr:
    input:
        vcf = "./vcf/combined.SNP.clean.vcf.gz"
    output:
        vcf = "./vcf/combined.SNP.clean.chr.vcf"
    params:
        out_prefix = "./vcf/combined.SNP.clean.chr"
    shell:
        '''
        {plink} --vcf {input.vcf} \
               --chr 1-18,X,Y \
               --recode vcf-iid \
               --out {params.out_prefix}
        '''

# 2. phase
rule phase:
    input:
        vcf = "./vcf/combined.SNP.clean.chr.vcf"
    output:
        vcf = "./vcf/combined.SNP.clean.chr.phased",
        temp_vcf=temp("temp.vcf")
    params:
        out_prefix = "./vcf/combined.SNP.clean.chr.phased",
        seed = 9823
    threads: 50
    shell:
        '''
        awk '$1 ~ /^#/ || $4 != "*"' {input.vcf} > {output.temp_vcf}

        java -Xmx128g -jar {beagle} \
            gt={output.temp_vcf} \
            out={output.vcf} \
            seed={params.seed} \
            nthreads={threads}
        '''
# 3.binary
rule convert_binary:
    input:
        vcf = "./vcf/combined.SNP.clean.chr.phased.vcf.gz"
    output:
        "./vcf/snp.fam",
        "./vcf/snp.bim",
        "./vcf/snp.bed"
    shell:
        '''
        {plink} --vcf {input.vcf} \
               --make-bed \
               --out ./vcf/snp
        '''
