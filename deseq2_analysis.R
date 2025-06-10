
#加载必要包
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(EnhancedVolcano)
library(ggrepel)
library(pheatmap)



# 加载表达矩阵
counts <- read.csv("gene_count_matrix.csv", row.names = 1)
colnames(counts) <- gsub("\\.\\d+$", "", colnames(counts))  # 去掉 Ensembl 编号中的版本号（可选）

# 构建样本信息表
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("A", 3), rep("B", 3))  # 根据样本数量和顺序填写
)

# 创建 DESeq2 数据对象
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)


# 运行 DESeq2 差异表达分析
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), "DEG_results.csv")
#============================ #
# 1. 统计差异基因数量
# ============================ #
padj_cutoff <- 0.05
logfc_cutoff <- 1

# 筛选出显著的差异表达基因
sig_DEGs <- subset(res, padj < padj_cutoff & abs(log2FoldChange) > logfc_cutoff)
write.csv(as.data.frame(sig_DEGs), file = "significant_DEGs.csv")

# 计算差异基因数量
total_DEGs <- sum(res$padj < padj_cutoff & abs(res$log2FoldChange) > logfc_cutoff, na.rm = TRUE)
up_DEGs <- sum(res$padj < padj_cutoff & res$log2FoldChange > logfc_cutoff, na.rm = TRUE)
down_DEGs <- sum(res$padj < padj_cutoff & res$log2FoldChange < -logfc_cutoff, na.rm = TRUE)

# 打印统计结果到终端
cat("Total significant DEGs (padj < 0.05 & |log2FC| > 1):", total_DEGs, "\n")
cat("Upregulated genes:", up_DEGs, "\n")
cat("Downregulated genes:", down_DEGs, "\n")


# ============================ #
# 2. PCA 图
# ============================ #
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot")
ggsave("PCA_plot.png", p, width=6, height=5)

# ============================ #
# 3. 火山图
# ============================ #
res_df <- as.data.frame(res)  # 将 DESeqResults 转成 data.frame

# 给显著性加标签
res_df$threshold <- with(res_df, ifelse(padj < padj_cutoff & abs(log2FoldChange) > logfc_cutoff,
                                        ifelse(log2FoldChange > logfc_cutoff, "Up", "Down"),
                                        "Not Sig"))

# 画火山图
volcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color=threshold)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="grey", "Up"="red")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black")

# 保存火山图
ggsave("volcano_plot.png", volcano, width=6, height=5)

# ============================ #
# 4. 韦恩图
# ============================ #
up_genes <- rownames(subset(res, padj < padj_cutoff & log2FoldChange > logfc_cutoff))
down_genes <- rownames(subset(res, padj < padj_cutoff & log2FoldChange < -logfc_cutoff))

venn <- venn.diagram(
  x = list(Up=up_genes, Down=down_genes),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.2,
  main = "Venn Diagram of DEGs"
)
png("Venn_plot.png", width=800, height=600)
grid.draw(venn)
dev.off()

# ============================ #
# 5. 热图（top 50 genes）
# ============================ #
topGenes <- head(order(res$padj), 50)  # 选择前50个padj最小的基因
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)  # 对每个基因进行去均值化处理
annotation_col <- data.frame(condition = colData(dds)$condition)
rownames(annotation_col) <- colnames(mat)

# 保存热图
pheatmap(mat, annotation_col = annotation_col,
         fontsize_row = 6, show_rownames = TRUE,
         filename = "Heatmap_top50.png")
