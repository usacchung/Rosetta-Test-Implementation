library(DESeq2)
library(airway)

data("airway")

dds <- DESeqDataSet(airway, design = ~ dex)

dds <- DESeq(dds)

res <- results(dds)

top5 <- head(res[order(res$padj), ], 5)

print("--- DESeq2 Easy Test Results (Top 5 Genes) ---")
print(top5)