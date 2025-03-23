
    library(DESeq2)
    df <- read.csv("            Abiotrophia.defectiva  Abiotrophia.sp.HMSC24B09  ...  Zhihengliuella.halotolerans  Zygoascus.ofunaensis
SampleID                                                     ...                                                   
N01_1_2                         0                         0  ...                            0                     0
N01_1_3                         0                         0  ...                            0                     0
N01_1_4                         0                         0  ...                            0                     0
N01_2_2                        89                       130  ...                            0                     0
N01_2_3                         0                         0  ...                            0                     0
...                           ...                       ...  ...                          ...                   ...
ZJH_N8_2_3                      0                         0  ...                            0                     0
ZJH_N8_2_4                      0                         0  ...                            0                     0
ZJH_N9_1_3                      0                         0  ...                            0                     0
ZJH_N9_2_2                      0                         0  ...                            0                     0
ZJH_N9_2_4                      0                         0  ...                            0                     0

[578 rows x 2047 columns]", row.names=1)
    dds <- DESeqDataSetFromMatrix(countData = df, colData = NULL, design = ~ 1)
    vsd <- varianceStabilizingTransformation(dds)
    write.csv(assay(vsd), "microbiome_abundance_vst.csv")
    