# Ontology and genenames
GO_2H <- read.table(file = "analysis_final/2H_GO_top2000.tsv", header = TRUE, sep = "\t", fill = TRUE )
GN_2H <- read.table(file = "analysis_final/genenames_2H.tsv", header = TRUE, sep = "\t", quote="", fill = TRUE)

# Genenames to Results
GN_Results_2H <- merge(x = results_2H, y = GN_2H, by.x = "mapped_id", all = TRUE)
GN_Results_8H <- merge(x = results_8H, y = GN_8H, by.x = "mapped_id", all = TRUE)

write.csv(GN_Results_2H, file = "", row.names = T)
write.csv(GN_Results_8H, file = "", row.names = T)

# Merge: gene names, specfic gene ontology
sGO <- read.table(file = "", header = TRUE, sep = "\t", fill = TRUE)
Gn_2H_p05 <- merge(x = sGO, y = results_2H, by.x = "mapped_id")
Gn_8H_p05 <- merge(x = sGO, y = results_8H, by.x = "mapped_id")

#making heatmaps
head(data)
data <- read.table(file = "", header = TRUE, sep = "\t")
data_subset <- as.matrix(data)
library(pheatmap)
pheatmap(data_subset,
         annotation_row = genenames,
         annotation_col = GO)
