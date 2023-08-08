setwd("/run/media/eric/analysis/RNAseq_DA/seqAnalysis/ClusterProfiler")
#### Libraries ####

library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.At.tair.db"
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)

# reading in data from EdgeR
df = read.csv("./edgeRout_2H_LRT0_all.csv", header=TRUE)
head(df)

# we want the log2 fold change
original_gene_list <- df$logFC

# name the vectorhttp://127.0.0.1:18303/graphics/plot_zoom_png?width=929&height=3484
names(original_gene_list) <- df$mapped_id

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             nPermSimple = 100000,
             keyType = "TAIR", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr",
             eps = 0)

#install.packages('ggridges')

gse2 <- simplify(gse, cutoff=0.65,
                by="p.adjust",
                select_fun=min,
                measure = "Wang")


ridgeplot(gse) + labs(x = "enrichment distribution")

require(DOSE)
dotplot(gse2, showCategory=200, split=".sign") + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)


# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="geneNum", foldChange=gene_list, showCategory = 3)


#### FOR KEGG ANALYSIS ####

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "TAIR", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("TAIR")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$mapped_id %in% dedup_ids$TAIR,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# create gseKEGG object

kegg_organism = "ath"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

ridgeplot(kk2) + labs(x = "enrichment distribution")
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

