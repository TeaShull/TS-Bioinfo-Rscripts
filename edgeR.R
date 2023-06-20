setwd("./")
library('edgeR')
library('RColorBrewer')
library("GenomicFeatures")
library("tximport")
library("readr")
library("tximport")
library("edgeR")

#### make TxDB for TAIR ####

#gtffile <- "./TAIR10_GFF3_genes.gff"
#file.exists(gtffile)
#txdb <- makeTxDbFromGFF(gtffile, format = "gff3", circ_seqs = character())
#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")
#write.csv(tx2gene, file = "./TAIR10tx2gene.gencode.v27.csv")
#head(tx2gene)

#File Variables, to avoid having to change the file names
TP = "2H"

#### Import data using tximport ####

#Read in tx2gene database
tx2gene <- read_csv("./TAIR10tx2gene.gencode.v27.csv")

#Read in sample list (file names and factors)
samples <- read.table(file = paste0(TP, "/samples.txt"), header = T)
head(samples)

#Retrieve file paths for salmon .sf files
files <- file.path(TP, paste0(samples$samples,"_quant.sf") )
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
head(files)

#Import salmon .sf files, and summarize to gene using tx2gene databas
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

#### Prep Data for EdgeR ####

cts <- txi$counts
head(cts)
seqDataGroups <- c(paste0(TP,"_0", TP, "_0", TP, "_0", TP,"_50", TP, "_50", TP, "_50"))
seqDataGroups

# use edgeR function DGEList to make read count list
d <- DGEList(counts=cts,group=factor(seqDataGroups))
head(d)

#### Data Filtering ####

dim(d)
d.full <- d # keep in case things go awry
head(d$counts)
head(cpm(d))

#total counts/sample
apply(d$counts,2,sum)

#Trim lowly-abundant genes. cpm number can be changed based on library size, and detected tags.
keep <- rowSums(cpm(d)>100) >= 2
keep
d <- d[keep,]
dim (d)

d$samples$lib.size <- colSums(d$counts)
d$samples

d <- calcNormFactors(d, logratioTrim = 0)
d

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))

d1 <- estimateCommonDisp(d, verbose=T)
names(d1)

d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")

# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
topTags(et12,n=100)


#de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
#for GSEA, uncomment above for standard analysis

de1 <- decideTestsDGE(et12, adjust.method="BH")
summary(de1)

# differentially expressed tags from the naive method in d1
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
design.mat
fit <- glmFit(d2, design.mat)
summary(de1)

DE <- topTags(et12,n=50000)
as.data.frame(DE)
write.csv(DE, paste0("./edgeRout_", TP, "_LRT0.csv"))

