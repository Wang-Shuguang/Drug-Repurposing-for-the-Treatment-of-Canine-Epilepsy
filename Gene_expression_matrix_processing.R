setwd("C:\\GSE7486")
library(affy)
dir_cels='C:\\GSE7486\\dir_files'
raw_data <- ReadAffy(celfile.path=dir_cels)
library(simpleaffy)
Data.qc <- qc(raw_data);
plot(Data.qc)
library(affyPLM);
library(RColorBrewer) ;
Pset <-fitPLM(raw_data);
colors <- brewer.pal(12, "Set3");
Mbox(Pset, ylim = c(-1, 1), col = colors, main = "RLE", las = 3) ;
boxplot(Pset, ylim = c(0.95, 1.22), col = colors, main = "NUSE", las = 3)
data <- rma(raw_data)
library(limma) ;
sampleNames(data) <-  gsub(".CEL.gz$", "", sampleNames(data));
data <- data[, -match(c("GSM181371", "GSM181372", "GSM181373", 
                        "GSM181374", "GSM181375", "GSM181376", 
                        "GSM181377", "GSM181378", "GSM181379", "GSM181380"), 
                      sampleNames(data))];
eset <- exprs(data) ;
design <- model.matrix(~ 0+factor(c(1,1,1,1,1,1,1,1,
                                    2,2,2,2,2,2,2,2,2,2,2,2)))
colnames(design) <- c("group1", "group2")
contrast.matrix <- makeContrasts (contrasts = "group1 - group2", levels = design);
fit <- lmFit(eset, design) ;
fit1 <- contrasts.fit(fit, contrast.matrix) ; 
fit2 <- eBayes(fit1) ;
dif <- topTable(fit2, coef = "group1 - group2",  n = nrow(fit2), lfc = log2(1.5)) ; 

library(annotate) ;
affydb <- annPkgName(raw_data@annotation, type = "db");
library(affydb, character.only = TRUE) ;

dif$symbols <- getSYMBOL(rownames(dif), affydb) ;
dif$EntrezID <- getEG(rownames(dif), affydb) ;
write.csv(dif, file = "Gene_expression_matrix.csv")


