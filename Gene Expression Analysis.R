file_name = "brca_tcga_pan_can_atlas_2018.tar.gz"
untar(file_name)
setwd(".\\brca_tcga_pan_can_atlas_2018")
getwd

#----import data and preprocess-------
clinical = read.delim("data_clinical_patient.txt",skip = 4, header = TRUE)

rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

#delete the genes for which there's more than one Hugo Symbol for simplicity
keep = !duplicated(rnaseq[,1])

rnaseq = rnaseq[keep,]
# delete Hugo_Symbol is "" row
rnaseq <- rnaseq[rnaseq$Hugo_Symbol != "", , drop = FALSE]

rownames(rnaseq)  = rnaseq[,1]

#load cna
cna = read.delim("data_cna.txt")
# find ERBB2 in cna

erbb2_indx = which(cna[,1] == 'ERBB2')

# Plot histogram to visualize explore the data.
par(mar = c(5, 4, 4, 2) + 0.1)

hist(as.numeric(cna[erbb2_indx,-c(1,2)]))

# match patients in rnaseq to patients in cna.
#which return index that is TRUE
rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases/patients which have cna data.(1068 rows)
rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna(total 1072 row,1068 true),all in there.
no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 

# sanity check.This will print an error if the result is not the same.
sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]
sanity_check

#prepare for survival_data
patientID<-colnames(rna_cna_sub)

#select columns
rownames(clinical) <-clinical[, 1]
survival_data<-clinical[,c("OS_STATUS","OS_MONTHS")]
survival_data$OS_STATUS<-factor(survival_data$OS_STATUS)

#select rows
converted_patientID <- gsub("\\.", "-", patientID)  # align patient id
converted_patientID <- gsub("\\-01$", "", converted_patientID)  # 

survival_data <- subset(survival_data, rownames(survival_data) %in%converted_patientID)
# If OS_STATUS is character, convert it to a numeric value
survival_data$OS_STATUS <- as.numeric(survival_data$OS_STATUS)-1

str(survival_data$OS_STATUS)
table(survival_data$OS_STATUS)
dim(survival_data)

# Pre-allocate memory for ERBB2+ cna
meta_erbb2 = matrix(0,length(rna_cna_id),1)
for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index(colnames) in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
  
}

hist(meta_erbb2)

#CHECK
# simple checks to make sure. 
col_i = colnames(rna_cna_sub)[1]
col_cna = which(colnames(cna)==col_i)

# sanity check
(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]
# see now if a positive meta_erbb2 is amplified.
pos_example = which(meta_erbb2==1)[1]
col_i = colnames(rna_cna_sub)[pos_example]
col_cna = which(colnames(cna)==col_i)

# sanity check
(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]
# both checks should print true.


# We will add a title to the metadata.
colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers
rna_cna_sub = round(rna_cna_sub)

# Install DESeq2.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install DeSeq2
BiocManager::install("DESeq2")
library(DESeq2)

#--------1. Differential Gene Expression：-------
dds <- DESeqDataSetFromMatrix(
  countData = rna_cna_sub,   # gene expression matrix
  colData = data.frame(erbb2_amplified =as.factor(meta_erbb2)),  # 包含元数据的数据框，这里是 ERBB2 放大信息
  design = ~ erbb2_amplified  # 设计矩阵，指定差异分析的模型，这里使用 ERBB2 放大信息作为差异的预测变量
)
# run DESeq2 analysis and normalize data  
dds <- DESeq(dds)

# Differentially expressed genes were obtained and compared between HER2 amplification and non-amplification groups
result <- results(dds)
result
class(result)
dim(result)#
# look at the differentially expressed genes
head(result)
result[1:5, ]
rownames(result[1:5,])
#genes with significant difference
# Count the number of significantly differentially expressed genes
num_significant_genes <- sum(!is.na(result$padj) & result$padj < 0.05)
num_significant_genes
# Print the result
cat("Number of significantly differentially expressed genes:", num_significant_genes, "\n")

#Volcano plot
library(ggplot2)
# Plot scatter plots of differentially expressed genes
ggplot(as.data.frame(result), aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not Significant"))) +
  labs(title = "Differential Expression Analysis",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_minimal()

#MA Plot:
#An MA plot (M stands for log-ratio, A stands for average abundance) is another option. It shows the relationship between the log2-fold change and the average expression level.
ggplot(as.data.frame(result), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not Significant"))) +
  labs(title = "MA Plot",
       x = "Average Expression",
       y = "log2(Fold Change)") +
  theme_minimal()

# select significantly differentially expressed genes with log2FoldChange
diff_genes <- as.data.frame(result[which(result$padj < 0.05), c("log2FoldChange")]) 
colnames(diff_genes)<-"log2FoldChange"
# Plot scatter plots of significantly differentially expressed genes
ggplot(diff_genes, aes(x = rownames(diff_genes), y = log2FoldChange, color = log2FoldChange > 0)) +
  geom_point(size = 2) +
  labs(title = "Differentially Expressed Genes",
       x = "Gene",
       y = "log2FoldChange") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#------2. Top 10 Differentially Expressed Genes Ranked by Fold Change：------
# Extract top 10 differentially expressed genes 
top10_genes <- head(order(abs(result$log2FoldChange), decreasing = TRUE), 10)
top10_genes_data <- result[top10_genes, ]
top10_genes_data

gene_color <- c("#FFB6C1", "#FFD700", "#FFA07A", "#00FA9A", "#40E0D0", "#9370DB", "#FF6347", "#00BFFF", "#98FB98", "#FF4500")
ggplot(as.data.frame(top10_genes_data), aes(x = log2FoldChange, y = -log10(pvalue),label = rownames(top10_genes_data))) +
  geom_point(size=3, aes(color = gene_color)) +
  geom_text(hjust = -0.1, vjust = 0) + 
  labs(title = "Top10 Genes",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_minimal()+
  guides(color = FALSE)
#------3. Pathway Enrichment：------
#A pathway enrichment analysis tool, such as Enrichr, GSEA and clusterProfiler, was used to analyze the differentially expressed genes. Here is an example of clusterProfiler:
  
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(clusterProfiler)

#Enrichment analysis using differential expression results  
  
# Format gene id
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# Use bitr with OrgDb argument,Convert Hugo Gene Symbols to Entrez Gene IDs把基因名Hugo转化成数字
gene_ids <- bitr(rownames(result), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform pathway enrichment analysis
enrich_result <- enrichKEGG(gene = gene_ids$ENTREZID, 
                            organism = "hsa",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH")
enrich_result
summary(as.data.frame(enrich_result) )

output_file <- "..//tables//Summary of enrichment result.xlsx"
write.xlsx(summary(as.data.frame(enrich_result) ), file = output_file, sheetName = "TopGenes", rowNames = TRUE, colNames = TRUE)
dotplot(enrich_result, showCategory=20)
cnetplot(enrich_result) 

#----4.Get the variance stabilized transformed expression values-----
#  variance stabilization transforms are primarily used for exploratory data analysis and visualization.
vst_data <- vst(dds)

# Extract the variance-stabilized expression values
vst_expression <- assay(vst_data)

# transpose the result, so that each row represent a gene, and each col represent a sample/patient
vst_expression <- t(vst_expression)
# Now, vst_expression contains expression values with stable variance, with each row corresponding to one gene and each column corresponding to one sample.

#-------5. PCA Plot：-----
#Preprocess
# Delete the constant column
vst_expression <- vst_expression[, apply(vst_expression, 2, function(x) length(unique(x)) > 1)]
#Delete all zero columns
vst_expression <- vst_expression[, colSums(vst_expression != 0) > 0]
# calculate PCA,  scale. = TRUE indicates normalization of the data before analysis
pca_result <- prcomp(vst_expression,center = TRUE, scale. = TRUE)
meta_erbb2<-factor(meta_erbb2)
summ<-summary(pca_result)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")

ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

library(ggplot2)
pca_plot <- ggplot(data = as.data.frame(pca_result$x), aes(x = PC1, y = PC2, color =meta_erbb2)) +
  geom_point(size = 3) +
  ggtitle("PCA Plot of VST Values") +
  labs(x=xlab,y=ylab)+
  theme_minimal()+
  stat_ellipse(level = 0.95, geom = "polygon", aes(fill = meta_erbb2), alpha = 0.2)
pca_plot

#--------6. Gene Expression Clustering:-----
#The purpose of  Gene expression clustering is to partition samples into groups with similar gene expression patterns. The commonly used clustering methods include hierarchical clustering, K-means clustering and so on. Here is an example of hierarchical clustering:
# Using differentially expressed gene data
gene_data <- assay(dds)
cat("There are", nrow(gene_data), "rows and", ncol(gene_data), "columns in gene_data.")
# Calculate the distance matrix of the gene
dist_matrix <- dist(t(gene_data))

# carry out hierarchical clustering 
clustering_result <- hclust(dist_matrix, method = "ward.D2")

print("The hierarchical clustering method is as below:") 
print(clustering_result)
print("Show the basic statistical inforamtion of hierarchical clustering result")
summary(clustering_result)

# plot tree
plot(clustering_result, main = "Hierarchical Clustering Tree")

# Cut the tree to get clusters
clusters <- cutree(clustering_result, k)

table(clusters)
cluster_counts <- table(clusters)  
  
for (i in 1:k) {
  cat("Cluster", i, "has", cluster_counts[i], "genes.\n")
}
# Perform PCA
pca_result <- prcomp(t(gene_data))

# Plot PCA with clusters
pca_plot <- ggplot(data = as.data.frame(pca_result$x), aes(x = PC1, y = PC2, color = factor(clusters))) +
  geom_point(size = 3) +
  ggtitle("PCA Plot with Gene Expression Clusters") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") +
  theme_minimal()

# Display PCA plot
print(pca_plot)  

#---------7. Cox survival regression model and differentially expressed genes：-------

#1.	With the vst values of the DE genes generate an overall survival model.-----
#The Cox proportional hazards model had to be positive for time to event, removing rows with zero survival time
survival_data <- survival_data[survival_data$OS_MONTHS > 0, ]
str(survival_data)
# Filters the vst_expression line
# String substitution using the gsub function
rownames(vst_expression) <- gsub("\\.", "-", rownames(vst_expression))  # 将点替换为破折号
rownames(vst_expression) <- gsub("\\-01$", "", rownames(vst_expression))  # 去掉末尾的 ".01"

vst_expression <- subset(vst_expression, rownames(vst_expression) %in%rownames(survival_data))

#Integrating gene expression data and survival data:
combined_data <- cbind(vst_expression, survival_data)

# Use survival analysis libraries (e.g. Survival)
library(survival)

# Create survival model
survival_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ERBB2 , data = combined_data)

# Evaluate the survival model and view the summary information for the survival model
summary(survival_model)


# Model evaluation
# How well the model fits:
# The overall fit of the data can be fitted by the model, such as Likelihood Ratio Test or Concordance Index.
# Likelihood ratio test
lr_test <- anova(survival_model, test = "Chisq")
print(lr_test)
#3. Model Application:
# Predicted survival probability:
# Using the fitted Cox model, the survival probability of an individual or group can be predicted.

# Make a prediction for a dataset
predicted_survival <- survfit(survival_model, newdata = new_data)

# Draw a survival curve:
# Kaplan-Meier curve
plot(survfit(Surv(OS_MONTHS, OS_STATUS) ~ ERBB2, data = combined_data))

# Concordance Index
concordance_index <- concordance.index(survival_model)
print(concordance_index)


#2.	Use lasso cross validation to find the set of genes which predict survival the best.---------

#install.packages("glmnet")

library(glmnet)

# Set up survival time and event status

survival_object <- Surv(time = survival_data$OS_MONTHS, event = survival_data$OS_STATUS)

# Use cv.glmnet for lasso cross-validation
cv_model <- cv.glmnet(vst_expression, survival_object, family = "cox")

# View summary information of the cross-validation results, providing information about the performance of the model at each regularization parameter value.
print(cv_model)

# Obtain the selected set of genes
best_lambda <- cv_model$lambda.min
lasso_model <- glmnet(x = as.matrix(vst_expression), 
                      y = Surv(survival_data$OS_MONTHS,survival_data$OS_STATUS), 
                      family = "cox", lambda = best_lambda)
summary(lasso_model)
selected_genes <- rownames(coef(lasso_model, s = best_lambda))[coef(lasso_model, s = best_lambda) != 0]
