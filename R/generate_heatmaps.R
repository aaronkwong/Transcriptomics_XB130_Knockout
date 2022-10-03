# R package management
library(renv)
renv::activate()

library("pheatmap")

# load expression and sample info
expr<-read.delim("./Data/2017-07-21, exported brain array RMA expression with annotation.tab",row.names=1)
exprs.eset.rma<-expr[,!colnames(expr) %in% c("entrez.id","gene.sym","chr","gene.name")]

# load the annotation data
load("./Data/2013-07-13, sample.info.rdata")
annotations<-expr[,c("entrez.id","gene.sym","chr","gene.name")]



# make rownames of sample.info the sample id
rownames(sample.info)<-sample.info$sample.ID

# clean up expression
datamat<-expr[,!colnames(expr) %in% c("entrez.id","gene.sym","chr","gene.name")]
if(!all(rownames(sample.info)==colnames(datamat))){
	print("Error. Expr and sample info columns arent in the right order.")
	quit()
}

#find and reorder for the groupings we want
sample_order<-c()
sample_order<-c(sample_order,which(sample.info$Genotype=="WT" & sample.info$Timepoint=="2 Weeks"))
sample_order<-c(sample_order,which(sample.info$Genotype=="KO" & sample.info$Timepoint=="2 Weeks"))
sample_order<-c(sample_order,which(sample.info$Genotype=="WT" & sample.info$Timepoint=="4 Weeks"))
sample_order<-c(sample_order,which(sample.info$Genotype=="KO" & sample.info$Timepoint=="4 Weeks"))
sample_order<-c(sample_order,which(sample.info$Genotype=="WT" & sample.info$Timepoint=="3 Months"))
sample_order<-c(sample_order,which(sample.info$Genotype=="KO" & sample.info$Timepoint=="3 Months"))

#rearrange the sample.info rows into the correct order
sample.info<-sample.info[sample_order,]
datamat<-datamat[,sample_order]

#read in differential expresion, and rank genes by fdr
x2_week_de<-read.delim("./results/limma_moderated_paired_t_test_2_Weeks.tab",row.names=1)
x2_week_de<-x2_week_de[order(x2_week_de$fdr,decreasing=FALSE),]
x4_week_de<-read.delim("./results/limma_moderated_paired_t_test_4_Weeks.tab",row.names=1)
x4_week_de<-x4_week_de[order(x4_week_de$fdr,decreasing=FALSE),]
x3_months_de<-read.delim("./results/limma_moderated_paired_t_test_3_Months.tab",row.names=1)
x3_months_de<-x3_months_de[order(x3_months_de$fdr,decreasing=FALSE),]

#create a text summary of the number of differentially expressed genes at cutoff FDR<0.05
sink("./results/de_gene_tabulation.txt")
fdr_cutoff<-0.05
print(paste0("Significant DE genes at fdr < ",fdr_cutoff))
print(paste0("Week 2 upregulated genes: ",sum(x2_week_de$fdr<fdr_cutoff & x2_week_de$fold.paired>0)))
print(paste0("Week 2 downregulated genes: ",sum(x2_week_de$fdr<fdr_cutoff & x2_week_de$fold.paired<0)))

print(paste0("Week 4 upregulated genes: ",sum(x4_week_de$fdr<fdr_cutoff & x4_week_de$fold.paired>0)))
print(paste0("Week 4 downregulated genes: ",sum(x4_week_de$fdr<fdr_cutoff & x4_week_de$fold.paired<0)))

print(paste0("Week 3 Months upregulated genes: ",sum(x3_months_de$fdr<fdr_cutoff & x3_months_de$fold.paired>0)))
print(paste0("Week 3 Months downregulated genes: ",sum(x3_months_de$fdr<fdr_cutoff & x3_months_de$fold.paired<0)))
sink()


row_fontsize<-15
png("./results/2_weeks_heatmap.png")
# get top 10 up and down regulated genes at 2 weeks between Mutant and WT
x2_week_de_pos<-x2_week_de[x2_week_de$fold.paired>0,]
x2_week_de_neg<-x2_week_de[x2_week_de$fold.paired<0,]
x2_week_heatmap_genes<-c(rownames(x2_week_de_pos)[1:10],rownames(x2_week_de_neg)[1:10])

# subset normalized expression for 2 weeks
datamat_2_weeks<-datamat[x2_week_heatmap_genes,sample.info$Timepoint=="2 Weeks"]
rownames(datamat_2_weeks)<-annotations[rownames(datamat_2_weeks),"gene.sym"]

# subset sample info for 2 weeks
x2_sample.info<-sample.info[sample.info$Timepoint=="2 Weeks",]
x2_sample.info$Genotype<-factor(x2_sample.info$Genotype)

#scale requires transpose, output from scale needs to be transposed back to original oritentation
pheatmap(t(scale(t(datamat_2_weeks))),cluster_rows=F,cluster_cols=F,annotation_col=x2_sample.info[,"Genotype",drop=F],main="2 Weeks",fontsize_row=row_fontsize)
dev.off()


png("./results/4_weeks_heatmap.png")
# get top 10 up and down regulated genes at 4 weeks between Mutant and WT
x4_week_de_pos<-x4_week_de[x4_week_de$fold.paired>0,]
x4_week_de_neg<-x4_week_de[x4_week_de$fold.paired<0,]
x4_week_heatmap_genes<-c(rownames(x4_week_de_pos)[1:10],rownames(x4_week_de_neg)[1:10])

# subset normalized expression for 4 weeks
datamat_4_weeks<-datamat[x4_week_heatmap_genes,sample.info$Timepoint=="4 Weeks"]
rownames(datamat_4_weeks)<-annotations[rownames(datamat_4_weeks),"gene.sym"]

# subset sample info for 4 weeks
x4_sample.info<-sample.info[sample.info$Timepoint=="4 Weeks",]
x4_sample.info$Genotype<-factor(x4_sample.info$Genotype)

#scale requires transpose, output from scale needs to be transposed back to original oritentation
pheatmap(t(scale(t(datamat_4_weeks))),cluster_rows=F,cluster_cols=F,annotation_col=x4_sample.info[,"Genotype",drop=F],main="4 Weeks",fontsize_row=row_fontsize)
dev.off()


png("./results/3_Months_heatmap.png")
# get top 10 up and down regulated genes at 3 Months between Mutant and WT
x3_Months_de_pos<-x3_months_de[x3_months_de$fold.paired>0,]
x3_Months_de_neg<-x3_months_de[x3_months_de$fold.paired<0,]
x3_Months_heatmap_genes<-c(rownames(x3_Months_de_pos)[1:10],rownames(x3_Months_de_neg)[1:10])

# subset normalized expression for 3 Months
datamat_3_Months<-datamat[x3_Months_heatmap_genes,sample.info$Timepoint=="3 Months"]
rownames(datamat_3_Months)<-annotations[rownames(datamat_3_Months),"gene.sym"]

# subset sample info for 3 Months
x3M_sample.info<-sample.info[sample.info$Timepoint=="3 Months",]
x3M_sample.info$Genotype<-factor(x3M_sample.info$Genotype)

#scale requires transpose, output from scale needs to be transposed back to original oritentation
pheatmap(t(scale(t(datamat_3_Months))),cluster_rows=F,cluster_cols=F,annotation_col=x3M_sample.info[,"Genotype",drop=F],main="3 Months",fontsize_row=row_fontsize)
dev.off()