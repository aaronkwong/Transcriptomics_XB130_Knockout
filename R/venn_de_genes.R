# R package management
library(renv)
renv::activate()

library("VennDiagram")

#load expression and sample info
expr<-read.delim("./Data/2017-07-21, exported brain array RMA expression with annotation.tab",row.names=1)
exprs.eset.rma<-expr[,!colnames(expr) %in% c("entrez.id","gene.sym","chr","gene.name")]

#load the sample info
load("./Data/2013-07-13, sample.info.rdata")
annotations<-expr[,c("entrez.id","gene.sym","chr","gene.name")]

#read in differential expresion, and rank genes by fdr
x2_week_de<-read.delim("./results/limma_moderated_paired_t_test_2_Weeks.tab",row.names=1)
x2_week_de<-x2_week_de[order(x2_week_de$fdr,decreasing=FALSE),]
x4_week_de<-read.delim("./results/limma_moderated_paired_t_test_4_Weeks.tab",row.names=1)
x4_week_de<-x4_week_de[order(x4_week_de$fdr,decreasing=FALSE),]
x3_months_de<-read.delim("./results/limma_moderated_paired_t_test_3_Months.tab",row.names=1)
x3_months_de<-x3_months_de[order(x3_months_de$fdr,decreasing=FALSE),]

# get list of differentially expressed genes
cutoff<-0.05
x2_week_genes<-rownames(x2_week_de)[x2_week_de$fdr<cutoff]
x4_week_genes<-rownames(x4_week_de)[x4_week_de$fdr<cutoff]
x3_months_genes<-rownames(x3_months_de)[x3_months_de$fdr<cutoff]

# Venn Diagram
venn.diagram(
  x = list(x2_week_genes, x4_week_genes, x3_months_genes),
  category.names = c("W2" , "W4" , "W12"),
  filename = './results/venn.png',
  output=TRUE,
  fill=rainbow(3),
  cex=4,
  cat.cex=0,
  imagetype="png"
)
