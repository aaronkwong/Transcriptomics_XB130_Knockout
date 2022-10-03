library(pca3d)

#define a plotting function
make_PCA_plot<-function(expression,sample.info,title,draw.legend){
	pca.result<-prcomp(t(expression),center = TRUE,scale=TRUE)

	#lets create a vector of colours to make WT blue and KO red
	p_colour<-gsub("WT","blue",sample.info$Genotype)
	p_colour<-gsub("KO","red",p_colour)

	#lets create a vector of shapes to make WT blue and KO red
	p_shape<-gsub("2 Weeks","sphere",sample.info$Timepoint)
	p_shape<-gsub("4 Weeks","tetrahedron",p_shape)
	p_shape<-gsub("3 Months","cube",p_shape)

	pca2d(pca.result,
		title=title,
		group=factor(sample.info$groupings),
		show.centroids = FALSE,
		show.plane = FALSE,
		legend=draw.legend,
		radius=2,
		col=p_colour,
		shape=p_shape
		)
	box(which = "plot", lty = "solid")
}


#load expression and sample info
load("./Data/2013-07-13, sample.info.rdata")
expr<-read.delim("./Data/2017-07-21, exported brain array RMA expression with annotation.tab",row.names=1)
exprs.eset.rma<-expr[,!colnames(expr) %in% c("entrez.id","gene.sym","chr","gene.name")]

#create new column in sample info for colouring points
sample.info$groupings<-paste0(sample.info$Genotype,sample.info$Timepoint)


#lets make sure the sample info match chip names match eset colnames
all(colnames(exprs.eset.rma)==sample.info$sampleID)

#lets make a PCA plot of only week 2
timepoint<-"2 Weeks"
exprs.eset.rma_2weeks<-exprs.eset.rma[,sample.info$Timepoint==timepoint]
sample.info_2weeks<-sample.info[sample.info$Timepoint==timepoint,]
all(colnames(exprs.eset.rma_2weeks)==sample.info_2weeks$sample.ID)
#TRUE

#lets make a PCA plot of only week 4
timepoint<-"4 Weeks"
exprs.eset.rma_4weeks<-exprs.eset.rma[,sample.info$Timepoint==timepoint]
sample.info_4weeks<-sample.info[sample.info$Timepoint==timepoint,]
all(colnames(exprs.eset.rma_4weeks)==sample.info_4weeks$sample.ID)
#TRUE

#lets make a PCA plot of only 3 months
timepoint<-"3 Months"
exprs.eset.rma_3months<-exprs.eset.rma[,sample.info$Timepoint==timepoint]
sample.info_3months<-sample.info[sample.info$Timepoint==timepoint,]
all(colnames(exprs.eset.rma_3months)==sample.info_3months$sample.ID)

#lets make plots without the legend
png("./results/XB130_pca_plot_all_no_legends.png")
make_PCA_plot(expression=exprs.eset.rma,sample.info=sample.info,title="All Groups: 3 timepoints x 2 groups",draw.legend=NULL)
dev.off()
png("./results/XB130_pca_plot_2W_no_legends.png")
make_PCA_plot(expression=exprs.eset.rma_2weeks,sample.info=sample.info_2weeks,title="2 Weeks KO and WT",draw.legend=NULL)
dev.off()
png("./results/XB130_pca_plot_4W_no_legends.png")
make_PCA_plot(expression=exprs.eset.rma_4weeks,sample.info=sample.info_4weeks,title="4 Weeks KO and WT",draw.legend=NULL)
dev.off()
png("./results/XB130_pca_plot_3M_no_legends.png")
make_PCA_plot(expression=exprs.eset.rma_3months,sample.info=sample.info_3months,title="3 Months KO and WT",draw.legend=NULL)
dev.off()

#lets make plots WITH the legend
png("./results/XB130_pca_plot_all_with_legends.png")
make_PCA_plot(expression=exprs.eset.rma,sample.info=sample.info,title="All Groups: 3 timepoints x 2 groups",draw.legend=TRUE)
dev.off()
png("./results/XB130_pca_plot_2W_with_legends.png")
make_PCA_plot(expression=exprs.eset.rma_2weeks,sample.info=sample.info_2weeks,title="2 Weeks KO and WT",draw.legend=TRUE)
dev.off()
png("./results/XB130_pca_plot_4W_with_legends.png")
make_PCA_plot(expression=exprs.eset.rma_4weeks,sample.info=sample.info_4weeks,title="4 Weeks KO and WT",draw.legend=TRUE)
dev.off()
png("./results/XB130_pca_plot_3M_with_legends.png")
make_PCA_plot(expression=exprs.eset.rma_3months,sample.info=sample.info_3months,title="3 Months KO and WT",draw.legend=TRUE)
dev.off()
