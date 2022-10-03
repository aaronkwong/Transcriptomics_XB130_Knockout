#Microarray pre processing already conducted by Dr. Ricardo Zamel

# R package management
library(renv)
renv::activate()

library(limma)

#source utilities file
source("./R/differential_gene_expression_utilities.R")

#create directory to write results to
dir.create("results")

#read in sample info
# loads a "sample_info" object
load("./data/2013-07-13, sample.info.rdata")

#read in expression data 
# loads exprs.eset.rma object
exprs.eset.rma_with_anno<-read.table("./data/2017-07-21, exported brain array RMA expression with annotation.tab", sep="\t", header=TRUE,row.names=1)
exprs.eset.rma<-exprs.eset.rma_with_anno[,sample.info$sample.ID]

#create table of annotations for later on
annotations<-exprs.eset.rma_with_anno[,c("entrez.id","gene.sym","chr","gene.name")]

#check to make sure sample info and exprs are in the same order
if(!all(sample.info$sample.ID == colnames(exprs.eset.rma))){
	cat("error. sample info and columns of gene expression don't match.")
		quit()
}


# this experiment contains gene expression from XB130 deficient mice at various timepoints
# this function  compares mutant versus wildtype gene expression of all samples from a given timepoitn
# sample.info_raw is the sample info for all samples
# exprs is the expression for all samples
# timepoint is the timepoint for which we want the comparison to be made. Passing "all" will include all samples.
compare_mutatant_to_wildtype<-function(sample.info_raw,exprs,annotations,timepoint){
	# check condition if we want to compare samples from all timepoints
	if (timepoint!="all"){
		# subset sample info rows from the timepoint we want
		sample.info<-sample.info_raw[sample.info_raw$Timepoint==timepoint,]
		# subset gene expression from the timepoint we want
		eset<-exprs[,sample.info$sample.ID]
	}else{
		sample.info<-sample.info_raw
		eset<-exprs
	}

	# create design matrix
	groups<-factor(sample.info$Genotype)
	design<-model.matrix(~groups)
	# fit model
	fit<-lmFit(eset,design)
	# emmpirical bayes
	fit<-eBayes(fit)

	# extract p values
	p.value<-fit$p.value[,"groupsWT"]

	#Fold Change
	log2.pair<-vector(mode="numeric", length=nrow(eset))
	names(log2.pair)<- rownames(eset)
	for (i in 1:nrow(eset)){
		x<-eset[i,sample.info$Genotype=="WT"]
		y<-eset[i,sample.info$Genotype=="KO"]
		log2<-y-x
		log2.pair[i]<-mean(as.numeric(log2))
	}

	# summary output
	summ<-data.frame(
		p.value=p.value,
		fdr=p.adjust(p.value, method="fdr")
		)

	# calculatet fold changes
	summ$log2ratio.paired<-log2.pair
	summ$fold.paired<-log2ratio2fold(log2.pair)

	summ_annotated<-cbind(summ,annotations)

	# write the output file
	write.table(summ, file =paste0("./results/limma_moderated_paired_t_test_",gsub(" ","_",timepoint),".tab"),row.names=TRUE,sep="\t",col.names=NA)
	write.table(summ_annotated, file =paste0("./results/limma_moderated_paired_t_test_",gsub(" ","_",timepoint),"_annotated.tab"),row.names=TRUE,sep="\t",col.names=NA)
}

#compute differential gene expression 2 weeks
compare_mutatant_to_wildtype(sample.info_raw=sample.info,exprs=exprs.eset.rma,annotations=annotations,timepoint="2 Weeks")
#compute differential gene expression 4 weeks
compare_mutatant_to_wildtype(sample.info_raw=sample.info,exprs=exprs.eset.rma,annotations=annotations,timepoint="4 Weeks")
#compute differential gene expression 4 Months
compare_mutatant_to_wildtype(sample.info_raw=sample.info,exprs=exprs.eset.rma,annotations=annotations,timepoint="3 Months")
#compute differential gene expression all samples
compare_mutatant_to_wildtype(sample.info_raw=sample.info,exprs=exprs.eset.rma,annotations=annotations,timepoint="all")