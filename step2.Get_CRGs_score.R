library(IOBR)

my_signature <- list(CRG_score = CRG_genes$Symbol)

ssgsea<-calculate_sig_score(eset            = exp_tmp, 
                            signature       = my_signature,
                            method          = "ssgsea",
                            mini_gene_count = 3)


ssgsea<-as.data.frame(ssgsea)
ssgsea <- column_to_rownames(ssgsea, var = "ID")
ssgsea <- ssgsea[rownames(meta_tmp),, drop = F]

ssgsea<-as.data.frame(t(ssgsea))
CRG_score<-as.matrix(ssgsea[1,])
CRG_score<-standarize.fun(CRG_score,halfwidth = 2)
CRG_score<-as.data.frame(t(CRG_score))

meta_tmp$CRG_score <- CRG_score$CRG_score

meta_tmp <- meta_tmp[order(meta_tmp$CRG_score),]
exp_tmp <- exp_tmp[,rownames(meta_tmp)]
