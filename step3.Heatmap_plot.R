library(ComplexHeatmap)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

## heatmap ------------------------------------------------------------------------------------
colormaps <- c(RColorBrewer::brewer.pal(name="Dark2", n = 8),RColorBrewer::brewer.pal(name="Paired", n = 12),RColorBrewer::brewer.pal(name="Set1", n = 9))
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")

# col_treatment
treatment_col <- colormaps[1:length(unique(meta_tmp$treatment))]
names(treatment_col) <- unique(meta_tmp$treatment)

# col_survival
surv_col <- c("#ED0000FF","grey")
names(surv_col) <- c("Alive", "Dead")

# col_response
response_col <- c("#EC7D21","#40548A")
names(response_col) <- c("RESPONDER", "NON RESPONDER")

# col_CRG_score
CRG_score_col <- colorRamp2(seq(-2, 2, length.out = 65), NMF:::ccRamp(x = heatmap.BlBkRd,n = 64))

color_all <- list(
  treatment_col = treatment_col,
  surv_col = surv_col,
  response_col = response_col,
  CRG_score_col = CRG_score_col
)

ha_immune = HeatmapAnnotation(Response = meta_tmp$response,
                              Survival = meta_tmp$Survival,
                              Treatment = meta_tmp$treatment,
                              CRG_score = meta_tmp$CRG_score,

                              annotation_legend_param=list(labels_gp = gpar(fontsize = 9),
                                                          title_gp = gpar(fontsize = 9, fontface = "bold"),
                                                          ncol=1),
                              
                              #gap=unit(c(1, rep(0, 4), 1, rep(0, 5)), "mm"),
                              col=list(Response = color_all$response_col,
                                      Survival = color_all$surv_col,
                                      Treatment = color_all$treatment_col,
                                      CRG_score = color_all$CRG_score_col
                              ),
                              #simple_anno_size = unit(1, "cm"), height = unit(6, "cm"),width = unit(15, "cm"),
                              show_annotation_name = TRUE,
                              annotation_name_side="right",
                              annotation_name_gp = gpar(fontsize = 15),
                              gap = unit(2, "mm"),
                              border = TRUE
)

## plot_data
exp_tmp <- exp_tmp[rownames(exp_tmp) %in% CRG_genes$Symbol,]
exp_tmp <- exp_tmp[rowSums(exp_tmp != 0) > 0, ]
plotdata <- exp_tmp
plotdata <- standarize.fun(plotdata, halfwidth = 2)
identical(rownames(meta_tmp), colnames(plotdata))

# plotdata <- plotdata[,rownames(id_match_final)]

hm <- ComplexHeatmap::pheatmap(plotdata,
                                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                                border_color = NA,
                                cluster_rows = T,
                                cluster_cols = T,
                                main = cancer_tmp,
                                row_title = "CRG-related genes",
                                #row_title_gp = gpar(fontsize = 15),
                                show_rownames = F,
                                show_colnames = F,
                                cellwidth = 10,
                                cellheight = 0.1, #5
                                top_annotation = ha_immune,
                                heatmap_legend_param = list(title = "CRG genes\n(Z-Score)",
                                                            title_gp = gpar(fontsize = 10,
                                                                            title_position = "topcenter",
                                                                            fontface = "bold")
                                )
)

pdf(file_name, width = 10, height = 12)
draw(hm, heatmap_legend_side = "right",
    annotation_legend_side = "bottom")
dev.off()
