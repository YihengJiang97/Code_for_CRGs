jyh_gs_pan_LR_forest = function (tpm, meta, geneset, geneset_alias, cancers, adjust = F) 
  {
    if (length(geneset) == 1) {
      if (geneset %in% names(msigb_genesets)) {
        genelist = msigb_genesets[geneset]
        names(genelist) = geneset_alias
      }
      else {
        stop("This geneset is not inluclded in the built in MSigDB geneset list.\n\n           Please use 'get_geneset()' function to get the geneset list,\n\n           or provide a character vector including user defined geneset.")
      }
    }
    else {
      genelist = data.frame(geneset_alias = rep(geneset_alias, 
        length(geneset)), geneset = geneset)
      genelist = split(genelist$geneset, genelist$geneset_alias)
    }
    if (adjust == T) {
      cox_results <- list()
      for (cancer in cancers) {
        exprSet = subset(tpm, Group == "Tumor" & Cancer == 
          cancer) %>% tibble::add_column(ID = stringr::str_sub(rownames(.), 
          1, 12), .before = "Cancer") %>% dplyr::filter(!duplicated(ID)) %>% 
          tibble::remove_rownames(.) %>% tibble::column_to_rownames("ID") %>% 
          dplyr::filter(rownames(.) %in% rownames(subset(meta, 
            Cancer == cancer)))
        exprSet = exprSet[, -(1:2)]
        exprSet = as.matrix(t(exprSet))
        gsvapar <- gsvaParam(exprData = exprSet, geneSets = genelist, 
          kcdf = "Gaussian")
        exprSet <- gsva(gsvapar)
        cl = meta[colnames(exprSet), ]
        cl$symbol = exprSet[geneset_alias, ]
        cl$symbol_bin <- ifelse(cl$symbol > median(cl$symbol), "High", "Low")
        cl$symbol_bin <- factor(cl$symbol_bin, levels = c("Low", "High"))
        # m = survival::coxph(survival::Surv(time, event) ~ symbol, data = cl)
        m = glm(flag ~ symbol_bin, data = cl, family = binomial)
        beta <- coef(m)
        se <- sqrt(diag(vcov(m)))
        HR <- exp(beta)
        HRse <- HR * se
        tmp <- round(cbind(coef = beta, se = se, z = beta/se, 
          p = 1 - pchisq((beta/se)^2, 1), HR = HR, HRse = HRse, 
          HRz = (HR - 1)/HRse, HRp = 1 - pchisq(((HR - 
            1)/HRse)^2, 1), HRCILL = exp(beta - qnorm(0.975, 
            0, 1) * se), HRCIUL = exp(beta + qnorm(0.975, 
            0, 1) * se)), 3)
        cox_results[[cancer]] = (tmp["symbol", ])
      }
      cox_results = do.call(rbind, cox_results)
      cox_results = as.data.frame(cox_results[, c(5, 9:10, 
        4)])
      cox_results <- cox_results[apply(cox_results, 1, function(x) all(is.finite(x))), ]
      np = paste0(cox_results$HR, " (", cox_results$HRCILL, 
        "-", cox_results$HRCIUL, ")")
      tabletext <- cbind(c("Cancer", rownames(cox_results)), 
        c("OR (95%CI)", np), c("P Value", cox_results$p))
      forestplot::forestplot(labeltext = tabletext, graph.pos = 3, 
        mean = c(NA, cox_results$HR), lower = c(NA, cox_results$HRCILL), 
        upper = c(NA, cox_results$HRCIUL), title = paste0("Odds Ratio Plot of ", 
          geneset_alias), hrzl_lines = list(`1` = grid::gpar(lwd = 2, 
          col = "black"), `2` = grid::gpar(lwd = 2, col = "black"), 
          `2` = grid::gpar(lwd = 2, col = "black")), 
        is.summary = c(TRUE, rep(FALSE, 33)), col = forestplot::fpColors(box = "#E64B35FF", 
          lines = "#E64B35FF", zero = "gray50"), zero = 1, 
        cex = 0.9, lineheight = "auto", colgap = unit(8, 
          "mm"), txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 1)), 
        boxsize = 0.5, ci.vertices = TRUE, ci.vertices.height = 0.3)
    }
    else {
      cox_results <- list()
      for (cancer in cancers) {
        exprSet = subset(tpm, Group == "Tumor" & Cancer == 
          cancer) %>% tibble::add_column(ID = stringr::str_sub(rownames(.), 
          1, 12), .before = "Cancer") %>% dplyr::filter(!duplicated(ID)) %>% 
          tibble::remove_rownames(.) %>% tibble::column_to_rownames("ID") %>% 
          dplyr::filter(rownames(.) %in% rownames(subset(meta, 
            Cancer == cancer)))
        exprSet = exprSet[, -(1:2)]
        exprSet = as.matrix(t(exprSet))
        gsvapar <- gsvaParam(exprData = exprSet, geneSets = genelist, 
          kcdf = "Gaussian")
        exprSet <- gsva(gsvapar)
        cl = meta[colnames(exprSet), ]
        cl$symbol = exprSet[geneset_alias, ]
        # m = survival::coxph(survival::Surv(time, event) ~ symbol, data = cl)
        m = glm(flag ~ symbol, data = cl, family = binomial)
        beta <- coef(m)
        se <- sqrt(diag(vcov(m)))
        HR <- exp(beta)
        HRse <- HR * se
        tmp <- round(cbind(coef = beta, se = se, z = beta/se, 
          p = 1 - pchisq((beta/se)^2, 1), HR = HR, HRse = HRse, 
          HRz = (HR - 1)/HRse, HRp = 1 - pchisq(((HR - 
            1)/HRse)^2, 1), HRCILL = exp(beta - qnorm(0.975, 
            0, 1) * se), HRCIUL = exp(beta + qnorm(0.975, 
            0, 1) * se)), 3)
        cox_results[[cancer]] = (tmp["symbol", ])
      }
      cox_results = do.call(rbind, cox_results)
      cox_results = as.data.frame(cox_results[, c(5, 9:10, 
        4)])
      cox_results <- cox_results[apply(cox_results, 1, function(x) all(is.finite(x))), ]
      np = paste0(cox_results$HR, " (", cox_results$HRCILL, 
        "-", cox_results$HRCIUL, ")")
      tabletext <- cbind(c("Cancer", rownames(cox_results)), 
        c("OR (95%CI)", np), c("P Value", cox_results$p))
      forestplot::forestplot(labeltext = tabletext, graph.pos = 3, 
        mean = c(NA, cox_results$HR), lower = c(NA, cox_results$HRCILL), 
        upper = c(NA, cox_results$HRCIUL), title = paste0("Odds Ratio Plot of ", 
          geneset_alias), hrzl_lines = list(`1` = grid::gpar(lwd = 2, 
          col = "black"), `2` = grid::gpar(lwd = 2, col = "black"), 
          `2` = grid::gpar(lwd = 2, col = "black")), 
        is.summary = c(TRUE, rep(FALSE, 33)), col = forestplot::fpColors(box = "#E64B35FF", 
          lines = "#E64B35FF", zero = "gray50"), zero = 1, 
        cex = 0.9, lineheight = "auto", colgap = unit(8, 
          "mm"), txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 1)), 
        boxsize = 0.5, ci.vertices = TRUE, ci.vertices.height = 0.3)
    }
  }
