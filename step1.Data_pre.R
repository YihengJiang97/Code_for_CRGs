for (i in cancer_types_input) {
  ### analysis
  cancer_tmp <- i
  meta_tmp <- get_meta_jyh(meta = meta_all_new, cancer = cancer_tmp)
  exp_tmp <- get_tpm_jyh(tpm = tpm_all_new, cancer = cancer_tmp)
  exp_tmp <- exp_tmp[which(exp_tmp$Group == "Tumor"),]

  refmap <- data.frame(
    patient_ID = substring(rownames(exp_tmp), 1, 12),
    sample_ID = rownames(exp_tmp)
  )

  # exp_tmp$sample_id <- rownames(exp_tmp)
  exp_tmp <- exp_tmp[,-c(1:2)]

  exp_tmp_uni <- anno_eset(
    exp_tmp,
    refmap,
    symbol = "patient_ID",
    probe = "sample_ID",
    method = "mean"
  )

  exp_tmp <- as.data.frame(t(exp_tmp_uni))

  meta_tmp$patient_ID <- rownames(meta_tmp)
  patient_keep <- intersect(meta_tmp$patient_ID, drug_all$patient_ID)
  meta_tmp <- meta_tmp[patient_keep,]
  exp_tmp <- exp_tmp[,patient_keep]

  meta_tmp <- left_join(meta_tmp, drug_all, by = "patient_ID")
  rownames(meta_tmp) <- meta_tmp$patient_ID
  meta_tmp$Survival <- ifelse(meta_tmp$event == 1, "Dead", "Alive")
}