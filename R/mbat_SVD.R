#' Function for mbat and mBAT-combo calculation with SVD
#'
#' @param gwas_noqc_final
#' @param eig
#' @param prop
#' @param fastBAT_call
#' @param z
#' @param chr_name
#' @param ID
#'
#' @return
#' @export
#'
#' @examples

mbat_SVD <- function(gwas_noqc_final, eig, prop, fastBAT_call, z, chr_name, ID, gene_annotate){
  posval <- eig$values[eig$values>1e-10]
  if (prop == 1) {
    k_prop <- length(posval)
  } else {
    k_prop <- which(cumsum(posval/sum(posval)) > as.numeric(prop))[1]
  }
  selected_prop <- 1:k_prop
  lambda_prop <- eig$values[selected_prop]
  U_prop <- eig$vectors[,selected_prop]
  UPz_prop <- crossprod(U_prop, z)
  chisq_prop <- crossprod(UPz_prop, 1/lambda_prop * UPz_prop)
  P_mbat_svd_prop <- pchisq(chisq_prop, df=k_prop, lower.tail=FALSE)
  res_mBAT <- data.frame(Gene = ID,
                        Chr = chr_name,
                        Start = as.numeric(gene_annotate%>%dplyr::filter(gene==ID,chr==chr_name)%>%dplyr::select(start)),
                        End = as.numeric(gene_annotate%>%dplyr::filter(gene==ID,chr==chr_name)%>%dplyr::select(end)),
                        No.SNPs = nrow(gwas_noqc_final),
                        SNP_start = gwas_noqc_final$snp[1],
                        SNP_end = gwas_noqc_final$snp[nrow(gwas_noqc_final)],
                        TopSNP = as.character(gwas_noqc_final%>%dplyr::filter(P==min(gwas_noqc_final$P))%>%dplyr::select(snp)),
                        TopSNP_Pvalue = min(gwas_noqc_final$P),
                        No.Eigenvalues = k_prop,
                        Chisq_mBAT = chisq_prop,
                        P_mBAT = P_mbat_svd_prop)
  P_fastbat <- fastBAT_call$Pvalue
  Chisq_fastBAT <- fastBAT_call$`Chisq(Obs)`
  P_combo <- ACATO(c(P_fastbat, P_mbat_svd_prop))
  res_mBAT$Chisq_fastBAT <- as.numeric(Chisq_fastBAT)
  res_mBAT$P_fastBAT <- as.numeric(P_fastbat)
  res_mBAT$P_mBATcombo <- as.numeric(P_combo)
  # reorder by column name
  res_mBAT <- res_mBAT[, c("Gene", "Chr", "Start", "End", "No.SNPs", "SNP_start", "SNP_end", "TopSNP", "TopSNP_Pvalue", "No.Eigenvalues", "Chisq_mBAT", "P_mBATcombo", "P_mBAT", "Chisq_fastBAT", "P_fastBAT")]
  return(res_mBAT)
}
