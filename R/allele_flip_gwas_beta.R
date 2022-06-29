#' functions for allele flip and dif_AF1 >= 0.2
#'
#' @param gwas
#' @param info_LD
#'
#' @return
#' @export
#'
#' @examples
allele_flip_gwas_beta <-function(gwas,info_LD) {
  #  PhysPos
  allele <- c("A1", "A2")
  gwas <- gwas %>% dplyr::mutate_at(allele, as.character)
  gwas$AF1 <- as.numeric(gwas$AF1)
  info_LD <- info_LD %>% dplyr::mutate_at(allele, as.character)
  info_LD$A1Freq <- as.numeric(info_LD$A1Freq)
  matched_id <- which(gwas$A1 == info_LD$A1 & gwas$A2 == info_LD$A2)

  gwas_tmp1 <- gwas[matched_id,]
  mismatched_id <- which(gwas$A1 == info_LD$A2 & gwas$A2 == info_LD$A1)
  gwas$beta[mismatched_id] <- -1 * gwas$beta[mismatched_id]
  gwas_tmp2 <- gwas[mismatched_id,]
  gwas_final <- rbind(gwas_tmp1,gwas_tmp2)
  gwas_final <- dplyr::inner_join(gwas_final,info_LD%>%dplyr::select(ID,PhysPos,A1Freq)%>%dplyr::rename(snp=ID))%>%dplyr::arrange(PhysPos)
  gwas_final$dif_AF <- abs(gwas_final$AF1 - gwas_final$A1Freq)
  gwas_final <- gwas_final%>%dplyr::filter(dif_AF<0.2)
  return(gwas_final)
}
