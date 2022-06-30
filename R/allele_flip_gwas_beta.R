#' functions for allele flip and dif_AF1 >= 0.2
#'
#' @param gwas
#' @param info_LD
#'
#' @return
#' @export
#'
#' @examples
allele_flip_gwas_beta <-function(gwas,info_LD,no_AF) {
  #  PhysPos
  allele <- c("A1", "A2")
  gwas <- gwas %>% mutate_at(allele, as.character)
  gwas$AF1 <- as.numeric(gwas$AF1)
  info_LD <- info_LD %>% mutate_at(allele, as.character)
  info_LD$A1Freq <- as.numeric(info_LD$A1Freq)
  matched_id <- which(gwas$A1 == info_LD$A1 & gwas$A2 == info_LD$A2)

  gwas_tmp1 <- gwas[matched_id,]
  mismatched_id <- which(gwas$A1 == info_LD$A2 & gwas$A2 == info_LD$A1)
  gwas$beta[mismatched_id] <- -1 * gwas$beta[mismatched_id]
  gwas$AF1[mismatched_id] <- 1 - gwas$AF1[mismatched_id]
  gwas_tmp2 <- gwas[mismatched_id,]
  gwas_final <- rbind(gwas_tmp1,gwas_tmp2)
  #badsnps(abs(gwas_AF1-LD_A1Freq)>=0.2)
  if(no_AF){
    gwas_final <- inner_join(gwas_final,info_LD%>%select(ID,PhysPos,A1Freq)%>%dplyr::rename(SNP=ID))%>%arrange(PhysPos)
  } else {
    gwas_final <- inner_join(gwas_final,info_LD%>%select(ID,PhysPos,A1Freq)%>%dplyr::rename(SNP=ID))%>%arrange(PhysPos)
    gwas_final$dif_AF <- abs(gwas_final$AF1 - gwas_final$A1Freq)
    gwas_final <- gwas_final%>%filter(dif_AF<0.2)
  }
  return(gwas_final)
}

