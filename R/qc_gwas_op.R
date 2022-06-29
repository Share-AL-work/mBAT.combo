#' Function for generating output of gwas after QC
#'
#' @param gwas
#' @param qc_snp_list
#'
#' @return
#' @export
#'
#' @examples
qc_gwas_op = function(gwas,qc_snp_list) {
  snplist = Reduce(intersect,list(gwas$SNP,qc_snp_list$SNP))
  gwas_qced = gwas[match(snplist, gwas$SNP),]
  return(gwas_qced)
}
