#' Function for generating output of LD after QC
#'
#' @param gwas
#' @param LD
#' @param qc_snp_list
#'
#' @return
#' @export
#'
#' @examples
qc_LD_op=function(gwas,LD,qc_snp_list){
  snplist = Reduce(intersect,list(gwas$SNP,qc_snp_list$SNP))
  #get the new order_position(id) of refLD snps (re-ordered according to common snps in gwas summary data).
  LD_qced_snp_id = colnames(LD)[match(snplist,colnames(LD))]
  #LD re-ordered according to the order of common SNPs between gwas and refLD.
  LD_qced = LD[LD_qced_snp_id,LD_qced_snp_id,drop=FALSE]
  return(LD_qced)
}
