
#' mBAT-combo, combining P-values of mBAT and fastBAT to gain the highest power.
#'
#' @param bim_file
#' @param map_file
#' @param assoc_file
#' @param inLD_prefix
#' @param result_path
#' @param fastBAT_output
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'


mBAT_combo <- function(bim_file,
                       map_file,
                       assoc_file,
                       LD_path,
                       result_path,
                       fastBAT_output,
                       prop,
                       gene_annotate) {
  gene_annotate <- data.table::fread(map_file,h=F)
  names(gene_annotate) <- c("chr","start","end","gene")
  #Read in original GWAS summary statistics in gene region.
  gwas_ori <- data.table::fread(assoc_file)
  #gwas_ori = read.table(assoc_file, sep = "", header = T, na.strings = "NA", stringsAsFactors = F)
  names(gwas_ori) <- c("snp","A1","A2","AF1","beta","se","P","N")

  bim <- data.table::fread(bim_file)
  names(bim) <- c("chr","snp","bp","pos","A1","A2")
  snp2gene_ma <- snp2gene(
    bim,
    gene_annotate,
    window_start = 50,
    window_end = 50,
    only_sets = FALSE
  )$map

  dt <- dplyr::inner_join(bim%>%dplyr::mutate(chr = as.character(chr)),snp2gene_ma%>%dplyr::select("snp","chr","pos","gene"),by= c("chr", "snp", "pos"))%>%dplyr::select(snp,gene,chr)
  dt <- na.omit(dt)
  dt_ma <- dplyr::inner_join(dt,gwas_ori,by="snp")
  # read in fastBAT pruning results including gene names
  fastBAT_file <- data.table::fread(fastBAT_output)
  fastBAT_file$Chr <- as.numeric(fastBAT_file$Chr)
  fastBAT_file$Gene <- as.character(fastBAT_file$Gene)

  gene_list <- data.frame(Gene=unique(dt_ma$gene))
  #gene_list = fastBAT_file%>%filter(chr==22)%>%select(Gene)

  all_res=NULL
  for(i in 1:nrow(gene_list)){
    ##Read in LD file
    inLD_prefix <- paste0(LD_path,"/Chr",dt_ma$chr[1],"_",gene_list[i,],".txt.ldm.full")
    LD_in <- ReadLDMbin(inLD_prefix)
    info_in <- read.table(paste0(inLD_prefix,".info"), sep = "", header = T, na.strings = "NA", stringsAsFactors = F)
    colnames(LD_in) <- info_in$ID
    rownames(LD_in) <- info_in$ID

    # Match gwas snplist with refLD
    #inLD: inLD_prefix
    gene <- basename(inLD_prefix) #chr1_KANK4.txt.ldm.full
    #print(gene)
    name=substr(gene,1,(nchar(gene)-13))
    ID=as.character(matrix(unlist(strsplit(name,'_')),ncol=2)[,2])
    Chr=as.character(matrix(unlist(strsplit(name,'_')),ncol=2)[,1])
    chr_name=as.numeric(substr(Chr,4,nchar(Chr)))
    gwas_ori = dt_ma%>%filter(chr==chr_name,gene==ID)
    common_snp_id=Reduce(intersect, list(gwas_ori$snp, colnames(LD_in))) #find-common-elements-from-refLD and gwas
    LD_common_snp_id = colnames(LD_in)[match(common_snp_id,colnames(LD_in))] #get the new order_position(id) of refLD snps (re-ordered according to common snps in gwas summary data).

    info_check=info_in[match(LD_common_snp_id,info_in$ID),] ## get MAF information of final refLD.
    gwas_noqc = gwas_ori[match(common_snp_id, gwas_ori$snp),] #gwas re-ordered according to the order of common SNPs in gwas.
    #allele flip for gwas summary statistics
    gwas_noqc_final=allele_flip_gwas_beta(gwas_noqc,info_check,no_AF=FALSE)
    if(nrow(gwas_noqc_final)==0){
      print(paste0("No snp left in ", gene))
    }else{
      ###############   remove SNPs that are not matched after allele flip
      snplist_alf=Reduce(intersect, list(gwas_noqc_final$snp, colnames(LD_in)))
      common_snp_id_02=Reduce(intersect, list(snplist_alf, common_snp_id))
      LD_common_snp_id_02 = colnames(LD_in)[match(common_snp_id_02,colnames(LD_in))]
      LD = LD_in[LD_common_snp_id_02,LD_common_snp_id_02]

      z = as.matrix(gwas_noqc_final$beta/gwas_noqc_final$se)
      R = as.matrix(LD)
      eig = eigen(R, symmetric=TRUE)

      fastBAT_call=fastBAT_file%>%dplyr::filter(Chr==chr_name,Gene==ID)
      if(nrow(fastBAT_call)==0){
        print(paste0("Chr",chr_name,"_",ID," fastBAT has no result"))
      }else{
        res_pre_prop = mbat_SVD(gwas_noqc_final, eig, prop, fastBAT_call, z, chr_name, ID,gene_annotate)
        all_res=rbind(all_res,res_pre_prop)
      }
    }
  }

  setwd(result_path)
  print(getwd())
  write.table(all_res, file=paste0(result_path,"/mBAT_R_Chr",chr_name,".txt"), sep = "\t", row.names= F, col.names = T, quote = F)
  return(all_res)
}

