
mBAT_combo <- function(bim_file,
                       map_file,
                       assoc_file,
                       inLD_prefix,
                       result_path,
                       fastBAT_output){
  gene_annotate=fread(map_file,h=F)
  names(gene_annotate)= c("chr","start","end","gene")
  #Read in original GWAS summary statistics in gene region.
  gwas_ori = fread(assoc_file)
  #gwas_ori = read.table(assoc_file, sep = "", header = T, na.strings = "NA", stringsAsFactors = F)
  names(gwas_ori)=c("snp","A1","A2","AF1","beta","SE","P","N")

  bim=fread(bim_file)
  names(bim)=c("chr","snp","bp","pos","A1","A2")
  snp2gene_ma = snp2gene(
    bim,
    gene_annotate,
    window_start = 50,
    window_end = 50,
    only_sets = FALSE
  )$map

  dt=inner_join(bim%>%mutate(chr = as.character(chr)),snp2gene_ma%>%select("snp","chr","pos","gene"))%>%select(snp,gene,chr)
  dt=na.omit(dt)
  dt_ma=inner_join(dt,gwas_ori,by="snp")
}
