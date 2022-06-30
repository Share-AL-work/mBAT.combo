#' Map SNPs to genes function by Jaehyun Joo, with changes
#' Annotate SNPs onto their arbitrary genomic regions to perform set-based association tests.
#'
#' @param info_snp
#' @param info_gene
#' @param window_start
#' @param window_end
#' @param only_sets
#'
#' @return
#' @export
#'
#' @examples
snp2gene <- function(info_snp, info_gene,
                     window_start = 50L,
                     window_end = 50L,
                     only_sets=FALSE) {

  info_snp_name <- deparse(substitute(info_snp))

  df_check(info_snp)
  df_check(info_gene)
  columns_check(info_snp, c("snp", "chr", "pos"))
  columns_check(info_gene, c("gene", "chr", "start", "end"))
  nonneg_num_check(window_start, "window_start")
  nonneg_num_check(window_end, "window_end")

  if (anyDuplicated(info_snp$snp) > 0L) {
    stop("SNP IDs in ", "'", info_snp_name, "'", " must be unique.",
         call. = FALSE)
  }

  snpdat <- data.table(
    snp = as.character(info_snp$snp),
    chr = as.character(info_snp$chr),
    start = as.integer(info_snp$pos),
    end = as.integer(info_snp$pos),
    key = c("chr", "start", "end")
  )

  genedat <- data.table(
    gene = as.character(info_gene$gene),
    chr = as.character(info_gene$chr),
    start = as.integer(
      ifelse(info_gene$start - window_start * 1000 < 1L,
             1,
             info_gene$start - window_start * 1000)
    ),
    end = as.integer(info_gene$end + window_end * 1000),
    gene.start = as.integer(info_gene$start),
    gene.end = as.integer(info_gene$end),
    key = c("chr", "start", "end")
  )

  mapped <- foverlaps(snpdat, genedat, type = "within")

  snp_sets <- lapply(split(mapped[!is.na(gene), .(snp, gene)],
                           by = "gene", keep.by = FALSE),
                     function(x) unique(unname(unlist(x))))

  if (only_sets) {
    list(sets = snp_sets)
  } else {
    setnames(mapped,
             old = c("start", "end", "i.start"),
             new = c("gene.adj.start", "gene.adj.end", "pos"))

    vkeep <- c("snp", "chr", "pos",
               "gene", "gene.start", "gene.end",
               "gene.adj.start", "gene.adj.end")

    list(sets = snp_sets,
         map = setDF(setorder(mapped[, vkeep, with = FALSE],
                              chr, gene.start, gene.end, pos,
                              na.last = TRUE)))
  }

}
