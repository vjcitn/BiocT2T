
.osn_bucket_to_cache <- function(
    entity, folder = "BiocScviR",
    prefix = "https://mghp.osn.xsede.org/bir190004-bucket01/",
    ca = BiocFileCache::BiocFileCache()) {
  pa <- BiocFileCache::bfcquery(ca, entity)
  if (nrow(pa) > 1) {
    stop(sprintf(
      "%s has multiple instances in cache, please inspect.",
      entity
    ))
  } else if (nrow(pa) == 1) {
    return(pa$rpath)
  }
  target <- paste0(prefix, folder, "/", entity)
  tf <- tempfile(entity) # for metadata
  download.file(target, tf)
  BiocFileCache::bfcrpath(ca, tf, action = "copy")
}

#' obtain a summarized experiment with GEUVADIS (Coriell) RNA-seq quantifications
#' @examples
#' requireNamespace("SummarizedExperiment")
#' gg = getGeuvRNA()
#' gg
#' SummarizedExperiment::assay(gg[1:4,1:5])
#' @export
getGeuvRNA = function() {
 ca_entry = .osn_bucket_to_cache( "urecog.rda", folder = "BiocGeuvadis" )
 get(load(ca_entry))
}

#' install early version of T2T TxDb, do nothing if already available
#' @export
install_early_t2t_txdb = function() {
  if (!("TxDb.Hsapiens.NCBI.CHM13v2" %in% rownames(installed.packages())))
    install.packages(
     "https://mghp.osn.xsede.org/bir190004-bucket01/BiocGeuvadis/TxDb.Hsapiens.NCBI.CHM13v2_0.99.0.tar.gz",
       type="source", repos=NULL)
  else message("TxDb.Hsapiens.NCBI.CHM13v2 already installed.")
}


