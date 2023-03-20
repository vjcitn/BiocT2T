library(GenomicFeatures)
mm = structure(list(name = c("method", "source", "Resource URL", "Genome"), value = c("Liftoff_v4", 
"JHU", "https://ccb.jhu.edu/T2T.shtml", "T2T-CHM13v2.0")), class = "data.frame", row.names = c(NA, 
-4L))
TxDb.Hsapiens.NCBI.CHM13v2 = makeTxDbFromGFF("chm13v2.0_RefSeq_Liftoff_v4.gff3", metadata=mm,
  organism="Homo sapiens", taxonomyId=9606)

makeTxDbPackage(TxDb.Hsapiens.NCBI.CHM13v2, version="0.99.0", person(given="Vincent", family="Carey", email="stvjc@channing.harvard.edu"), person(given="Vincent", family="Carey", email="stvjc@channing.harvard.edu"), pkgname="TxDb.Hsapiens.NCBI.CHM13v2")
