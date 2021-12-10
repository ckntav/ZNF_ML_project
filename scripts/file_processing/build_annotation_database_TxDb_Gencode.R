setwd("/Users/chris/Desktop/ZNF_ML_project")

library(GenomicFeatures)

##### Build new annotation database v38
gff_file_v38 <- "input/gencode/gencode.v38.annotation.gff3"
gencode38_basic_chr_TxDb <- makeTxDbFromGFF(file = gff_file_v38, format = c("gff3"),
                                            dataSource = "Gencode_v38",
                                            organism = "Homo sapiens",
                                            taxonomyId = 9606)

saveDb(gencode38_basic_chr_TxDb, file = "input/gencode/txdb.gencode38BasicChr.sqlite")

##### Build new annotation database
# gff_file_v22 <- "input/gencode/gencode.v22.annotation.gff3"
# gencode22_basic_chr_TxDb <- makeTxDbFromGFF(file = gff_file_v22, format = c("gff3"),
#                                             dataSource = "Gencode_v22",
#                                             organism = "Homo sapiens",
#                                             taxonomyId = 9606)

# saveDb(gencode22_basic_chr_TxDb, file = "input/gencode/txdb.gencode22BasicChr.sqlite")

#####
# raw_promoters <- promoters(gencode22_basic_chr_TxDb, upstream = 3000, downstream = 3000, columns = "gene_id")
# ensg_without_version <- str_replace(unlist(raw_promoters$gene_id), "\\.[0-9]+$", "")
# raw_promoters$gene_id <- ensg_without_version
# names(raw_promoters) <- raw_promoters$gene_id
