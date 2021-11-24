library(dplyr)
library("biomaRt") 

listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
# Organism: Mus musculus
datasets[,1] %>% {.[grep("musculus", .)]}
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)

# Platforms: [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array
AttAffy <- listAttributes(ensembl)[,1] %>% {.[grep("affy_mouse430_",.)]}
listAttributes(ensembl)[,1] %>% {.[grep("gene_name",.)]}

Identifiers <- c(AttAffy, "external_gene_name")

AffyToSymbol <- 
  getBM(attributes=Identifiers, 
        filters = AttAffy, 
        values = "1416519_at",
        mart = ensembl)

AffyToSymbol
