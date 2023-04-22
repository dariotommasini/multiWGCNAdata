### =========================================================================
### multiWGCNAdata metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c(paste0("Regional expression profiles of autistic brain"),
            paste0("Metadata for regional expression profiles of autistic brain"),
            paste0("Regional astrocyte expression profiles in experimental autoimmune encephalomyelitis"),
            paste0("Metadata for regional astrocyte expression profiles in experimental autoimmune encephalomyelitis"),
            paste0("Astrocyte weighted gene coexpression networks")),
  Description = c(paste0("The autism dataset from Voineagu et al. 2011, extracted from the Gene Expression Omnibus (accession number GSE28521)."),
                  paste0("The metadata for the autism dataset (Voineagu et al. 2011)."),
                  paste0("The astrocyte dataset from Itoh et al. PNAS. 2018, extracted from the supplementary materials."),
                  paste0("The metadata for the astrocyte Ribotag dataset (Itoh et al. PNAS. 2018)."),
                  paste0("The astrocyte networks from Tommasini and Fogel, BMC Bioinformatics, 2023. derived from the astrocyte Ribotag data from Itoh et al. PNAS. 2018.")),
  BiocVersion = "3.12",
  Genome = c(rep(NA, 2),
             rep("mm10", 3)),
  SourceType = rep("TXT", 5),
  SourceUrl = c(rep("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28521", 2),
                rep("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100329", 3)),
  SourceVersion = "Jan 28 2015",
  Species = c(rep("Homo sapiens", 2),
              rep("Mus musculus", 3)),
  TaxonomyId = c(rep(9606, 2),
                 rep(10090, 3)),
  Coordinate_1_based = NA,
  DataProvider = "GEO",
  Maintainer = "Dario Tommasini <dtommasini0@gmail.com>",
  RDataClass = "list",
  DispatchClass = c(rep("Rda",5)),
  RDataPath = paste0("multiWGCNAdata/",
                       c("astrocyte_data.rda",
                         "astrocyte_metadata.rda",
                         "astrocyte_networks.rda",
                         "autism_data.rda",
                         "autism_metadata.rda")),
  Tags = ""
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
