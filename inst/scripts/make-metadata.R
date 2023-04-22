### =========================================================================
### multiWGCNAdata metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c(paste0("Regional expression profiles of autistic brain"),
            paste0("Metadata for regional expression profiles of autistic brain"),
            paste0("Regional astrocyte expression profiles in experimental autoimmune encephalomyelitis"),
            paste0("Metadata for regional astrocyte expression profiles in experimental autoimmune encephalomyelitis"),
            paste0("Astrocyte weighted gene coexpression networks"),
            paste0("Timecourse expression profiles of tau pathology in entorhinal cortex"),
            paste0("Metadata for timecourse expression profiles of tau pathology in entorhinal cortex"),
            paste0("Tau weighted gene coexpression networks")),
  Description = c(paste0("The autism dataset from Voineagu et al. 2011, extracted from the Gene Expression Omnibus (accession number GSE28521)."),
                  paste0("The metadata for the autism dataset (Voineagu et al. 2011)."),
                  paste0("The astrocyte dataset from Itoh et al. PNAS. 2018, extracted from the supplementary materials."),
                  paste0("The metadata for the astrocyte Ribotag dataset (Itoh et al. PNAS. 2018)."),
                  paste0("The astrocyte networks from Tommasini and Fogel, BMC Bioinformatics, 2023. derived from the astrocyte Ribotag data from Itoh et al. PNAS. 2018."),
                  paste0("The tau dataset from Castanho et al. Cell Rep. 2020, extracted from GEO (GSE125957) processed data file."),
                  paste0("The metadata for the tau timecourse dataset (Castanho et al. Cell Rep. 2020)."),
                  paste0("The tau networks from Tommasini and Fogel, BMC Bioinformatics, 2023. derived from the tau timecourse data from Castanho et al. Cell Rep. 2020.")),
  BiocVersion = "3.12",
  Genome = c(rep(NA, 2),
             rep("mm10", 3)),
  SourceType = rep("TXT", 5),
  SourceUrl = c(rep("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28521", 2),
                rep("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100329", 3),
                rep("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125957", 3)),
  SourceVersion = "Jan 28 2015",
  Species = c(rep("Homo sapiens", 2),
              rep("Mus musculus", 6)),
  TaxonomyId = c(rep(9606, 2),
                 rep(10090, 6)),
  Coordinate_1_based = NA,
  DataProvider = "GEO",
  Maintainer = "Dario Tommasini <dtommasini0@gmail.com>",
  RDataClass = "list",
  DispatchClass = c(rep("Rda", 8)),
  RDataPath = paste0("multiWGCNAdata/",
                       c("autism_data.rda",
                         "autism_metadata.rda",
                         "astrocyte_data.rda",
                         "astrocyte_metadata.rda",
                         "astrocyte_networks.rda",
                         "tau_data.rda",
                         "tau_metadata.rda",
                         "tau_networks.rda")),
  Tags = ""
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
