# Deena M.A. Gendoo
# Created on January 10, 2017
# Generat3 the MVA-DNF matrix, using set of MVA-defined genes (genes identified by Dr. Linda Penn's lab)
########################################################################
########################################################################
source("Preprocessing/Dependencies.R")
source("Preprocessing/FunctionsBank.R")
source("Preprocessing/Logger.R")
load("Data/PrePros.RData")

DNF_logger <- get_logger("DNF.log")
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

gmt_path <- "Data/GMT/c2.cp.kegg.v7.4.symbols.gmt"

gmt <- read_gmt(gmt_path = gmt_path, logger = DNF_logger)

i <- 1
pathway_name <- gmt$geneset.names[i]
pathway_genes <- as.character(sort(gmt$genesets[[i]]))
# mva <- read.csv("Data/MVA_genes.csv")
# mvagenes <- as.character(sort(mva$Gene.name))

pert_genes <- rownames(pertData)
common_genes <- Reduce(intersect, list(pert_genes, pathway_genes))
num_common_genes <- length(common_genes)
pertData <- pertData[common_genes, ]

log4r::info(DNF_logger, paste(
    "{{", num_common_genes, "}} common genes between pertData and pathway {{",
    pathway_name, "}}: {{", common_genes, "}}"
))
# dim(pertData)
# rownames(pertData)

# Reduce All Matrices to lowest common set of drugs
# Get 237 drugs now in the reduced sets

pert_drugs <- if (is(pertData, "matrix")) sort(colnames(pertData)) else sort(names(pertData))
strc_drugs <- names(strcData)
senes_drugs <- sort(colnames(sensData))

commonDrugs <- Reduce(intersect, list(pert_drugs, strc_drugs, senes_drugs))

pertData <- if (is(pertData, "matrix")) pertData[, commonDrugs] else pertData[commonDrugs]
strcData <- strcData[commonDrugs]
sensData <- sensData[, commonDrugs]

# Sanity Checks
drug_sanity_check(pertData, sensData, strcData, DNF_logger)
# if (num_common_genes == 1) {
#     if (ncol(sensData) != length(pertData)) stop(sprintf("error!"))
#     if (all(names(pertData) != colnames(sensData))) stop(sprintf("error!"))
#     if (all(names(pertData) != names(strcData))) stop(sprintf("error!"))
# } else {
#     if (ncol(sensData) != ncol(pertData)) stop(sprintf("error!"))
#     if (all(colnames(pertData) != colnames(sensData))) stop(sprintf("error!"))
#     if (all(colnames(pertData) != names(strcData))) stop(sprintf("error!"))
# }
# if (ncol(sensData) != length(strcData)) stop(sprintf("error!"))


## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayer(sensData)
pertAffMat <- constPerturbationLayer(pertData)
integrtStrctSensPert <- integrateStrctSensPert(sensAffMat, strcAffMat, pertAffMat)

# Sanity Check - should all have the same dimensions: 237 X 237 --> 238...
dim(strcAffMat)
dim(sensAffMat)
dim(pertAffMat)
dim(integrtStrctSensPert)

# Save an RData Object with all matrices: MVA-DNF and single layer taxonomies
save(integrtStrctSensPert, strcAffMat, sensAffMat, pertAffMat, file = "MVA_DNF.RData")