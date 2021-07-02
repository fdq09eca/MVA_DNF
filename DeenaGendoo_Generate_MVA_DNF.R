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
pathway <- as.character(sort(gmt$genesets[[i]]))
mva <- read.csv("Data/MVA_genes.csv")
mvagenes <- as.character(sort(mva$Gene.name))
num_gene_match_with_pertData <- length(which(rownames(pertData) %in% pathway))

pertData <- pertData[(which(rownames(pertData) %in% mva$Gene.name)), ]
dim(pertData)
rownames(pertData)

# Reduce All Matrices to lowest common set of drugs
# Get 237 drugs now in the reduced sets
commonDrugs <- Reduce(intersect, list(
    sort(names(strcData)), sort(colnames(sensData)),
    sort(colnames(pertData))
))
strcData <- strcData[commonDrugs]
sensData <- sensData[, commonDrugs]
pertData <- pertData[, commonDrugs]

# Sanity Checks
if (ncol(sensData) != ncol(pertData)) stop(sprintf("error!"))
if (ncol(sensData) != length(strcData)) stop(sprintf("error!"))
if (all(colnames(pertData) != colnames(sensData))) stop(sprintf("error!"))
if (all(colnames(pertData) != names(strcData))) stop(sprintf("error!"))

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