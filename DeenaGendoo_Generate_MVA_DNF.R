# Deena M.A. Gendoo
# Created on January 10, 2017
# Generat3 the MVA-DNF matrix, using set of MVA-defined genes (genes identified by Dr. Linda Penn's lab)
########################################################################
########################################################################
source("Preprocessing/Logger.R")
source("Preprocessing/Dependencies.R")
source("Preprocessing/FunctionsBank.R")
load("Data/PrePros.RData")

bad_chars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

get_dnf <- function(gmt, gmt_file, idx,
                    pertData, sensData, strcData,
                    min_num_common_genes = 2, logger = get_logger("DNF.log")) {
    pathway_name <- gmt$geneset.names[idx]
    pathway_genes <- as.character(sort(gmt$genesets[[idx]]))

    pert_genes <- rownames(pertData)
    common_genes <- Reduce(intersect, list(pert_genes, pathway_genes))
    num_common_genes <- length(common_genes)

    log4r::info(logger, paste(
        "** Matching gene betweem `pertData` and ", "{", pathway_name, "} **\n",
        "\t - {", num_common_genes, "} genes matched.\n",
        "\t - Matched genes: {", toString(common_genes), "}"
    ))

    pass_sanity_check <- num_common_genes_sanity_check(
        num_common_genes = num_common_genes,
        pathway_name = pathway_name,
        min_num_common_genes = min_num_common_genes
    )

    if (!pass_sanity_check) {
        log4r::debug(
            logger = get_logger("DNF_unconsidered_pathway_min_gene_2_.log"),
            message = paste(
                "gmt: ", gmt_file, "[", idx, "]:", pathway_name,
                "matched genes", "{", num_common_genes, "}",
                "< min_num_common_genes {", min_num_common_genes, "}"
            )
        )
        return(NULL)
    }


    pertData <- pertData[common_genes, ]

    # Reduce All Matrices to lowest common set of drugs
    # Get 237 drugs now in the reduced sets
    pert_drugs <- sort(colnames(pertData))
    strc_drugs <- names(strcData)
    senes_drugs <- sort(colnames(sensData))
    commonDrugs <- Reduce(intersect, list(pert_drugs, strc_drugs, senes_drugs))

    pertData <- pertData[, commonDrugs]
    strcData <- strcData[commonDrugs]
    sensData <- sensData[, commonDrugs]

    # Sanity Checks
    drug_sanity_check(pertData, sensData, strcData)

    ## network layer construction and integration by SNF
    pertAffMat <- constPerturbationLayer(pertData)
    sensAffMat <- constSensitivityLayer(sensData)
    strcAffMat <- constStructureLayer(strcData)
    integrtStrctSensPert <- integrateStrctSensPert(
        sensAffMat, strcAffMat, pertAffMat
    )

    # Sanity Check - should all have the same dimensions: 237 X 237 --> 238...

    dimension_sanity_check(
        pert_aff_mat = pertAffMat,
        sens_aff_mat = sensAffMat,
        strc_aff_mat = strcAffMat,
        integrt_strct_sens_pert = integrtStrctSensPert
    )

    dnf <- list(
        "pathway" = pathway_name,
        "common_genes" = common_genes,
        "strc_layer" = strcAffMat,
        "sens_layer" = sensAffMat,
        "pert_layer" = pertAffMat,
        "integreated_network" = integrtStrctSensPert
    )
    return(dnf)
}

gmt_paths <- get_gmt_paths()
gmt_path <- "Data/GMT/c2.cp.kegg.v7.4.symbols.gmt"
gmt_file <- basename(gmt_path)
gmt <- read_gmt(gmt_path = gmt_path, logger = get_logger("DNF_gmt_reading.log"))

num_pathways <- length(gmt$geneset.names)
dnfs <- list()

for (idx in 1:num_pathways) {
    dnf <- get_dnf(
        gmt = gmt,
        gmt_file = gmt_file,
        idx = idx,
        pertData = pertData,
        sensData = sensData,
        strcData = strcData,
        min_num_common_genes = 2,
        logger = get_logger("DNF_loop.log")
    )
    if (is.null(dnf)) next()
    dnfs[[length(dnfs) + 1]] <- dnf
}

print(length(dnfs))
# Save an RData Object with all matrices: MVA-DNF and single layer taxonomies
# save(integrtStrctSensPert, strcAffMat, sensAffMat, pertAffMat, file = "MVA_DNF.RData")