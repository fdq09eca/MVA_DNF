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

LOG_LEVEL <- "INFO"
MIN_NUM_COMMON_GENES <- 2

DNFs <- list()

gmt_paths <- get_gmt_paths()
num_skipped_pathway <- 0
for (gmt_path in gmt_paths) {
    gmt_file <- basename(gmt_path)
    gmt <- GSA.read.gmt(gmt_path)
    gmt <- read_gmt(
        gmt_path = gmt_path,
        logger = get_logger(
            log_file = "DNF_gmt_reading.log",
            log_lv = LOG_LEVEL
        )
    )
    for (idx in 1:length(gmt$geneset.names)) {
        pathway_name <- gmt$geneset.names[idx]
        pathway_genes <- as.character(sort(gmt$genesets[[idx]]))
        dnf <- get_dnf(
            pathway_name = pathway_name,
            pathway_genes = pathway_genes,
            pertData = pertData,
            sensData = sensData,
            strcData = strcData,
            min_num_common_genes = MIN_NUM_COMMON_GENES,
            logger = get_logger("DNF_loop.log", log_lv = LOG_LEVEL)
        )
        if (is.null(dnf)) {
            num_skipped_pathway <- num_skipped_pathway + 1
            log4r::info(
                logger = get_logger(
                    log_file = "DNF_unconsidered_pathway_min_gene_2_.log",
                    log_lv = LOG_LEVEL
                ),
                message = paste(
                    "gmt: ", gmt_file, "[", idx, "]:", pathway_name,
                    "< min_num_common_genes {", MIN_NUM_COMMON_GENES, "}"
                )
            )
            next()
        }
        DNFs[[length(DNFs) + 1]] <- dnf
    }
}

logger_DNFs_report(
    gmt_files = get_gmt_files(),
    DNFs = DNFs,
    logger = get_logger("DNFs_generation_report.log"),
    num_skipped_pathway = num_skipped_pathway
)
# Save an RData Object with all matrices: MVA-DNF and single layer taxonomies
# save(integrtStrctSensPert, strcAffMat, sensAffMat, pertAffMat, file = "MVA_DNF.RData")