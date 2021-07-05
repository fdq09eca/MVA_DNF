# LIST OF FUNCTIONS
# Functions developed by Gendoo, Ghoraie, and El-Hachem for DNF publication (Cancer Research 2017)
source("Preprocessing/Logger.R")
##### constPert #####

constPerturbationLayer <- function(pertDat) {
  # Correlation for Perturbation
  pertCor <- cor(pertDat, method = "pearson", use = "pairwise.complete.obs")
  ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
  pertAff <- SNFtool::affinityMatrix(1 - pertCor, 20, 0.5)
  return(pertAff)
}

##### constSens #####

constSensitivityLayer <- function(sensDat) {

  # Correlation for Sensivity
  sensCor <- cor(sensDat, method = "pearson", use = "pairwise.complete.obs")
  ## if NA remaining in cor matrix, replace with 0s, not very clean but no other choices for now
  sensCor <- apply(sensCor, 1, function(x) ifelse(is.na(x), 0, x))
  ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
  sensAff <- SNFtool::affinityMatrix(1 - sensCor, 20, 0.5)

  return(sensAff)
}

##### constStruct #####

constStructureLayer <- function(targFps) {
  ## Correlation for Structure (Tanimoto metric)
  fpSim <- fingerprint::fp.sim.matrix(targFps, method = "tanimoto")
  rownames(fpSim) <- names(targFps)
  colnames(fpSim) <- names(targFps)
  ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
  fpAff <- SNFtool::affinityMatrix(1 - fpSim, 20, 0.5)
  return(fpAff)
}

##### Merger  #####

integrateStrctSensPert <- function(sensAff, strcAff, pertAff) {
  integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff))
  colnames(integration) <- rownames(integration) <- colnames(strcAff)
  return(integration)
}

read_gmt <- function(gmt_path, logger = get_logger("DNF.log")) {
  #' read_gmt_file with error checking
  curr_func <- match.call()[[1]]

  logger_init_msg(
    logger = logger, curr_func = curr_func,
    msg = paste("loading ", gmt_path)
  )

  gmt <- GSA.read.gmt(gmt_path)

  sanity_check(
    condition = length(gmt$genesets) == length(gmt$geneset.names),
    case = "length(gmt$genesets) != length(gmt$geneset.names)",
    logger = logger,
    caller = curr_func
  )

  log4r::info(logger, paste(
    "*** [", curr_func, "] {", basename(gmt_path), "} loaded successfully. ***\n",
    "\t- number of pathways: {", length(gmt$genesets), "}"
  ))
  return(gmt)
}


drug_sanity_check <- function(pertData, sensData, strcData,
                              logger = get_logger("DNF.log")) {
  #' Sanity check for data dimension after reducing to the common drugs.
  #'
  #' @description pertData, sensData, strcData should share the followings:
  #'  - a same number of drugs
  #'  - a same name of all drugs
  #' Errors will be thrown if any failure happens.
  #' @param pertData matrix. Perturbation data. a matrix with gene name rows and drug columns
  #' @param sensData matrix. Sensitivity data. a matrix with cell line rows and drug columns
  #' @param strcData list. Drug Strcture data, a list of fingerprint object.
  #' @param logger logger.
  #' @note This sanity check should be logged by ("DNF.log"), while optional.

  curr_func <- match.call()[[1]]
  logger_init_msg(logger = logger, curr_func = curr_func)

  sanity_check(
    condition = ncol(sensData) == ncol(pertData),
    case = "ncol(sensData) != ncol(pertData)",
    logger = logger,
    caller = curr_func
  )

  sanity_check(
    condition = all(colnames(pertData) == colnames(sensData)),
    case = "all(colnames(pertData) != colnames(sensData))",
    logger = logger,
    caller = curr_func
  )

  sanity_check(
    condition = all(colnames(pertData) == names(strcData)),
    case = "all(colnames(pertData) != names(strcData))",
    logger = logger,
    caller = curr_func
  )

  sanity_check(
    condition = ncol(sensData) == length(strcData),
    case = "ncol(sensData) != length(strcData)",
    logger = logger,
    caller = curr_func
  )

  logger_complete_msg(logger = logger, curr_func = curr_func)
}

get_gmt_files <- function(gmts_dir) {
  return(list.files(gmts_dir, pattern = "*.gmt$"))
}

get_gmt_paths <- function(gmts_dir) {
  files <- get_gmt_files(gmts_dir)
  paths <- list()
  for (i in 1:length(files)) {
    paths[i] <- file.path(gmts_dir, files[i])
  }
  return(paths)
}

num_common_genes_sanity_check <- function(num_common_genes, pathway_name,
                                          logger = get_logger("DNF.log"),
                                          min_num_common_genes = 2) {
  #' Sanity check for number of common gene between
  #' pertData and the pathway_genes.
  #'
  #' @description num_common_gene must be >= 2.
  #' Otherwise, there is no perason correlation.
  #' The pathway, therefore, will not be considered.
  #' @note Realistically, num_common_gene should be >= 5
  #' @note This sanity check should be logged by ("DNF.log"), while optional.

  curr_func <- match.call()[[1]]
  logger_init_msg(logger = logger, curr_func = curr_func)

  result <- sanity_check(
    condition = num_common_genes >= min_num_common_genes,
    case = paste("num_common_genes <", min_num_common_genes),
    logger = logger,
    caller = curr_func,
    warning = TRUE
  )

  if (result == FALSE) {
    log4r::warn(
      logger = logger,
      message = paste(
        ":::", "{", pathway_name, "} has NOT enough (", num_common_genes,
        ") common genes with `pertData` for consideration. :::"
      )
    )
    logger_complete_msg(logger = logger, curr_func = curr_func)
    return(FALSE)
  }

  log4r::info(
    logger = logger,
    message = paste(
      "{", pathway_name, "} has enough (", num_common_genes,
      ") common genes with `pertData` for consideration."
    )
  )

  logger_complete_msg(logger = logger, curr_func = curr_func)
  return(result)
}

dimension_sanity_check <- function(pert_aff_mat, sens_aff_mat, strc_aff_mat,
                                   integrt_strct_sens_pert,
                                   logger = get_logger("DNF.log")) {
  curr_func <- match.call()[[1]]
  logger_init_msg(logger = logger, curr_func = curr_func)

  sanity_check(
    condition = all(
      dim(pert_aff_mat) == dim(strc_aff_mat),
      dim(strc_aff_mat) == dim(sens_aff_mat),
      dim(sens_aff_mat) == dim(integrt_strct_sens_pert)
    ),
    case = "Layers do not share the same dimensions.",
    logger = logger,
    caller = curr_func
  )

  logger_complete_msg(logger = logger, curr_func = curr_func)
}

get_dnf <- function(pathway_name, pathway_genes,
                    pertData, sensData, strcData,
                    min_num_common_genes = 2, logger = get_logger("DNF.log", log_lv = "DEBUG")) {
  

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
    min_num_common_genes = min_num_common_genes,
    logger = logger
  )

  if (!pass_sanity_check) {
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
  drug_sanity_check(pertData, sensData, strcData, logger = logger)

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