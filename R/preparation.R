#----------------------------------------------------------------------------------------------------

#' Build augmented structure comprising CM and CP values from coexpression, ME and PE values.
#'
#' @export
#' @param ortholog_ME_PE_structure Data frame comprising columns: `regulator`, `target_candidate`, `ortholog_module_status` (1 if orthologous, 0 otherwise), `ME` (1 if high affinity for regulator motif, 0 otherwise) and `PE` (1 if orthologous interaction, 0 otherwise).
#' @param coexpression Numerical symmetric matrix providing gene coexpression values.
#' @param operons Data frame in DOOR format comprising operon annotations for all target candidates in `ortholog_ME_PE_structure`.
#' @param delta_CM Numerical cutoff point for coexpression module for CM: any values above `threshold` are considered to form a coexpression module with `target_candidate`; defaults to "auto" which represents the 95th percentile of coexpression scores involving `target_candidate`.
#' @param delta_CP Numerical cutoff point for coexpression module for CM: any values above `threshold` are considered to form a coexpression module with `target_candidate`; defaults to "auto" which represents the 95th percentile of coexpression scores involving `target_candidate`.
#' @return Data frame implementing augmented structure comprising columns: `regulator`, `target_candidate`, `ortholog_module_status` (1 if orthologous, 0 otherwise), `ME` (1 if high affinity for regulator motif, 0 otherwise), `PE` (1 if orthologous interaction, 0 otherwise), `CM` and `CP`.
#'
build_proxy_structure <- function(ortholog_ME_PE_structure, coexpression, operons, delta_CM="auto", delta_CP="auto"){
  regulator <- unique(ortholog_ME_PE_structure$regulator)
  target_candidates <- ortholog_ME_PE_structure$target_candidate
  
  ortholog_structure <- ortholog_ME_PE_structure[, c("target_candidate", "ortholog_module_status")]
  ME_structure <- ortholog_ME_PE_structure[, c("target_candidate", "ME")]
  PE_structure <- ortholog_ME_PE_structure[, c("target_candidate", "PE")]
  
  ortholog_module <- ortholog_ME_PE_structure$target_candidate[ortholog_ME_PE_structure$ortholog_module_status == 1]
  
  ME_module <- ortholog_ME_PE_structure$target_candidate[ortholog_ME_PE_structure$ME == 1]
  CM_structure <- compute_CMs(target_candidates, ME_module, coexpression, operons, delta_CM)
  
  PE_module <- ortholog_ME_PE_structure$target_candidate[ortholog_ME_PE_structure$PE == 1]
  CP_structure <- compute_CPs(target_candidates, PE_module, ortholog_module, coexpression, operons, delta_CP)
  
  proxy_structure <- plyr::join_all(list(ortholog_structure, ME_structure, PE_structure, CM_structure, CP_structure), by="target_candidate", type="full")
  proxy_structure$regulator <- regulator
  return(proxy_structure)
}

#----------------------------------------------------------------------------------------------------

#' Prepare data for use in BINDER model.
#'
#' @export
#' @param proxy_structure Data frame returned by `build_proxy_structure` implementing augmented structure comprising columns: `regulator`, `target_candidate`, `ortholog_module_status` (1 if orthologous, 0 otherwise), `ME` (1 if high affinity for regulator motif, 0 otherwise), `PE` (1 if orthologous interaction, 0 otherwise), `CM` and `CP`.
#' @return List comprising `regulator`, `target_candidate`, `N` (number of regulator-target observations), `K` (number of variables in the auxiliary stratum := 2 (i.e. ME and PE)), `M` (number of variables in the primary stratum := 2 (i.e. CM and CP)), `X` (matrix comprising the variables in the auxiliary stratum (ME, PE)), `Y` (matrix comprising the variables in the primary stratum (CM, CP)).
#'
prepare_data <- function(proxy_structure){
  regulator <- unique(proxy_structure$regulator)
  target_candidate <- proxy_structure$target_candidate
  
  X <- as.matrix(proxy_structure[, c("ME", "PE")])
  rownames(X) <- target_candidate
  Y <- as.matrix(proxy_structure[, c("CM", "CP")])
  Y <- ifelse(Y == 1, 0.9999, Y); Y <- ifelse(Y == 0, 0.0001, Y)
  rownames(Y) <- target_candidate
  
  N <- nrow(X)
  K <- ncol(X)
  M <- ncol(Y)
  
  prepared_data <- list(regulator=regulator, target_candidate=target_candidate, N=N, K=K, M=M, X=X, Y=Y)
  return(prepared_data)
}

#----------------------------------------------------------------------------------------------------
