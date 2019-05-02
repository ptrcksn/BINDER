#----------------------------------------------------------------------------------------------------

#' Generate coexpression module for target candidate.
#'
#' @export
#' @param target_candidate Target candidate for which coexpression module is generated.
#' @param coexpression Numerical symmetric matrix providing gene coexpression values.
#' @param operon Character vector comprising genes of operon to which `target_candidate` belongs.
#' @param threshold Numerical cutoff point: any values above `threshold` are considered to form a coexpression module with `target_candidate`; defaults to "auto" which represents the 95th percentile of coexpression scores involving `target_candidate`.
#' @return Character vector containing genes sharing a coexpression with `target_candidate` greater than `threshold`.
#'
find_coexpression_module <- function(target_candidate, coexpression, operon, threshold="auto"){
  threshold <- (if(identical(threshold, "auto")) quantile(coexpression[target_candidate, ], 0.95, na.rm=TRUE) else threshold)
  coexpression_module <- setdiff(colnames(coexpression)[coexpression[target_candidate, ] > threshold], operon)
  return(coexpression_module)
}

#----------------------------------------------------------------------------------------------------


#' Compute CM value for target candidate.
#'
#' @export
#' @param target_candidate Target candidate for which CM value is computed.
#' @param target_candidates Character vector containing all target candidates for regulator.
#' @param ME_module Character vector comprising all target candidates with a high affinity for regulator motif.
#' @param coexpression Numerical symmetric matrix providing gene coexpression values.
#' @param operon Character vector comprising genes of operon to which `target_candidate` belongs.
#' @param threshold Numerical cutoff point: any values above `threshold` are considered to form a coexpression module with `target_candidate`; defaults to "auto" which represents the 95th percentile of coexpression scores involving `target_candidate`.
#' @return CM value for `target_candidate`.
#'
compute_CM <- function(target_candidate, target_candidates, ME_module, coexpression, operon, threshold="auto"){
  coexpression_module <- find_coexpression_module(target_candidate, coexpression, threshold)
  
  n_in_coexpression_module <- length(coexpression_module)
  n_not_in_coexpression_module <- length(target_candidates) - n_in_coexpression_module - sum(ifelse(operon %in% target_candidates, 1, 0))
  n_in_ME_module <- length(ME_module) - sum(ifelse(operon %in% ME_module, 1, 0))
  n_in_ME_module_in_coexpression_module <- length(intersect(ME_module, coexpression_module))
  
  p_value <- phyper(n_in_ME_module_in_coexpression_module, m=n_in_coexpression_module, n=n_not_in_coexpression_module, k=n_in_ME_module, lower.tail=FALSE)
  p_value <- if(n_in_coexpression_module == 0 | n_in_ME_module_in_coexpression_module == 0 | n_in_ME_module == 0) 1 else p_value
  CM <- 1 - p_value
  return(CM)
}

#' Compute CM values for all target candidates.
#'
#' @export
#' @param target_candidates Character vector containing all target candidates for regulator.
#' @param ME_module Character vector comprising all target candidates with a high affinity for regulator motif.
#' @param coexpression Numerical symmetric matrix providing gene coexpression values.
#' @param operons Data frame in DOOR format comprising operon annotations for all genes in `target_candidates`.
#' @param threshold Numerical cutoff point: any values above `threshold` are considered to form a coexpression module with `target_candidate`; defaults to "auto" which represents the 95th percentile of coexpression scores involving `target_candidate`.
#' @return Numerical vector containing CM values for `target_candidates`.
#'
compute_CMs <- function(target_candidates, ME_module, coexpression, operons, threshold="auto"){
  CMs <- as.numeric(sapply(target_candidates, function(x){operonID <- operons$OperonID[which(operons$Synonym == x)]; operon <- operons$Synonym[operons$OperonID == operonID]; compute_CM(x, target_candidates, ME_module, coexpression, operon, threshold)}))
  CMs <- 1-p.adjust((1-CMs), method="fdr") # FDR-based adjustment.
  CM_structure <- data.frame(target_candidate=target_candidates, CM=CMs)
  return(CM_structure)
}

#--------------------------------------------------

#' Compute CP value for target candidate.
#'
#' @export
#' @param target_candidate Target candidate for which CP value is computed.
#' @param target_candidates Character vector containing all target candidates for regulator.
#' @param PE_module Character vector comprising all target candidates with precedent orthologous interaction with regulator.
#' @param ortholog_module Character vector comprising all orthologous target candidates.
#' @param coexpression Numerical symmetric matrix providing gene coexpression values.
#' @param operon Character vector comprising genes of operon to which `target_candidate` belongs.
#' @param threshold Numerical cutoff point: any values above `threshold` are considered to form a coexpression module with `target_candidate`; defaults to "auto" which represents the 95th percentile of coexpression scores involving `target_candidate`.
#' @return CP value for `target_candidate`.
#'
compute_CP <- function(target_candidate, target_candidates, PE_module, ortholog_module, coexpression, operon, threshold="auto"){
  coexpression_module <- find_coexpression_module(target_candidate, coexpression, threshold)
  
  n_in_coexpression_module <- length(intersect(coexpression_module, ortholog_module)) #length(coexpression_module)
  n_not_in_coexpression_module <- length(ortholog_module) - n_in_coexpression_module - sum(ifelse(operon %in% ortholog_module, 1, 0)) #length(target_candidates) - n_in_coexpression_module - sum(ifelse(operon %in% target_candidates, 1, 0))
  n_in_PE_module <- length(PE_module) - sum(ifelse(operon %in% PE_module, 1, 0))
  n_in_PE_module_in_coexpression_module <- length(intersect(PE_module, coexpression_module))
  
  p_value <- phyper(n_in_PE_module_in_coexpression_module, m=n_in_coexpression_module, n=n_not_in_coexpression_module, k=n_in_PE_module, lower.tail=FALSE)
  p_value <- if(n_in_coexpression_module == 0 | n_in_PE_module_in_coexpression_module == 0 | n_in_PE_module == 0) 1 else p_value
  CP <- 1 - p_value
  return(CP)
}

#' Compute CP values for all target candidates.
#'
#' @export
#' @param target_candidates Character vector containing all target candidates for regulator.
#' @param PE_module Character vector comprising all target candidates with precedent orthologous interaction with regulator.
#' @param coexpression Numerical symmetric matrix providing gene coexpression values.
#' @param operons Data frame in DOOR format comprising operon annotations for all genes in `target_candidates`.
#' @param threshold Numerical cutoff point: any values above `threshold` are considered to form a coexpression module with `target_candidate`; defaults to "auto" which represents the 95th percentile of coexpression scores involving `target_candidate`.
#' @return Numerical vector containing CP values for `target_candidates`.
#'
compute_CPs <- function(target_candidates, PE_module, ortholog_module, coexpression, operons, threshold="auto"){
  CPs <- as.numeric(sapply(target_candidates, function(x){operonID <- operons$OperonID[which(operons$Synonym == x)]; operon <- operons$Synonym[operons$OperonID == operonID]; compute_CP(x, target_candidates, PE_module, ortholog_module, coexpression, operon, threshold)}))
  CPs <- 1-p.adjust((1-CPs), method="fdr") # FDR-based adjustment.
  CP_structure <- data.frame(target_candidate=target_candidates, CP=CPs)
  return(CP_structure)
}

#----------------------------------------------------------------------------------------------------

#' Derive symmetric coexpression matrix from expression matrix with genes represented on rows and experimental conditions represented on columns.
#'
#' @export
#' @param expression Numerical expression matrix with genes represented on rows and experimental conditions represented on columns.
#' @param target_candidates Character vector containing all target candidates for regulator.
#' @return Derived numerical symmetric coexpression matrix.
#'
compute_coexpression <- function(expression, target_candidates){
  target_candidates <- sort(intersect(target_candidates,rownames(expression)))
  expression <- expression[target_candidates, ]
  coexpression <- abs(cor(t(expression), use="pairwise.complete.obs"))
  return(coexpression)
}

#----------------------------------------------------------------------------------------------------
