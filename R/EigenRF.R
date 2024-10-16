library(parallel)
library(randomForest)

#' Title Metabolites Normalization Method:Eigen_RF
#'
#' @param peak the table of peak values, each column represents a sample, and each row represents a metabolite
#' @param groups the vector of groups
#' @param metabolites the table or vector of metabolites
#'
#' @return the normalized peak table
#' @export
#'
#' @importFrom parallel makeCluster detectCores parSapply stopCluster
#' @importFrom randomForest randomForest
#'
#' @examples Eigen_RF(peak, groups, metabolites)
#' \donttest {
#' Eigen_RF(peak, groups, metabolites)
#' }
rf_norm = function(peak, groups, metabolites) {
  cat(crayon::green("Getting bias with EigenMS......\n"))
  m_ints_eig1 = eig_norm1(m = peak, treatment = groups, prot.info = metabolites)
  Norm_EigenMS = eig_norm2(rv = m_ints_eig1)
  bias_EigenMS <- data.frame(Norm_EigenMS[["bias"]])
  bias_EigenMS <- data.frame(bias_EigenMS[,-1])
  row.names(bias_EigenMS) <- metabolites$metabolites
  bias_EigenMS <- data.frame(t(bias_EigenMS))
  rownames <- row.names(bias_EigenMS)
  bias_EigenMS <- data.frame(apply(bias_EigenMS, 2, as.numeric))
  row.names(bias_EigenMS) <- rownames
  cat(crayon::bgRed("Got bias!\n"))
  cat(crayon::green("Eigen_RF normalization......\n"))
  cl = makeCluster(detectCores() - 5)
  norm = parSapply(cl, X = 1:nrow(metabolites), function(j, peak, randomForest, groups, bias_EigenMS) {
    set.seed(42)
    df = data.frame(y = c(t(peak[j, ])), t(peak[-j, ]),groups)
    model = randomForest(y~., data = df,importance = F)
    newdata = data.frame(t(peak[-j, ]),groups)
    rf= data.frame(predict(model, newdata))
    new = peak[j, ] - bias_EigenMS[, j] + median(scale(rf, center = TRUE, scale = FALSE))
    return(new)
  }, peak, randomForest, groups, bias_EigenMS)
  norm = data.frame(apply(t(norm), 2, as.numeric))
  norm[is.na(norm)] <- 0
  for(i in 1:nrow(norm)) {
    norm[i, norm[i, ] < 0] = .5 * min(norm[i, norm[i, ] > 0])
  }
  row.names(norm) <- metabolites$metabolites
  cat(crayon::bgRed("All done!\n"))
  return(list(normalized = norm))
  stopCluster(cl)
}


