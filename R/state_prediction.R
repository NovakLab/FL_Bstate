dat_dir <- "./data"

#Data alignment
align_data <- function(signature, exp_df){
  id_match <- apply(signature[,1:2], 2, function(x) sum(rownames(exp_df) %in% x))
  if(sum(id_match > 75) == 0) stop("Too few overlapping genes!!")

  gene_id_chosen <- names(which(id_match > 75))
  signature <- signature[,3:5] %>% set_rownames(signature[,gene_id_chosen])

  matching_genes <- intersect(rownames(signature), rownames(exp_df))

  if((lapply(list(1:50, 51:100, 101:150), function(x) which(match(matching_genes, rownames(signature)) %in% x) %>% length()) %>%
      unlist() %>% min()) < 25) stop("One or more states are unsuported by the overlapping genes!")

  return(signature[matching_genes,])
}


#This is the SVR itself
svr_coef <- function(sample_num, expression_mat, basis_mat, nu, norm = FALSE){
  #SVM
  model <- svm(x = basis_mat[,], y = expression_mat[rownames(basis_mat), sample_num],
               type = "nu-regression", nu = nu, kernel = "linear")
  coef <- (t(model$coefs) %*% model$SV)

  if(norm){
    coef[coef < 0] <- 0
    coef <- coef/sum(coef)
  }
  return(coef)
}



#SVR model
#' Predict FL B cell States
#' @import dplyr
#' @import magrittr
#' @import e1071
#' @importFrom stats model.matrix
#' @importFrom utils data
#'
#' @param Expression_df A dataframe containing log2TPM+1 values with samples as columns and genes as rows. Rownames should either be HUGO symbols or
#' or ensembl gene ID's without decimals. Non-log scale will be autodetected and corrected.
#' @param nu_val The nu parameter for support vector regression. Default is nu = 0.1. Do not change unless you have good reason to.
#' @param batch_correct Logical. If TRUE, data is realigned before an additional SVR run. Use only if data comes from alternate platforms than RNAseq.
#' Default = FALSE.
#'
#' @return A data frame of predicted, normalized FL B cell state values.
#' @export
#'
#' @examples
#' \donttest{
#' data("test_exp", envir = environment())
#' Coef_predictor(test_exp, nu = 0.1, batch_correct = FALSE)
#' }
Coef_predictor <- function(Expression_df, nu_val = 0.1, batch_correct = FALSE){
  data("signature_mat", envir = environment())

  #Get a signature matrix with the overlapping genes from expression set
  signature_mat <- align_data(signature_mat, Expression_df)

  #Is the data in the right space
  if(max(Expression_df) > 50){
    Expression_df <- log2(Expression_df[rownames(signature_mat),]+1)
  }else{
    Expression_df <- Expression_df[rownames(signature_mat),]
  }

  #Batch correction or no?
  if(batch_correct){
    if(!requireNamespace("sva", quietly = TRUE)) {
      stop("Package \"sva\" must be installed to use batch correction!", call. = FALSE)
    }

    coef_pred <- lapply(1:ncol(Expression_df), function(x) svr_coef(x, Expression_df, log2(signature_mat+1), nu = nu_val, norm = FALSE)) %>% do.call("rbind", .) %>% set_rownames(colnames(Expression_df))

    new_exp <- as.matrix(log2(signature_mat+1))%*%as.matrix(t(coef_pred))
    colnames(new_exp) <- paste0(colnames(new_exp), "_new")

    comb_exp <- cbind(Expression_df, new_exp)
    batch <- c(rep("1", ncol(Expression_df)), rep("2", ncol(Expression_df)))
    names(batch) <- colnames(comb_exp)

    group <- factor(rep(colnames(Expression_df), 2))
    mod_mat <- model.matrix(~ group)

    adjusted <- sva::ComBat(as.matrix(comb_exp), batch=batch, mod = mod_mat, par.prior = TRUE, mean.only = TRUE)
    new_exp <- adjusted[,colnames(Expression_df)]

    coef_pred <- lapply(1:ncol(new_exp), function(x) svr_coef(x, new_exp, log2(signature_mat+1), nu = nu_val, norm = TRUE)) %>% do.call("rbind", .) %>% set_rownames(colnames(new_exp))

    return(coef_pred)
  }else{
    coef_pred <- lapply(1:ncol(Expression_df), function(x) svr_coef(x, Expression_df, log2(signature_mat+1), nu = nu_val, norm = TRUE)) %>% do.call("rbind", .) %>% set_rownames(colnames(Expression_df))

    return(coef_pred)
  }
}


#Class Prediction
#' Predict FL B cell Group Assignments
#' @description Predicts the absolute class/group designation of samples based on their normalized state values.
#'
#' @import dplyr
#'
#' @param state_coef_df A dataframe containing B cell state coefficients, predicted or otherwise generated, which are organized into samples
#' as rows and coefficients as columns. The columns must have the column names of c("INFM", "PDZ", "CMI"), in that order.
#'
#' @return A named list of class assignments.
#' @export
#'
#' @examples
#' \donttest{
#' data("test_exp", envir = environment())
#' state_coef_df <- Coef_predictor(test_exp, nu = 0.1, batch_correct = FALSE)
#' predict_class(state_coef_df)
#' }
predict_class <- function(state_coef_df){
  if(any(!between(rowSums(state_coef_df), 0.99, 1.01))){
    stop("Coefficient numbers are not normalized to 1! These values will not work.", call. = FALSE)
  }
  if(ncol(state_coef_df) < 3) {
    stop("Input df is not in the correct orientation! Please check again!", call. = FALSE)
  }

  max_class <- apply(state_coef_df, 1, function(x) which.max(x[1:3])) %>%
    sub("1", "INFM", .) %>% sub("2", "PDZ", .) %>% sub("3", "CMI", .)

  new_class <-
    data.frame("Max" = max_class, state_coef_df[,1:3]) %>%
    apply(1, function(x) if(x[1] == "INFM" & x["CMI"] > 0.38){
      "CMI"
    }else if(x[1] == "PDZ" & x["INFM"] > 0.42){
      "INFM"
    }else if(x[1] == "CMI" & x["PDZ"] > 0.36){
      "PDZ"
    }else{
      x[1]
    })

  return(new_class)
}



