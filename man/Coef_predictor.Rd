% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state_prediction.R
\name{Coef_predictor}
\alias{Coef_predictor}
\title{Predict FL B cell States}
\usage{
Coef_predictor(Expression_df, nu_val = 0.1, batch_correct = FALSE)
}
\arguments{
\item{Expression_df}{A dataframe containing log2TPM+1 values with samples as columns and genes as rows. Rownames should either be HUGO symbols or
or ensembl gene ID's without decimals. Non-log scale will be autodetected and corrected.}

\item{nu_val}{The nu parameter for support vector regression. Default is nu = 0.1. Do not change unless you have good reason to.}

\item{batch_correct}{Logical. If TRUE, data is realigned before an additional SVR run. Use only if data comes from alternate platforms than RNAseq.
Default = FALSE.}
}
\value{
A data frame of predicted, normalized FL B cell state values.
}
\description{
Predict FL B cell States
}
\examples{
\donttest{
data("test_exp", envir = environment())
Coef_predictor(test_exp, nu = 0.1, batch_correct = FALSE)
}
}
