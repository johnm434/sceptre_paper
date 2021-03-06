% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_process_data_functions.R
\name{load_fbm}
\alias{load_fbm}
\title{Load FBM}
\usage{
load_fbm(fbm_metadata)
}
\arguments{
\item{fbm_metadata}{a "metadata" list returned by the create_fbm_from_mtx function.}
}
\value{
an FBM
}
\description{
Load FBM
}
\examples{
fbm_metadata <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/expression_FBM_metadata.rds")
load_fbm(fbm_metadata)
}
