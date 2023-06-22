#' Aquaporin Dataset
#'
#' Dataset with molecular dynamics simulations for the yeast aquaporin (Aqy1) - the gated water channel of the yeast Pichi pastoris.
#' The dataset contains only the diameter Y of the channel which is used in the data analysis in (Klockmann and Krivobokova, 2023).
#' The diameter Y is measured by the distance between two centers of mass of certain residues of the protein.
#' The dataset includes a 100 nanosecond time frame, split into 20000 equidistant observations.
#' The full dataset, including the Euclidean coordinates of all 783 atoms, is available from the authors.
#' For more details see (Klockmann, Krivobokova; 2023).
#'
#' @format A data frame with 20000 rows and 1 variable:
#' \itemize{
#'     \item{\code{Y}:  }{the diameter of the channel}
#'}
#'@source   see (Klockmann, Krivobokova; 2023).
#'@examples
#'data(aquaporin)
"aquaporin"
