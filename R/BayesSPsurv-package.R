#' BayesSPsurv: Bayesian Spatial Split-Population Survival Models
#'
#' Fits Bayesian parametric split-population (SP) survival models that
#' account for structural and spatial heterogeneity among "at risk" and
#' "immune" populations. Implements Weibull, Exponential, and Log-Logistic
#' forms for the duration component.
#'
#' @section Models:
#' Three estimation functions are provided:
#' \itemize{
#'   \item \code{\link{pooledSPsurv}}: SP survival model with no frailties.
#'   \item \code{\link{exchangeSPsurv}}: SP survival model with exchangeable (i.i.d.) frailties.
#'   \item \code{\link{spatialSPsurv}}: SP survival model with spatially correlated frailties (CAR prior).
#' }
#'
#' @section Data format:
#' All models require a \code{data.frame} with the following structure:
#' \describe{
#'   \item{Duration variable}{The dependent variable in the \code{duration} formula.
#'     Measures time from the start of the current spell to failure or censoring.}
#'   \item{Immunity indicator}{The dependent variable in the \code{immune} formula.
#'     A binary variable: 1 if the unit experienced the event (is "at risk"), 0 if immune/cured.}
#'   \item{Y0}{Elapsed time since inception until the beginning of the current time period (t-1).
#'     For the first observation of each unit, this is typically 0.}
#'   \item{LY}{Last-year indicator. Coded as 1 if this is the unit's last observation
#'     (due to censoring or failure), 0 otherwise.}
#' }
#'
#' For models with frailties (\code{exchangeSPsurv}, \code{spatialSPsurv}), an additional
#' column \code{S} is needed containing a spatial/unit ID (integer) for each observation.
#'
#' For the spatial model (\code{spatialSPsurv}), a square spatial weights matrix \code{A}
#' is also required, where rows/columns correspond to the unique values of \code{S}.
#'
#' @section Preparing data without spduration:
#' While the examples use \code{spduration::add_duration()} for convenience, you can
#' prepare the data manually. The key variables are:
#'
#' \preformatted{
#' # Example: manual data preparation
#' # Suppose 'mydata' has columns: unit_id, year, event (0/1), and covariates
#'
#' # Sort by unit and time
#' mydata <- mydata[order(mydata$unit_id, mydata$year), ]
#'
#' # Duration: time since start of current spell
#' # (count consecutive 0s in 'event' within each unit, reset on event)
#'
#' # Y0: duration at the start of the period (lag of duration, 0 for first obs)
#'
#' # LY: 1 if last observation for the unit, 0 otherwise
#'
#' # cured/atrisk: 1 if unit ever experienced event, 0 if always censored
#' }
#'
#' @section Spatial weights:
#' The \code{\link{spatial_SA}} function creates spatial weights matrices using
#' Gleditsch-Ward country codes and capital-city distances. For non-country data
#' (e.g., municipalities, districts), construct the matrix \code{A} directly:
#'
#' \preformatted{
#' # Example: custom spatial weights matrix
#' # A = 1 if units are neighbors, 0 otherwise
#' n_units <- length(unique(mydata$unit_id))
#' A <- matrix(0, n_units, n_units)
#' # Fill in adjacencies:
#' A[1, 2] <- A[2, 1] <- 1  # units 1 and 2 are neighbors
#' # etc.
#' }
#'
"_PACKAGE"
