#' @title spatial_SA
#' @description Generates a spatial weights matrix (A) and spatial unit ID (sp_id/S) for
#'     country-level data using Gleditsch-Ward (G&W) country codes and capital-city distances.
#'     This function only works with G&W country codes. For subnational data (municipalities,
#'     districts, etc.), construct the spatial weights matrix \code{A} manually.
#'     See \code{\link{BayesSPsurv-package}} for an example.
#'
#' @param data data.frame containing the observations.
#' @param var_ccode name of the variable containing Gleditsch-Ward country codes.
#' @param threshold distance in kilometers defining spatial adjacency. By default is 800.
#'
#' @return list. Contains database with variable sp_id (S) and matrix A.
#'
#' @examples
#'
#' walter <- spduration::add_duration(Walter_2015_JCR,"renewed_war",
#'                                    unitID = "ccode",
#'                                    tID = "year",
#'                                    freq = "year",
#'                                    ongoing = FALSE)
#'
#' walter <- spatial_SA(data = walter,
#'                      var_ccode = "ccode",
#'                      threshold = 800L)
#' @export
#'
spatial_SA <- function(data, var_ccode, threshold = 800L) {

    names(data)[which(names(data) == var_ccode)] <- 'ccode'
    ccodes    <- unique(data$ccode)
    # internal spatialSPsurv
    distcc    <- dplyr::distinct(capdist,
                                 numa,
                                 numb,
                                 kmdist,
                                 midist) ## Gleditsch and Ward Distance data
    distcc$km <- ifelse(distcc$kmdist > threshold, 0, 1)
    distcc_a  <- distcc[which(distcc$numa %in% ccodes),]
    distcc_b  <- distcc_a[which(distcc_a$numb %in% ccodes),]
    # matrix A
    matA      <- reshape2::acast(distcc_b, distcc_b$numa ~ distcc_b$numb, value.var = "km")
    diag(matA)<- 0L
    # S
    buildS    <- data.frame(
        ccode = rownames(matA),
        sp_id = 1:nrow(matA)
    )
    outdata   <- merge(data, buildS, by = "ccode")
    return(list(data_sp = outdata, matrixA = matA))
}

