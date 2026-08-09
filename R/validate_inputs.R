# Input validation for BayesSPsurv model functions

validate_inputs <- function(duration, immune, Y0, LY, data, N, burn, thin,
                            w, m, form, S = NULL, A = NULL,
                            prop.varV = NULL, prop.varW = NULL) {

    # Check data is a data.frame
    if (!is.data.frame(data)) {
        stop("'data' must be a data.frame.", call. = FALSE)
    }

    # Check N, burn, thin are valid
    if (!is.numeric(N) || length(N) != 1 || N < 1) {
        stop("'N' must be a positive integer (number of MCMC iterations).", call. = FALSE)
    }
    if (!is.numeric(burn) || length(burn) != 1 || burn < 0) {
        stop("'burn' must be a non-negative integer.", call. = FALSE)
    }
    if (burn >= N) {
        stop("'burn' (", burn, ") must be less than 'N' (", N, ").", call. = FALSE)
    }
    if (!is.numeric(thin) || length(thin) != 1 || thin < 1) {
        stop("'thin' must be a positive integer.", call. = FALSE)
    }
    if (thin >= (N - burn)) {
        stop("'thin' (", thin, ") must be less than N - burn (", N - burn, ").", call. = FALSE)
    }

    # Check Y0 and LY columns exist in data
    if (!Y0 %in% names(data)) {
        stop("Column '", Y0, "' (Y0) not found in data. Available columns: ",
             paste(names(data), collapse = ", "), call. = FALSE)
    }
    if (!LY %in% names(data)) {
        stop("Column '", LY, "' (LY) not found in data. Available columns: ",
             paste(names(data), collapse = ", "), call. = FALSE)
    }

    # Check S column if provided
    if (!is.null(S) && !S %in% names(data)) {
        stop("Column '", S, "' (S) not found in data. Available columns: ",
             paste(names(data), collapse = ", "), call. = FALSE)
    }

    # Check spatial matrix A dimensions
    if (!is.null(A) && !is.null(S)) {
        n_units <- length(unique(data[[S]]))
        if (nrow(A) != ncol(A)) {
            stop("Spatial weights matrix 'A' must be square. Got ",
                 nrow(A), " x ", ncol(A), ".", call. = FALSE)
        }
        if (nrow(A) != n_units) {
            stop("Spatial weights matrix 'A' has ", nrow(A),
                 " rows but there are ", n_units,
                 " unique spatial units in '", S, "'.", call. = FALSE)
        }
    }

    # Check w vector
    if (!is.numeric(w) || length(w) != 3) {
        stop("'w' must be a numeric vector of length 3, e.g. c(1, 1, 1).", call. = FALSE)
    }

    # Check proposal variances for spatial/frailty models
    if (!is.null(prop.varV) && (!is.numeric(prop.varV) || prop.varV <= 0)) {
        stop("'prop.varV' must be a positive number.", call. = FALSE)
    }
    if (!is.null(prop.varW) && (!is.numeric(prop.varW) || prop.varW <= 0)) {
        stop("'prop.varW' must be a positive number.", call. = FALSE)
    }

    invisible(TRUE)
}
