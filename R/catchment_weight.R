#' @rdname catchment_ratio
#' @examples
#'
#' # gives weight only to costs under 2
#' catchment_weight(matrix(c(1, 2, 1, 3, 1, 2), 3), 2)
#' @return \code{catchment_weight}: A sparse matrix of weights.
#' @export

catchment_weight <- function(cost, weight = NULL, max_cost = NULL, adjust_zeros = 1e-6, scale = 2,
                             normalize_weight = FALSE, verbose = FALSE) {
  if (is.null(dim(cost))) cost <- matrix(cost, ncol = 1)
  if (is.data.frame(cost)) cost <- as.matrix(cost)
  if (!is.null(weight) && is.numeric(adjust_zeros)) {
    zero_costs <- !is.na(cost) & cost == 0
    if (any(zero_costs)) {
      if (verbose) cli_alert_info(paste("setting non-NA 0s to", adjust_zeros))
      cost[zero_costs] <- adjust_zeros
      rm(zero_costs)
    }
  }
  if (anyNA(cost)) cost[is.na(cost)] <- 0
  if (is.numeric(weight) && is.null(dim(weight)) && length(weight) > 1) {
    weight <- matrix(weight, ncol = ncol(cost))
  }
  if (is.null(weight) && !is.null(max_cost)) {
    weight <- max_cost
    max_cost <- NULL
  }
  if (is.null(weight)) {
    if (verbose) cli_alert_info("weight: cost")
    w <- as(cost, "CsparseMatrix")
  } else if (is.null(dim(weight))) {
    if (is.numeric(weight)) {
      # single buffer value means a uniformly weighted catchment area (original)
      if (verbose) cli_alert_info("weight: cost over {.field 0} and under {.field {weight[[1]]}}")
      w <- as(cost > 0 & cost < weight[[1]], "CsparseMatrix") * 1
    } else if (is.list(weight)) {
      # list of steps for roughly graded weightings (enhanced)
      if (verbose) cli_alert_info("weight: cost over {.field 0} and under steps of {.arg weight}")
      weight <- weight[order(-vapply(weight, "[[", 1, 1))]
      w <- as((cost <= weight[[1]][1]) * weight[[1]][2], "CsparseMatrix")
      for (s in weight[-1]) w[cost <= s[1]] <- s[2]
      w[cost <= 0] <- 0
    } else if (is.character(weight)) {
      # name of a weight function (kernel density)
      if (verbose) cli_alert_info("weight: {weight} weight function")
      weight <- tolower(weight[[1]])
      if ("units" %in% class(cost)) {
        cost <- matrix(
          as.numeric(cost), nrow(cost), ncol(cost),
          dimnames = dimnames(cost)
        )
      }
      if (grepl("^(?:gr|n)", weight)) {
        # gravity / normal kernel
        w <- as(sqrt(1 / cost^scale), "CsparseMatrix")
      } else if (grepl("^e", weight)) {
        # exponential kernel
        w <- as(exp(-cost * scale), "CsparseMatrix")
      } else if (grepl("^loga", weight)) {
        # logarithmic kernel
        w <- as(1 / (1 + log(cost, scale)), "CsparseMatrix")
      } else if (grepl("^li", weight)) {
        # linear kernel
        w <- as((scale - cost) / scale, "CsparseMatrix")
        w[w < 0] <- 0
      } else if (grepl("^l", weight)) {
        # logistic kernel
        w <- as(1 / (1 + exp(scale * cost)), "CsparseMatrix")
      } else if (grepl("^ga", weight)) {
        # Gaussian kernel
        w <- as(exp(-cost^2 / (2 * scale^2)), "CsparseMatrix")
      } else if (grepl("^d", weight) && exists(weight)) {
        # assumed to be a density function like dnorm
        w <- as(cost, "CsparseMatrix")
        w@x <- do.call(weight, list(w@x, 0, scale))
      } else if (grepl("^p", weight) && exists(weight)) {
        # assumed to be a distribution function like pnorm
        w <- as(cost, "CsparseMatrix")
        w@x <- do.call(weight, list(-w@x, 0, scale))
      } else {
        cli_abort("{.arg weight} not recognized")
      }
      w[cost <= 0 | !is.finite(w)] <- 0
    } else if (is.function(weight)) {
      if (verbose) cli_alert_info("weight: custom weight function")
      w <- as(weight(cost), "CsparseMatrix")
      w[cost <= 0 | !is.finite(w)] <- 0
    }
  } else {
    if (verbose) cli_alert_info("weight: {.arg weight}")
    w <- as(if (is.data.frame(weight)) as.matrix(weight) else weight, "CsparseMatrix")
  }
  if (!is.null(max_cost)) {
    if (verbose) cli_alert_info("zerod-out cost over {max_cost}")
    w[cost > max_cost] <- 0
  }
  if (normalize_weight) {
    wr <- rowSums(w)
    wr[wr == 0] <- 1
    w <- w * sweep(w, 1, wr, "/")
  }
  w
}
