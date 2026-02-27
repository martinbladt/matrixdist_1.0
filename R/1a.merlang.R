#' Mixture of Erlang (continuous) constructor
#'
#' Builds a continuous phase-type model with block-diagonal Erlang components
#' (a.k.a. hyper-Erlang with fixed orders).
#'
#' @param block_sizes Integer vector with Erlang orders for each component.
#' @param probs Optional mixture weights. If omitted, random normalized weights
#'  are used.
#' @param rates Optional positive rates for each component. If omitted, random
#'  rates are used.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
merlang_ph <- function(block_sizes, probs = NULL, rates = NULL) {
  blocks <- .merlang_validate_blocks(block_sizes)
  m <- length(blocks)
  p <- sum(blocks)

  if (is.null(probs)) {
    probs <- stats::runif(m)
    probs <- probs / sum(probs)
  } else {
    if (length(probs) != m) stop("probs should have same length as block_sizes")
    if (any(probs < 0)) stop("probs should be non-negative")
    s <- sum(probs)
    if (s <= 0) stop("probs should sum to a positive number")
    probs <- probs / s
  }

  if (is.null(rates)) {
    rates <- stats::runif(m)
  } else {
    if (length(rates) != m) stop("rates should have same length as block_sizes")
    if (any(rates <= 0)) stop("rates should be positive")
  }

  alpha <- rep(0, p)
  S <- matrix(0, p, p)
  starts <- cumsum(c(1L, blocks))[seq_len(m)]

  for (j in seq_len(m)) {
    s <- starts[j]
    e <- s + blocks[j] - 1L
    alpha[s] <- probs[j]
    diag(S)[s:e] <- -rates[j]
    if (e > s) {
      S[cbind(s:(e - 1L), (s + 1L):e)] <- rates[j]
    }
  }

  obj <- ph(alpha = alpha, S = S)
  obj@name <- paste0("merlang ph(", p, ")")
  S_obj <- obj@pars$S
  attr(S_obj, "merlang_blocks") <- as.integer(blocks)
  obj@pars$S <- S_obj
  obj
}

#' Mixture of Erlang (discrete) constructor
#'
#' Builds a discrete phase-type model with block-diagonal Erlang components.
#'
#' @param block_sizes Integer vector with Erlang orders for each component.
#' @param probs Optional mixture weights. If omitted, random normalized weights
#'  are used.
#' @param q Optional forward probabilities in (0,1) for each component. If
#'  omitted, random probabilities are used.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
merlang_dph <- function(block_sizes, probs = NULL, q = NULL) {
  blocks <- .merlang_validate_blocks(block_sizes)
  m <- length(blocks)
  p <- sum(blocks)

  if (is.null(probs)) {
    probs <- stats::runif(m)
    probs <- probs / sum(probs)
  } else {
    if (length(probs) != m) stop("probs should have same length as block_sizes")
    if (any(probs < 0)) stop("probs should be non-negative")
    s <- sum(probs)
    if (s <= 0) stop("probs should sum to a positive number")
    probs <- probs / s
  }

  if (is.null(q)) {
    q <- stats::runif(m)
  } else {
    if (length(q) != m) stop("q should have same length as block_sizes")
    if (any(q <= 0 | q >= 1)) stop("q should be in (0, 1)")
  }

  alpha <- rep(0, p)
  S <- matrix(0, p, p)
  starts <- cumsum(c(1L, blocks))[seq_len(m)]

  for (j in seq_len(m)) {
    s <- starts[j]
    e <- s + blocks[j] - 1L
    alpha[s] <- probs[j]
    diag(S)[s:e] <- 1 - q[j]
    if (e > s) {
      S[cbind(s:(e - 1L), (s + 1L):e)] <- q[j]
    }
  }

  obj <- dph(alpha = alpha, S = S)
  obj@name <- paste0("merlang dph(", p, ")")
  S_obj <- obj@pars$S
  attr(S_obj, "merlang_blocks") <- as.integer(blocks)
  obj@pars$S <- S_obj
  obj
}

.merlang_validate_blocks <- function(block_sizes) {
  if (length(block_sizes) == 0) stop("block_sizes should be non-empty")
  blocks <- as.integer(block_sizes)
  if (any(is.na(blocks))) stop("block_sizes should be integer")
  if (any(blocks <= 0)) stop("block_sizes should be positive")
  blocks
}

.merlang_fit_blocks <- function(x, merlang_blocks = NULL) {
  blocks <- merlang_blocks
  if (is.null(blocks)) {
    blocks <- attr(x@pars$S, "merlang_blocks")
  }
  if (is.null(blocks)) {
    return(integer(0))
  }
  blocks <- .merlang_validate_blocks(blocks)
  p <- length(x@pars$alpha)
  if (sum(blocks) != p) {
    stop("sum(merlang_blocks) should equal model dimension")
  }
  blocks
}
