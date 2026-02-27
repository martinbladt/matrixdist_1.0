% Manually maintained
\name{fit_multiple}
\alias{fit_multiple}
\title{Fit Multiple Competing Phase-Type Models}
\usage{
fit_multiple(
  models,
  y,
  weight = numeric(0),
  rcen = numeric(0),
  rcenweight = numeric(0),
  stepsEM = 1000,
  methods = c("RK", "RK"),
  rkstep = NA,
  uni_epsilon = NA,
  maxit = 100,
  reltol = 1e-08,
  optim_method = "BFGS",
  r = 1,
  burn = 0,
  erlang = NULL,
  merlang_blocks = NULL,
  plot_progress = TRUE,
  plot_every = 100,
  plot_breaks = 50,
  plot_pause = 0.01,
  parallel = TRUE,
  cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  colors = NULL
)
}
\arguments{
\item{models}{List of competing models inheriting from class \code{ph}.}

\item{y}{Vector of observations.}

\item{weight}{Optional observation weights.}

\item{rcen}{Optional right-censored observations (vector) or interval-censored observations (matrix).}

\item{rcenweight}{Optional weights for censored observations.}

\item{stepsEM}{Total EM iterations per model.}

\item{methods}{Methods for matrix exponential calculation.}

\item{rkstep}{Runge-Kutta step size.}

\item{uni_epsilon}{Epsilon parameter for uniformization.}

\item{maxit}{Maximum iterations in transformation-parameter optimization.}

\item{reltol}{Relative tolerance in transformation-parameter optimization.}

\item{optim_method}{Optimization method passed to \code{optim}.}

\item{r}{Sub-sampling proportion.}

\item{burn}{Burn percentage for histogram midpoint fitting (uncensored only).}

\item{erlang}{Optional logical flag for exact Erlang updates.}

\item{merlang_blocks}{Optional Erlang block sizes for exact mixed-Erlang updates.}

\item{plot_progress}{Logical. If \code{TRUE}, a two-panel progress display is updated.}

\item{plot_every}{Positive integer. EM chunk size per asynchronous worker update.}

\item{plot_breaks}{Number of histogram breaks for uncensored plotting.}

\item{plot_pause}{Pause in seconds after each plot refresh.}

\item{parallel}{Logical. If \code{TRUE}, uses asynchronous parallel chunks when available.}

\item{cores}{Number of worker processes in asynchronous mode.}

\item{colors}{Optional vector of colors for competing models.}
}
\value{
A list with fitted models and trajectories.
}
\description{
Runs several competing \code{ph}/\code{iph} models from user-provided starting values.
Progress is updated asynchronously when a model finishes its current EM chunk,
so faster models refresh earlier without waiting for slower models.
}
\examples{
\donttest{
set.seed(1)
truth <- ph(structure = "merlang", block_sizes = c(5, 5), probs = c(0.5, 0.5), rates = c(0.1, 0.02))
y <- sim(truth, n = 1500)

cands <- list(
  merlang = ph(structure = "merlang", block_sizes = c(5, 5)),
  gcoxian = ph(structure = "gcoxian", dimension = 10),
  general = ph(structure = "general", dimension = 10)
)

fit_multiple(cands, y, stepsEM = 120, plot_progress = TRUE, plot_every = 20)

set.seed(2)
truth2 <- ph(structure = "merlang", block_sizes = c(4, 4), probs = c(0.35, 0.65), rates = c(0.08, 0.02))
y2 <- sim(truth2, n = 1800)
cutoff <- stats::quantile(y2, 0.8)
y_obs <- y2[y2 <= cutoff]
y_rc <- y2[y2 > cutoff]

cands2 <- list(
  merlang = ph(structure = "merlang", block_sizes = c(4, 4)),
  erlang = ph(structure = "erlang", dimension = 8),
  coxian = ph(structure = "coxian", dimension = 8)
)

fit_multiple(cands2, y_obs, rcen = y_rc, stepsEM = 120, plot_progress = TRUE, plot_every = 20)
}
}
