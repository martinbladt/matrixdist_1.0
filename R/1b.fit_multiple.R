#' Fit Multiple Competing Phase-Type Models
#'
#' Runs several competing `ph`/`iph` models from user-provided starting values.
#' Progress is updated asynchronously when a model finishes its current EM chunk,
#' so faster models refresh earlier without waiting for slower models.
#'
#' @param models List of competing models (objects inheriting from `ph`).
#' @param y Vector of observations.
#' @param weight Optional observation weights.
#' @param rcen Optional right-censored observations (vector) or interval-censored
#'  observations (matrix).
#' @param rcenweight Optional weights for censored observations.
#' @param stepsEM Total EM iterations per model.
#' @param methods Methods for matrix exponential calculation.
#' @param rkstep Runge-Kutta step size.
#' @param uni_epsilon Epsilon parameter for uniformization.
#' @param maxit Maximum iterations in transformation-parameter optimization.
#' @param reltol Relative tolerance in transformation-parameter optimization.
#' @param optim_method Optimization method passed to `optim`.
#' @param r Sub-sampling proportion.
#' @param burn Burn percentage for histogram midpoint fitting (uncensored only).
#' @param erlang Optional logical flag for exact Erlang updates.
#' @param merlang_blocks Optional Erlang block sizes for exact mixed-Erlang updates.
#' @param plot_progress Logical. If `TRUE`, a two-panel progress display is updated.
#' @param plot_every Positive integer. EM chunk size per asynchronous worker update.
#' @param plot_breaks Number of histogram breaks for uncensored plotting.
#' @param plot_pause Pause in seconds after each plot refresh.
#' @param parallel Logical. If `TRUE`, uses asynchronous parallel chunks when available.
#' @param cores Number of worker processes in asynchronous mode.
#' @param colors Optional vector of colors for competing models.
#'
#' @return A list with fitted models and trajectories.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' truth <- ph(
#'   structure = "merlang",
#'   block_sizes = c(3, 4),
#'   probs = c(0.4, 0.6),
#'   rates = c(0.12, 0.03)
#' )
#' y <- sim(truth, n = 400)
#' cands <- list(
#'   merlang = ph(structure = "merlang", block_sizes = c(3, 4)),
#'   gcoxian = ph(structure = "gcoxian", dimension = 7),
#'   general = ph(structure = "general", dimension = 7)
#' )
#' fm <- fit_multiple(cands, y, stepsEM = 40, plot_progress = TRUE, plot_every = 10)
#' }
fit_multiple <- function(models,
                         y,
                         weight = numeric(0),
                         rcen = numeric(0),
                         rcenweight = numeric(0),
                         stepsEM = 1000,
                         methods = c("RK", "RK"),
                         rkstep = NA,
                         uni_epsilon = NA,
                         maxit = 100,
                         reltol = 1e-8,
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
                         colors = NULL) {
  methods <- .normalize_em_methods(methods)
  if (!is.list(models) || length(models) < 1) {
    stop("models should be a non-empty list")
  }
  if (!all(vapply(models, methods::is, logical(1), class2 = "ph"))) {
    stop("all elements in models should inherit from class ph")
  }
  if (length(plot_every) != 1 || !is.finite(plot_every) || plot_every < 1) {
    stop("plot_every should be a positive integer")
  }
  if (length(plot_breaks) != 1 || !is.finite(plot_breaks) || plot_breaks < 1) {
    stop("plot_breaks should be a positive integer")
  }
  if (length(plot_pause) != 1 || !is.finite(plot_pause) || plot_pause < 0) {
    stop("plot_pause should be a non-negative number")
  }
  if (length(stepsEM) != 1 || !is.finite(stepsEM) || stepsEM < 1) {
    stop("stepsEM should be a positive integer")
  }
  if (length(cores) != 1) {
    stop("cores should be a positive integer")
  }
  if (!is.finite(cores) || cores < 1) {
    cores <- 1
  }
  cores <- as.integer(cores)

  n_models <- length(models)
  model_names <- names(models)
  if (is.null(model_names) || any(model_names == "")) {
    model_names <- paste0("model_", seq_len(n_models))
  }

  if (is.null(colors)) {
    base_colors <- c(
      "#0B6E4F", "#C44536", "#2D7DD2", "#F4A259", "#7D4E57",
      "#1D4E89", "#3CAEA3", "#9A031E", "#5F0F40", "#8D99AE"
    )
    colors <- rep(base_colors, length.out = n_models)
  }
  if (length(colors) < n_models) {
    stop("colors should have length at least length(models)")
  }
  colors <- colors[seq_len(n_models)]

  stepsEM <- as.integer(stepsEM)
  chunk_size <- as.integer(plot_every)
  ll_start <- max(1L, as.integer(plot_every))
  plot_breaks <- as.integer(plot_breaks)
  plot_pause <- as.numeric(plot_pause)

  right_censored <- is.vector(rcen) && length(rcen) > 0
  km_plot <- NULL
  y_plot <- NULL
  x_data <- NULL
  if (isTRUE(plot_progress)) {
    if (!is.vector(y) || is.matrix(y) || !is.numeric(y)) {
      stop("plot_progress requires y to be a numeric vector")
    }
    y_plot <- as.numeric(y)
    y_plot <- y_plot[is.finite(y_plot)]
    if (length(y_plot) < 2) {
      stop("plot_progress requires at least two finite observations")
    }

    if (right_censored) {
      rcen_plot <- as.numeric(rcen)
      rcen_plot <- rcen_plot[is.finite(rcen_plot)]
      w_plot <- if (length(weight) == 0) rep(1, length(y_plot)) else as.numeric(weight)
      rcw_plot <- if (length(rcenweight) == 0) rep(1, length(rcen_plot)) else as.numeric(rcenweight)
      if (length(w_plot) != length(y_plot) || length(rcw_plot) != length(rcen_plot)) {
        stop("incompatible plotting inputs for right-censored data")
      }
      km_raw <- data.frame(
        time = c(y_plot, rcen_plot),
        weight = c(w_plot, rcw_plot),
        event = c(rep(1, length(y_plot)), rep(0, length(rcen_plot)))
      )
      km_total <- stats::aggregate(weight ~ time, data = km_raw, FUN = sum)
      km_event <- stats::aggregate(weight ~ time, data = km_raw[km_raw$event == 1, , drop = FALSE], FUN = sum)
      event_weight <- numeric(nrow(km_total))
      if (nrow(km_event) > 0) {
        event_weight[match(km_event$time, km_total$time)] <- km_event$weight
      }
      risk <- sum(km_total$weight)
      surv <- numeric(nrow(km_total))
      s_now <- 1
      for (i in seq_len(nrow(km_total))) {
        if (risk > 0 && event_weight[i] > 0) {
          s_now <- s_now * (1 - event_weight[i] / risk)
        }
        surv[i] <- s_now
        risk <- risk - km_total$weight[i]
      }
      km_plot <- list(
        time = c(0, km_total$time),
        surv = c(1, pmax(0, pmin(1, surv)))
      )
      x_data <- c(y_plot, rcen_plot)
    } else {
      x_data <- y_plot
    }
  }

  if (isTRUE(plot_progress)) {
    x_left <- min(x_data, na.rm = TRUE)
    x_right <- max(x_data, na.rm = TRUE)
    if (!is.finite(x_left) || !is.finite(x_right)) {
      stop("plot_progress requires finite plotting range")
    }
    if (right_censored) {
      x_left <- 0
    }
    if (x_right <= x_left) {
      pad <- max(abs(x_left) * 0.05, 1e-8)
      x_left <- x_left - pad
      x_right <- x_right + pad
    }
    grid_plot <- seq(x_left, x_right, length.out = 300)
  } else {
    x_left <- 0
    x_right <- 1
    grid_plot <- 0
  }

  states <- models
  done_steps <- integer(n_models)
  remaining <- rep(stepsEM, n_models)
  track <- matrix(NA_real_, nrow = n_models, ncol = stepsEM)
  can_rewrite <- isTRUE(interactive()) &&
    isTRUE(tryCatch(base::isatty(stdout()), error = function(e) FALSE)) &&
    sink.number(type = "output") == 0L &&
    sink.number(type = "message") == 0L &&
    nzchar(Sys.getenv("TERM")) &&
    Sys.getenv("TERM") != "dumb"
  status_initialized <- FALSE
  status_line_count <- n_models + 1L
  status_width <- 0L
  plot_state <- function() {
    if (!isTRUE(plot_progress)) {
      return(invisible(NULL))
    }
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = c(1, 2))

    if (right_censored) {
      graphics::plot(
        km_plot$time,
        km_plot$surv,
        type = "s",
        col = "black",
        xlim = c(x_left, x_right),
        ylim = c(0, 1),
        xlab = "",
        ylab = "",
        main = "",
        xaxs = "i"
      )
      for (i in seq_len(n_models)) {
        if (done_steps[i] < 1) {
          next
        }
        surv_fit <- tryCatch(
          cdf(states[[i]], grid_plot, lower.tail = FALSE),
          error = function(e) rep(NA_real_, length(grid_plot))
        )
        if (all(is.finite(surv_fit))) {
          graphics::lines(grid_plot, surv_fit, col = colors[i], lwd = 2)
        }
      }
    } else {
      hobj <- graphics::hist(
        y_plot,
        breaks = plot_breaks,
        plot = FALSE
      )
      dfit_list <- vector("list", n_models)
      y_top <- max(hobj$density, na.rm = TRUE)
      for (i in seq_len(n_models)) {
        if (done_steps[i] < 1) {
          next
        }
        dfit <- tryCatch(
          dens(states[[i]], grid_plot),
          error = function(e) rep(NA_real_, length(grid_plot))
        )
        if (all(is.finite(dfit))) {
          dfit_list[[i]] <- dfit
          y_top <- max(y_top, dfit, na.rm = TRUE)
        }
      }
      if (!is.finite(y_top) || y_top <= 0) {
        y_top <- 1
      }
      y_top <- 1.05 * y_top
      graphics::hist(
        y_plot,
        probability = TRUE,
        breaks = plot_breaks,
        col = "grey90",
        border = "white",
        ylim = c(0, y_top),
        xlab = "",
        ylab = "",
        main = ""
      )
      for (i in seq_len(n_models)) {
        if (!is.null(dfit_list[[i]])) {
          graphics::lines(grid_plot, dfit_list[[i]], col = colors[i], lwd = 2)
        }
      }
    }
    if (n_models <= 8) {
      graphics::legend("topright", legend = model_names, col = colors, lwd = 2, bty = "n", cex = 0.8)
    }

    if (ncol(track) >= ll_start) {
      finite_track <- track[, ll_start:ncol(track), drop = FALSE]
      finite_track <- finite_track[is.finite(finite_track)]
    } else {
      finite_track <- numeric(0)
    }
    if (length(finite_track) == 0) {
      graphics::plot.new()
      graphics::box()
    } else {
      yr <- range(finite_track)
      if (!is.finite(yr[1]) || !is.finite(yr[2]) || yr[1] == yr[2]) {
        yr <- yr + c(-1, 1)
      }
      graphics::plot(
        c(ll_start, max(ll_start, stepsEM)),
        yr,
        type = "n",
        xlab = "",
        ylab = "",
        main = ""
      )
      for (i in seq_len(n_models)) {
        if (done_steps[i] < ll_start) {
          next
        }
        idx <- seq.int(ll_start, done_steps[i])
        graphics::lines(idx, track[i, idx], col = colors[i], lwd = 2)
      }
    }

    utils::flush.console()
    if (plot_pause > 0) {
      Sys.sleep(plot_pause)
    }
  }

  fit_chunk <- function(model, steps_chunk) {
    fit_obj <- NULL
    utils::capture.output({
      fit_obj <- fit(
        model,
        y,
        weight = weight,
        rcen = rcen,
        rcenweight = rcenweight,
        stepsEM = steps_chunk,
        methods = methods,
        rkstep = rkstep,
        uni_epsilon = uni_epsilon,
        maxit = maxit,
        reltol = reltol,
        every = steps_chunk + 1L,
        optim_method = optim_method,
        r = r,
        burn = burn,
        erlang = erlang,
        merlang_blocks = merlang_blocks,
        plot_progress = FALSE
      )
    }, type = "output")
    ll_hist <- as.numeric(fit_obj@fit$logLikHist)
    if (length(ll_hist) == 0) {
      ll_hist <- rep(as.numeric(logLik(fit_obj)), steps_chunk)
    }
    if (length(ll_hist) < steps_chunk) {
      ll_hist <- c(ll_hist, rep(tail(ll_hist, 1), steps_chunk - length(ll_hist)))
    }
    list(model = fit_obj, ll = ll_hist[seq_len(steps_chunk)])
  }

  update_model <- function(i, chunk_result) {
    start <- done_steps[i] + 1L
    end <- min(stepsEM, start + length(chunk_result$ll) - 1L)
    states[[i]] <<- chunk_result$model
    track[i, start:end] <<- chunk_result$ll[seq_len(end - start + 1L)]
    done_steps[i] <<- end
    remaining[i] <<- stepsEM - end
  }
  
  current_ll <- function(i) {
    if (done_steps[i] < 1) {
      return(NA_real_)
    }
    as.numeric(track[i, done_steps[i]])
  }
  
  render_status <- function(start_elapsed) {
    elapsed_now <- unname(proc.time()["elapsed"]) - start_elapsed
    lines <- character(n_models + 1L)
    lines[1] <- paste0("elapsed ", formatC(elapsed_now, format = "f", digits = 1), "s")
    for (i in seq_len(n_models)) {
      ll_i <- current_ll(i)
      ll_txt <- if (is.finite(ll_i)) formatC(ll_i, format = "f", digits = 4) else "--"
      state_txt <- if (done_steps[i] >= stepsEM) "done" else "running"
      lines[i + 1L] <- paste0(
        model_names[i],
        " ",
        done_steps[i],
        "/",
        stepsEM,
        " ll ",
        ll_txt,
        " ",
        state_txt
      )
    }
    if (can_rewrite) {
      if (status_initialized) {
        cat("\033[", status_line_count, "A", sep = "")
      }
      for (j in seq_along(lines)) {
        cat("\r\033[2K", lines[j], "\n", sep = "")
      }
      status_initialized <<- TRUE
    } else {
      line <- paste(lines, collapse = " | ")
      line_width <- nchar(line, type = "width")
      pad <- if (status_width > line_width) strrep(" ", status_width - line_width) else ""
      cat("\r", line, pad, sep = "")
      status_width <<- max(status_width, line_width)
    }
    utils::flush.console()
  }

  options(digits.secs = 4)
  cat(format(Sys.time(), format = "%H:%M:%OS"), ": multiple EM started", sep = "")
  cat("\n", sep = "")

  start_elapsed <- unname(proc.time()["elapsed"])
  render_status(start_elapsed)
  plot_state()

  parallel_mode <- isTRUE(parallel) && n_models > 1L && cores > 1L && .Platform$OS.type != "windows"
  if (isTRUE(parallel) && !parallel_mode && .Platform$OS.type == "windows" && n_models > 1L) {
    warning("asynchronous parallel mode is not available on Windows; running sequentially")
  }

  if (parallel_mode) {
    cores <- as.integer(max(1L, min(cores, n_models)))
    active <- list()
    ready <- seq_len(n_models)

    launch_ready <- function() {
      while (length(active) < cores && length(ready) > 0) {
        idx <- ready[1]
        ready <<- ready[-1]
        if (remaining[idx] <= 0) {
          next
        }
        steps_chunk <- min(chunk_size, remaining[idx])
        job <- parallel::mcparallel(fit_chunk(states[[idx]], steps_chunk), silent = TRUE, mc.set.seed = TRUE)
        active[[as.character(job$pid)]] <<- list(job = job, idx = idx)
      }
    }

    on.exit({
      if (exists("active", inherits = FALSE) && length(active) > 0) {
        for (pid in suppressWarnings(as.integer(names(active)))) {
          if (is.finite(pid)) {
            try(tools::pskill(pid), silent = TRUE)
          }
        }
      }
    }, add = TRUE)

    while (length(active) > 0 || length(ready) > 0) {
      launch_ready()
      if (length(active) == 0) {
        next
      }
      collected <- parallel::mccollect(lapply(active, function(z) z$job), wait = FALSE)
      if (is.null(collected) || length(collected) == 0) {
        Sys.sleep(0.01)
        next
      }
      for (pid in names(collected)) {
        meta <- active[[pid]]
        active[[pid]] <- NULL
        res <- collected[[pid]]
        if (inherits(res, "try-error")) {
          stop(sprintf("model '%s' failed in asynchronous fitting", model_names[meta$idx]))
        }
        update_model(meta$idx, res)
        if (remaining[meta$idx] > 0) {
          ready <- c(ready, meta$idx)
        }
        render_status(start_elapsed)
        plot_state()
      }
    }
  } else {
    while (any(remaining > 0)) {
      for (i in seq_len(n_models)) {
        if (remaining[i] <= 0) {
          next
        }
        steps_chunk <- min(chunk_size, remaining[i])
        update_model(i, fit_chunk(states[[i]], steps_chunk))
        render_status(start_elapsed)
        plot_state()
      }
    }
  }

  elapsed <- unname(proc.time()["elapsed"]) - start_elapsed
  final_ll <- vapply(states, function(obj) as.numeric(logLik(obj)), numeric(1))
  best <- which.max(final_ll)

  if (!can_rewrite) {
    cat("\n", sep = "")
  }
  cat(format(Sys.time(), format = "%H:%M:%OS"), ": multiple EM finalized", sep = "")
  cat("\n", sep = "")

  list(
    models = states,
    model_names = model_names,
    logLikHist = track,
    final_logLik = final_ll,
    best = best,
    best_name = model_names[best],
    elapsed = elapsed
  )
}
