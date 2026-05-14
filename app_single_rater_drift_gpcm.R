# app.R - Single-rater Drift GPCM
#
# Shiny application for estimating time-varying rater severity drift
# in a single-rater setting using a GPCM/PCM-style ordinal model with
# a non-centered random-walk prior for beta[t].
#
# Main features:
#   - CSV/TSV/TXT/XLSX input
#   - Flexible column selection
#   - Optional time rebinning to a target number of time points
#   - Optional score shift for 0-based scores
#   - Weak or JSON-based prior for sigma_rw
#   - WAIC, PSIS-LOO, PPP, RMSE, and prior-vs-posterior diagnostics
#
# Expected input columns, or equivalent columns selected in the UI:
#   StudentID, TimeID, Score
# Optional:
#   ItemID (read but not used in this simplified single-rater model)

APP_VERSION <- "1.0.0"

suppressPackageStartupMessages({
  library(shiny)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(rstan)
  library(loo)
  library(jsonlite)
})

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(shiny.maxRequestSize = 200 * 1024^2)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

pkg_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

read_file_as_text <- function(path) {
  path <- as.character(path[1])
  if (is.na(path) || !nzchar(path)) stop("Empty file path.")
  if (!file.exists(path)) stop("File not found: ", path)

  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  txt <- tryCatch(rawToChar(raw, multiple = TRUE), error = function(e) NULL)
  if (is.null(txt)) stop("Failed to decode file as text.")
  paste(txt, collapse = "")
}

parse_prior_json <- function(path) {
  path <- as.character(path[1])
  if (is.na(path) || !nzchar(path)) stop("Empty prior JSON path.")
  if (!file.exists(path)) stop("Prior JSON file not found: ", path)

  txt <- read_file_as_text(path)
  first <- sub("^\\s+", "", sub("^\\ufeff", "", txt))
  if (!startsWith(first, "{")) {
    stop(
      "Uploaded prior file does not look like JSON. ",
      "The file should start with '{' and contain keys mu_log_sigma_rw and sd_log_sigma_rw."
    )
  }

  obj <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) NULL)
  if (!is.list(obj)) stop("Failed to parse prior JSON.")

  mu <- suppressWarnings(as.numeric(obj$mu_log_sigma_rw))
  sd <- suppressWarnings(as.numeric(obj$sd_log_sigma_rw))
  T_used <- suppressWarnings(as.integer(obj$T_used))

  if (!is.finite(mu)) stop("Prior JSON is missing a finite numeric key: mu_log_sigma_rw")
  if (!is.finite(sd) || sd <= 0) stop("Prior JSON is missing a positive numeric key: sd_log_sigma_rw")

  list(
    mu_log_sigma_rw = mu,
    sd_log_sigma_rw = sd,
    T_used = T_used,
    method = if (!is.null(obj$method)) as.character(obj$method) else NA_character_
  )
}

read_csv_robust <- function(path) {
  path <- as.character(path[1])

  try_read <- function(expr) tryCatch(expr, error = function(e) NULL)

  out <- try_read(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (!is.null(out)) return(out)

  out <- try_read(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8"))
  if (!is.null(out)) return(out)

  out <- try_read(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "CP932"))
  if (!is.null(out)) return(out)

  out <- try_read(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, sep = ";"))
  if (!is.null(out)) return(out)

  # Handle UTF-16 CSV files exported from some spreadsheet applications.
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  raw <- readBin(con, what = "raw", n = min(4000, file.info(path)$size))
  if (any(raw == as.raw(0x00))) {
    parse_lines <- function(lines) {
      tc <- textConnection(lines)
      on.exit(close(tc), add = TRUE)
      read.csv(tc, stringsAsFactors = FALSE, check.names = FALSE)
    }

    out <- try_read({
      lines <- readLines(file(path, open = "r", encoding = "UTF-16LE"), warn = FALSE)
      parse_lines(lines)
    })
    if (!is.null(out)) return(out)

    out <- try_read({
      lines <- readLines(file(path, open = "r", encoding = "UTF-16BE"), warn = FALSE)
      parse_lines(lines)
    })
    if (!is.null(out)) return(out)
  }

  stop("Failed to read CSV. Please export the file as UTF-8 CSV with headers.")
}

read_any_table <- function(datapath, filename, sheet = 1) {
  datapath <- as.character(datapath[1])
  filename <- as.character(filename[1])
  ext <- tolower(tools::file_ext(filename))

  if (ext %in% c("xlsx", "xls")) {
    return(as.data.frame(readxl::read_excel(datapath, sheet = sheet)))
  }

  if (ext %in% c("csv", "tsv", "txt")) {
    if (ext == "tsv") {
      out <- tryCatch(read.delim(datapath, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    return(read_csv_robust(datapath))
  }

  stop("Unsupported file extension: .", ext, ". Use .csv, .tsv, .txt, .xlsx, or .xls.")
}

rebin_time <- function(df, target_T, method = c("bin", "drop")) {
  method <- match.arg(method)

  if (is.null(target_T) || is.na(target_T) || target_T < 2) {
    return(list(data = df, mapping = NULL))
  }

  old_levels <- sort(unique(df$TimeID))
  old_T <- length(old_levels)

  if (old_T <= target_T) {
    mapping <- data.frame(old_time = old_levels, new_time = seq_along(old_levels))
    df$TimeID <- as.integer(match(df$TimeID, old_levels))
    return(list(data = df, mapping = mapping))
  }

  old_index <- as.integer(match(df$TimeID, old_levels))

  if (method == "bin") {
    new_time <- floor((old_index - 1) * target_T / old_T) + 1L
    mapping <- data.frame(
      old_time = old_levels,
      new_time = floor((seq_len(old_T) - 1) * target_T / old_T) + 1L
    )
    df$TimeID <- new_time
    return(list(data = df, mapping = mapping))
  }

  keep_levels <- old_levels[seq_len(target_T)]
  df <- df[df$TimeID %in% keep_levels, , drop = FALSE]
  new_levels <- sort(unique(df$TimeID))
  df$TimeID <- as.integer(match(df$TimeID, new_levels))
  mapping <- data.frame(old_time = new_levels, new_time = seq_along(new_levels))

  list(data = df, mapping = mapping)
}

make_column_selectors <- function(df) {
  names_original <- names(df)
  names_key <- tolower(gsub("\\s+", "", names_original))

  pick <- function(keys) {
    idx <- which(names_key %in% tolower(keys))
    if (length(idx) == 0) return("")
    names_original[idx[1]]
  }

  tagList(
    selectInput(
      "col_student", "Student/person ID column",
      choices = names_original,
      selected = pick(c("studentid", "personid", "examineeid", "student", "person", "sid", "id", "j_id", "j"))
    ),
    selectInput(
      "col_time", "Time column",
      choices = names_original,
      selected = pick(c("timeid", "time", "tid", "tau", "t"))
    ),
    selectInput(
      "col_item", "Item column (optional; not used by this model)",
      choices = c("", names_original),
      selected = pick(c("itemid", "item", "iid", "i_id", "i"))
    ),
    selectInput(
      "col_score", "Score column",
      choices = names_original,
      selected = pick(c("score", "response", "y", "category", "k"))
    )
  )
}

# -----------------------------------------------------------------------------
# Stan model
# -----------------------------------------------------------------------------

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> T;
  int<lower=2> K;

  array[N] int<lower=1, upper=J> jj;
  array[N] int<lower=1, upper=T> tt;
  array[N] int<lower=1, upper=K> y;

  int<lower=0> P;
  matrix[N, P] X;

  real<lower=0> D;

  real mu_log_sigma_rw;
  real<lower=0> sd_log_sigma_rw;
}
parameters {
  vector[J] z_theta;
  real mu_theta;
  real<lower=0> sigma_theta;

  vector[T - 1] z_beta;
  real<lower=0> sigma_rw;

  vector<lower=0>[K - 1] d_step;
  vector[P] gamma;
}
transformed parameters {
  vector[J] theta;
  vector[T] beta;
  vector[K] step_raw;
  vector[K] step;

  theta = mu_theta + sigma_theta * z_theta;

  beta[1] = 0;
  for (t in 2:T) {
    beta[t] = beta[t - 1] + sigma_rw * z_beta[t - 1];
  }

  step_raw[1] = 0;
  for (k in 2:K) {
    step_raw[k] = step_raw[k - 1] + d_step[k - 1];
  }

  {
    real m = mean(step_raw);
    for (k in 1:K) {
      step[k] = step_raw[k] - m;
    }
  }
}
model {
  z_theta ~ normal(0, 1);
  mu_theta ~ normal(0, 1);
  sigma_theta ~ student_t(3, 0, 1.5);

  z_beta ~ normal(0, 1);
  sigma_rw ~ lognormal(mu_log_sigma_rw, sd_log_sigma_rw);

  d_step ~ lognormal(-1.0, 0.6);
  gamma ~ normal(0, 1);

  for (n in 1:N) {
    real eta = theta[jj[n]] - beta[tt[n]];
    if (P > 0) {
      eta += dot_product(to_row_vector(X[n]), gamma);
    }

    vector[K] util;
    for (k in 1:K) {
      util[k] = D * (k * eta - step[k]);
    }

    y[n] ~ categorical_logit(util);
  }
}
generated quantities {
  vector[N] log_lik;
  real T_obs = 0;
  real T_rep = 0;
  real rmse = 0;
  real rmse_rep = 0;

  for (n in 1:N) {
    real eta = theta[jj[n]] - beta[tt[n]];
    if (P > 0) {
      eta += dot_product(to_row_vector(X[n]), gamma);
    }

    vector[K] util;
    for (k in 1:K) {
      util[k] = D * (k * eta - step[k]);
    }

    log_lik[n] = categorical_logit_lpmf(y[n] | util);

    {
      vector[K] p = softmax(util);
      int y_rep = categorical_rng(p);
      real mu = 0;
      real v = 0;

      for (k in 1:K) {
        mu += k * p[k];
      }
      for (k in 1:K) {
        real d = k - mu;
        v += square(d) * p[k];
      }

      T_obs += square(y[n] - mu) / (v + 1e-9);
      T_rep += square(y_rep - mu) / (v + 1e-9);
      rmse += square(y[n] - mu);
      rmse_rep += square(y_rep - y[n]);
    }
  }

  rmse = sqrt(rmse / N);
  rmse_rep = sqrt(rmse_rep / N);
}
"

# -----------------------------------------------------------------------------
# User interface
# -----------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel(paste("Single-rater Drift GPCM", APP_VERSION)),

  sidebarLayout(
    sidebarPanel(
      fileInput(
        "datafile", "Data file",
        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
        placeholder = "CSV/XLSX with StudentID, TimeID, and Score"
      ),
      numericInput("sheet", "Excel sheet index", value = 1, min = 1),
      uiOutput("col_selectors"),

      hr(),
      numericInput("time_target", "Target number of time points", value = 10, min = 2, step = 1),
      selectInput("time_method", "Time rebinning method", choices = c("bin", "drop"), selected = "bin"),
      checkboxInput("score_shift01", "Shift scores by +1 before category coding", value = FALSE),

      hr(),
      selectInput("prior_mode", "Prior for sigma_rw", choices = c("weak", "json"), selected = "weak"),
      fileInput("priorfile", "Prior JSON", accept = c(".json")),
      verbatimTextOutput("prior_used"),

      hr(),
      numericInput("Dscale", "Link scaling D", value = 1.7, min = 0.5, step = 0.1),
      numericInput("chains", "Chains", value = 4, min = 1),
      numericInput("iter", "Iterations", value = 4000, min = 1000, step = 500),
      numericInput("warmup", "Warmup", value = 2000, min = 500, step = 500),
      numericInput("adapt", "adapt_delta", value = 0.995, min = 0.8, max = 0.999, step = 0.001),
      numericInput("treedepth", "max_treedepth", value = 15, min = 10, max = 15, step = 1),
      numericInput("seed", "Seed", value = 123, min = 1),
      actionButton("run", "Run Stan", class = "btn-primary"),

      hr(),
      tags$details(
        tags$summary("Session information"),
        verbatimTextOutput("session_info")
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Beta drift", plotOutput("beta_plot", height = "420px"), tableOutput("beta_table")),
        tabPanel("WAIC / LOO", verbatimTextOutput("waicloo_text")),
        tabPanel("PPP / RMSE", verbatimTextOutput("ppp_rmse_text")),
        tabPanel("Prior vs posterior", plotOutput("sigma_prior_post_plot", height = "340px"), tableOutput("prior_post_table")),
        tabPanel("Diagnostics", verbatimTextOutput("diag_text")),
        tabPanel("Data preview", tableOutput("data_preview"))
      )
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {
  stan_model_cache <- reactiveVal(NULL)

  output$session_info <- renderPrint({
    cat("App version:", APP_VERSION, "\n")
    cat("Working directory:", getwd(), "\n")
    cat("R:", R.version.string, "\n")
    cat("shiny:", pkg_version("shiny"), "\n")
    cat("rstan:", pkg_version("rstan"), "\n")
    cat("loo:", pkg_version("loo"), "\n")
    cat("jsonlite:", pkg_version("jsonlite"), "\n")
  })

  raw_data <- reactive({
    req(input$datafile)
    datapath <- as.character(input$datafile$datapath[1])
    filename <- as.character(input$datafile$name[1])

    validate(need(file.exists(datapath), "Uploaded file path not found."))

    if (tolower(tools::file_ext(filename)) %in% c("xlsx", "xls")) {
      validate(need(is.numeric(input$sheet) && input$sheet >= 1, "Excel sheet index must be >= 1."))
    }

    tryCatch(
      read_any_table(datapath, filename, sheet = input$sheet),
      error = function(e) stop("Failed to read data file: ", conditionMessage(e))
    )
  })

  output$col_selectors <- renderUI({
    req(raw_data())
    tryCatch(
      make_column_selectors(raw_data()),
      error = function(e) tags$div(style = "color:#b00020;", paste("Column selector failed:", conditionMessage(e)))
    )
  })

  current_prior <- reactive({
    if (input$prior_mode == "weak") {
      return(list(mu_log_sigma_rw = -1.0, sd_log_sigma_rw = 0.6, T_used = NA_integer_, method = "weak"))
    }

    req(input$priorfile)
    parse_prior_json(input$priorfile$datapath[1])
  })

  output$prior_used <- renderPrint({
    tryCatch({
      pr <- current_prior()
      cat(sprintf("sigma_rw ~ LogNormal(%.6f, %.6f)\n", pr$mu_log_sigma_rw, pr$sd_log_sigma_rw))
      cat("Prior mode:", input$prior_mode, "\n")
      if (is.finite(pr$T_used)) cat("T used in prior:", pr$T_used, "\n")
      if (!is.na(pr$method)) cat("Prior method:", pr$method, "\n")
    }, error = function(e) {
      cat("Prior parsing failed:\n")
      cat(conditionMessage(e), "\n")
    })
  })

  observeEvent(input$run, {
    df0 <- raw_data()

    validate(need(!is.null(input$col_student) && nzchar(input$col_student), "Select a student/person ID column."))
    validate(need(!is.null(input$col_time) && nzchar(input$col_time), "Select a time column."))
    validate(need(!is.null(input$col_score) && nzchar(input$col_score), "Select a score column."))

    df <- data.frame(
      StudentID = df0[[input$col_student]],
      TimeID = df0[[input$col_time]],
      ItemID = if (!is.null(input$col_item) && nzchar(input$col_item)) df0[[input$col_item]] else NA,
      Score = df0[[input$col_score]],
      stringsAsFactors = FALSE
    )

    df <- df %>%
      mutate(
        StudentID = suppressWarnings(as.numeric(StudentID)),
        TimeID = suppressWarnings(as.numeric(TimeID)),
        ItemID = suppressWarnings(as.numeric(ItemID)),
        Score = suppressWarnings(as.numeric(Score))
      ) %>%
      filter(!is.na(StudentID), !is.na(TimeID), !is.na(Score))

    validate(need(nrow(df) > 0, "No valid rows after numeric conversion."))

    if (isTRUE(input$score_shift01)) {
      df$Score <- df$Score + 1
    }

    rb <- rebin_time(df, target_T = input$time_target, method = input$time_method)
    df <- rb$data

    df <- df %>%
      mutate(
        j = as.integer(factor(StudentID, levels = sort(unique(StudentID)))),
        t = as.integer(factor(TimeID, levels = sort(unique(TimeID))))
      )

    score_levels <- sort(unique(df$Score))
    validate(need(length(score_levels) >= 2, "Score must have at least two distinct values."))

    score_map <- setNames(seq_along(score_levels), score_levels)
    df$y <- as.integer(score_map[as.character(df$Score)])
    K <- length(score_levels)

    pr <- tryCatch(current_prior(), error = function(e) e)
    if (inherits(pr, "error")) {
      showNotification(paste("Prior parsing failed:", conditionMessage(pr)), type = "error", duration = NULL)
      return(NULL)
    }

    stan_data <- list(
      N = nrow(df),
      J = length(unique(df$j)),
      T = length(unique(df$t)),
      K = as.integer(K),
      jj = as.integer(df$j),
      tt = as.integer(df$t),
      y = as.integer(df$y),
      P = 0L,
      X = matrix(0, nrow = nrow(df), ncol = 0),
      D = as.numeric(input$Dscale),
      mu_log_sigma_rw = pr$mu_log_sigma_rw,
      sd_log_sigma_rw = pr$sd_log_sigma_rw
    )

    mu0 <- stan_data$mu_log_sigma_rw
    sd0 <- stan_data$sd_log_sigma_rw

    output$data_preview <- renderTable({
      head(df %>% select(StudentID, TimeID, Score, j, t, y), 20)
    })

    if (is.null(stan_model_cache())) {
      stan_model_cache(rstan::stan_model(model_code = stan_code))
    }

    fit <- tryCatch({
      rstan::sampling(
        object = stan_model_cache(),
        data = stan_data,
        chains = input$chains,
        iter = input$iter,
        warmup = input$warmup,
        seed = input$seed,
        control = list(adapt_delta = input$adapt, max_treedepth = input$treedepth),
        refresh = max(1, input$iter %/% 10)
      )
    }, error = function(e) e)

    if (inherits(fit, "error")) {
      showNotification(paste("Stan sampling failed:", conditionMessage(fit)), type = "error", duration = NULL)
      output$diag_text <- renderPrint({ print(fit) })
      return(NULL)
    }

    # Beta drift summary -------------------------------------------------------
    beta_mat <- as.matrix(fit, pars = "beta")
    beta_df <- data.frame(
      time = seq_len(ncol(beta_mat)),
      mean = colMeans(beta_mat),
      q05 = apply(beta_mat, 2, stats::quantile, 0.05),
      median = apply(beta_mat, 2, stats::median),
      q95 = apply(beta_mat, 2, stats::quantile, 0.95)
    )

    output$beta_plot <- renderPlot({
      ggplot(beta_df, aes(x = time, y = mean)) +
        geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
        geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        scale_x_continuous(breaks = beta_df$time) +
        labs(
          x = "Time point",
          y = "beta[t]",
          title = "Rater severity drift: posterior mean and 90% interval",
          subtitle = "Random-walk beta[t] with anchor beta[1] = 0"
        ) +
        theme_minimal(base_size = 13)
    })

    output$beta_table <- renderTable({ beta_df }, digits = 3)

    # WAIC / LOO ---------------------------------------------------------------
    output$waicloo_text <- renderPrint({
      log_lik <- as.matrix(fit, pars = "log_lik")

      waic_res <- tryCatch(loo::waic(log_lik), error = function(e) e)
      loo_res <- tryCatch(loo::loo(log_lik), error = function(e) e)

      cat("=== WAIC ===\n")
      if (inherits(waic_res, "error")) {
        cat("WAIC failed: ", conditionMessage(waic_res), "\n\n", sep = "")
      } else {
        print(waic_res$estimates, digits = 6)
        cat("\nIf p_waic warnings appear, prefer PSIS-LOO.\n\n")
      }

      cat("=== PSIS-LOO ===\n")
      if (inherits(loo_res, "error")) {
        cat("LOO failed: ", conditionMessage(loo_res), "\n", sep = "")
      } else {
        print(loo_res$estimates, digits = 6)

        k <- loo_res$diagnostics$pareto_k
        if (!is.null(k)) {
          k_cat <- cut(
            k,
            breaks = c(-Inf, 0.5, 0.7, 1, Inf),
            labels = c("good (<=0.5)", "ok (0.5-0.7)", "bad (0.7-1)", "very bad (>1)"),
            right = TRUE
          )
          cat("\nPareto-k counts:\n")
          print(table(k_cat, useNA = "ifany"))
          cat("\nPareto-k summary:\n")
          print(summary(k))
        }
      }
    })

    # PPP / RMSE ---------------------------------------------------------------
    output$ppp_rmse_text <- renderPrint({
      T_obs <- as.vector(as.matrix(fit, pars = "T_obs"))
      T_rep <- as.vector(as.matrix(fit, pars = "T_rep"))
      ppp <- mean(T_rep > T_obs)

      rmse <- as.vector(as.matrix(fit, pars = "rmse"))
      rmse_rep <- as.vector(as.matrix(fit, pars = "rmse_rep"))

      cat(sprintf("Posterior predictive p-value: %.6f\n\n", ppp))

      cat("RMSE: expected category vs observed category\n")
      print(stats::quantile(rmse, c(0.025, 0.5, 0.975)), digits = 6)
      cat(sprintf("Mean: %.6f\n\n", mean(rmse)))

      cat("RMSE: replicated category vs observed category\n")
      print(stats::quantile(rmse_rep, c(0.025, 0.5, 0.975)), digits = 6)
      cat(sprintf("Mean: %.6f\n", mean(rmse_rep)))
    })

    # Prior vs posterior for sigma_rw -----------------------------------------
    output$sigma_prior_post_plot <- renderPlot({
      posterior <- as.vector(as.matrix(fit, pars = "sigma_rw"))
      x_max <- max(
        as.numeric(stats::quantile(posterior, 0.995)),
        stats::qlnorm(0.995, meanlog = mu0, sdlog = sd0)
      )

      posterior_density <- density(posterior, from = 0, to = x_max, n = 400)
      plot_data <- data.frame(
        sigma_rw = rep(posterior_density$x, 2),
        density = c(
          stats::dlnorm(posterior_density$x, meanlog = mu0, sdlog = sd0),
          posterior_density$y
        ),
        distribution = rep(c("Prior", "Posterior"), each = length(posterior_density$x))
      )

      ggplot(plot_data, aes(x = sigma_rw, y = density, linetype = distribution)) +
        geom_line(linewidth = 1) +
        labs(
          x = "sigma_rw",
          y = "Density",
          title = "sigma_rw: prior vs posterior",
          subtitle = sprintf("Prior: LogNormal(%.4f, %.4f)", mu0, sd0)
        ) +
        theme_minimal(base_size = 13)
    })

    output$prior_post_table <- renderTable({
      posterior <- as.vector(as.matrix(fit, pars = "sigma_rw"))
      prior_q <- stats::qlnorm(c(0.05, 0.5, 0.95), meanlog = mu0, sdlog = sd0)
      posterior_q <- stats::quantile(posterior, c(0.05, 0.5, 0.95))

      data.frame(
        quantity = c("mean", "sd", "q05", "median", "q95"),
        prior = c(
          exp(mu0 + sd0^2 / 2),
          sqrt((exp(sd0^2) - 1) * exp(2 * mu0 + sd0^2)),
          prior_q[1], prior_q[2], prior_q[3]
        ),
        posterior = c(
          mean(posterior),
          stats::sd(posterior),
          posterior_q[[1]], posterior_q[[2]], posterior_q[[3]]
        )
      )
    }, digits = 6)

    # Diagnostics --------------------------------------------------------------
    output$diag_text <- renderPrint({
      cat("=== Prior hyperparameters passed to Stan ===\n")
      cat(sprintf("prior_mode = %s\n", input$prior_mode))
      cat(sprintf("mu_log_sigma_rw = %.6f\n", mu0))
      cat(sprintf("sd_log_sigma_rw = %.6f\n\n", sd0))

      print(fit, pars = c("mu_theta", "sigma_theta", "sigma_rw", "beta"), probs = c(0.05, 0.5, 0.95))

      sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      divergences <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
      treedepth_hits <- sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= input$treedepth))

      cat("\nDivergences per chain: ", paste(divergences, collapse = ", "),
          " (total: ", sum(divergences), ")\n", sep = "")
      cat("Max-treedepth hits per chain: ", paste(treedepth_hits, collapse = ", "),
          " (total: ", sum(treedepth_hits), ")\n", sep = "")
    })
  })
}

shinyApp(ui, server)
