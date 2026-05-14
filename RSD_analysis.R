# app.R
# Single-Rater Drift GPCM
#
# This Shiny application fits a single-rater drift model using RStan.
# The time-varying rater severity parameter beta[t] follows a non-centered
# random walk with beta[1] fixed to 0.
#
# Main outputs:
#   - Posterior trajectory of beta[t]
#   - WAIC and PSIS-LOO based on pointwise log-likelihood
#   - Posterior predictive p-value (PPP)
#   - RMSE summaries
#   - Prior vs posterior comparison for sigma_rw
#
# Required input data columns:
#   - StudentID
#   - TimeID
#   - Score
#
# Optional input column:
#   - ItemID
#
# Prior JSON format, if used:
# {
#   "mu_log_sigma_rw": -0.290109,
#   "sd_log_sigma_rw": 0.600000,
#   "T_used": 10,
#   "method": "learned_prior_revised"
# }

APP_VERSION <- "1.0.1"

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
options(shiny.sanitize.errors = FALSE)
options(shiny.fullstacktrace = TRUE)

# -------------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------------

package_version_safe <- function(pkg) {
  tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
}

is_blank <- function(x) {
  is.null(x) || length(x) == 0 || is.na(x[1]) || !nzchar(as.character(x[1]))
}

looks_like_xlsx <- function(path) {
  path <- as.character(path[1])
  if (is_blank(path) || !file.exists(path)) return(FALSE)

  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)

  raw4 <- readBin(con, what = "raw", n = 4)
  length(raw4) >= 4 &&
    identical(raw4[1:4], as.raw(c(0x50, 0x4B, 0x03, 0x04))) # ZIP / XLSX signature: PK\003\004
}

looks_like_xls <- function(path) {
  path <- as.character(path[1])
  if (is_blank(path) || !file.exists(path)) return(FALSE)

  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)

  raw8 <- readBin(con, what = "raw", n = 8)
  length(raw8) >= 8 &&
    identical(raw8[1:8], as.raw(c(0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1))) # legacy XLS signature
}

read_file_as_text <- function(path) {
  path <- as.character(path[1])

  if (is_blank(path)) {
    stop("Empty file path.")
  }

  if (!file.exists(path)) {
    stop(paste0("File not found: ", path))
  }

  raw <- readBin(path, what = "raw", n = file.info(path)$size)

  if (length(raw) >= 4 && identical(raw[1:4], as.raw(c(0x50, 0x4B, 0x03, 0x04)))) {
    stop(
      "The uploaded file appears to be an Excel/ZIP file, not a plain text file. ",
      "For the prior, upload a valid JSON file. For the data, upload a valid CSV/XLSX file."
    )
  }

  txt <- tryCatch(
    rawToChar(raw, multiple = FALSE),
    error = function(e) {
      stop(paste0("Failed to decode the file as text: ", conditionMessage(e)))
    }
  )

  txt <- enc2utf8(txt)

  if (!is.character(txt) || length(txt) != 1L) {
    stop("Internal read error: decoded text is not a single character string.")
  }

  txt
}

parse_prior_json <- function(path) {
  path <- as.character(path[1])

  if (is_blank(path)) {
    stop("Empty prior JSON path.")
  }

  if (!file.exists(path)) {
    stop(paste0("Prior JSON file not found: ", path))
  }

  txt <- read_file_as_text(path)
  txt <- sub("^\ufeff", "", txt)

  obj <- tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = TRUE),
    error = function(e) NULL
  )

  if (is.list(obj)) {
    mu <- suppressWarnings(as.numeric(obj$mu_log_sigma_rw))
    sd <- suppressWarnings(as.numeric(obj$sd_log_sigma_rw))
    T_used <- suppressWarnings(as.integer(obj$T_used))
    method <- if (!is.null(obj$method)) as.character(obj$method) else NA_character_

    if (is.finite(mu) && is.finite(sd) && sd > 0) {
      return(list(
        mu_log_sigma_rw = mu,
        sd_log_sigma_rw = sd,
        T_used = T_used,
        method = method,
        parse_method = "jsonlite_text"
      ))
    }
  }

  # Fallback for simple JSON-like files.
  grab_num <- function(key) {
    pattern <- paste0('"', key, '"\\s*:\\s*([-+0-9.eE]+)')
    m <- regexec(pattern, txt, perl = TRUE)
    r <- regmatches(txt, m)
    if (length(r) == 0 || length(r[[1]]) < 2) return(NA_real_)
    suppressWarnings(as.numeric(r[[1]][2]))
  }

  mu <- grab_num("mu_log_sigma_rw")
  sd <- grab_num("sd_log_sigma_rw")
  T_used <- suppressWarnings(as.integer(grab_num("T_used")))

  if (!is.finite(mu) || !is.finite(sd) || sd <= 0) {
    stop(
      "Could not extract valid prior hyperparameters from JSON. ",
      "The JSON must contain numeric keys: mu_log_sigma_rw and sd_log_sigma_rw."
    )
  }

  list(
    mu_log_sigma_rw = mu,
    sd_log_sigma_rw = sd,
    T_used = T_used,
    method = NA_character_,
    parse_method = "regex"
  )
}

read_csv_robust <- function(path) {
  path <- as.character(path[1])

  if (looks_like_xlsx(path) || looks_like_xls(path)) {
    stop(
      "The uploaded file has an Excel workbook signature. ",
      "Please upload it as .xlsx/.xls, or export it as a true CSV file."
    )
  }

  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  raw <- readBin(con, what = "raw", n = 4000)
  has_nul <- any(raw == as.raw(0x00))

  parse_lines_as_csv <- function(lines) {
    tc <- textConnection(lines)
    on.exit(close(tc), add = TRUE)
    read.csv(tc, stringsAsFactors = FALSE, check.names = FALSE)
  }

  out <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  out <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8"),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  if (has_nul) {
    out <- tryCatch({
      lines <- readLines(file(path, open = "r", encoding = "UTF-16LE"), warn = FALSE)
      parse_lines_as_csv(lines)
    }, error = function(e) NULL)
    if (!is.null(out)) return(out)

    out <- tryCatch({
      lines <- readLines(file(path, open = "r", encoding = "UTF-16BE"), warn = FALSE)
      parse_lines_as_csv(lines)
    }, error = function(e) NULL)
    if (!is.null(out)) return(out)
  }

  out <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "CP932"),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  out <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, sep = ";"),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  stop("Failed to read the CSV file. Please export the data as UTF-8 CSV with headers.")
}

read_excel_safe <- function(path, sheet = 1, ext = "xlsx") {
  path <- as.character(path[1])

  file_ext <- if (tolower(ext) %in% c("xls", "xlsx")) paste0(".", tolower(ext)) else ".xlsx"
  tmp <- tempfile(fileext = file_ext)
  ok <- file.copy(path, tmp, overwrite = TRUE)

  if (!isTRUE(ok)) {
    stop("Failed to copy the uploaded Excel file to a temporary path.")
  }

  df <- readxl::read_excel(tmp, sheet = sheet)
  as.data.frame(df)
}

read_any_table <- function(datapath, filename, sheet = 1) {
  datapath <- as.character(datapath[1])
  filename <- as.character(filename[1])
  ext <- tolower(tools::file_ext(filename))

  # Detect file type by signature first.
  if (looks_like_xlsx(datapath)) {
    return(read_excel_safe(datapath, sheet = sheet, ext = "xlsx"))
  }

  if (looks_like_xls(datapath)) {
    return(read_excel_safe(datapath, sheet = sheet, ext = "xls"))
  }

  if (ext %in% c("xlsx", "xls")) {
    return(read_excel_safe(datapath, sheet = sheet, ext = ext))
  }

  if (ext %in% c("csv", "tsv", "txt")) {
    if (ext == "tsv") {
      df <- tryCatch(
        read.delim(datapath, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
      )
      if (!is.null(df)) return(df)
    }

    return(read_csv_robust(datapath))
  }

  stop(paste0(
    "Unsupported file extension: .", ext,
    ". Please upload a CSV, TSV, TXT, XLSX, or XLS file."
  ))
}

rebin_time <- function(df, T_target, method = c("bin", "drop")) {
  method <- match.arg(method)

  if (is.null(T_target) || is.na(T_target) || T_target < 2) {
    return(list(df = df, mapping = NULL))
  }

  levels_old <- sort(unique(df$TimeID))
  T_old <- length(levels_old)

  if (T_old <= T_target) {
    return(list(
      df = df,
      mapping = data.frame(old = levels_old, new = seq_along(levels_old))
    ))
  }

  idx_map <- setNames(seq_along(levels_old), levels_old)
  t_old <- as.integer(idx_map[as.character(df$TimeID)])

  if (method == "bin") {
    t_new <- floor((t_old - 1) * T_target / T_old) + 1L
    mapping <- data.frame(
      old = levels_old,
      new = floor((seq_len(T_old) - 1) * T_target / T_old) + 1L
    )
    df$TimeID <- t_new
    return(list(df = df, mapping = mapping))
  }

  keep_levels <- levels_old[seq_len(T_target)]
  df2 <- df[df$TimeID %in% keep_levels, , drop = FALSE]
  idx_map2 <- setNames(seq_along(sort(unique(df2$TimeID))), sort(unique(df2$TimeID)))
  df2$TimeID <- as.integer(idx_map2[as.character(df2$TimeID)])

  list(
    df = df2,
    mapping = data.frame(old = keep_levels, new = seq_len(length(keep_levels)))
  )
}

# -------------------------------------------------------------------------
# Stan model
# -------------------------------------------------------------------------

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> T;
  int<lower=2> K;

  int<lower=1, upper=J> jj[N];
  int<lower=1, upper=T> tt[N];
  int<lower=1, upper=K> y[N];

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

    {
      vector[K] util;
      for (k in 1:K) {
        util[k] = D * (k * eta - step[k]);
      }
      y[n] ~ categorical_logit(util);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  real T_obs;
  real T_rep;
  real rmse;
  real rmse_rep;

  T_obs = 0;
  T_rep = 0;
  rmse = 0;
  rmse_rep = 0;

  for (n in 1:N) {
    real eta = theta[jj[n]] - beta[tt[n]];

    if (P > 0) {
      eta += dot_product(to_row_vector(X[n]), gamma);
    }

    {
      vector[K] util;
      vector[K] p;
      int y_rep;
      real mu_y;
      real var_y;

      for (k in 1:K) {
        util[k] = D * (k * eta - step[k]);
      }

      log_lik[n] = categorical_logit_lpmf(y[n] | util);

      p = softmax(util);
      y_rep = categorical_rng(p);

      mu_y = 0;
      var_y = 0;

      for (k in 1:K) {
        mu_y += k * p[k];
      }

      for (k in 1:K) {
        real d = k - mu_y;
        var_y += d * d * p[k];
      }

      T_obs += square(y[n] - mu_y) / (var_y + 1e-9);
      T_rep += square(y_rep - mu_y) / (var_y + 1e-9);

      rmse += square(y[n] - mu_y);
      rmse_rep += square(y_rep - y[n]);
    }
  }

  rmse = sqrt(rmse / N);
  rmse_rep = sqrt(rmse_rep / N);
}
"

# -------------------------------------------------------------------------
# Shiny UI
# -------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel(paste0("Single-Rater Drift GPCM ", APP_VERSION)),

  sidebarLayout(
    sidebarPanel(
      fileInput(
        "datafile",
        "Data file",
        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
      ),
      helpText("Required columns: StudentID, TimeID, Score. ItemID is optional and not used by this model."),

      numericInput("sheet", "Excel sheet index", value = 1, min = 1, step = 1),
      uiOutput("col_selectors"),

      hr(),

      numericInput("time_target", "Target number of time points", value = 10, min = 2, step = 1),
      selectInput("time_method", "Time rebinning method", choices = c("bin", "drop"), selected = "bin"),
      checkboxInput("score_shift01", "Shift scores by +1 before category recoding", value = FALSE),

      hr(),

      selectInput("prior_mode", "Prior for sigma_rw", choices = c("weak", "json"), selected = "weak"),
      fileInput("priorfile", "Prior JSON", accept = c(".json", ".txt")),
      verbatimTextOutput("prior_used"),

      hr(),

      numericInput("Dscale", "Link scaling constant D", value = 1.7, min = 0.5, step = 0.1),
      numericInput("chains", "Chains", value = 4, min = 1, step = 1),
      numericInput("iter", "Iterations", value = 4000, min = 1000, step = 500),
      numericInput("warmup", "Warmup", value = 2000, min = 500, step = 500),
      numericInput("adapt", "adapt_delta", value = 0.995, min = 0.8, max = 0.999, step = 0.001),
      numericInput("treedepth", "max_treedepth", value = 15, min = 10, max = 15, step = 1),
      numericInput("seed", "Seed", value = 123, min = 1, step = 1),

      actionButton("run", "Run Stan", class = "btn-primary"),

      hr(),

      tags$details(
        tags$summary("Session information"),
        verbatimTextOutput("session_info", placeholder = TRUE)
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Severity drift", plotOutput("beta_plot", height = "420px"), tableOutput("beta_table")),
        tabPanel("WAIC / LOO", verbatimTextOutput("waicloo_text")),
        tabPanel("PPP / RMSE", verbatimTextOutput("ppp_rmse_text")),
        tabPanel("Prior vs posterior", plotOutput("sigma_prior_post_plot", height = "340px"), tableOutput("prior_post_table")),
        tabPanel("Diagnostics", verbatimTextOutput("diag_text")),
        tabPanel("Cleaned data", tableOutput("data_preview"))
      )
    )
  )
)

# -------------------------------------------------------------------------
# Shiny server
# -------------------------------------------------------------------------

server <- function(input, output, session) {
  stan_model_cached <- reactiveVal(NULL)

  output$session_info <- renderPrint({
    cat("APP_VERSION:", APP_VERSION, "\n")
    cat("Working directory:", getwd(), "\n")
    cat("R:", R.version.string, "\n")
    cat("shiny:", package_version_safe("shiny"), "\n")
    cat("rstan:", package_version_safe("rstan"), "\n")
    cat("loo:", package_version_safe("loo"), "\n")
    cat("jsonlite:", package_version_safe("jsonlite"), "\n")
    cat("readxl:", package_version_safe("readxl"), "\n")
  })

  raw_df <- reactive({
    req(input$datafile)

    dp <- as.character(input$datafile$datapath[1])
    fn <- as.character(input$datafile$name[1])

    shiny::validate(
      shiny::need(file.exists(dp), paste0("Uploaded file path not found: ", dp))
    )

    sh <- input$sheet
    if (looks_like_xlsx(dp) || looks_like_xls(dp) || tolower(tools::file_ext(fn)) %in% c("xlsx", "xls")) {
      shiny::validate(shiny::need(is.numeric(sh) && sh >= 1, "Sheet index must be >= 1."))
    }

    tryCatch(
      {
        df0 <- read_any_table(dp, fn, sheet = sh)
        shiny::validate(
          shiny::need(is.data.frame(df0), "The uploaded file could not be converted to a data frame."),
          shiny::need(nrow(df0) > 0, "The uploaded file has no rows."),
          shiny::need(ncol(df0) > 0, "The uploaded file has no columns.")
        )
        df0
      },
      error = function(e) {
        stop(paste0(
          "Data file read error.\n",
          "Filename: ", fn, "\n",
          "Detected XLSX signature: ", looks_like_xlsx(dp), "\n",
          "Detected XLS signature: ", looks_like_xls(dp), "\n",
          "Details: ", conditionMessage(e)
        ))
      }
    )
  })

  output$col_selectors <- renderUI({
    # Do not call raw_df() before a file is uploaded.
    # req() raises a silent Shiny condition, but tryCatch() can accidentally
    # convert it into a visible "Column selector error" message.
    if (is.null(input$datafile)) {
      return(tags$div(
        style = "color: #666; font-size: 0.9em;",
        "Upload a data file to select columns."
      ))
    }

    df0 <- tryCatch(
      raw_df(),
      error = function(e) {
        return(structure(list(message = conditionMessage(e)), class = "column_selector_error"))
      }
    )

    if (inherits(df0, "column_selector_error")) {
      return(tags$div(
        style = "color: #b00020; white-space: pre-wrap;",
        paste0("Column selector error:\n", df0$message)
      ))
    }

    nm0 <- names(df0)
    nm <- tolower(gsub("\\s+", "", nm0))

    pick <- function(keys) {
      idx <- which(nm %in% tolower(keys))
      if (length(idx) == 0) return("")
      nm0[idx[1]]
    }

    tagList(
      selectInput(
        "col_student",
        "StudentID column",
        choices = nm0,
        selected = pick(c("studentid", "personid", "examineeid", "student", "person", "sid", "id", "j_id", "j"))
      ),
      selectInput(
        "col_time",
        "TimeID column",
        choices = nm0,
        selected = pick(c("timeid", "time", "tid", "tau", "t"))
      ),
      selectInput(
        "col_item",
        "ItemID column (optional; not used)",
        choices = c("", nm0),
        selected = pick(c("itemid", "item", "iid", "i_id", "i"))
      ),
      selectInput(
        "col_score",
        "Score column",
        choices = nm0,
        selected = pick(c("score", "response", "y", "category", "k"))
      )
    )
  })

  output$prior_used <- renderPrint({
    if (input$prior_mode == "weak") {
      cat("Using weak prior:\n")
      cat("  sigma_rw ~ LogNormal(-1.0000, 0.6000)\n")
      return(invisible(NULL))
    }

    if (is.null(input$priorfile)) {
      cat("Prior mode is JSON, but no JSON file has been uploaded.\n")
      return(invisible(NULL))
    }

    pr <- tryCatch(
      parse_prior_json(input$priorfile$datapath[1]),
      error = function(e) e
    )

    if (inherits(pr, "error")) {
      cat("Prior JSON parse error:\n")
      cat(conditionMessage(pr), "\n")
      return(invisible(NULL))
    }

    cat("Using JSON prior:\n")
    cat(sprintf("  sigma_rw ~ LogNormal(%.6f, %.6f)\n", pr$mu_log_sigma_rw, pr$sd_log_sigma_rw))
    cat(sprintf("  Parsed by: %s\n", pr$parse_method))

    if (is.finite(pr$T_used)) {
      cat(sprintf("  T_used: %s\n", pr$T_used))
    }

    if (!is.na(pr$method)) {
      cat(sprintf("  Method: %s\n", pr$method))
    }
  })

  observeEvent(input$run, {
    req(input$datafile)

    df0 <- raw_df()

    shiny::validate(
      shiny::need(!is_blank(input$col_student), "Select the StudentID column."),
      shiny::need(!is_blank(input$col_time), "Select the TimeID column."),
      shiny::need(!is_blank(input$col_score), "Select the Score column.")
    )

    df <- data.frame(
      StudentID = df0[[input$col_student]],
      TimeID = df0[[input$col_time]],
      ItemID = if (!is_blank(input$col_item)) df0[[input$col_item]] else NA,
      Score = df0[[input$col_score]]
    )

    df <- df %>%
      mutate(
        StudentID = suppressWarnings(as.numeric(StudentID)),
        TimeID = suppressWarnings(as.numeric(TimeID)),
        ItemID = suppressWarnings(as.numeric(ItemID)),
        Score = suppressWarnings(as.numeric(Score))
      ) %>%
      filter(!is.na(StudentID), !is.na(TimeID), !is.na(Score))

    shiny::validate(
      shiny::need(nrow(df) > 0, "No valid rows remain after numeric conversion and missing-value filtering.")
    )

    if (isTRUE(input$score_shift01)) {
      df$Score <- df$Score + 1
    }

    rb <- rebin_time(df, T_target = input$time_target, method = input$time_method)
    df2 <- rb$df

    df2 <- df2 %>%
      mutate(
        j = as.integer(factor(StudentID, levels = sort(unique(StudentID)))),
        t = as.integer(factor(TimeID, levels = sort(unique(TimeID))))
      )

    score_levels <- sort(unique(df2$Score))
    shiny::validate(
      shiny::need(length(score_levels) >= 2, "Score must have at least two distinct values.")
    )

    score_map <- setNames(seq_along(score_levels), score_levels)
    df2$y <- as.integer(score_map[as.character(df2$Score)])
    K <- length(score_levels)

    stan_data <- list(
      N = nrow(df2),
      J = length(unique(df2$j)),
      T = length(unique(df2$t)),
      K = as.integer(K),
      jj = as.integer(df2$j),
      tt = as.integer(df2$t),
      y = as.integer(df2$y),
      P = 0L,
      X = matrix(0, nrow = nrow(df2), ncol = 0),
      D = as.numeric(input$Dscale),
      mu_log_sigma_rw = -1.0,
      sd_log_sigma_rw = 0.6
    )

    if (input$prior_mode == "json") {
      if (is.null(input$priorfile) || is_blank(input$priorfile$datapath[1])) {
        showNotification("Prior mode is JSON, but no JSON file has been uploaded.", type = "error", duration = NULL)
        return(NULL)
      }

      pr <- tryCatch(
        parse_prior_json(as.character(input$priorfile$datapath[1])),
        error = function(e) e
      )

      if (inherits(pr, "error")) {
        showNotification(paste0("Prior JSON parse error: ", conditionMessage(pr)), type = "error", duration = NULL)
        output$diag_text <- renderPrint({
          cat("Prior JSON parse error:\n")
          cat(conditionMessage(pr), "\n")
        })
        return(NULL)
      }

      stan_data$mu_log_sigma_rw <- pr$mu_log_sigma_rw
      stan_data$sd_log_sigma_rw <- pr$sd_log_sigma_rw

      if (!is.finite(stan_data$mu_log_sigma_rw) ||
          !is.finite(stan_data$sd_log_sigma_rw) ||
          stan_data$sd_log_sigma_rw <= 0) {
        showNotification("Invalid prior hyperparameters in JSON.", type = "error", duration = NULL)
        output$diag_text <- renderPrint({
          cat("Invalid prior hyperparameters in JSON.\n")
          cat("mu_log_sigma_rw:", stan_data$mu_log_sigma_rw, "\n")
          cat("sd_log_sigma_rw:", stan_data$sd_log_sigma_rw, "\n")
        })
        return(NULL)
      }

      showNotification(
        sprintf(
          "Using JSON prior: sigma_rw ~ LogNormal(%.4f, %.4f)",
          stan_data$mu_log_sigma_rw,
          stan_data$sd_log_sigma_rw
        ),
        type = "message",
        duration = 6
      )
    }

    mu0 <- stan_data$mu_log_sigma_rw
    sd0 <- stan_data$sd_log_sigma_rw

    output$data_preview <- renderTable({
      head(df2 %>% select(StudentID, TimeID, Score, j, t, y), 20)
    }, digits = 3)

    if (is.null(stan_model_cached())) {
      stan_model_cached(rstan::stan_model(model_code = stan_code))
    }

    fit <- tryCatch(
      rstan::sampling(
        stan_model_cached(),
        data = stan_data,
        chains = as.integer(input$chains),
        iter = as.integer(input$iter),
        warmup = as.integer(input$warmup),
        seed = as.integer(input$seed),
        control = list(
          adapt_delta = as.numeric(input$adapt),
          max_treedepth = as.integer(input$treedepth)
        ),
        refresh = max(1, as.integer(input$iter) %/% 10)
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      showNotification(paste0("Stan sampling failed: ", conditionMessage(fit)), type = "error", duration = NULL)
      output$diag_text <- renderPrint({
        cat("Stan sampling failed.\n")
        cat(conditionMessage(fit), "\n")
      })
      return(NULL)
    }

    beta_mat <- as.matrix(fit, pars = "beta")
    Tn <- ncol(beta_mat)

    beta_df <- data.frame(
      t = seq_len(Tn),
      mean = colMeans(beta_mat),
      q05 = apply(beta_mat, 2, quantile, 0.05),
      q50 = apply(beta_mat, 2, quantile, 0.50),
      q95 = apply(beta_mat, 2, quantile, 0.95)
    )

    output$beta_plot <- renderPlot({
      ggplot(beta_df, aes(x = t, y = mean)) +
        geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
        geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        scale_x_continuous(breaks = beta_df$t) +
        labs(
          x = "Time point",
          y = "beta[t]",
          title = "Rater severity drift",
          subtitle = "Posterior mean with 90% credible interval; beta[1] is fixed to 0"
        ) +
        theme_minimal(base_size = 13)
    })

    output$beta_table <- renderTable({
      beta_df
    }, digits = 3)

    output$waicloo_text <- renderPrint({
      ll_mat <- as.matrix(fit, pars = "log_lik")

      waic_res <- tryCatch(loo::waic(ll_mat), error = function(e) e)
      loo_res <- tryCatch(loo::loo(ll_mat), error = function(e) e)

      cat("=== WAIC ===\n")
      if (inherits(waic_res, "error")) {
        cat("WAIC failed: ", conditionMessage(waic_res), "\n\n", sep = "")
      } else {
        print(waic_res$estimates, digits = 6)
        cat("\nNote: If p_waic warnings appear, prefer PSIS-LOO.\n\n")
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

    output$ppp_rmse_text <- renderPrint({
      Tobs <- as.vector(as.matrix(fit, pars = "T_obs"))
      Trep <- as.vector(as.matrix(fit, pars = "T_rep"))
      ppp <- mean(Trep > Tobs)

      rmse <- as.vector(as.matrix(fit, pars = "rmse"))
      rmse_rep <- as.vector(as.matrix(fit, pars = "rmse_rep"))

      cat(sprintf("Posterior predictive p-value: PPP = %.6f\n\n", ppp))

      cat("RMSE: E[y | parameters] vs observed y\n")
      print(quantile(rmse, probs = c(0.025, 0.5, 0.975)), digits = 6)
      cat(sprintf("mean = %.6f\n\n", mean(rmse)))

      cat("RMSE: replicated y_rep vs observed y\n")
      print(quantile(rmse_rep, probs = c(0.025, 0.5, 0.975)), digits = 6)
      cat(sprintf("mean = %.6f\n", mean(rmse_rep)))
    })

    output$sigma_prior_post_plot <- renderPlot({
      post <- as.vector(as.matrix(fit, pars = "sigma_rw"))

      x_max <- max(
        as.numeric(quantile(post, 0.995)),
        qlnorm(0.995, meanlog = mu0, sdlog = sd0)
      )

      x_max <- max(x_max, 1e-6)
      dpost <- density(post, from = 0, to = x_max, n = 400)

      df_plot <- data.frame(
        x = rep(dpost$x, 2),
        density = c(
          dlnorm(dpost$x, meanlog = mu0, sdlog = sd0),
          dpost$y
        ),
        distribution = rep(c("Prior", "Posterior"), each = length(dpost$x))
      )

      ggplot(df_plot, aes(x = x, y = density, linetype = distribution)) +
        geom_line(linewidth = 1) +
        labs(
          x = "sigma_rw",
          y = "Density",
          title = "Prior vs posterior for sigma_rw",
          subtitle = sprintf("Prior: LogNormal(%.4f, %.4f)", mu0, sd0)
        ) +
        theme_minimal(base_size = 13)
    })

    output$prior_post_table <- renderTable({
      post <- as.vector(as.matrix(fit, pars = "sigma_rw"))

      prior_mean <- exp(mu0 + (sd0^2) / 2)
      prior_sd <- sqrt((exp(sd0^2) - 1) * exp(2 * mu0 + sd0^2))
      prior_q <- qlnorm(c(0.05, 0.5, 0.95), meanlog = mu0, sdlog = sd0)
      post_q <- quantile(post, probs = c(0.05, 0.5, 0.95))

      data.frame(
        quantity = c("mean", "sd", "q05", "median", "q95"),
        prior = c(prior_mean, prior_sd, prior_q[1], prior_q[2], prior_q[3]),
        posterior = c(mean(post), sd(post), post_q[[1]], post_q[[2]], post_q[[3]])
      )
    }, digits = 6)

    output$diag_text <- renderPrint({
      cat("=== Prior hyperparameters passed to Stan ===\n")
      cat(sprintf("prior_mode = %s\n", input$prior_mode))
      cat(sprintf("mu_log_sigma_rw = %.6f\n", mu0))
      cat(sprintf("sd_log_sigma_rw = %.6f\n\n", sd0))

      cat("=== Parameter summary ===\n")
      print(fit, pars = c("mu_theta", "sigma_theta", "sigma_rw", "beta"), probs = c(0.05, 0.5, 0.95))

      sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      div <- sapply(sp, function(x) sum(x[, "divergent__"]))
      td <- sapply(sp, function(x) sum(x[, "treedepth__"] >= as.integer(input$treedepth)))

      cat("\nDivergences per chain: ", paste(div, collapse = ", "), " (total: ", sum(div), ")\n", sep = "")
      cat("Max-treedepth hits per chain: ", paste(td, collapse = ", "), " (total: ", sum(td), ")\n", sep = "")
    })
  })
}

shinyApp(ui, server)
