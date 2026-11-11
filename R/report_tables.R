# Presentation only: scientific values, inclusion and row order stay upstream.
# Quarto owns each numbered caption through its tbl-* label and tbl-cap.
report_column_labels <- function(columns) {
  labels <- gsub('_', ' ', columns, fixed = TRUE)
  labels <- gsub('R2', 'R²', labels, fixed = TRUE)
  labels <- gsub('participant days', 'participant-days', labels, fixed = TRUE)
  labels <- gsub('half hours', 'valid half-hours', labels, fixed = TRUE)
  labels <- gsub('extremity', 'chronotype strength', labels, fixed = TRUE)
  labels[columns == 'direction'] <- 'chronotype direction'
  labels <- sub('^[.]', '', labels)
  labels <- paste0(toupper(substr(labels, 1, 1)), substring(labels, 2))
  stats::setNames(labels, columns)
}

show_table <- function(x, digits = 3, labels = NULL, stub = NULL,
                       integer_cols = NULL, missing_text = 'NA') {
  stopifnot(is.data.frame(x), !anyDuplicated(names(x)),
            length(digits) == 1L, is.numeric(digits), is.finite(digits),
            digits >= 0, digits == as.integer(digits),
            length(missing_text) == 1L, is.character(missing_text),
            !is.na(missing_text))
  if (!is.null(stub)) {
    stopifnot(is.character(stub), length(stub) == 1L, stub %in% names(x))
  }
  display_labels <- report_column_labels(names(x))
  if (!is.null(labels)) {
    stopifnot(is.character(labels), !is.null(names(labels)),
              !anyNA(labels), !anyDuplicated(names(labels)),
              all(nzchar(names(labels))), all(names(labels) %in% names(x)))
    display_labels[names(labels)] <- labels
  }
  numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
  # Read.csv and the analysis count helpers retain integer counts. Explicit
  # integer columns cover double-storage counts without guessing from values.
  if (is.null(integer_cols)) {
    integer_cols <- names(x)[vapply(x, is.integer, logical(1))]
  } else {
    stopifnot(is.character(integer_cols), all(integer_cols %in% numeric_cols))
    stopifnot(all(vapply(x[integer_cols], function(z) {
      all(is.na(z) | (is.finite(z) & z == trunc(z)))
    }, logical(1))))
  }
  # Disable inherited grouped-data behavior so gt never gathers/reorders rows.
  out <- gt::gt(x, rowname_col = stub, groupname_col = NULL,
                auto_align = FALSE, locale = 'en') |>
    gt::cols_label(.list = as.list(display_labels)) |>
    gt::cols_align(align = 'left', columns = tidyselect::everything())
  if (length(numeric_cols)) {
    out <- out |>
      gt::fmt_number(columns = tidyselect::all_of(numeric_cols),
                     decimals = digits, use_seps = TRUE,
                     drop_trailing_zeros = FALSE, locale = 'en') |>
      gt::cols_align(align = 'right', columns = tidyselect::all_of(numeric_cols))
  }
  if (length(integer_cols)) {
    out <- out |>
      gt::fmt_integer(columns = tidyselect::all_of(integer_cols),
                      use_seps = TRUE, locale = 'en')
  }
  out <- out |>
    gt::sub_missing(columns = tidyselect::everything(), missing_text = missing_text)
  # Reader-facing names only; raw values and exported identifiers stay intact.
  text_cols <- names(x)[vapply(x, function(z) is.character(z) || is.factor(z), logical(1))]
  if (length(text_cols)) {
    out <- gt::text_case_match(out,
      "primary" ~ "Near-eye", "chest" ~ "Chest",
      "direction" ~ "Chronotype direction", "extremity" ~ "Chronotype strength",
      "Direction by extremity" ~ "Chronotype direction by chronotype strength",
      "Direction by day type" ~ "Chronotype direction by day type",
      "Extremity by day type" ~ "Chronotype strength by day type",
      "Direction by extremity by day type" ~ "Chronotype direction by chronotype strength by day type",
      .locations = gt::cells_body(columns = tidyselect::all_of(text_cols)))
  }
  if (length(text_cols)) {
    out <- gt::text_transform(out,
      locations = gt::cells_body(columns = tidyselect::all_of(text_cols)),
      fn = function(z) gsub("dMSFsc", "MSF<sub>sc(diary)</sub>", z, fixed = TRUE))
  }
  if (!is.null(stub)) {
    out <- gt::tab_stubhead(out, label = unname(display_labels[stub]))
  }
  out |>
    gt::tab_options(
      table.width = gt::pct(100), table.layout = 'auto', table.align = 'left',
      table.font.names = 'inherit', table.font.size = '0.9rem',
      data_row.padding = gt::px(6), column_labels.font.weight = '600',
      container.width = gt::pct(100), container.overflow.x = 'auto',
      container.overflow.y = 'visible',
      quarto.disable_processing = FALSE, quarto.use_bootstrap = FALSE
    )
}
