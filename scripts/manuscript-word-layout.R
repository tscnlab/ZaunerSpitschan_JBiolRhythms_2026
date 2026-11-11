# Word-only layout matching the main MeLiDos manuscript.
# This changes presentation properties, never manuscript text or figure content.
source(".Rprofile")
source("scripts/word-xml-normalize.R")

format_manuscript_word <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  work <- tempfile("chronotype-word-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = work)
  document_path <- file.path(work, "word/document.xml")
  document <- xml2::read_xml(document_path)
  # Consolidate Quarto's duplicate caption properties before applying styles.
  normalize_word_xml(document)
  ns <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")
  find <- function(x, xpath) xml2::xml_find_all(x, xpath, ns)
  first <- function(x, xpath) xml2::xml_find_first(x, xpath, ns)
  missing <- function(x) inherits(x, "xml_missing")
  node <- function(name, attributes = "") xml2::xml_root(xml2::read_xml(paste0(
    '<w:', name, ' xmlns:w="', ns[["w"]], '" ', attributes, '/>'
  )))
  ensure <- function(parent, name) {
    out <- first(parent, paste0("./w:", name))
    if (missing(out)) out <- xml2::xml_add_child(parent, node(name), .where = 0)
    out
  }
  property <- function(parent, name, attributes = "") {
    xml2::xml_remove(find(parent, paste0("./w:", name)))
    xml2::xml_add_child(parent, node(name, attributes))
  }
  before_text <- xml2::xml_text(find(document, "//w:t"))
  before_pictures <- length(find(document, "//w:drawing"))
  body <- first(document, "//w:body")

  # Match the main paper's section starts and heading spacing.
  major <- c("Abstract", "Introduction", "Methods", "Results", "Discussion",
             "Acknowledgements", "References")
  title_page <- TRUE
  abstract_page <- FALSE
  for (p in find(body, "./w:p")) {
    label <- trimws(xml2::xml_text(p))
    if (label == "Abstract") {
      title_page <- FALSE
      abstract_page <- TRUE
    }
    if (label == "Introduction") abstract_page <- FALSE
    style <- xml2::xml_attr(first(p, "./w:pPr/w:pStyle"), "val")
    pp <- ensure(p, "pPr")
    if (!is.na(style) && grepl("^(Heading|berschrift)", style)) {
      property(pp, "keepNext")
      property(pp, "spacing", 'w:before="240" w:after="160" w:line="240" w:lineRule="auto"')
      if (label %in% major) property(pp, "pageBreakBefore")
    }
    if (abstract_page && label != "Abstract") {
      property(pp, "spacing", 'w:before="0" w:after="120" w:line="300" w:lineRule="auto"')
    }
    if (title_page) {
      property(pp, "spacing", 'w:before="0" w:after="80" w:line="300" w:lineRule="auto"')
      property(pp, "suppressLineNumbers")
      size <- if (!is.na(style) && style %in% c("Title", "Titel")) 36L else
        if (!is.na(style) && style %in% c("Author", "FirstParagraph")) 24L else 22L
      for (run in find(p, ".//w:r")) {
        property(ensure(run, "rPr"), "sz", paste0('w:val="', size, '"'))
      }
    }
  }

  # Preserve image proportions; large figure captions may continue on a new page.
  drawing_index <- 0L
  for (drawing in find(document, "//w:drawing")) {
    drawing_index <- drawing_index + 1L
    extents <- xml2::xml_find_all(drawing, ".//*[local-name()='extent' or local-name()='ext'][@cx and @cy]")
    if (!length(extents)) next
    width <- as.numeric(xml2::xml_attr(extents[[1L]], "cx"))
    height <- as.numeric(xml2::xml_attr(extents[[1L]], "cy"))
    # Figures 3–5 fill the usable line width. Word's inline rendering adds
    # an inset of about 9 pt; reserve 10 pt so the rightmost labels remain visible.
    scale <- if (drawing_index >= 3L) ((9021 - 200) / 1440) * 914400 / width else
      min(1, 6.05 * 914400 / width, 6.2 * 914400 / height)
    xml2::xml_set_attr(extents, "cx", as.character(round(width * scale)))
    xml2::xml_set_attr(extents, "cy", as.character(round(height * scale)))
  }
  for (table in find(body, "./w:tbl[.//w:drawing]")) {
    # Match the reference document's 9021-twip text area explicitly. Quarto's
    # default 7920-twip grid can clip wide images despite a 100% table width.
    properties <- ensure(table, "tblPr")
    property(properties, "tblW", 'w:type="dxa" w:w="9021"')
    property(properties, "tblInd", 'w:type="dxa" w:w="0"')
    margins <- ensure(properties, "tblCellMar")
    for (side in c("top", "left", "bottom", "right")) {
      property(margins, side, 'w:type="dxa" w:w="0"')
    }
    xml2::xml_remove(find(table, "./w:tblGrid"))
    grid <- xml2::xml_add_sibling(properties, node("tblGrid"), .where = "after")
    xml2::xml_add_child(grid, node("gridCol", 'w:w="9021"'))
    for (cell in find(table, "./w:tr/w:tc")) {
      cell_properties <- ensure(cell, "tcPr")
      property(cell_properties, "tcW", 'w:type="dxa" w:w="9021"')
      cell_margins <- ensure(cell_properties, "tcMar")
      for (side in c("top", "left", "bottom", "right")) {
        property(cell_margins, side, 'w:type="dxa" w:w="0"')
      }
    }
    for (row in find(table, "./w:tr")) property(ensure(row, "trPr"), "cantSplit")
    for (p in find(table, ".//w:p")) {
      pp <- ensure(p, "pPr")
      property(pp, "spacing", 'w:before="0" w:after="0" w:line="240" w:lineRule="auto"')
      property(pp, "suppressLineNumbers")
    }
  }
  for (p in find(document, "//w:p[.//w:drawing]")) {
    pp <- ensure(p, "pPr")
    xml2::xml_remove(find(pp, "./w:keepNext"))
    property(pp, "jc", 'w:val="center"')
    property(pp, "ind", 'w:left="0" w:right="0" w:firstLine="0"')
    property(pp, "suppressLineNumbers")
  }
  for (p in find(document, '//w:p[w:pPr/w:pStyle[@w:val="ImageCaption"]]')) {
    pp <- ensure(p, "pPr")
    xml2::xml_remove(find(pp, "./w:keepLines"))
    property(pp, "suppressLineNumbers")
    property(pp, "spacing", 'w:before="80" w:after="120" w:line="240" w:lineRule="auto"')
    for (run in find(p, ".//w:r")) property(ensure(run, "rPr"), "sz", 'w:val="20"')
  }

  styles_path <- file.path(work, "word/styles.xml")
  styles <- xml2::read_xml(styles_path)
  style_nodes <- find(styles, "//w:style")
  ids <- xml2::xml_attr(style_nodes, "styleId")
  names <- vapply(style_nodes, function(s) xml2::xml_attr(first(s, "./w:name"), "val"), character(1))
  for (reference in find(document, "//w:pStyle | //w:rStyle | //w:tblStyle")) {
    id <- xml2::xml_attr(reference, "val")
    if (id %in% ids) next
    position <- match(tolower(id), tolower(names))
    if (!is.na(position)) xml2::xml_set_attr(reference, "w:val", ids[[position]])
  }
  for (id in c("Author", "Date")) if (!id %in% ids) {
    style <- xml2::xml_add_child(styles, node("style", paste0('w:type="paragraph" w:customStyle="1" w:styleId="', id, '"')))
    property(style, "name", paste0('w:val="', id, '"'))
    normal <- xml2::xml_attr(first(styles, '//w:style[@w:type="paragraph" and @w:default="1"]'), "styleId")
    property(style, "basedOn", paste0('w:val="', normal, '"'))
  }
  for (footer_path in list.files(file.path(work, "word"), pattern = "^footer[0-9]+[.]xml$", full.names = TRUE)) {
    footer <- xml2::read_xml(footer_path)
    for (p in find(footer, "//w:p")) property(ensure(p, "pPr"), "suppressLineNumbers")
    normalize_word_xml(footer)
    xml2::write_xml(footer, footer_path)
  }
  normalize_word_xml(document)
  normalize_word_xml(styles)
  stopifnot(identical(before_text, xml2::xml_text(find(document, "//w:t"))),
            before_pictures == length(find(document, "//w:drawing")))
  xml2::write_xml(document, document_path)
  xml2::write_xml(styles, styles_path)
  output <- tempfile("chronotype-manuscript-", tmpdir = normalizePath(dirname(path)), fileext = ".docx")
  on.exit(unlink(output), add = TRUE)
  zip::zipr(output, list.files(work, all.files = TRUE, no.. = TRUE), root = work,
            include_directories = FALSE)
  stopifnot(file.rename(output, path))
  message("MeLiDos Word style applied; text and ", before_pictures, " figures preserved.")
  invisible(TRUE)
}

paths <- commandArgs(trailingOnly = TRUE)
if (!length(paths)) paths <- strsplit(Sys.getenv("QUARTO_PROJECT_OUTPUT_FILES"), "\n", fixed = TRUE)[[1L]]
for (path in paths[grepl("[.]docx$", paths, ignore.case = TRUE)]) format_manuscript_word(path)
