# Normalize the WordprocessingML properties that HTML/gt conversion emits.
# ECMA-376 defines a fixed child order, even when Word/LibreOffice tolerate
# some variations. Keep this independent of presentation and of model results.
normalize_word_xml <- function(document) {
  ns <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")
  # Quarto 1.9.37 can add an inline caption pPr after Pandoc's table-cell pPr.
  # Word allows only one pPr per paragraph. Combine them before ordering and
  # deduplication below, retaining later properties and the caption style.
  for (paragraph in xml2::xml_find_all(document, '//w:p[count(w:pPr) > 1]', ns)) {
    properties <- xml2::xml_find_all(paragraph, './w:pPr', ns)
    first <- properties[[1L]]
    for (extra in properties[-1L]) {
      for (child in xml2::xml_children(extra)) {
        xml2::xml_add_child(first, child, .copy = FALSE)
      }
      xml2::xml_remove(extra)
    }
  }
  orders <- list(
    pPr = "pStyle keepNext keepLines pageBreakBefore framePr widowControl numPr suppressLineNumbers pBdr shd tabs suppressAutoHyphens kinsoku wordWrap overflowPunct topLinePunct autoSpaceDE autoSpaceDN bidi adjustRightInd snapToGrid spacing ind contextualSpacing mirrorIndents suppressOverlap jc textDirection textAlignment textboxTightWrap outlineLvl divId cnfStyle rPr sectPr pPrChange",
    rPr = "rStyle rFonts b bCs i iCs caps smallCaps strike dstrike outline shadow emboss imprint noProof snapToGrid vanish webHidden color spacing w kern position sz szCs highlight u effect bdr shd fitText vertAlign rtl cs em lang eastAsianLayout specVanish oMath rPrChange",
    sectPr = "headerReference footerReference footnotePr endnotePr type pgSz pgMar paperSrc pgBorders lnNumType pgNumType cols formProt vAlign noEndnote titlePg textDirection bidi rtlGutter docGrid printerSettings sectPrChange",
    tblPr = "tblStyle tblpPr tblOverlap bidiVisual tblStyleRowBandSize tblStyleColBandSize tblW jc tblCellSpacing tblInd tblBorders shd tblLayout tblCellMar tblLook tblCaption tblDescription tblPrChange",
    tcPr = "cnfStyle tcW gridSpan hMerge vMerge tcBorders shd noWrap tcMar textDirection tcFitText vAlign hideMark headers cellIns cellDel cellMerge tcPrChange",
    tblBorders = "top start left bottom end right insideH insideV",
    tcBorders = "top start left bottom end right insideH insideV tl2br tr2bl",
    tblCellMar = "top start left bottom end right",
    tcMar = "top start left bottom end right",
    style = "name aliases basedOn next link autoRedefine hidden uiPriority semiHidden unhideWhenUsed qFormat locked personal personalCompose personalReply rsid pPr rPr tblPr trPr tcPr tblStylePr"
  )
  # gt's inline count notation can leave RTF control groups outside Word runs.
  # Translate their formatting while retaining every existing text run.
  for (paragraph in xml2::xml_find_all(document, '//w:p[text()[normalize-space()]]', ns)) {
    format <- NULL
    for (child in xml2::xml_contents(paragraph)) {
      if (xml2::xml_type(child) == "text") {
        value <- xml2::xml_text(child)
        if (!nzchar(trimws(value))) next
        format <- switch(value, "{\\i " = "italic", "}{\\sub " = "subscript", "}" = "plain", NULL)
        if (is.null(format)) stop("Unexpected bare text in a Word paragraph: ", value)
        xml2::xml_remove(child)
      } else if (xml2::xml_name(child) == "r" && !is.null(format) && format != "plain") {
        rp <- xml2::xml_find_first(child, './w:rPr', ns)
        if (inherits(rp, 'xml_missing')) rp <- xml2::xml_add_child(child, 'w:rPr', .where = 0)
        if (format == "italic") xml2::xml_add_child(rp, 'w:i')
        if (format == "subscript") xml2::xml_add_child(rp, 'w:vertAlign', 'w:val' = 'subscript')
      }
    }
  }
  # gt uses missing/zero border widths to mean no visible border.
  # Encode that intent explicitly instead of leaving an ambiguous single border.
  for (border in xml2::xml_find_all(document, '//w:tcBorders/* | //w:tblBorders/*', ns)) {
    width <- suppressWarnings(as.numeric(xml2::xml_attr(border, 'sz')))
    value <- xml2::xml_attr(border, 'val')
    if (!is.na(value) && value == 'single' && (is.na(width) || width <= 0)) {
      xml2::xml_set_attr(border, 'w:val', 'nil')
      xml2::xml_set_attr(border, 'w:sz', NULL, ns = ns)
    }
  }
  for (border in xml2::xml_find_all(document,
      '//w:tcBorders/*[@w:val="nil" or @w:val="none"] | //w:tblBorders/*[@w:val="nil" or @w:val="none"]', ns)) {
    for (attribute in c('sz', 'space', 'color')) xml2::xml_set_attr(border, paste0('w:', attribute), NULL, ns = ns)
  }
  # Give empty table cells an explicit paragraph for Word compatibility.
  # gt can omit this otherwise conventional cell content.
  for (cell in xml2::xml_find_all(document, '//w:tc[not(w:p) and not(w:tbl)]', ns)) {
    xml2::xml_add_child(cell, 'w:p')
  }
  # Imported table styles in the reference template need distinct display names.
  for (style in xml2::xml_find_all(document, '//w:style[not(w:name)]', ns)) {
    xml2::xml_add_child(style, 'w:name',
      'w:val' = paste('Table', xml2::xml_attr(style, 'styleId')), .where = 0)
  }
  for (tag in names(orders)) {
    ordering <- strsplit(orders[[tag]], ' ', fixed = TRUE)[[1L]]
    for (parent in xml2::xml_find_all(document, paste0('//w:', tag), ns)) {
      children <- xml2::xml_children(parent)
      repeated <- if (tag == "sectPr") c("headerReference", "footerReference") else if (tag == "style") "tblStylePr" else character()
      duplicate <- duplicated(xml2::xml_name(children), fromLast = TRUE) & !xml2::xml_name(children) %in% repeated
      xml2::xml_remove(children[duplicate])
      children <- xml2::xml_children(parent)
      index <- order(match(xml2::xml_name(children), ordering), na.last = TRUE)
      for (child in children[index]) xml2::xml_add_child(parent, child, .copy = FALSE)
    }
  }
  invisible(document)
}

