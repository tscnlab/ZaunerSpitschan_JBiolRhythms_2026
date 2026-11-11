#!/bin/sh
# Only Word output needs the local R-based typography and OOXML normalization.
case "${QUARTO_PROJECT_OUTPUT_FILES:-}" in
  *chronotype-manuscript.docx*) Rscript scripts/manuscript-word-layout.R ;;
esac
