# Manuscript citation style

`american-medical-association.csl` is based on the **AMA Manual of Style 11th edition** style from the [Citation Style Language repository](https://github.com/citation-style-language/styles/blob/master/american-medical-association.csl), retrieved 29 October 2026. A local adjustment added on 11 November 2026 displays an explicitly supplied publication genre, such as `Preprint`, instead of `Published online` for unpaginated articles. The file retains its contributor credits and Creative Commons Attribution-ShareAlike 3.0 license.

The [current JBR instructions](https://journals.sagepub.com/author-instructions/jbr) specify AMA style. This CSL is applied to `index.qmd`; the analytical supplement retains its existing author-date citations.

# Word export

`melidos-reference.docx` is an unchanged copy of `assets/reference.docx` from the main MeLiDos paper (`ZaunerEtAl_bioRxiv_2026`), copied on 7 November 2026. Its styles, A4 page size, margins, line numbering and page-number footer are used for the manuscript Word export. The authors-block Lua helpers come from the same project and retain their license notices. Source paths and the template checksum are recorded in `melidos-word-provenance.json`.

Run `R_PROFILE_USER=.Rprofile quarto render index.qmd --to docx` from the project root to create `report/chronotype-manuscript.docx`. The Word-only filter positions the existing abstract and keywords after the affiliations. The post-render layout script follows the main paper's heading, title-page and caption conventions, keeps images within portrait pages and verifies that it preserves all text and figures. The manuscript retains its JBR citation style.
