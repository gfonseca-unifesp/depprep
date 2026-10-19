# depprep — project instructions

## What it is
Preprocessing that corrects spatial, temporal and spatio-temporal autocorrelation in data
sets with many variables of different dependence regimes, while preserving the
correlations between variables. Output: a corrected matrix for PCA, MDS, k-means and SOM.

Reference documents (Portuguese, private): interno/docs/projeto_revisado.md,
interno/docs/guia_execucao.md, interno/docs/parecer_projeto_original.md.

## Primary estimand
Correlation between the residual components of the variables after conditioning on the
large-scale structure of their group. Core rule: variables in the same group use the SAME
basis (projection P_g); pairs across groups use the union basis; a group with no
dependence has no basis.

## Rules
- Never modify data/raw. Every analysis comes from analysis/NN_*.R and writes to results/.
- Every simulation calls set.seed() explicitly and records the seed in its output.
- Dependence signatures use only marginal information of each variable
  (variogram, ACF), never cross-correlations.
- Every new function has a test in tests/testthat before the commit. Run devtools::test().
- Do not change simulation parameters frozen in docs/protocol.md without asking.
- Small commits, messages in English in the imperative mood.

## Two repositories and language
- This folder is the PUBLIC repository (github.com/gfonseca-unifesp/depprep). Everything
  in it is in English: code, function/file/folder names, comments, tests, error messages,
  roxygen docs, README, DESCRIPTION, analysis/ scripts, outputs and figures in results/,
  docs/protocol.md, docs/glossary.md, docs/literature.csv and commit messages.
- interno/ is ignored here and is a separate PRIVATE repository
  (github.com/gfonseca-unifesp/depprep_int), in Portuguese, for the maintainer's own
  follow-up: the progress report, the decision log and the planning documents.
  Commit there with `git -C interno ...`; messages there may be in Portuguese.
- Never copy Portuguese material from interno/ into the public repository.

## Progress report (interno/, in Portuguese)
- Source: interno/report/relatorio.qmd -> interno/report/relatorio.html (embed-resources: true).
- Visual and content model: interno/report/relatorio_v0.html.
- Numbers and figures ALWAYS come from results/ through code; never type results by hand.
- Each stage has: objective, definitions, didactic example, results, interpretation,
  pending decision (options + consequences) and status.
- Validated decisions go to interno/report/decisoes.md.
- Render with: quarto render interno/report/relatorio.qmd

## Commands
- Tests: Rscript -e 'devtools::test()'
- Check: Rscript -e 'devtools::check()'
- Note (this machine): Rscript and quarto are not on PATH. Use the full paths:
  - Rscript: "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" (R 4.5.1 also exists; use 4.6.1)
  - quarto (bundled with RStudio): "C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe"
  - quarto needs QUARTO_R="C:/Program Files/R/R-4.6.1/bin" to find R.
