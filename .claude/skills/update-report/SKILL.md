---
name: update-report
description: Updates interno/report/relatorio.qmd (three layers: panel, current results by question, history) with the results of one stage of the depprep project and renders the HTML. Use when a stage is finished or when the user asks to update the report (e.g. "/update-report etapa 7").
---

The report is the maintainer's private follow-up document and is written in Portuguese.
It lives in interno/ (a separate private repository); never copy it into the public one.
Format approved by the maintainer on 2026-10-23 (preview in
interno/report/formato_aprovado/): three layers, in this order.

1. **Panel (Camada 1, chunk `abertura`).** What the reader sees first:
   - the stage strip, from the `etapas` data frame in the `setup` chunk;
   - key figures of the current stage;
   - "O que mudou desde a última versão": automatic, from git, since the last commit that
     changed report/relatorio.html;
   - the open decisions, as cards linked to their decision cards;
   - the hypothesis scoreboard (placar), one row per pre-registered hypothesis;
   - one decision card per open decision, built with `cartao()`.

   A card has these fields: what the data say, options with consequences (one marked
   recommended), recommendation, "mudaria se", and the source files.
2. **Current results by question (Camada 2, chunk `camada2`).** One block per question
   (H1, H2, H3 with R and X, H4, Q3, then descriptive comparisons and robustness). Each
   block opens with a "Resposta curta" box that says:
   - the answer;
   - whether it holds for every correction and technique.

   Then come:
   - a "Como ler" guide;
   - the figure or table;
   - an "origem:" line naming the results/ files.

   Only the current version of a result appears here. Amended results (e.g.
   protocol-v1.1) appear beside the pre-registered ones, with the label; they never
   replace them.
3. **History (Camada 3).** Everything closed in `<details class="hist">` blocks:
   - the timeline (`linha`);
   - the fundamentals;
   - one block per stage narrative (`abre_etapa(n)` … `html('</details>')`);
   - the superseded runs (`subst`);
   - the decision register;
   - the glossary.

   A link into a closed block opens it (assets/relatorio-js.html).

To add or close a stage:

1. Read interno/docs/protocol.md, interno/report/decisoes.md and the stage's files in
   results/.
2. Add a `calcN` chunk that reads every number from results/ (never type a result in the
   text; build sentences with the computed values).
3. Update the panel:
   - the stage status in `etapas`;
   - the key figures;
   - the scoreboard;
   - the cards for the stage's open decisions.

   A new pending decision first gets a row in decisoes.md with an empty choice and
   "Pendente".
4. Update the results layer:
   - replace the affected question blocks with the new version;
   - move what they replace to "Rodadas substituídas" in the history (folder, what it
     was, why it changed).
5. Add the stage narrative to the history with the stage template:
   - Objetivo, Passo a passo with numbered steps (N.1, N.2…), Definições;
   - one `.sub` per step, linking to the question blocks for the results;
   - Checklist de validação.

   Add the steps to `linha`.
6. Add new terms to the glossary (term, group, definition, example).
7. Write in clear Portuguese for a reader from a neighbouring field. Short sentences.
   Every figure and table has a "Como ler" note and an "origem:" line.
8. Render from the project root with quarto (paths and QUARTO_R in CLAUDE.md). Fix errors
   until it renders, then check the page in a browser.
9. Set the stage status to "Aguardando validação" and list the pending decisions for the
   user. Do not mark any decision as taken.
10. Commit in the private repository: `git -C interno commit`, message
    "Atualiza relatório: etapa N".
