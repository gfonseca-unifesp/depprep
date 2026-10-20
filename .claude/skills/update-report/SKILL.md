---
name: update-report
description: Updates interno/report/relatorio.qmd with the section of one stage of the depprep project and renders the HTML. Use when a stage is finished or when the user asks to update the report (e.g. "/update-report etapa 3").
---

The report is the maintainer's private follow-up document and is written in Portuguese.
It lives in interno/ (a separate private repository); never copy it into the public one.

1. Read interno/docs/protocol.md, interno/report/decisoes.md and the files in results/ for the
   given stage.
2. Update ONLY that stage's section of interno/report/relatorio.qmd, following the template:
   Status, Objetivo, Definições (term + short definition + concrete example),
   Exemplo didático (simple figure), Resultados (only from results/, through code),
   Interpretação (what it means and what is uncertain), Decisão pendente
   (question, 2-3 options, consequence of each, recommendation), Checklist de validação.
3. Add new terms to the Glossário and new rows to the Registro de decisões.
4. Write in clear Portuguese for a reader from a neighbouring field. Short sentences.
   Every figure has a caption that says what to look at.
5. Never type result numbers in the text: use inline R (`r valor`).
6. Run: quarto render interno/report/relatorio.qmd (see CLAUDE.md for the quarto path
   and QUARTO_R). Fix errors until it renders.
7. Set the stage status to "Aguardando validação" and list the pending decisions for the
   user. Do not mark any decision as taken.
8. Commit in the private repository: git -C interno commit, message
   "Atualiza relatório: etapa N".
