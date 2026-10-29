# Deployable copy of the depprep app (for example on a Shiny server); locally use
# depprep::run_app().
shiny::shinyApp(depprep:::.app_ui(), depprep:::.app_server)
