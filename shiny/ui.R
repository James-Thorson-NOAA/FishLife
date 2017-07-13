library(shiny)

library(FishLife)
Estimate_database = Load_previous_results()

# Load stuff
# Page for user interface
fluidPage(

  titlePanel("Visualize fish traits"),

  # Panel for settings
  sidebarPanel(
    # Display useful info
    h1("Background"),
    h4("This page visualizes predictions of life history traits"),
    br(),

    # Time series plots
    h1("Select a taxon"),
    # Choose region
    #textInput(inputId="Class", label="Taxonomic class", value = "Predictive"),
    #textInput(inputId="Order", label="Taxonomic order", value = "Predictive"),
    #textInput(inputId="Family", label="Taxonomic family", value = "Predictive"),
    #textInput(inputId="Genus", label="Taxonomic genus", value = "Predictive"),
    #textInput(inputId="Species", label="Taxonomic species", value = "Predictive")
    selectInput(inputId="Class", label="Taxonomic class", choices=sort(unique(Estimate_database$Z_ik[,'Class'])), multiple=FALSE, selected="Actinopterygii"),
    uiOutput("orderSelex"),
    uiOutput("familySelex"),
    uiOutput("genusSelex"),
    uiOutput("speciesSelex"),
    checkboxInput( inputId="plotAncestors", label="Plot ancestors for taxon?", value=TRUE)
  ),

  # Configuration of plotting tabs
  mainPanel(
    tabsetPanel(
      # Bivariate traits tab
      tabPanel("Traits: bivariate",
        #textOutput('debug_text1'),
        #textOutput('debug_text2'),
        plotOutput('plot1', height="800px")    #
      )
    )
  )
)
