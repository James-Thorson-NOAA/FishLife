library(shiny)

library(FishLife)
library(rfishbase)

# Function containing things to display
function(input, output, session){

  #### Dynamic user inputs
  # The following reactive function returns orders from a selected class
  order_subset <- reactive({ sort(unique(fishbase$Order[which(fishbase$Class==input$Class)])) })
  output$orderSelex <- renderUI({
    selectInput(inputId="Order", label="Taxonomic order", choices=order_subset(), multiple=FALSE, selected=order_subset()[1])
  })
  # The following reactive function returns family from a selected order
  #family_subset <- reactive({ sort(unique(Estimate_database$Z_ik[which(Estimate_database$Z_ik[,'Class']==input$Class & Estimate_database$Z_ik[,'Order']==input$Order),'Family'])) })
  family_subset <- reactive({ sort(unique(fishbase$Family[which(fishbase$Class==input$Class & fishbase$Order==input$Order)])) })
  output$familySelex <- renderUI({
    selectInput(inputId="Family", label="Family", choices=family_subset(), multiple=FALSE, selected=family_subset()[1])
  })
  # The following reactive function returns genus from a selected fanily
  #genus_subset <- reactive({ sort(unique(Estimate_database$Z_ik[which(Estimate_database$Z_ik[,'Order']==input$Order & Estimate_database$Z_ik[,'Family']==input$Family),'Genus'])) })
  genus_subset <- reactive({ sort(unique(fishbase$Genus[which(fishbase$Order==input$Order & fishbase$Family==input$Family)])) })
  output$genusSelex <- renderUI({
    selectInput(inputId="Genus", label="Genus", choices=genus_subset(), multiple=FALSE, selected=genus_subset()[1])
  })
  # The following reactive function returns species from a selected genus
  #species_subset <- reactive({ sort(unique(Estimate_database$Z_ik[which(Estimate_database$Z_ik[,'Family']==input$Family & Estimate_database$Z_ik[,'Genus']==input$Genus),'Species'])) })
  species_subset <- reactive({ sort(unique(fishbase$Species[which(fishbase$Family==input$Family & fishbase$Genus==input$Genus)])) })
  output$speciesSelex <- renderUI({
    selectInput(inputId="Species", label="Species", choices=species_subset(), multiple=FALSE, selected=species_subset()[1])
  })

  # Match species from text
  Match_taxonomy <- reactive({
    Search_species(Class=input$Class, Order=input$Order, Family=input$Family, Genus=input$Genus, Species=input$Species, add_ancestors=input$plotAncestors, ParentChild_gz=Estimate_database$ParentChild_gz)$match_taxonomy
  })

  # Disply species
  output$debug_text1 <- renderPrint({
    paste0( c(input$Class,input$Order,input$Family,input$Genus,input$Species), collapse="_" )
  })
  output$debug_text2 <- renderPrint({
    Match_taxonomy()
  })

  #### Plots
  # Plot taxa
  output$plot1 <- renderPlot({
    #input$activate
    #isolate({
    Plot_taxa( Taxa=Match_taxonomy(), Cov_gjj=Estimate_database$Cov_gjj, Mean_gj=Estimate_database$ParHat$beta_gj, ParentChild_gz=Estimate_database$ParentChild_gz, Y_ij=Estimate_database$Y_ij )
    #})
  })

}
