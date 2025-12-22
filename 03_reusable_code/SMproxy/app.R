# Load required packages

library(shiny)
library(here)
library(tidyverse)
library(ggplot2)

# Load external scripts

source("../Thesis scripts/02_utils.R") # Load utility functions 
source("../Thesis scripts/05_generate_population.R") # Load SIM function
source("../Thesis scripts/06_perform_mc.R") # Load Monte Carlo function

# Define the input parameters

inpar <- function(input){
  
  # These parameters are not directly required for the simulation code
  # Nevertheless they are included because the user of the app is able to modify these parameters
  
  n = as.integer(input$n)
  percExt = as.numeric(input$percExternal)
  
  list(
    
    # These parameters are directly required for the simulation code
    
    nA = n,
    nB = n,
    nE = round(n * percExt),
    
    PercOverlap = as.numeric(input$percOverlap),
    tranmat_diag = as.numeric(input$p_diag),
    tranmat_sym = (input$w == "Symmetrical"),
    cia =  as.numeric(input$CIA),
    
    N = 1e6,
    zpx = 0.4,
    Xsource = input$Xsource
  )
}

ui <- fluidPage(
  
  # Application title
  
  titlePanel("Exploration of Statistical Matching with a Proxy Variable"),
  
  # Create action button
  
  actionButton("run", "Run the simulation"),
  
  # Create slider input for the user to modify the number of Monte Carlo runs
  
  sliderInput("MCnum", "Number of Monte Carlo runs",
              min = 1,
              max = 10,
              value = 3,
              step = 1
              ),
  
  tableOutput("tab"), # Output results table
  plotOutput("rmse_plot") # Output RMSE plot

)

server <- function(input, output) {
  
  # Retrieve the results
  
  res <- eventReactive(input$run, {   # eventReactive to run the simulation only when user presses the action button
    
    withProgress(message = "Running the simulation...", value = 0, { # Initalize progression bar
      
      incProgress(0.1, detail = "Inserting the parameters")
      param <- inpar(input)
      
      incProgress(0.4, detail = "Generating the population")
      simout <- do.call(SIM, param)
      
      incProgress(0.4, detail = "Performing the Monte Carlo simulation")
      mcout <- MC(
        SIMout = simout,
        MCnum = as.integer(input$MCnum), # Override the default number of Monte Carlo runs
        seed = 1L
      )
      
      incProgress(0.1, detail = "Finishing touch")
      list(
        sim = simout,
        mc = mcout
      )
      
    })
    
  }, ignoreInit = TRUE)   # ignoreInit prevents the function from running at app initialization
  
  # Construct results table
  
  output$tab <- renderTable({
    
    req(res()) # Ensure the results to be available
    mc <- res()$mc # Retrieve Monte Carlo results
    
    data.frame(
      Estimator = c("DRE", # Doubly robust estimator
                "IPF", # Iterative proportional fitting estimator
                "EXT" # External model based estimator
                ), 
      RMSE = c(mc$DRE$rmse_mean, # Mean RMSE value per estimator
               mc$IPF$rmse_mean,
               mc$EXT$rmse_mean
              )
      )
    }, digits = 4)
  
  # Construct RMSE plot
  
  output$rmseplot <- renderPlot({
    
    req(res()) # Ensure the results to be available
    mc <- res()$mc # Retrieve Monte Carlo results
    
    rmsedata <- data.frame(
      
      Iteration = rep(1:input$MCnum, 3), # Iteration numbers for all three estimators
      
      Estimator = factor(rep(c("DRE", "IPF", "EXT"), 
                             levels = c("DRE", "IPF", "EXT"))), 
      RMSEmean = c(
        mc$DRE$rmse,
        mc$IPF$rmse,
        mc$EXT$rmse
      ),
      
      # Generate the plot
      
      ggplot(rmsedata, aes(x = Estimator, y = RMSE)) +
        geom_col(fill = "pink") +
        labs(
          title = "Average RMSE value per estimator",
          x = "Estimator",
          y = "RMSE"
        ) +
        theme_minimal()
      
    )
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
