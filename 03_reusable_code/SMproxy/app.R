# Load required packages

library(shiny)
library(here)
library(ggplot2)
library(DT) # Add colour to output table

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
  
  # Side bar 
  
  sidebarLayout(
    
    sidebarPanel(
      
      # Create slider input for the user to modify the number of Monte Carlo runs
      
      sliderInput("MCnum", "Number of Monte Carlo runs",
                  min = 100,
                  max = 1000,
                  value = 200,
                  step = 1
      ),
      
      tags$hr(), # Horizontal line
      
      # Additional parameter input
      
      sliderInput("n", "Size of samples A and B",
                   min = 1000,
                   max = 10000,
                   value = 2000,
                   step = 100
      ),
      
      
      sliderInput("percExternal", "External sample proportion",
                   min = 0.1,
                   max = 0.9,
                   value = 0.3,
                   step = 0.1
      ),
      
      sliderInput("percOverlap", "Overlap proportion",
                   min = 0.1,
                   max = 0.5,
                   value = 0.2,
                   step = 0.1
      ),
      
      sliderInput("p_diag", "Proxy association strength",
                   min = 0.1,
                   max = 0.9,
                   value = 0.5,
                   step = 0.1
      ),
      
      selectInput("w", "Misclassification probabilities",
                  choices = c("Symmetrical", "Asymmetrical"),
                  selected = "Symmetrical"
      ),
      
      sliderInput("CIA", "Conditional independence assumption violation",
                   min = 0.0,
                   max = 1,
                   value = 0.3,
                   step = 0.1
      ),
      
      selectInput("Xsource", "Source of the marginal distribution of X",
                  choices = c("pop", "B"),
                  selected = "B"
                  ),
      
      tags$hr(), 
      
      # Create the action button
      
      actionButton("run", "Run the simulation")
      
    ),
    
    # Main panel
    
    mainPanel(
      
      wellPanel( # Ensure separation of output sections, when there is no simulation run yet
        h4("Output table; mean RMSE value per estimator"),
        DTOutput("tab"), # Output results table
      ),
      
      tags$hr(),
      
      wellPanel(
        h4("Output boxplots; RMSE distribution across the Monte Carlo runs"),
        plotOutput("rmseplot") # Output RMSE plot
      )
      
    )
    
  )
  
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
        mc = mcout,
        para = list(
          MCnum = as.integer(input$MCnum),
          n = as.integer(input$n),
          percExt = as.numeric(input$percExternal),
          percOverlap = as.numeric(input$percOverlap),
          p_diag = as.numeric(input$p_diag),
          w = input$w,
          CIA = as.numeric(input$CIA),
          Xsource = input$Xsource
      )
      )
    })
    
  }, ignoreInit = TRUE)   # ignoreInit prevents the function from running at app initialization
  
  # Construct results table
  
  output$tab <- renderDT({
    
    req(res()) # Ensure the results to be available
    mc <- res()$mc # Retrieve Monte Carlo results
    
    tabdat <- data.frame(
      Estimator = c("DRE", # Doubly robust estimator
                "IPF", # Iterative proportional fitting estimator
                "EXT" # External model based estimator
                ), 
      RMSE = c(mc$DRE$rmse_mean, # Mean RMSE value per estimator
               mc$IPF$rmse_mean,
               mc$EXT$rmse_mean
              )
      )
    
    # Rank estimator performance based on RMSE values
    
    tabdat$rank <- rank(tabdat$RMSE, ties.method = "first") 
    
    # Verbal performance indicator to enhance color blind accessibility
    
    tabdat$Performance <- c("Best", "Intermediate", "Worst")[tabdat$rank]
    
    datatable(
      tabdat[, c("Estimator", "RMSE", "Performance")],
      rownames = FALSE,
      options = list(
        dom = "t", # Show only the table
        ordering = FALSE # Disable column ordering
      )
    ) %>%
      formatRound("RMSE", 4) %>% # Round the RMSE values
      formatStyle("Performance",
        backgroundColor = styleEqual(
          c("Best", "Intermediate", "Worst"),
          c("green", "orange", "red") 
        )
      )
  }, options = list(dom = "t"))
  
  # Construct RMSE plot
  
  output$rmseplot <- renderPlot({
    
    req(res()) # Ensure the results to be available
    mc <- res()$mc # Retrieve Monte Carlo results
    m <- res()$para$MCnum # Retrieve number of Monte Carlo runs
    
    dre <- mc$DRE$rmse
    ipf <- mc$IPF$rmse
    ext <- mc$EXT$rmse
    
    m_use <- min(m, length(dre), length(ipf), length(ext)) # Ensure equal lengths for all estimators
    
    rmsedata <- data.frame(
      
      Iteration = c(seq_len(m_use), seq_len(m_use), seq_len(m_use)), # Iteration numbers for all three estimators
      
      Estimator = factor(
        c(rep("DRE", m_use), rep("IPF", m_use), rep("EXT", m_use)),
        levels = c("DRE", "IPF", "EXT")),
      
      RMSE = c(
        dre[seq_len(m_use)],
        ipf[seq_len(m_use)],
        ext[seq_len(m_use)]
      )
    )
      
      # Generate the plot
      
      ggplot(rmsedata, aes(x = Estimator, y = RMSE, fill = Estimator)) +
        geom_boxplot(fill = "pink") +
        labs(
          x = "Estimator",
          y = "RMSE"
        ) +
        theme_bw()
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
