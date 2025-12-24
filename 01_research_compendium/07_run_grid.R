# The RUN function executes the simulation study over a grid of parameter settings
# The SIM and MC functions are requisiite for this function to work

RUN <- function(fungrid, MCnum = 1000) { 
  
  outlist <- vector("list", nrow(fungrid))  # Create an empty list to store output per row
  pb <- txtProgressBar(min = 0, max = nrow(fungrid), style = 3)  # Initialize a progress bar
  
  for (i in seq_len(nrow(fungrid))) {  
    
    g <- fungrid[i, , drop = FALSE]  # Select row of the parameter grid
    
    # Run SIM 
    
    S <- SIM( 
      N = g$N, # Population size
      nA = g$nA, # Sample A size          
      nB = g$nB, # Sample B size     
      nE = g$nE, # External sample size      
      PercOverlap = g$PercOverlap, # Percentage overlap between samples A and B
      tranmat_diag = g$tranmat_diag, # Diagonal elements of the transition matrix
      tranmat_sym = g$tranmat_sym, # Whether the transition matrix is symmetric # WAS VOORHEEN isTRUE
      seed = g$ID, # Random seed for reproducibility       
      zpx = g$zpx, # Dependency of Z on P and X        
      cia = g$cia, # Conditional independence assumption parameter    
      
      Xsource = if ("Xsource" %in% names(g)) as.character(g$Xsource) else "B", # Source of X variable
      mechanism = if ("mechanism" %in% names(g)) as.character(pull(g$mechanism)) else "MAR", # Selectivity mechanism
      pattern = if ("pattern" %in% names(g)) as.character(pull(g$pattern)) else NA_character_, # Selectivity pattern
      degree = if ("degree" %in% names(g)) as.character(pull(g$degree)) else NA_character_ # Selectivity degree
    )
    
    # Run Monte Carlo
    
    MCsim <- MC(S, MCnum = MCnum, seed = g$ID) # Monte Carlo simulation 
    
    # Store results
    
    outlist[[i]] <- data.frame(
      
      ID = as.integer(pull(g$ID)), # Simulation ID
      N = as.integer(pull(g$N)), # Population size
      nA = as.integer(pull(g$nA)), # Sample A size
      nB = as.integer(pull(g$nB)), # Sample B size
      nE = as.integer(pull(g$nE)), # External sample size
      PercOverlap = as.numeric(pull(g$PercOverlap)), # Percentage overlap between samples A and B
      tranmat_diag = as.numeric(pull(g$tranmat_diag)), # Diagonal elements of the transition matrix
      tranmat_sym = g$tranmat_sym, # Whether the transition matrix is symmetric
      cia = as.numeric(pull(g$cia)), # Conditional independence assumption parameter
      zpx = as.numeric(pull(g$zpx)), # Dependency of Z on P and X
      
      Xsource = if ("Xsource" %in% names(g)) as.character(pull(g$Xsource)) else "B", # Source of X variable
      mechanism = if ("mechanism" %in% names(g)) as.character(pull(g$mechanism)) else "MAR", # Selectivity mechanism
      pattern = if ("pattern" %in% names(g)) as.character(pull(g$pattern)) else NA_character_, # Selectivity pattern
      degree = if ("degree" %in% names(g)) as.character(pull(g$degree)) else NA_character_, # Selectivity degree
      
      # Doubly robust estimator results
      DRE_totbias_mean = as.numeric(MCsim$DRE$totbias_mean)[1], # Mean of total bias corresponding to the DRE estimate
      DRE_totbias_var = as.numeric(MCsim$DRE$totbias_var)[1], # Variance of total bias corresponding to the DRE estimate
      DRE_meanbias_mean = as.numeric(MCsim$DRE$meanbias_mean)[1], # Mean of mean bias corresponding to the DRE estimate
      DRE_meanbias_var = as.numeric(MCsim$DRE$meanbias_var)[1], # Variance of mean bias corresponding to the DRE estimate
      DRE_rmse_mean = as.numeric(MCsim$DRE$rmse_mean)[1], # Mean of RMSE corresponding to the DRE estimate
      DRE_rmse_var = as.numeric(MCsim$DRE$rmse_var)[1], # Variance of RMSE corresponding to the DRE estimate
      DRE_bias2_mean = MCsim$DRE$bias2_mean[1], # Mean of bias squared corresponding to the DRE estimate
      DRE_var_mean = MCsim$DRE$var_mean[1], # Mean of variance corresponding to the DRE estimate
      
      # Iterative proportional fitting results
      IPF_totbias_mean = as.numeric(MCsim$IPF$totbias_mean)[1], # Mean of total bias corresponding to the IPF estimate
      IPF_totbias_var = as.numeric(MCsim$IPF$totbias_var)[1], # Variance of total bias corresponding to the IPF estimate
      IPF_meanbias_mean = as.numeric(MCsim$IPF$meanbias_mean)[1], # Mean of mean bias corresponding to the IPF estimate
      IPF_meanbias_var = as.numeric(MCsim$IPF$meanbias_var)[1], # Variance of mean bias corresponding to the IPF estimate
      IPF_rmse_mean = as.numeric(MCsim$IPF$rmse_mean)[1], # Mean of RMSE corresponding to the IPF estimate
      IPF_rmse_var= as.numeric(MCsim$IPF$rmse_var)[1], # Variance of RMSE corresponding to the IPF estimate
      IPF_bias2_mean = MCsim$IPF$bias2_mean[1], # Mean of bias squared corresponding to the IPF estimate
      IPF_var_mean = MCsim$IPF$var_mean[1], # Mean of variance corresponding to the IPF estimate
      
      # External estimator results
      EXT_totbias_mean = as.numeric(MCsim$EXT$totbias_mean)[1], # Mean of total bias corresponding to the EXT estimate
      EXT_totbias_var = as.numeric(MCsim$EXT$totbias_var)[1], # Variance of total bias corresponding to the EXT estimate
      EXT_meanbias_mean = as.numeric(MCsim$EXT$meanbias_mean)[1], # Mean of mean bias corresponding to the EXT estimate
      EXT_meanbias_var = as.numeric(MCsim$EXT$meanbias_var)[1], # Variance of mean bias corresponding to the EXT estimate
      EXT_rmse_mean = as.numeric(MCsim$EXT$rmse_mean)[1], # Mean of RMSE corresponding to the EXT estimate
      EXT_rmse_var = as.numeric(MCsim$EXT$rmse_var)[1], # Variance of RMSE corresponding to the EXT estimate
      EXT_bias2_mean = MCsim$EXT$bias2_mean[1], # Mean of bias squared corresponding to the EXT estimate
      EXT_var_mean = MCsim$EXT$var_mean[1], # Mean of variance corresponding to the EXT estimate
      
      # Relative improvement results
      REL_mean_IPFminDRE = as.numeric(MCsim$REL_IPF$mean_IPFminusDRE)[1], # Mean relative improvement of IPF over DRE
      REL_var_IPFminDRE  = as.numeric(MCsim$REL_IPF$var_IPFminusDRE)[1], # Variance of relative improvement of IPF over DRE
      
      REL_mean_EXTminDRE = as.numeric(MCsim$REL_EXT$mean_EXTminusDRE)[1], # Mean relative improvement of EXT over DRE
      REL_var_EXTminDRE  = as.numeric(MCsim$REL_EXT$var_EXTminusDRE)[1] # Variance of relative improvement of EXT over DRE
    )
    
    setTxtProgressBar(pb, i)  # Update progress bar
    
  }
  
  close(pb)  # Close progress bar
  do.call(rbind, outlist)  # Combine all list elements
  
}