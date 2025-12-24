# The MC function performs a Monte Carlo simulation to evaluate the performance of the three different estimators 
# It resamples the input data multiple times, computes the estimators for each resample, and calculates bias and root mean square error (RMSE) 
# The SIM function is requisite for this function to work properly

MC <- function(SIMout, # SIM output
               MCnum = 1000, # Number of simulation sets
               seed = 1) {
  
  set.seed(seed) # Ensure reproducibility
  
  A0 <- SIMout$samples$A # Extract sample A from the SIM output
  B0 <- SIMout$samples$B # Extract sample B from the SIM output
  E0 <- SIMout$samples$E # Extract sample E from the SIM output
  Overlap0 <- SIMout$samples$overlap # Extract the overlap sample from the SIM output
  
  Ylvl <- SIMout$poplvl$Ylvl # Extract population Y levels from the SIM output
  Zlvl <- SIMout$poplvl$Zlvl # Extract population Z levels from the SIM output
  Xlvl <- SIMout$poplvl$Xlvl # Extract population X levels from the SIM output 
  Plvl <- SIMout$poplvl$Plvl # Extract population P levels from the SIM output
  
  YZpop <- SIMout$popjoint$YZpop[Ylvl, Zlvl, drop = FALSE] # Extract population joint distribution of Y and Z from the SIM output
  
  Xsource <- SIMout$meta$Xsource # Extract source of X from the SIM output
  
  if (is.null(Xsource)) Xsource <- "B" # Default source of X is sample B
  Xsource <- match.arg(as.character(Xsource), c("B", "pop")) # Match the source for the marginal distribution of X
  
  mean_theta_DRE <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory for mean bias per cell for the DRE estimator
  M2_theta_DRE <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory for sum of squares per cell for the DRE estimator
  tb_DRE  <- numeric(MCnum) # Memory for total bias for the doubly robust estimator
  mb_DRE  <- numeric(MCnum) # Memory for mean bias for the doubly robust estimator
  rmse_DRE  <- numeric(MCnum) # Memory for RMSE for the doubly robust estimator
  
  mean_theta_IPF <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory for mean bias per cell for the IPF estimator
  M2_theta_IPF <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory for sum of squares per cell for the IPF estimator    
  tb_IPF  <- numeric(MCnum) # Memory for total bias for the IPF estimator
  mb_IPF  <- numeric(MCnum) # Memory for mean bias for the IPF estimator
  rmse_IPF  <- numeric(MCnum) # Memory for RMSE for the IPF estimator
  
  mean_theta_EXT <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory for mean bias per cell for the external estimator
  M2_theta_EXT <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory for sum of squares per cell for the external estimator     
  tb_EXT   <- numeric(MCnum) # Memory for total bias for the external estimator
  mb_EXT   <- numeric(MCnum) # Memory for mean bias for the external estimator
  rmse_EXT <- numeric(MCnum) # Memory for RMSE for the external estimator
  
  # Relative improvements showing how much better DRE is compared to IPF and EXT  
  # Positive values denote improvement of DRE over the other referenced methods
  
  mean_rel_IPF <- numeric(MCnum) # Memory for relative improvement of DRE over IPF  
  mean_rel_EXT <- numeric(MCnum) # Memory for relative improvement of DRE over EXT
  
  for (b in seq_len(MCnum)) { # Monte Carlo loop
    
    Ab <- A0[sample.int(nrow(A0), nrow(A0), replace = TRUE), , drop = FALSE] # Resample A with replacement
    Bb <- B0[sample.int(nrow(B0), nrow(B0), replace = TRUE), , drop = FALSE] # Resample B with replacement
    Eb <- E0[sample.int(nrow(E0), nrow(E0), replace = TRUE), , drop = FALSE] # Resample E with replacement
    Ob <- Overlap0[sample.int(nrow(Overlap0), nrow(Overlap0), replace = TRUE), , drop = FALSE] # Resample overlap with replacement
    
    # Conditional probability of Y given P from sample E
    YconP_E <- conSM(Eb$Y, Eb$P, given = "col",
                     rlv = Ylvl, clv = Plvl, a = 0.5)
    
    # Conditional probability of P given X from sample B
    PconX_B <- conSM(Bb$P, Bb$X, given = "col",
                     rlv = Plvl, clv = Xlvl, a = 0.5)
    
    marXb <- getXmar(
      source = Xsource,
      B = Bb,
      pop = SIMout$popdat,
      Xlvl = Xlvl,
      a = 0.5
    )
    
    marXb <- marXb / sum(marXb) # Normalize the marginal distribution of X
    
    ZconPX_B <- setNames(vector("list", length(Xlvl)), Xlvl) # Create memory
    
    for (x in Xlvl) { # Loop per level of X
      
      idx <- Bb$X == x # Filter per level of X
      ZconPX_B[[x]] <- conSM(Bb$Z[idx], Bb$P[idx], given = "col", # Conditional probability of Z given P and X per level of X
                             rlv = Zlvl, clv = Plvl, a = 0.5)
    }
    
    aligned  <- alignP(YconP_E, PconX_B, ZconPX_B, Plvl) # Apply alignment function
    
    YconP_EB <- aligned$YconP # Store aligned objects
    PconX_B  <- aligned$PconX 
    ZconPX_B <- aligned$ZconPX
    
    # Compute double robust estimator
    # This part assumes that the external sample reflects the population structure correctly
    # Compute YZ_EB
    
    YZ_EB <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory matrix
    
    for (x in Xlvl) { # Loop per level of X
      
      ZconPx <- ZconPX_B[[x]] # Extract per level of X
      wPx <- as.numeric(PconX_B[, x]) * as.numeric(marXb[x]) # Weights
      YZ_EB <- YZ_EB + tcrossprod(YconP_EB, sweep(ZconPx, 2, wPx, "*")) # Update YZ_EB
      
    }
    
    YZ_EB <- YZ_EB / sum(YZ_EB) # Normalize the resulting matrix
    jointEXT <- orderYZ(YZ_EB) # Store the external estimator joint distribution of Y and Z
    
    # Compute YZ_EO 
    # This part assumes that the overlap sample reflects the population structure correctly
    
    if (nrow(Ob) == 0) { # If there is no unit overlap
      
      YZ_EO <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), 
                      dimnames = list(Ylvl, Zlvl)) # Memory matrix             
      
      YZoverlap <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl),
                          dimnames = list(Ylvl, Zlvl)) # Memory matrix           
      
    } else { # If there is unit overlap
      
      # Conditional probability of P given X (overlap) with smoothing
      PconX_O <- conSM(Ob$P, Ob$X, given = "col",
                       rlv = Plvl, clv = Xlvl, a = 0.5)
      
      # Marginal of X in overlap with smoothing
      X_O <- marSM(Ob$X, levels = Xlvl, a = 0.5)
      
      # Conditional probability of Z given P and X (per X-level) with smoothing
      ZconPX_O <- setNames(vector("list", length(Xlvl)), Xlvl)
      for (x in Xlvl) {
        idx <- Ob$X == x
        ZconPX_O[[x]] <- conSM(Ob$Z[idx], Ob$P[idx], given = "col",
                               rlv = Zlvl, clv = Plvl, a = 0.5)
      }
      
      # Align P levels across Y|P (from E), P|X (overlap) and Z|P,X (overlap)
      aligned_O <- alignP(YconP_E, PconX_O, ZconPX_O, Plvl)
      YconP_O   <- aligned_O$YconP
      PconX_O   <- aligned_O$PconX
      ZconPX_O  <- aligned_O$ZconPX
      
      YZ_EO <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl),
                      dimnames = list(Ylvl, Zlvl)) # Memory matrix
      
      for (x in Xlvl) { # Loop per level of X
        
        ZconPxO <- ZconPX_O[[x]]  # Extract per level of X
        wPxO    <- as.numeric(PconX_O[, x]) * as.numeric(X_O[x]) # Weights
        YZ_EO   <- YZ_EO + tcrossprod(YconP_O, sweep(ZconPxO, 2, wPxO, "*")) # Update YZ_EO
      }
      
      YZ_EO <- YZ_EO / sum(YZ_EO) # Normalize the resulting matrix
      YZoverlap <- prop.table(table(Ob$Y, Ob$Z)) # Joint distribution of Y and Z in overlap sample
      
    }
    
    jointDRE <- (orderYZ(YZoverlap) - orderYZ(YZ_EO)) + orderYZ(YZ_EB) # Doubly robust estimator joint distribution of Y and Z
    jointDRE <- jointDRE / sum(jointDRE) # Normalize the resulting matrix
    
    # Iterative Proportional Fitting
    
    jointIPFx <- array(0, dim = c(length(Ylvl), length(Zlvl), length(Xlvl)),
                       dimnames = list(Ylvl, Zlvl, Xlvl)) # Memory array
    
    # Integrate smoothing into the conditional probabilities
    YconX_A  <- conSM(Ab$Y, Ab$X, given = "col", rlv = Ylvl, clv = Xlvl, a = 0.5)
    ZconX_Bb <- conSM(Bb$Z, Bb$X, given = "col", rlv = Zlvl, clv = Xlvl, a = 0.5)
    
    for (x in Xlvl) { # Loop per level of X
      
      rowtar <- as.numeric(YconX_A[Ylvl, x, drop = FALSE]) # Row targets
      coltar <- as.numeric(ZconX_Bb[Zlvl, x, drop = FALSE]) # Column targets
      
      rowtar[is.na(rowtar)] <- 0 # Insert zeroes in case of missings in the row target
      coltar[is.na(coltar)] <- 0 # Insert zeroes in case of missings in the column target
      
      if (sum(rowtar) == 0) rowtar[] <- 1/length(rowtar) # Avoid zeroes in de the row target
      if (sum(coltar) == 0) coltar[] <- 1/length(coltar) # Avoid zeroes in de the column target
      
      idxx <- Ob$X == x # Per level of X
      
      if (any(idxx)) { # Seed if per indicated level of X is available
        
        tabO <- table(Ob$Y[idxx], Ob$Z[idxx]) # Joint distribution in overlap per level of X
        seedIPF <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory matrix
        seedIPF[rownames(tabO), colnames(tabO)] <- as.numeric(tabO) # Insert joint distribution values
        
      } else { # If no overlapping units are available per indicated level of X
        
        seedIPF <- outer(rowtar, coltar) # use independence as seed via the outer product
        dimnames(seedIPF) <- list(Ylvl, Zlvl) # Assign the correct dimension names
        
      }
      
      jointIPFx[, , x] <- ipfSafe(seedIPF, rowtar, coltar)$matrix # Perform IPF and store the resulting matrix
      
    }
    
    PX <- marXb # Integrate smoothing into marginal of X from sample B
    jointIPF <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory matrix
    
    for (x in Xlvl) { # Loop per level of X
      
      jointIPF <- jointIPF + jointIPFx[, , x] * as.numeric(PX[x]) # Update jointIPF by weighing per level of X
      
    }
    
    jointIPF <- jointIPF / sum(jointIPF) # Normalize the resulting matrix
    
    # Update running means and variances per cell for each estimator using Welford's method
    
    delta_DRE <- jointDRE - mean_theta_DRE # Difference of current estimate to running mean for the DRE estimator
    mean_theta_DRE <- mean_theta_DRE + delta_DRE / b # Update running mean for the DRE estimator
    M2_theta_DRE <- M2_theta_DRE + delta_DRE * (jointDRE - mean_theta_DRE) # Update sum of squares for the DRE estimator
    
    delta_IPF <- jointIPF - mean_theta_IPF # Difference of current estimate to running mean for the IPF estimator
    mean_theta_IPF <- mean_theta_IPF + delta_IPF / b # Update running mean for the IPF estimator
    M2_theta_IPF <- M2_theta_IPF + delta_IPF * (jointIPF - mean_theta_IPF) # Update sum of squares for the IPF estimator
    
    delta_EXT <- jointEXT - mean_theta_EXT # Difference of current estimate to running mean for the external estimator
    mean_theta_EXT <- mean_theta_EXT + delta_EXT / b # Update running mean for the external estimator
    M2_theta_EXT <- M2_theta_EXT + delta_EXT * (jointEXT - mean_theta_EXT) # Update sum of squares for the external estimator   
    
    # Compute bias and rmse
    
    bDRE <- abs(jointDRE - YZpop) # Bias per cell for the doubly robust estimator
    bIPF <- abs(jointIPF - YZpop) # Bias per cell for the IPF estimator
    bEXT <- abs(jointEXT - YZpop) # Bias per cell for the external estimator
    
    # DRE
    tb_DRE[b] <- sum(bDRE) # Total bias for the doubly robust estimator
    mb_DRE[b] <- mean(bDRE) # Mean bias for the doubly robust estimator
    rmse_DRE[b] <- sqrt(mean((jointDRE - YZpop)^2)) # RMSE for the doubly robust estimator
    
    # IPF
    tb_IPF[b] <- sum(bIPF) # Total bias for the IPF estimator
    mb_IPF[b] <- mean(bIPF) # Mean bias for the IPF estimator
    rmse_IPF[b] <- sqrt(mean((jointIPF - YZpop)^2)) # RMSE for the IPF estimator
    
    # EXT
    tb_EXT[b] <- sum(bEXT) # Total bias for the external estimator
    mb_EXT[b] <- mean(bEXT) # Mean bias for the external estimator
    rmse_EXT[b] <- sqrt(mean((jointEXT - YZpop)^2)) # RMSE for the external estimator
    
    # Relative improvements 
    mean_rel_IPF[b] <- mean(bIPF - bDRE) # Relative improvement of DRE over IPF
    mean_rel_EXT[b] <- mean(bEXT - bDRE) # Relative improvement of DRE over EXT
  }
  
  # Bias-variance decomposition per cell and aggregated per estimator
  
  var_theta_DRE <- M2_theta_DRE / (MCnum - 1) # Variance per cell for the DRE estimator
  var_theta_IPF <- M2_theta_IPF / (MCnum - 1) # Variance per cell for the IPF estimator
  var_theta_EXT <- M2_theta_EXT / (MCnum - 1) # Variance per cell for the external estimator
  
  bias_cell_DRE <- mean_theta_DRE - YZpop # Bias per cell for the DRE estimator
  bias_cell_IPF <- mean_theta_IPF - YZpop # Bias per cell for the IPF estimator
  bias_cell_EXT <- mean_theta_EXT - YZpop # Bias per cell for the external estimator
  
  DRE_bias2_mean <- mean(bias_cell_DRE^2) # Mean squared bias per cell for the DRE estimator
  IPF_bias2_mean <- mean(bias_cell_IPF^2) # Mean squared bias per cell for the IPF estimator
  EXT_bias2_mean <- mean(bias_cell_EXT^2) # Mean squared bias per cell for the external estimator
  
  DRE_var_mean <- mean(var_theta_DRE) # Mean variance per cell for the DRE estimator
  IPF_var_mean <- mean(var_theta_IPF) # Mean variance per cell for the IPF estimator
  EXT_var_mean <- mean(var_theta_EXT) # Mean variance per cell for the external estimator
  
  # Output list
  
  MCout <- list( # Output list
    MCnum = MCnum, # Number of Monte Carlo runs
    DRE = list( # Doubly robust estimator results
      totbias = tb_DRE, meanbias = mb_DRE, rmse = rmse_DRE, #
      totbias_mean = mean(tb_DRE), totbias_var = var(tb_DRE),
      meanbias_mean = mean(mb_DRE), meanbias_var = var(mb_DRE),
      rmse_mean = mean(rmse_DRE), rmse_var = var(rmse_DRE),
      bias_cell = bias_cell_DRE,var_cell = var_theta_DRE,
      bias2_mean = DRE_bias2_mean,var_mean = DRE_var_mean
    ),
    IPF = list( # IPF estimator results
      totbias = tb_IPF, meanbias = mb_IPF, rmse = rmse_IPF,
      totbias_mean = mean(tb_IPF), totbias_var = var(tb_IPF),
      meanbias_mean = mean(mb_IPF), meanbias_var = var(mb_IPF),
      rmse_mean = mean(rmse_IPF), rmse_var = var(rmse_IPF),
      bias_cell = bias_cell_IPF, var_cell = var_theta_IPF,
      bias2_mean = IPF_bias2_mean, var_mean = IPF_var_mean
    ),
    EXT = list( # External estimator results
      totbias = tb_EXT, meanbias = mb_EXT, rmse = rmse_EXT,
      totbias_mean = mean(tb_EXT),  totbias_var = var(tb_EXT),
      meanbias_mean = mean(mb_EXT), meanbias_var = var(mb_EXT),
      rmse_mean = mean(rmse_EXT),   rmse_var = var(rmse_EXT),
      bias_cell = bias_cell_EXT, var_cell = var_theta_EXT,
      bias2_mean = EXT_bias2_mean, var_mean = EXT_var_mean
    ),
    REL_IPF = list( # Relative improvement of DRE over IPF
      mean_IPFminusDRE = mean(mean_rel_IPF),
      var_IPFminusDRE  = var(mean_rel_IPF),
      series = mean_rel_IPF
    ),
    REL_EXT = list( # Relative improvement of DRE over EXT
      mean_EXTminusDRE = mean(mean_rel_EXT),
      var_EXTminusDRE  = var(mean_rel_EXT),
      series = mean_rel_EXT
    )
    
  )
  
  return(MCout)
  
}

