# The SIM function generates a synthetic population 
# Including samples A, B, the overlap and external samples

SIM <- function(
    nY = 3, # Number of Y categories
    nZ = 3, # Number of Z categories
    nP = 3, # Number of P categories
    nX = 6, # Number of X categories
    
    Ypopmar = rep(1/3,3), # Vector of population margin Y
    Zpopmar = rep(1/3,3), # Vector of population margin Z
    Xpopmar = rep(1/6,6), # Vector of population margin Z
    
    N = 1000000, # Population size
    
    tranmat_diag, # Strength relation Y and P (vector)
    tranmat_sym = TRUE, # Symmetrical property of transition matrix
    
    tol = 1e-6, # Tolerance value
    seed = 0,
    
    nA, # Sample A size
    nB, # Sample B size
    nE, # External sample size
    PercOverlap, # Sample size overlap percentage
    
    zpx = 0.4, # Dependency Z on P and X
    cia = 0, # If 0, CIA is valid
    
    Xsource = c("B", "pop"), # Source of the marginal distribution of X
    
    mechanism = c("MAR", "MNAR"), # Selectivity mechanism
    pattern, # Selectivity pattern when mechanism is MAR
    degree # Selectivity degree when mechanism is MAR
    
) {
  
  set.seed(seed) # Ensure reproducibility
  
  mechanism <- match.arg(as.character(mechanism), c("MAR", "MNAR")) # Match selectivity mechanism
  
  if (any(Ypopmar < 0) || any(Zpopmar < 0) || any(Xpopmar < 0))
    stop("Population margins contain negative values.") # Check whether population margins are positive
  
  if (abs(sum(Ypopmar) - 1) > tol) # Check whether population margins add up to 1
    stop("Sum of Ypopmar does not add up to 1.")
  if (abs(sum(Zpopmar) - 1) > tol)
    stop("Sum of Zpopmar does not add up to 1.")
  if (abs(sum(Xpopmar) - 1) > tol)
    stop("Sum of Xpopmar does not add up to 1.")
  
  Ylvl <- paste0("Y", seq_len(nY)) # Labels of Y levels
  Zlvl <- paste0("Z", seq_len(nZ)) # Labels of Z levels
  Xlvl <- paste0("X", seq_len(nX)) # Labels of X levels
  Plvl <- paste0("P", seq_len(nP)) # Labels of P levels
  
  if (length(Ypopmar) != nY) # Check whether population margins align with the number of categories
    stop("Length of Ypopmar does not match nY.")
  if (length(Zpopmar) != nZ) 
    stop("Length of Zpopmar does not match nZ.")
  if (length(Xpopmar) != nX) 
    stop("Length of Xpopmar does not match nX.")
  
  Ypop <- factor( # Generate population Y variable
    sample(Ylvl, N, replace = TRUE, prob = Ypopmar),
    levels = Ylvl
  )
  
  Xpop <- factor( # Generate population X variable
    sample(Xlvl, N, replace = TRUE, prob = Xpopmar),
    levels = Xlvl
  )
  
  if (length(tranmat_diag) == 1L) tranmat_diag <- rep(tranmat_diag, nY) 
  # Ensure that the diagonal of transition matrix same length as nY
  
  # Create memory for transition matrix
  tranmat <- matrix(0, nrow = nY, ncol = nP, dimnames = list(Y = Ylvl, P = Plvl))
  
  if (tranmat_sym) { # If transition matrix is symmetric
    
    if (nY != nP) # Stopping condition in case of differing number of categories
      stop("Number of Categories across Y and P differ")
    
    # Construct empty transition matrix
    tranmat <- matrix(0, nrow = nY, ncol = nP,
                      dimnames = list(Y = paste0("Y", seq_len(nY)),
                                      P = paste0("P", seq_len(nP))))
    
    off_diag <- (1 - tranmat_diag) / (nY - 1) # Off-diagonal elements
    
    for (i in seq_len(nY)) { # Define the transition matrix
      
      tranmat[i, ] <- off_diag[i]   # Insert off-diagonal elements
      tranmat[i, i] <- tranmat_diag[i] # Insert diagonal elements
      
    }
    
  } else { # If transition matrix is not symmetric
    
    for (i in seq_len(nY)) {
      
      d <- tranmat_diag[i] # Define diagonal element
      massoff <- 1 - d # Off-diagonal probability mass                
      
      j <- i # Connect Y and P levels                    
      offcol <- setdiff(seq_len(nP), j) # All columns except the diagonal one
      
      if (length(offcol) == 0L) { # In case nP == 1
        
        tranmat[i, j] <- d # Only diagonal element
        
      } else { # In case nP > 1
        
        # Random off-diagonal elements scaled to off-diagonal mass
        
        off <- randoff(length(offcol)) * massoff
        
        tranmat[i, offcol] <- off   # Insert off-diagonal elements
        tranmat[i, j]      <- d     # Insert diagonal element
        
        # In this case, a gamma(1,1) distribution is used to ensure positive values
        # a gamma(1,1) distribution is shape wise comparable to an uniform distribution
        # A gamma distribution is used instead of an uniform distribution to avoid zero values and is more flexible
        
      }
    }
  }  
  
  Ppop <- character(length(Ypop)) # Memory vector
  
  for (y in Ylvl) { # Loop across the Y categories
    
    idx <- which(Ypop == y) # Select indices for Y level
    
    if (length(idx)) { # Use transition matrix to generate Ppop
      
      Ppop[idx] <- sample(Plvl, size = length(idx), replace = TRUE, prob = tranmat[y, ])
      
    }
  }
  
  Ppop <- factor(Ppop, levels = Plvl) # Impose factorial structure
  
  Zpop <- ciadag(Ypop, Ppop, Xpop, Ylvl, Zlvl, Plvl, zpx = zpx, cia = cia) # Generate population Z variable
  
  popdat <- data.frame( # Create population data frame
    ID = seq_len(N), # Identification variable
    Y = Ypop, 
    Z = Zpop, 
    X = Xpop,
    P = Ppop)  
  
  # Sample A
  
  A <- popdat |> # Draw sample A
    sample_n(nA, replace = FALSE) |> # Sample without replacement
    select(ID, Y, X) # Select relevant variables
  
  # Sample B
  
  Poverlap <- PercOverlap # Percentage of unit overlap
  Noverlap <- round(Poverlap * nB) # Number overlapping units
  IDoverlap <- sample(A$ID, Noverlap, replace = FALSE) # Sample overlap from sample A
  
  PopNonA <- setdiff(popdat$ID, A$ID) # Create subset of population data excluding units in sample A
  restB <- nB - Noverlap # Compute the number of units that needs to be sampled from the upper subset
  
  # Stopping condition in case not enough units are available for sampling
  
  if (restB > length(PopNonA)) stop("Not enough population units outside A for requested overlap.") 
  
  subB <- sample(PopNonA, restB, replace = FALSE) # Sample additional units for sample B
  B_ID <- c(IDoverlap, subB) # The ID's of the units in sample B
  
  B <- popdat |> # Draw sample B
    filter(ID %in% B_ID) |> # Filter population data for sample B ID's
    select(ID, Z, X, P) # Select relevant variables
  
  # Overlap sample
  
  overlap <- inner_join(A, B, by = c("ID", "X")) |> # Create overlap sample
    select(ID, Y, X, Z, P) # Select relevant variables
  
  # Sample E
  
  if (mechanism  == "MAR") {
    
    # Retrieve pattern and degree selectivity weights for X
    
    wX <- mar_patdeg[[pattern]][[degree]] 
    
    # Fallback to original if there is no degree specified
    
    if (is.null(wX)) wX <- mar_patdeg[[pattern]][["original"]] 
    
    # Convey the weights to the population data
    
    wE <- wX[as.character(popdat$X)] 
    
  } else {
    
    # Retrieve the YZ matrix for the specified pattern
    
    matYZ <- mnar_YZ[[pattern]] 
    
    # Convert factor variables into characters for safety reasons
    
    Yc <- as.character(popdat$Y)
    Zc <- as.character(popdat$Z)
    
    # Convey the weights to the population data
    
    wE <- mapply(function(y, z) matYZ[y, z], Yc, Zc)
    
  }
  
  # Ensure weights are all positive
  
  wE <- pmax(wE, 0)
  
  # SRS as fallback in case all weights are zero
  
  if (all(wE == 0)) wE[] <- 1
  
  # Retrieve normalized probabilities
  
  probE <- wE / sum(wE)
  
  # Sample E based on indices drawn according to the weights
  
  idxE <- sample( 
    seq_len(nrow(popdat)),
    size = nE,
    replace = FALSE,
    prob = probE
  )
  
  E <- popdat |> # Draw sample E
    slice(idxE) |> # Select sampled indices
    select(ID, Y, P) # Select relevant variables
  
  
  # Get the marginal distribution of X
  
  Xsource <- match.arg(as.character(Xsource), c("B", "pop")) # Match argument for the source of the marginal distribution of X
  marX <- getXmar( # Get marginal distribution of X
    source = Xsource,
    B = B,
    pop = popdat,
    Xlvl = Xlvl,
    a = 0.5
  )
  
  # True joint distribution of Y and Z in population
  
  YZpop <- prop.table(table(popdat$Y, popdat$Z)) # True joint distribution 
  YZpop <- YZpop[Ylvl, Zlvl, drop = FALSE] # Order YZ levels
  
  SIMout <- list( # Output list
    meta = list(
      Xsource = Xsource,
      mechanism = mechanism,
      pattern = pattern,
      degree = degree), # Source of the marginal distribution of X
    popdat = popdat, # Population data
    poplvl = list( # List variable levels
      Ylvl = Ylvl, 
      Zlvl = Zlvl,
      Xlvl = Xlvl,
      Plvl = Plvl
    ),
    popdat_summary = summary(popdat), # Summary statistics population data
    tranmat = tranmat, # Transition Matrix
    popjoint = list(
      YZpop = YZpop, # True joint Y and Z in population
      YPpop = prop.table(table(popdat$Y, popdat$P)), # True joint Y and P in population
      PconY = prop.table(table(popdat$Y, popdat$P), margin = 1)), # Conditional of P given Y
    samples = list( # Samples A, B, E and the overlap
      A = A,
      B = B,
      E = E, 
      overlap = overlap)
    
  )
  
  return(SIMout) 
  
}
