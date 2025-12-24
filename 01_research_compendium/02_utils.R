# This script contains all helper functions required to perform the analysis of this research

# Helper function to ensure the correct factor level structure
# Variable x is ensured to have the factor levels, specified in lvls

fixLvl <- function(x, lvls) { # Ensure correct factor structure 
  
  factor(x, levels = lvls)
  
}

# Helper function to draw latent variable Z given Y, P, X 
# Integrating the conditional independence assumption dependency parameter cia

ciadag <- function(Ypop, Ppop, Xpop, Ylvl, Zlvl, Plvl, 
                   zpx = 0.4, # Probability Z depends on P and X
                   cia = 0) { # If zero, the conditional independence assumption holds
  
  nZ <- length(Zlvl) # Number of categories of latent variable Z
  nP <- length(Plvl) # Number of categories of proxy variable P
  
  tranmatZP <- matrix((1 - zpx) / (nP - 1), # Transition matrix Z given P
                      nrow = nZ, ncol = nP, 
                      dimnames = list(Zlvl, Plvl)) 
  
  # In case nZ does not equal nP
  # zpx probability is inserted on the diagonal
  
  for (i in seq_len(min(nZ, nP))) { 
    
    tranmatZP[i, i] <- zpx # Insert dependency Z on P and X
    
  }
  
  tranmatZP <- sweep(tranmatZP, 2, colSums(tranmatZP), "/") # Normalize the transition matrix
  
  nY <- length(Ylvl) # Number of categories of variable Y
  
  tranmatZY <- matrix((1 - zpx) / (nY - 1), nrow = nZ, ncol = nY, # Transition matrix Z given Y
                      dimnames = list(Zlvl, Ylvl))
  
  
  # In case nZ does not equal nY
  # zpx probability is inserted on the diagonal
  
  for (i in seq_len(min(nZ, nY))) {
    
    tranmatZY[i, i] <- zpx # Insert dependency Z on P and X
    
  }
  
  tranmatZY <- sweep(tranmatZY, 2, colSums(tranmatZY), "/") # Normalize the transition matrix
  
  idxP <- as.integer(factor(Ppop, levels = Plvl)) # column index for P
  idxY <- as.integer(factor(Ypop, levels = Ylvl)) # column index for Y
  
  ZPXmat <- tranmatZP[, idxP, drop = FALSE] # Transition matrix Z given P and X
  ZYmat  <- tranmatZY[, idxY, drop = FALSE] # Transition matrix Z given Y
  
  Zmat <- (1 - cia) * ZPXmat + cia * ZYmat # Integrate cia parameter 
  Zmat <- sweep(Zmat, 2, colSums(Zmat), "/") # Normalize the resulting matrix
  
  Zdraw <- apply(Zmat, 2, function(p) which(rmultinom(1, 1, p) == 1)) # Draw Z given Y, P, X
  Zpop <- factor(Zlvl[Zdraw], levels = Zlvl) # Ensure correct factor structure
  
  factor(Zpop, levels = Zlvl) # Return variable Z
  
}

# Helper function to create random off-diagonal elements for the transition matrix

randoff <- function(r) { # Function to create random off-diagonal elements
  
  w <- rexp(r) # Draw randomly values from a gamma(1,1) distribution
  w / sum(w) # Normalize such that row sums equal 1
  
  # In this case, a gamma(1,1) distribution is used to ensure positive values
  # a gamma(1,1) distribution is shape wise comparable to an uniform distribution
  # A gamma distribution is used instead of an uniform distribution to avoid zero values and is more flexible
  
}

# Helper function to compute the marginal probability of variable VAR 
# First the variable is ensured to be a factor variable 
# If factor levels are specified, these are imposed as well
# Then, the marginal probability is computed using the table()-function
# A smoothing factor is integrated to avoid cells containing zero counts

marSM <- function(VAR, levels = NULL, a = 0.5) {  
  
  if (is.null(levels)) { # If there are no factor levels specified
    
    VAR <- factor(VAR) # VAR in ensured to be a factor variable          
    
  } else { # If there are factor levels specified
    
    VAR <- fixLvl(VAR, levels) # VAR is ensured to be a factor variable 
    
  }
  
  tab <- table(VAR) # Compute the marginal probabilities
  (tab + a) / sum(tab + a) # Integrate the smoothing factor and normalize
  
}

# Helper function to compute the conditional probability of a variable given a second variable
# The given argument is used to condition either on the column or row variable
# If factor levels are specified, they can be imposed on the column and row variables
# A smoothing factor is integrated to avoid cells containing zero counts

conSM <- function(RVAR, CVAR, given = c("col","row"),
                  rlv = NULL, 
                  clv = NULL, 
                  a = 0.5) { 
  
  given <- match.arg(given) # Condition on either the column or row variable
  
  if (is.null(rlv)) { # If there are no factor levels specified for the row variable
    
    RVAR <- factor(RVAR) # RVAR is ensured to be a factor variable          
    
  } else { # If there are factor levels specified for the row variable
    
    RVAR <- fixLvl(RVAR, rlv) # Enforce factor levels to the row variable if present        
    
  }
  
  if (is.null(clv)) { # If there are no factor levels specified for the column variable
    
    CVAR <- factor(CVAR) # CVAR is ensured to be a factor variable
    
  } else { # If there are factor levels specified for the column variable
    
    CVAR <- fixLvl(CVAR, clv) # Enforce factor levels to the column variable if present
    
  }
  
  # Construct memory matrix with correct number of dimensions and factor levels if specified
  
  J <- table(RVAR, CVAR)        
  J <- (J + a) / sum(J + a) # Integrate the smoothing factor and normalize
  
  if (given == "col") { # In case the distribution is conditioned on the column variable
    
    conprob <- sweep(J, 2, colSums(J), "/") # Condition on the column variable
    
  } else { # In case the distribution is condition on the row variable
    
    conprob <- sweep(J, 1, rowSums(J), "/") # Condition on the row variable
  }
  
  return(conprob)  # Return the conditional probability
}

# Helper function to order the YZ matrix correctly

orderYZ <- function(M, Ylvl = NULL, Zlvl = NULL) { 
  
  if (is.null(Ylvl)) { # If no Y levels are specified
    
    # Extract Y levels from row names of M or from the table structure of M
    Ylvl <- if (!is.null(rownames(M))) rownames(M) else sort(unique(rownames(as.table(M)))) 
    
  }
  
  if (is.null(Zlvl)) { # If no Z levels are specified
    
    # Extract Z levels from column names of M or from the table structure of M
    Zlvl <- if (!is.null(colnames(M))) colnames(M) else sort(unique(colnames(as.table(M))))
    
  }
  
  M[Ylvl, Zlvl, drop = FALSE]# Return the ordered matrix
  
}

# Helper function to align all proxy variable levels

alignP <- function(YconP, PconX, ZconPXlist, Plvl) {  
  
  # Construct a list containing all proxy variable level vectors
  
  Psets <- c(list(colnames(YconP), rownames(PconX)), lapply(ZconPXlist, colnames)) 
  P <- Reduce(intersect, Psets) # Extract the intersection of the proxy variable level vectors
  
  if (length(P) == 0L) P <- Plvl # If the intersection is empty, the default proxy factor levels are inserted
  
  list( # Enforce the intersected proxy variable levels
    YconP  = YconP[, P, drop = FALSE],
    PconX  = PconX[P, , drop = FALSE],
    ZconPX = lapply(ZconPXlist, function(M) M[, P, drop = FALSE]),
    Pused  = P
  )
}

# Helper function to perform Iterative Proportional Fitting (IPF)

ipfSafe <- function(seed, # Seed matrix, the unit overlap of samples A and B in principle
                    rowtar, coltar, # Target marginal distributions of the row and column variables
                    maxit = 5000, # Maximal number of iterations
                    tol = 1e-10, # Tolerance factor to specify the break condition
                    eps = 1e-10) { # Smoothing factor to avoid cells containing zero counts
  
  M <- seed # Seed matrix
  M <- M + eps # Integrate additive smoothing factor
  M <- M / sum(M) # Normalize the resulting seed matrix
  
  r <- as.numeric(rowtar) # Target marginal distribution of the row variable 
  c <- as.numeric(coltar) # Target marginal distribution of the column variable 
  
  r[is.na(r)] <- 0 # Replace NAs in specified target row marginal by 0
  c[is.na(c)] <- 0 # Replace NAs in specified target column marginal by 0
  
  if (sum(r) == 0) { # If all elements in the target row marginal are zero
    
    r[] <- 1 / length(r) # The target row marginal is uniformly distributed
    
  }
  
  if (sum(c) == 0) { # If all elements in the target column marginal are zero
    
    c[] <- 1 / length(c) # The target column marginal is uniformly distributed
    
  }
  
  i <- 0 # Initialize iteration counter 
  
  for (it in 1:maxit) {
    
    i <- it # Iteration count
    
    rs <- rowSums(M) # Compute the (old) row totals
    rs[rs == 0] <- NA_real_ # Replace zeroes avoiding zero division
    M <- sweep(M, 1, r / rs, "*") # Rescale rows using the row target
    
    cs <- colSums(M) # Compute the (old) column totals
    cs[cs == 0] <- NA_real_  # Replace zeroes avoiding zero division
    M <- sweep(M, 2, c / cs, "*") # Rescale columns using the column target
    
    row_err <- abs(rowSums(M) - r) # Compute the difference between the row totals and the target
    col_err <- abs(colSums(M) - c) # Compute the difference between the column totals and the target
    
    if (max(row_err, na.rm = TRUE) < tol && max(col_err, na.rm = TRUE) < tol) { 
      break # Break condition
    }
  }
  
  M <- M / sum(M) # Normalize resulting joint distribution
  
  out <- list(
    matrix = M,                
    iterations = i              
  )
  
  return(out)
}


# Helper function to pull the first element of a vector or NA if the vector is empty

pull <- function(x) {
  
  if (length(x) == 0L) return(NA)
  x[1]
  
}

# Helper function to retrieve the marginal of X

getXmar <- function(source = c("B", "pop"), B, pop, Xlvl, a = 0.5) {
  
  source <- match.arg(source) # Match the source argument
  
  if (source == "B") { # If the marginal of X is retrieved from sample B
    
    marX <- marSM(B$X, levels = Xlvl, a = a) # Compute the marginal of X in sample B
    
  } else { # If the marginal of X is retrieved from the population
    
    marX <- marSM(pop$X, levels = Xlvl, a = a) # Compute the marginal of X in the population
    
  }
  
  return(marX) # Return the marginal of X
  
}
