# This script contains all selectivity patterns and degrees assuming Missing at Random
# Also, this script contains all main and interaction scenarios assuming Missing Not at Random
# The derived patterns and scenarios are retrieved from Sojka (2025)

# Missing at Random
# Selectivity patterns and degrees

mar_patdeg <- list(
  MAR_LinDec = list( # Linear decrease pattern
    reduced = c(X1 = 3, X2 = 2.6, X3 = 2.2, X4 = 1.8, X5 = 1.4, X6 =1), # Reduced selectivity degree
    original = c(X1 = 6, X2 = 5, X3 = 4, X4 = 3, X5= 2, X6 = 1), # Original selectivity degree
    increased = c(X1 = 12, X2 = 9.8, X3 = 7.6, X4 = 5.4, X5 = 3.2, X6 = 1) # Increased selectivity degree
  ),
  MAR_Ushape = list( # U-shaped pattern
    original = c(X1 = 8, X2 = 2, X3 = 1, X4 = 1, X5 = 2, X6 = 8) # Original selectivity degree
  ),
  MAR_Step = list( # Step function pattern
    original = c(X1 = 4, X2 = 4, X3 = 4, X4 = 1, X5 = 1, X6 = 1) # Original selectivity degree
  ),
  MAR_ExtInc = list( # Extreme increase pattern
    original = c(X1 = 1, X2 = 2, X3 = 3, X4 = 5, X5 = 10, X6 = 40) # Original selectivity degree
  )
)

# Missing Not at Random

mnar_YZ <- list(
  
  # Main scenarios
  
  MNAR_ClasInc = matrix( # Classic increase main scenario
    c(
      1, 2, 4,
      2, 4, 8,
      3, 6, 12
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_NonMono = matrix( # Non-monotonic main scenario
    c(
      6, 2, 12,
      15, 5, 30,
      3, 1, 6
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_Yonly = matrix( # Y-only main scenario
    c(
      4, 4, 4,
      1, 1, 1,
      7, 7, 7
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  # Interaction scenarios
  
  MNAR_ClasInc_WeakInt = matrix( # Classic increase with weak interaction
    c(
      1, 1.7, 3.3,  
      1.7, 4.3, 6.7,  
      2.5, 5, 11  
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_ClasInc_ModInt = matrix( # Classic increase with moderate interaction
    c(
      1, 1, 2,
      1, 5, 4,
      1.5, 3, 10.8
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_ClasInc_StrongInt = matrix( # Classic increase with strong interaction
    c(
      1.9, 1, 1.7,
      1, 8.8, 4,
      1.3, 3, 18.8
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_ClasInc_ExtInt = matrix( # Classic increase with extreme interaction
    c(
      10, 1, 1,
      1, 10, 1,
      1, 1, 10
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  )
)
