# This script can be used to load and prepare the NHANES data for the simulation study

# Get the data

d1 <- read_xpt(here("Data","DEMO_I.xpt"))
d2 <- read_xpt(here("Data","BPX_I.xpt"))
d3 <- read_xpt(here("Data","BMX_I.xpt"))
d4 <- read_xpt(here("Data","GHB_I.xpt"))
d5 <- read_xpt(here("Data","TCHOL_I.xpt"))

# Convert all variable names to lowercase 

d1 <- rename_with(d1, tolower)
d2 <- rename_with(d2, tolower)
d3 <- rename_with(d3, tolower)
d4 <- rename_with(d4, tolower)
d5 <- rename_with(d5, tolower)


# Subset the data

d1.t <- dplyr::select(d1, seqn, riagendr, ridageyr)
d2.t <- dplyr::select(d2, seqn, bpxsy1)
d3.t <- dplyr::select(d3, seqn, bmxbmi)
d4.t <- dplyr::select(d4, seqn, lbxgh)
d5.t <- dplyr::select(d5, seqn, lbdtcsi)

# Merge the data

d <- list(d1.t, d2.t, d3.t, d4.t, d5.t) |>
  purrr::reduce(dplyr::left_join, by = "seqn")