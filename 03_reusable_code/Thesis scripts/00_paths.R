# This code defines helper functions to create file paths for different project directories

# In case the "here" package is not installed yet
# install.packages("here") 

library(here) # Load the here package

# Define helper functions for different project directories

path_data <- function(...) here::here("data", ...)
path_outputs <- function(...) here::here("outputs", ...)
path_scripts <- function(...) here::here("scripts", ...)
path_R <- function(...) here::here("R", ...)
path_figs <- function(...) here::here("outputs", "figures", ...)
path_tabs <- function(...) here::here("outputs", "tables", ...)