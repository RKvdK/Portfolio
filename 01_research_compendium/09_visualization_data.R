# This script visualizes the simulation results stored in 'mlarpis_dat.rds'
# It is important to note that this figure is purely for illustrative purposes given this is a subset
# In my thesis project, I intend to create more comprehensive visualizations based on the full dataset

rmse_summary <- mlarpis_dat %>%
  select(ends_with("_rmse_mean")) %>% # Select the variables with RMSE means
  summarise(across(everything(), mean)) %>% # Average across all subsetted parameter combinations
  pivot_longer(cols = everything(), # Convert to long format
               names_to = "estimator",
               values_to = "avg_rmse") %>%
  mutate(estimator = sub("_rmse_mean$", "", estimator))

barplot <- ggplot(rmse_summary, aes(x = estimator, y = avg_rmse)) +
  geom_col(fill = "pink2") +
  labs(
    title = "Average RMSE value across all subsetted parameter conbinations",
    x = NULL,
    y = "Average RMSE"
  ) +
  theme_bw()

print(barplot)

ggsave("rmse_barplot.png", plot = barplot)