# Rename the variables

d$age <- d$ridageyr # Age in years at screening
d$sex <- d$riagendr # Gender
d$bp <- d$bpxsy1 # Systolic blood pressure
d$bmi <- d$bmxbmi # Body Mass Index
d$HbA1C <- d$lbxgh # Glycohemoglobin
d$chol <- d$lbdtcsi # Total Cholesterol

d$age[d$age<18] <- NA # Recode all minors to missing

dc <- cc(subset(d,select=c("age","sex","bmi","HbA1C","bp"))) # Select the complete cases

# Data analysis 

summary(lm(bp ~ HbA1C + age + as.factor(sex), data=dc)) # Summary of linear model without BMI
confint(lm(bp ~ HbA1C + age + as.factor(sex), data=dc)) # Confidence intervals of model without BMI

summary(lm(bp ~ HbA1C + bmi + age + as.factor(sex), data=dc)) # Summary of linear model with BMI
confint(lm(bp ~ HbA1C + bmi + age + as.factor(sex), data=dc)) # Confidence intervals of model with BMI

# Simulation of the measurement error

ref <- lm(bp ~ HbA1C + bmi + age + as.factor(sex), data=dc)$coef[2] # Reference value of the coefficient of the without measurement error

n.sim <- 1e3 # Number of simulations

perc.me.exp <- seq(0,.5,.1) # Percentage of measurement error variance for exposure
perc.me.conf<- seq(0,.5,.1) # Percentage of measurement error variance for confounder
scenarios <- expand.grid(perc.me.exp,perc.me.conf) # Create all scenarios for measurement error via simulation grid

var.exp <- var(dc$HbA1C) # Variance of exposure
var.conf <- var(dc$bmi) # Variance of confounder

n <- dim(dc)[1] # Sample size
beta.hat <- matrix(ncol=dim(scenarios)[1], nrow=n.sim) # Memory matrix to store the simulation results

for (k in 1:n.sim){ # Start of simulation
  
  print(k) # Print the simulation number
  set.seed(k) # Set the seed for reproducibility, varying by simulation number
  
  for (i in 1:dim(scenarios)[1]){ # Loop over all measurement error scenarios
    
    var.me.exp <- var.exp*scenarios[i,1]/(1-scenarios[i,1]) # Calculate the measurement error variance for exposure
    var.me.conf <- var.conf*scenarios[i,2]/(1-scenarios[i,2]) # Calculate the measurement error variance for confounder
    
    dc$HbA1C.me <- dc$HbA1C + rnorm(dim(dc)[1], 0, sqrt(var.me.exp)) # Simulate measurement error for exposure
    dc$bmi.me <- dc$bmi + rnorm(dim(dc)[1], 0, sqrt(var.me.conf)) # Simulate measurement error for confounder
    
    beta.hat[k,i] <- lm(bp ~ HbA1C.me + age + bmi.me + as.factor(sex), data=dc)$coef[2] # Fit the model with measurement error and store the coefficient estimate
    
    } # End of loop over all measurement error scenarios
  } # End of simulation

# Data visualization

tot.mat <- cbind(100*scenarios,apply(beta.hat,2,mean)) # Create a matrix with measurement error percentages and mean estimates
colnames(tot.mat) <- c("me.exp","me.conf","estimate") # Rename the columns

FIGURE <- ggplot(tot.mat, aes(me.exp, me.conf)) + #
  geom_tile(color="white",aes(fill = estimate)) +
  geom_text(aes(label = round(estimate, 2))) +
  scale_fill_gradient2(low="#D55E00",mid="white",high = "#56B4E9", midpoint=ref) +
  labs(x=paste("% of total variance of HbA1c due to measurement error"),
       y=paste("% of total variance of BMI due to measurement error")) +
  coord_equal()+
  scale_y_continuous(breaks=unique(tot.mat[,1]))+
  scale_x_continuous(breaks=unique(tot.mat[,1]))+
  theme(panel.background = element_rect(fill='white', colour='grey'),
        plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10))

print(FIGURE)

out_dir <- here::here("Outputs") # Define output directory
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE) # Create output directory if it doesn't exist

ggsave( # Save the figure
  filename = file.path(out_dir, "Figure_STRATOS.tif"),
  plot     = FIGURE,
  dpi      = 300,
  width    = 160,   
  height   = 120,  
  units    = "mm"
)
