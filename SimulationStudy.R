library(mvtnorm)
library(dplyr)
library(car)
library(reshape2)
library(ggplot2)

n<- 800 # sample size per simulation
n_pfas<- 2 # number of simulated PFAS
nsim<- 10000 # number of simulations

# correlation strength 
approx_correlation = 0.9

sigma_sim<- sqrt(approx_correlation) + diag(n_pfas)
diag(sigma_sim)<- 1

# Create coefficients for the PFAS effects
set.seed(678020) # change this to do different values of coefficients in simulation
#I assume all PFAS are negatively associated
pfas_betas<- -runif(n_pfas, 0, 0.015)
if(n_pfas>=3){
  pfas_betas[3]<- 0 # setting PFAS3 conditional effect ot 0
}

X_PFAS = rmvnorm(n,mean=rep(0,n_pfas), sigma=sigma_sim)

# Create simulated PFAS values with a realistic range/distribution of values
PFAS_MeanSimulated<- list()
for(j in 1:n_pfas){
  PFAS_MeanSimulated[[j]] = 8+exp(rnorm(n, 1 + X_PFAS[,j],0.25) ) 
}
X_PFAS_Simulated<- do.call(cbind, PFAS_MeanSimulated)

cor(X_PFAS_Simulated) # correlations between the variables range from 0.7 to 0.8
cor(X_PFAS_Simulated, method = "spearman") # correlations between the variables range from 0.7 to 0.8

X_PFAS_Simulated %>% hist(main='Histogram of Simulated PFAS Values', 20, xlim=c(0,100))

# Create a simple log-linear regression model for the response, as in paper
# create 2 uncorrelated X variables, for example, gender and booster_type
x_sex<- rbinom(n,1,  0.5)
x_boostertype<- rbinom(n, 1, 0.5)

betas<- c(0.045, 0.015, 0.015,  pfas_betas)
X<- cbind(1, x_sex, x_boostertype, X_PFAS_Simulated) # Design-matrix, intercept, sex, boosttype, then PFAS (correlated)

#y<- exp(X %*% betas + rnorm(n, 0, 0.25))  # generate y values for this simulation
#hist(y) # range of values seems reasonable 

Xpfas<- as.data.frame(X[,4:(3+n_pfas)])
colnames(Xpfas)<- paste0("PFAS",1:n_pfas)

simdata<- data.frame(y=y, sex=X[,2], boostertype=X[,3]
)
simdata<- cbind(simdata, Xpfas) # design matrix for simulations

# Do repeated simulations

BMD_Results<- list()
for(i in 1:nsim){
  y<- exp(X %*% betas + rnorm(n, 0, 0.25)) 

  Xpfas<- as.data.frame(X[,4:(3+n_pfas)])
  colnames(Xpfas)<- paste0("PFAS",1:n_pfas)
  
  simdata<- data.frame(y=y, sex=X[,2], boostertype=X[,3]
  )
  simdata<- cbind(simdata, Xpfas)
  
  lm_mod1<- lm(log(y)~sex+boostertype + PFAS1 + PFAS2, data=simdata) #"True" Model, still has issues with multicollinearity
  lm_mod2<- lm(log(y)~sex+boostertype + PFAS1 + log(PFAS2), data=simdata) # Proposed Model, has issues with multicollinearity
  lm_mod_single<- lm(log(y)~sex+boostertype + PFAS1, data=simdata) # Unadjusted model, marginalizes the relationship among predictors

  tmp<- cbind(
    # calculate BMD
  log(0.95)/coef(lm_mod_single)[4],
  log(0.95) / coef(lm_mod1)[4],
  log(0.95) / coef(lm_mod2)[4],

  # calculate BMDL
  log(0.95)/ (coef(lm_mod_single)[4] - 1.96*summary(lm_mod_single)$coefficients["PFAS1","Std. Error"]),
  log(0.95) / (coef(lm_mod1)[4] - 1.96*summary(lm_mod1)$coefficients["PFAS1","Std. Error"]),
  log(0.95) / (coef(lm_mod2)[4] - 1.96*summary(lm_mod1)$coefficients["PFAS1","Std. Error"])
  )
  
  BMD_Results[[i]]<- tmp
}

BMD_Results<- do.call(rbind, BMD_Results)
colnames(BMD_Results)<- c("BMD_Unadjusted","BMD_Adjusted","BMD_GeneratingModel",
                          "BMDL_Unadjusted", "BMDL_Adjusted", "BMDL_GeneratingModel")
BMD_Results<- as.data.frame(BMD_Results)

par(mfrow=c(3,2))
hist(BMD_Results$BMD_Unadjusted)
hist(BMD_Results$BMDL_Unadjusted)
hist(BMD_Results$BMD_Adjusted)
hist(BMD_Results$BMDL_Adjusted)
hist(BMD_Results$BMD_GeneratingModel)
hist(BMD_Results$BMDL_GeneratingModel)

colMeans(BMD_Results) # Summarise the average value estimated from the 3 models and 2 quantities

#####################################################################################################################
# Results seem to confirm that the proposed method is likely inappropriate. Make a higher-quality graphic for appendix #
#####################################################################################################################

bmd_df<- melt(BMD_Results) %>%
  mutate(
    model = case_when(
      grepl("_Unadjusted", variable) ~ "Unadjusted",
      grepl("_Adjusted", variable) ~ "Adjusted",
      TRUE ~ '"True" Model'
    ),
    estimate_type = if_else(grepl("BMD_", variable), "BMD","BMDL"))

# Histogram of BMD estimates, capped at +20 or -20
bmd_df %>% filter(estimate_type=="BMD") %>%
  mutate(value = 
  if_else(value > 20,
          20,
          if_else(value< (-20),
                  -20,
                  value))
) %>% 
  ggplot(aes(x=value)) +
  geom_histogram() + 
  facet_grid(.~model, scales='free_x')  #+ xlim(0,11)

BMD_est<- bmd_df %>% filter(estimate_type=="BMD") %>%
  mutate(value =
           if_else(value > 50,
                   50,
                   if_else(value< (-50),
                           -50,
                           value))
  ) %>%
  ggplot(aes(y=value, x=model, fill=model)) +
  geom_boxplot() + xlab("Log-Linear Model Fit") +
  ylab("Estimated Benchmark Dose (Simulated Data)") + 
  ggtitle("Benchmark Dose Estimate\n(10,000 Simulations)") + 
  scale_fill_brewer("Model Type",type='qual', palette = 'Dark2')
BMD_est

ggsave("BMD_Estimate.png", BMD_est, dpi=300,width=7,height=6)

BMDL_est <- bmd_df %>% filter(estimate_type=="BMDL") %>%
  mutate(value = 
           if_else(value > 20,
                   20,
                   if_else(value< (-20),
                           -20,
                           value))
  ) %>% 
  ggplot(aes(y=value, x=model, fill=model)) +
  geom_boxplot() + xlab("Log-Linear Model Fit") +
  ylab("Estimated Benchmark Dose Lower Confidence Limit") + 
  ggtitle("Benchmark Dose Lower Confidence Interval Limit Estimate\n(10,000 Simulations)") + 
  scale_fill_brewer("Model Type",type='qual', palette = 'Dark2')

BMDL_est

ggsave("BMDL_Estimate.png", BMDL_est, dpi=300,width=7,height=6)
# no limits

BMDL_est <- bmd_df %>% filter(estimate_type=="BMDL") %>%
  # mutate(value = 
  #          if_else(value > 20,
  #                  20,
  #                  if_else(value< (-20),
  #                          -20,
  #                          value))
  # ) %>% 
  ggplot(aes(y=value, x=model, fill=model)) +
  geom_boxplot() + xlab("Log-Linear Model Fit") +
  ylab("Estimated Benchmark Dose Lower Confidence Limit") + 
  ggtitle("Benchmark Dose Lower Confidence Interval Limit Estimate\n(10,000 Simulations)") + 
  scale_fill_brewer("Model Type",type='qual', palette = 'Dark2')
BMDL_est
ggsave("BMDL_Estimate_NoLimit.png", BMDL_est, dpi=300,width=7,height=6)



BMD_est<- bmd_df %>% filter(estimate_type=="BMD") %>%
  # mutate(value =
  #          if_else(value > 50,
  #                  50,
  #                  if_else(value< (-50),
  #                          -50,
  #                          value))
  # ) %>%
  ggplot(aes(y=value, x=model, fill=model)) +
  geom_boxplot() + xlab("Log-Linear Model Fit") +
  ylab("Estimated Benchmark Dose Lower Confidence Limit") + 
  ggtitle("Benchmark Dose Estimate\n(10,000 Simulations)") + 
  scale_fill_brewer("Model Type",type='qual', palette = 'Dark2')
BMD_est

ggsave("BMD_Estimate_NoLimit.png", BMD_est, dpi=300,width=7,height=6)

# All plots together, capped at +-50
BMD_both_est<- bmd_df %>% 
  mutate(value =
           if_else(value > 50,
                   50,
                   if_else(value< (-50),
                           -50,
                           value))
  ) %>%
  ggplot(aes(y=value, x=model, fill=model)) +
  geom_boxplot() + xlab("Log-Linear Model Fit") +
  ylab("Estimated Benchmark Dose Lower Confidence Limit") + 
  ggtitle("Benchmark Dose Lower Confidence Interval Limit Estimate\n(10,000 Simulations)") + 
  scale_fill_brewer("Model Type",type='qual', palette = 'Dark2') + 
facet_grid(estimate_type~., scales='free')  #+ xlim(0,11)

BMD_both_est
ggsave("BMDL_Estimate_NoLimit.png", BMD_est, dpi=300,width=7,height=6)

# plot means with +-1sd as barchart
bmd_agg<- bmd_df %>% 
  group_by(estimate_type, model) %>%
  summarise(meanvalue=mean(value), 
            sd=sd(value)) 

bmd_agg %>% filter(estimate_type=="BMDL") %>%
  ggplot(aes(y=meanvalue, x=model, fill=model)) +
  geom_bar(stat='identity') + 
  geom_errorbar(aes(ymin=meanvalue-sd, ymax=meanvalue+sd), width=.2,
                position=position_dodge(.9)) +
  xlab("Log-Linear Model Fit") + ylab("Estimated Benchmark Dose Lower Confidence Limit") + 
  ggtitle("Benchmark Dose Lower Confidence Interval Limit Estimate\n(100,000 Simulations)") + 
  scale_fill_brewer("Model Type",type='qual', palette = 'Dark2')

