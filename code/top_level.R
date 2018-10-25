require(readr)
require(readxl)
require(plyr)
require(dplyr)
require(tidyr)
require(ggplot2)

# --------------------------------------------------------------------
# Data Preparation
# --------------------------------------------------------------------
# Load in main data source & rename columns
happiness <- read_excel("data/read/online-data-chapter-2-whr-2017.xlsx", 
                        sheet = "Data behind Table 2.1 WHR 2017")
happiness_names <- c("wp5_country","country","year","life_ladder",
                     "log_gdp_percapita","social_support",
                     "healthy_life_expectancy","freedom_of_choice","generosity",
                     "corruption_perception","positive affect",
                     "negative affect","govt_confidence","democratic_quality",
                     "delivery_quality","sd_lifeladder",
                     "sd_over_mean_lifeladder","gini_ind","gini_ind_avg",
                     "gine_house_income","trust_gallup","trust_81to84",
                     "trust89to93","trust94to98","trust99to04",
                     "trust05to09","trust10to14")
happiness <- as_data_frame(happiness)
colnames(happiness) <- happiness_names

# Load in data source with regional labels
regions <- read_delim("data/read/regions.csv", delim=",")
colnames(regions) <- c("region", "country")

# Join the main data set with the regional labels,
# Extract another data set with just 2016
hap_df <- left_join(happiness, regions, by='country')
hap2016 <- filter(hap_df, year==2016)
happiest2016 <- arrange(hap2016,desc(life_ladder))[1:10,]
happiest2016_countries <- as.matrix(happiest2016[,"country"])
happiestAll <- hap_df[,c("country","year","life_ladder","region")] %>%
  filter(country %in% happiest2016_countries)


# --------------------------------------------------------------------
# Visualization
# --------------------------------------------------------------------
# Print out wide format, uncleaned data
happiestAll_wide <- happiestAll %>% 
  spread(key=year,value=life_ladder) %>%
  arrange(desc(`2016`)) %>%
  dplyr::select(Country=country,Region=region,everything())
write_csv(format(happiestAll_wide, digits=3, nsmall=2),
          "data/write/happiest.csv")

# Plot all longitudinal data
ggplot(data=hap_df, mapping=aes(x=year,y=life_ladder)) +
  geom_point(aes(color=wp5_country),show.legend=FALSE) +
  geom_line(aes(color=wp5_country),show.legend=FALSE) +
  facet_wrap(~ region, nrow=2)

# Plot all longitudinal data, highlight happiest
ggplot(data=hap_df, mapping=aes(x=year,y=life_ladder,group=country)) +
  geom_point(color='darkgrey',alpha=0.7,show.legend=FALSE) +
  geom_line(color='darkgrey',alpha=0.7,show.legend=FALSE) +
  geom_point(data=happiestAll,aes(x=year,y=life_ladder,color=country)) +
  geom_line(data=happiestAll,aes(x=year,y=life_ladder,color=country))

# Plot all longitudinal data, highlight USA
usa <- hap_df %>% filter(country=="United States")
ggplot(data=hap_df, mapping=aes(x=year,y=life_ladder,group=country)) +
  geom_point(color='darkgrey',alpha=0.7,show.legend=FALSE) +
  geom_line(color='darkgrey',alpha=0.7,show.legend=FALSE) +
  geom_point(data=usa,aes(x=year,y=life_ladder),color='blue') +
  geom_line(data=usa,aes(x=year,y=life_ladder),color='blue')

# Plot all longitudinal data
nearestvalue <- matrix(c(2008,2009,7.524521,7.524521),ncol=2)
colnames(nearestvalue) <- c("year","life_ladder")
nearestvalue <- cbind(country="Switzerland",as.data.frame(nearestvalue))
interpolated <- matrix(c(2013,2013,2014,7.634506,7.561374,7.499732),ncol=2)
colnames(interpolated) <- c("year","life_ladder")
interpolated <- cbind(country=c("Switzerland","Norway","Iceland"),
                      as.data.frame(interpolated))
happiestAll %>% filter(year >= 2008) %>%
  filter(country %in% c("Switzerland","Iceland","Norway")) %>%
  ggplot(mapping=aes(x=year,y=life_ladder)) +
  geom_point(aes(color=country),show.legend=FALSE) +
  geom_line(aes(color=country)) +
  geom_point(data=nearestvalue[nearestvalue$year==2008,], 
             aes(x=year,y=life_ladder,color=country),
             shape=17,size=3) + 
  geom_line(data=nearestvalue, aes(x=year,y=life_ladder,color=country),
            linetype="dashed") +
  geom_point(data=interpolated,aes(x=year,y=life_ladder,color=country),
             shape=17,size=3)
  

# Plot pairs of predictors and response variables
pairs(hap2016[,c(4:7)])
pairs(hap2016[,c(4,8:10)])
pairs(hap2016[,c(4,11:13)])


# --------------------------------------------------------------------
# Regression Work
# --------------------------------------------------------------------
# Subset data
hap_2016_r <- hap2016[,c('life_ladder','log_gdp_percapita','social_support',
                         'healthy_life_expectancy','freedom_of_choice',
                         'generosity','corruption_perception')]
# Set initial parameters
B <- 20000

hap_2016_r <- hap_2016_r[complete.cases(hap_2016_r),]
cor(hap_2016_r)

y <- as.matrix(hap_2016_r[,1])
X <- hap_2016_r[,-1]

results <- controller(X,y,B)
results$summary

# Burn-in is done in the regress function output
# Plot convergence
par(mfrow=c(1,1))
matplot(results$converge,type='l',ylim=c(-0.1,0.1),
        ylab='Running Mean')

# Plot regression parameter distributions
par(mfrow=c(3,2))
lnames <- names(results$densities)
for (i in 2:7) {
  plot_dist(lnames[i],results$summary[3:4,i],results$densities[[i]])
}

# Compute DIC
dic <- dic_norm(y,results)
dic


# --------------------------------------------------------------------
# Longitudinal Mixed Model Work
# --------------------------------------------------------------------
# REQUIRES: cleaned data set: happiestAll
require(MCMCpack)
require(mvtnorm)
require(Matrix)

# Function to summarize based on column vectors
summary_coeff_row <- function(coeff_mat) {
  m1 <- apply(coeff_mat,1,mean)
  m2 <- apply(coeff_mat,1,median)
  ci <- apply(coeff_mat,1,quantile,probs=c(0.025,0.975))
  s <- rbind(mean=m1, median=m2, ci)
  return(s)
}

# First we reshape the data to ensure observations in
# the same years across countries

# This interpolates the life_ladder data from 2012, 2014
# for Norway and Switzerland to populate a value for 2013
i2013 <- happiestAll %>%
  filter(country %in% c("Norway","Switzerland")) %>%
  filter(year %in% c(2012,2014)) %>%
  group_by(country) %>%
  dplyr::summarize(life_ladder=mean(life_ladder)) %>%
  cbind(year=rep(2013,2)) %>%
  cbind(region=rep("Western Europe",2)) %>%
  dplyr::select(country,year,life_ladder,region)

# This interpolates the life_ladder data from 2013, 2015
# for Iceland to populate a value for 2014
i2014 <- happiestAll %>%
  filter(country == "Iceland") %>%
  filter(year %in% c(2013,2015)) %>%
  group_by(country) %>%
  dplyr::summarize(life_ladder = mean(life_ladder)) %>%
  cbind(year=2014) %>%
  cbind(region="Western Europe") %>%
  dplyr::select(country,year,life_ladder,region)
  
# This binds the observed data with the generated
# observations from above and with 
# taking the 2009 Switzerland value
# for 2008, to complete the data set
happiestAug <- happiestAll %>% 
  rbind(c(country="Switzerland",year=2008,
          happiestAll[(happiestAll$country=="Switzerland" 
                       & happiestAll$year==2009),"life_ladder"],
          region="Western Europe")) %>%
  rbind(i2013) %>%
  rbind(i2014) %>%
  filter(year %in% c(2008,2012,2013,2014,2015,2016)) %>%
  arrange(country, year)

# Output the cleaned data set as wide format CSV
happiestAug_wide <- happiestAug %>% 
  spread(key=year,value=life_ladder) %>%
  arrange(desc(`2016`)) %>%
  dplyr::select(Country=country,Region=region,everything())
write_csv(format(happiestAug_wide, digits=3, nsmall=2),
          "data/write/happiest_clean.csv")

# Pull out the response variable
Y <- as.matrix(happiestAug$life_ladder)

N	<- dim(Y)[1] #Total number of observations: N = nT
n	<- length(unique(happiestAug$country)) # Num of countries
# T is number of years

# Generate the Z matrix by:
# adding intercept column to predictor data
# split final design matrix based on country
# use country blocks in diagonal matrix
ones <- rep(1,N)
X <- happiestAug %>%
  cbind(intercept=ones) %>%
  dplyr::select(country,intercept,year)
Xeach <- split(X,X$country)
Xeach <- lapply(Xeach,"[", i= ,c("intercept","year"))
Xeach <- lapply(Xeach,as.matrix)

Z <- bdiag(Xeach)
ztz	<- t(Z)%*%Z

X <- as.matrix(X[,c("intercept","year")])
xtxi <- solve(t(X)%*%X)

# Initialize storage
B <- 200
beta	<- matrix(0, nrow = 2, ncol = B)
rownames(beta) <- c("intercept","slope")
U	<- matrix(0, nrow = 2*n, ncol = B)
tau2_e	<- vector('numeric', length = B)
tau2_u0	<- vector('numeric', length = B)
tau2_u1	<- vector('numeric', length = B)

# Matrices to predict 2017
hap2017 <- matrix(NA, nrow=n, ncol=B)
rownames(hap2017) <- names(Xeach)
X2017i <- matrix(c(1,2017),nrow=1)
X2017 <- matrix(rep(X2017i,n),ncol=2,byrow=TRUE)
Z2017 <- bdiag(rep(list(X2017i),n))

# Calculate initial values
beta[,1] <- xtxi%*%t(X)%*%Y
U[,1]	<- rep(1, 2*n)
tau2_e[1]	<- 1
tau2_u0[1] <- 1
tau2_u1[1] <- 1
hap2017[,1] <- as.matrix(X2017%*%beta[,1] + Z2017%*%U[,1])

set.seed(1895)
for(t in 2:B){
  #print(t)
  # Generate reused terms
  Xb <- X%*%beta[,t-1]
  ZU <- Z%*%U[,t-1]
  
  ## update beta ##
  vbeta	<- xtxi/tau2_e[t-1]
  mbeta	<- xtxi%*%t(X)%*%(Y - ZU)
  beta[,t]	<- rmvnorm(1, mbeta, vbeta)
  #print(beta[,t])
  
  ## update U ##
  tau_little_vec <- c(tau2_u0[t-1],tau2_u1[t-1])
  tau_full_vec <- rep(tau_little_vec,n)
  #print(tau_full_vec)
  Sigmainv <- diag(tau_full_vec)
  vU <- as.matrix(solve(tau2_e[t-1]*ztz+Sigmainv))
  mU <- as.matrix(tau2_e[t-1]*vU%*%t(Z)%*%(Y-Xb))
  U[,t]	<- rmvnorm(1, mU, vU)
  #print(U[,t])
  
  ## update tau2_e ##
  rate_e <- as.numeric(t(Y - Xb - ZU)%*%(Y - Xb - ZU))
  tau2_e[t]	<- rgamma(1, N/2, rate_e/2)
  #print(tau2_e[t])
  
  ## update tau2_u0 ##
  idx <- seq(1,2*n,by=2)
  u0 <- U[idx,t-1]
  rate_u0 <- sum(u0^2)
  tau2_u0[t] <- rgamma(1, n/2, rate_u0/2)
  #print(tau2_u0[t])
  
  ## update tau2_u1 ##
  u1 <- U[idx+1,t-1]
  rate_u1 <- sum(u1^2)
  tau2_u1[t]	<- rgamma(1, n/2, rate_u1/2)
  #print(tau2_u1[t])
  
  ## predict 2017 ##
  hap2017[,t] <- as.numeric(X2017%*%beta[,t] + Z2017%*%U[,t])
}

# Compute Summaries
beta_s <- summary_coeff_row(beta[,floor(B/2):B])
U_s <- summary_coeff_row(U[,floor(B/2):B])
Y_pred <- X%*%matrix(beta_s[1,],ncol=1) + Z%*%U_s[1,]
happiestAug2 <- cbind(happiestAug,happy_fit=as.matrix(Y_pred))

# Assess Convergence
U_conv <- apply(U,1,function(x) (cumsum(x)/(1:length(x)))-mean(x))
matplot(U_conv,type='l',ylim=c(-0.1,0.1),
        ylab='Running Mean')

beta_conv <- apply(beta,1,function(x) (cumsum(x)/(1:length(x)))-mean(x))
matplot(beta_conv,type='l',ylim=c(-200,200), ylab='Running Mean')

# Generate output data for report
full_beta <- do.call(rbind, replicate(n, beta, simplify=FALSE))
combined_form <- full_beta + U
full_s <- summary_coeff_row(combined_form[,floor(B/2):B])
output_table <- cbind("Intercept"=full_s[1,seq(1,19,2)], 
                      "Slope"=full_s[1,seq(2,20,2)],
                      "Intercept Low"=full_s[3,seq(1,19,2)], 
                      "Intercept High"=full_s[4,seq(1,19,2)],
                      "Slope Low"=full_s[3,seq(2,20,2)], 
                      "Slope High"=full_s[4,seq(2,20,2)])
output_table <- cbind(Country=names(Xeach),as.data.frame(output_table))
write_csv(format(output_table, digits=3, nsmall=2),
          "data/write/output_long_params.csv")

# Visualize each fitted line
ggplot(data=happiestAug2, mapping=aes(x=year,y=life_ladder)) +
  geom_point(show.legend = FALSE,color="darkgrey") +
  geom_line(show.legend = FALSE,color="darkgrey") + 
  geom_line(mapping=aes(x=year,y=happy_fit),
            linetype="dashed",
            show.legend = FALSE) +
  facet_wrap(~country,nrow=2) +
  theme(panel.spacing.x=unit(0.75, "lines"))

# Posterior Inference
hap2017_burn <- hap2017[,floor(B/2):B]

# (1) The probability that Finland continues 
# to be the happiest country
argmax <- apply(hap2017_burn,2,which.max)
finland_win <- sum(argmax==4)/length(argmax)

# (2) The probability that Finland drops
# to the least happiest country
argmin <- apply(hap2017_burn,2,which.min)
finland_lose <- sum(argmin==4)/length(argmin)

# (3) The probability that New Zealand overtakes Canada
nz_over_canada <- sum(hap2017_burn[7,]>hap2017_burn[2,])/dim(hap2017_burn)[2]

# (4) The probability that Denmark will be the happiest
denmark_win <- sum(argmax==3)/length(argmax)
