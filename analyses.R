# Modeling Air Quality and COVID-19 Deaths in South Korea
 
#### Table of Contents
# 1. Background (#bg)
# 2. Data (#dat)
#    i. Data Preparation (#datprep)
# 3. Model (with interaction terms) (#model1)
#    i. Model Code (#code1)
#    ii. Define Model Parameters (#paramdef1)
#    iii. Run the Model (#run1)
# 4. Results (with interaction terms) (#result1)
#    i. Trace plots (#plot1)
#    ii. 95% Highest Density Intervals (HDI's) (#hdi1)
# 5. Model (without interaction terms) (#model2)
#    i. Model Code (#code2)
#    ii. Define Model Parameters (#paramdef2)
#    iii. Run the Model (#run2)
# 6. Results (without interaction terms) (#result2)
#    i. Trace plots (#plot2)
#    ii. 95% Highest Density Intervals (HDI's) (#hdi2)

##=============================================================================
## Set up                                                 
library(nimble)
library(coda)
library(dplyr)
library(HDInterval)
library(ggplot2)
library(kableExtra)

##-----------------------------------------------------------------------------
## Data {#dat}
path = 'data/cleaned/'
covid = read.csv(paste0(path, 'covid.csv'))
covid = covid[, c('province', 'deceased')]
pop = read.csv(paste0(path, 'SKpopulation.csv'))
colnames(pop) <- c('province', 'population')
suiAQ = read.csv(paste0(path, 'suiAQ.csv'))
elderly = read.csv(paste0(path, 'elderly.csv'))
elderlyratio = merge(elderly, pop, by = 'province')
elderlyratio = mutate(elderlyratio, elderlyratio = elderly / population)
elderlyratio = elderlyratio[, c('province', 'elderlyratio')]
socio = read.csv(paste0(path, 'socio.csv'))
socio = socio[, c('province', 'university_count')]
data = merge(covid, pop, by = 'province') %>% merge(suiAQ) %>% 
  merge(elderlyratio) %>% merge(socio)


### Data Preparation {#datprep}
## define the data
y <- data[, 'deceased'] # response
x <- data[, c('co', 'pm10', 'no2', 'ozone', 'so2', 
              'elderlyratio', 'university_count')] # predictors
population <- data[, 'population']
nX <- ncol(x) # number of covariates
nTotal <- nrow(x) # number of observations

## standardize the data
xMean <- rep(0, nX)
xSD <- rep(0, nX)
zX <- matrix(0, nrow = nTotal, ncol = nX) # standardized x
for (j in 1:nX) {
  xMean[j] <- mean(x[, j])
  xSD[j]   <-   sd(x[, j])
  for (i in 1:nTotal) {
    zX[i, j] <- (x[i, j] - xMean[j]) / xSD[j]
  }
}

## find the interacting terms
interact <- abs(cor(zX)) > 0.6 # indicator matrix of highly correlated covariates 
ind <- which(interact == TRUE, arr.ind = T); rownames(ind) <- 
  NULL # index  pairs of highly correlated covariates
uniq <- unique(data.frame(t(apply(ind, 1, sort)))) # unique index pairs
interactInd <- uniq[uniq$X1 != uniq$X2, ] # excluding interactions with oneselves
nInt <- nrow(interactInd) # number of interacting terms
zX1 <- zX[, interactInd$X1]  
zX2 <- zX[, interactInd$X2]

## define the initial values
f <- 'y~'
for (i in 1:nX) { f <- paste0(f,'zX[,',as.character(i),']+') }
for (r in 1:nInt) { f <- paste0(f,'zX[,',as.character(interactInd[r,1]),
                                ']*zX[,',as.character(interactInd[r,2]),']+') }
f <- paste0(f,'offset(log(population))') # formula
fit <- glm(formula = f, family = poisson, control = list(maxit = 100)) 
beta0Init <- coef(fit)[1] # intercept
betaInits <- coef(fit)[2:(nX+1)] # beta coefficients of the regular terms
betaIInits <- coef(fit)[(nX+2):(nX+nInt+1)] # beta coefficients of the Interacting terms
lambdaInits <- rep(mean(y), nTotal)


##-----------------------------------------------------------------------------
## Model (with interaction terms) {#model1}
### Model Code {#code1}
interactCode <- nimbleCode({
  # likelihoods
  for (i in 1:nTotal) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + sum(beta[1:nX] * zX[i, 1:nX]) + 
                              sum(betaI[1:nInt] * (zX1[i, 1:nInt] * zX2[i, 1:nInt])) + 
                              log(population[i])
  }
  # priors
  beta0 ~ dnorm(0, 1) # intercept
  for (j in 1:nX) {
    beta[j] ~ dnorm(0, 1) # beta coefficients of the regular terms
  }
  for (r in 1:nInt) {
    betaI[r] ~ dnorm(0, 1) # beta coefficients of the Interacting terms
  }
})


### Define Model Parameters {#paramdef1}
## data
interactData <- list(zX = zX, zX1 = zX1, zX2 = zX2, 
                     y = y, population = population, xMean = xMean, xSD = xSD)

## constants of the model 
interactConsts <- list(nTotal = nTotal, nX = nX, nInt = nInt, 
                       interactInd = interactInd)

## values to initialize the algorithm
interactInits <- list(beta0 = beta0Init, beta = betaInits, 
                      betaI = betaIInits, lambda = lambdaInits)

### Run the Model {#run1}
outInteract <- nimbleMCMC( code = interactCode, constants = interactConsts,
                           data = interactData, inits = interactInits,
                           niter = 22000, nburnin = 2000,
                           thin = 20, nchains = 1,
                           monitors = c("beta0", "beta", "betaI"),
                           samplesAsCodaMCMC = TRUE)


##-----------------------------------------------------------------------------
## Results (with interaction terms) {#result1}
summary(outInteract)


### Trace plots {#plot1}
plot(outInteract)


### 95% Highest Density Intervals (HDI's) {#hdi1}
## calclulate the HDI's
intHDI <- data.frame(t(apply(outInteract, 2, hdi))) 
intHDI <- intHDI[-8,] # without the HDI of beta0
intHDI <- cbind(X = rownames(intHDI), intHDI) 

## name the covariates to be plotted 
shrtColnames <- c(colnames(x)[1:3], 'ozon', 'so2', 'elder','univ')
predictors <- c(shrtColnames, paste0(shrtColnames[interactInd$X1],'*',
                                     shrtColnames[interactInd$X2]))

## plot the HDI's
ggplot(intHDI, aes(x=X)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,size=1, color='black') +
  geom_hline(yintercept = 0, linetype="dashed", size=0.8) +
  coord_cartesian(ylim = c(-4, 7), xlim = c(0.5, 11.5), expand = FALSE, clip = "off") +
  ylab("beta") +
  annotate(geom = "text", x = seq_len(nrow(intHDI)), y = -4.2, 
           label = predictors, size = 3) +
  theme_bw() +
  ggtitle('95% Highest Posterior Density Intervals') + 
  theme(plot.margin = unit(c(1, 1, 1, 0), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        legend.position = "none")
  
  
##-----------------------------------------------------------------------------
## Model (without interaction terms) {#model2}
### Model Code {#code2}
woInteractCode <- nimbleCode({
  # likelihoods
  for (i in 1:nTotal) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + sum(beta[1:nX] * zX[i, 1:nX]) + log(population[i])
  }
  # priors
  beta0 ~ dnorm(0, 1) # intercept
  for (j in 1:nX) {
    beta[j] ~ dnorm(0, 1) # beta coefficients of the regular terms
  }
})


### Define Model Parameters {#paramdef2}
## data
woInteractData <- list(zX = zX, y = y, population = population, 
                       xMean = xMean, xSD = xSD)

## constants of the model 
woInteractConsts <- list(nTotal = nTotal, nX = nX)

## values to initialize the algorithm
woInteractInits <- list(beta0 = beta0Init, beta = betaInits, lambda = lambdaInits)


### Run the Model {#run2}
outWoInteract <- nimbleMCMC( code = woInteractCode, constants = woInteractConsts,
                             data = woInteractData, inits = woInteractInits,
                             niter = 22000, nburnin = 2000,
                             thin = 20, nchains = 1,
                             monitors = c("beta0", "beta"),
                             samplesAsCodaMCMC = TRUE)


##-----------------------------------------------------------------------------
## Results (without interaction terms) {#result2} 
summary(outWoInteract)


### Trace plots {#plot2}
plot(outWoInteract)


### 95% Highest Density Intervals (HDI's) {#hdi2}
## calclulate the HDI's
woIntHDI <- data.frame(t(apply(outWoInteract, 2, hdi))) 
woIntHDI <- woIntHDI[-8,] # without the HDI of beta0
woIntHDI <- cbind(X = rownames(woIntHDI), woIntHDI) 

## name the covariates to be plotted 
shrtColnames <- c(colnames(x)[1:3], 'ozon', 'so2', 'elder','univ')

## plot the HDI's
ggplot(woIntHDI, aes(x=X)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,size=1, color='black') +
  geom_hline(yintercept = 0, linetype="dashed", size=0.8) +
  coord_cartesian(ylim = c(-4, 7), xlim = c(0.5, 7.5), expand = FALSE, clip = "off") +
  ylab("beta") +
  annotate(geom = "text", x = seq_len(nrow(woIntHDI)), y = -4.2, 
           label = shrtColnames, size = 3) +
  theme_bw() +
  ggtitle('95% Highest Posterior Density Intervals') + 
  theme(plot.margin = unit(c(1, 1, 1, 0), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
        legend.position = "none")