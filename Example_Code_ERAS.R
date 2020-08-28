#### Example: R code  to estimate the latent class 2PL model  proposed in the paper titled: ####
# Studying Enhanced Recovery After Surgery (ERAS) core items in colorectal surgery:
# A causal model with latent variables
# Authors: Gemma, M., Pennoni, F. and Braga, M. #
#######################################################

load("sim_dataERAS.Rdata")

# sim_data contains: 
# - Pre-treatment covariates: V2,V3,V4,V5,V6
# (e.g. a dummy for the year 1 of the surgery,
# weight loss, age, diabete, gender)
# - Preoperative-interventions: V1, V7 (e.g. counseling, premedication)
# - Intraoperative-interventions:  V8 (e.g. fluid administration, surgical drainage)

# intervensions
lapply(sim_data[,c(1,7,8)],table)
# pre-tretments covariates
lapply(sim_data[,c(2:3,5,6)],table); summary(sim_data$V4)

#### Estimate of weights considering  the first two blocks of interventions ####
# (Figure 1)  #

### weights: first block ###
require(nnet)
res1 <- multinom(V1 ~ V2 + V3 + V4 + V5, sim_data)
res2 <- multinom(V7 ~   V2 + V3 + V4 + V5, sim_data)
P1 <- res1$fitted.values; summary(P1)
P2 <- res2$fitted.values; summary(P2)
w1 <- 1/((sim_data$V1+1*P1)*(sim_data$V7+1*P2))
### weights: second block ### 
# example with one variable e.g. surgical drainage (V7)
res21 <- multinom(V8 ~ V2 + V3 + V4 + V5 + V6 + V1 + V7, sim_data)            
P21 <- res21$fitted.values
w2 <- 1/(sim_data$V8+1*P21)
# trim weights
w2 <- pmin(w2,5)
# combine the weights 
w <- as.vector(w1+w2)


## Define 2 dimensions for the 6 items  ###
# dimension 1:  by items 1,2,3 (e.g. complications,SSI, medical com.)
# dimension 2:  by items 4,5 (e.g. time ready for discharge, actual discharge)

multi1 <- rbind(c(1,2,3),c(4,5,0)); multi1

#### Estimate a two-parameters logistic and select the best number of latent states ####
require(MultiLCIRT)
# S1 contains the observed items
lapply(data.frame(S1),table)
# X1 contains the observed treatments 

#### Determine the suitable number of latent classes ####
ksearch1 <- search.model(S=S1,
                         yv=w, 
                         kv=c(1:3), 
                         X = NULL, 
                         fort = TRUE, 
                         tol = 10^-10, 
                         multi = multi1)

#### LC-2PL causal model with k = 2 (simple example) latent classes ####
X <-cbind(sim_data$V1,sim_data$V7,sim_data$V8)
modpr <- est_multi_poly(S = S1,
                          yv = w,
                          X = X, 
                          multi=multi1,
                          link=2, 
                          disc=1,
                          difl=0,
                          tol = 10^-7,
                          k = 2,
                          start = 0,
                          fort = TRUE,
                          out_se = TRUE,
                          output = TRUE)
# print the results
summary(modpr)
# point on the latent continuum  
modpr$Th
# item difficulty
modpr$Bec
# discriminating power of the item 
modpr$gac
# averaged weights of each LC
modpr$piv
# conditional response probabilities
round(modpr$Phi,3)
# estimated  effects
round(modpr$De,3)


