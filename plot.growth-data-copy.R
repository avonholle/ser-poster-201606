# plot.growth-data-manuscript.R

# Note: adapted from plot.growth-data-copy.R in ~Dropbox\unc.grad.school\presentations\Epi-congress-of-americas-2016\ser-poster-201606

# Make data to plot the growth curves for each of the countries
# and accompanying interaction term change
# and different group differences at baseline (birth)

# Simulate a Reed growth curve as per page 211 of 
# Methods in Human Growth Research, isbn 9780511542411

# .....................................................

library(reshape2)
library(ggplot2)
library(car)
library(lattice)
library(nlme)
library(sitar)
require(data.table)
library(plyr)

## @knitr plotg-data

# Note: this macro comes from http://www.who.int/childgrowth/software/en/
# will set up the igrowup_standard-mod.r function
source("run-who.R")

# 1) Parameters for model
# ............................................
# Get parameters for growth curves
# ................................................................

# Parameters from following papers: http://doi.wiley.com/10.1111/rssa.12020 (On modelling early life weight trajectories, Table 2)
# and
# http://doi.wiley.com/10.1002/ajhb.22575 (Postnatal Growth Velocity and Overweight in Early Adolescents: A Comparison of Rural and Urban African Boys and Girls, Table 2)

param.vals.boys = 
  t(matrix(
    c( 12.52,  10.95,  19.80,
      0.85,  0.66,	3.29,
      -0.08,	1.48,	-10.47,
      -10.25,	-8.28,	-19.74),
    nrow=4,
    ncol=3, byrow=T))
param.vals.boys

param.vals.girls=
  t(matrix(
    c(9.72,	6.76,	16.87,
      0.79,	-0.41,	3.20,
      1.33,	5.87,	-8.96,
      -7.49,	-3.22,	-16.77),
    nrow=4,
    ncol=3, byrow=T))

p.df = cbind(param.vals.boys,param.vals.girls)

p.df = as.data.frame(p.df)
colnames(p.df) = c("a.1", "b.1", "c.1", "d.1", "a.2", "b.2", "c.2", "d.2")
p.df$study=rownames(p.df)
p.df

sigmas = c(0.75)
gends = 2 # 1=male, 2=female

# Function to determine beta4 paramter given the interaction between slopes, beta5
# ------------------------------------------------------------


# original function to generate outcome for girls in the Pizzi paper
# yij = a.2 + b.2*t.star + c.2*log(t.star) + d.2/t.star + e*x.1 + f.1*t.star*x.1 + eij, # Pizzi

# Some definitions:           
#   e = beta4, baseline difference term
#   f.1 = beta5, interaction term
#   t.star = as.numeric((time+9)/9) 
#   f.1.star = f.1 /( (1+9)/9 - 1) # Note: a one unit change in months is equivalent to a 0.11 change in this scale -- so 1/0.11*0.1 = 0.9 should convert to 0.1 change on time scale

# If you want the mean intercept to be the same for each x.1 group at baseline, t=0 and t.star=1, then 
# should have e*x.1 + f.1*t.star*x.1 = 0
# Then e = -f.1*t.star

# if f.1 = 0.5 and t.star = (0+9)/9 = 10/9 = 1
# so e = -0.5 * 1 = -0.5

get.int = function(b5) {
  e = -b5*1
  return(e)
}

get.int(0.50) # -0.50
get.int(-0.50) # 0.50

# combine the differences at baseline in the first col and
# the interaction term in the second col
# Note: for the latter models based on the Pizzi paper you have to adjust the baseline value so that 
# you get the values at birth to be equal. Since baseline is approximately time of conception in these models
# you need to change the baseline to a different one than a model with a different interaction term.

base.int = matrix(c(0, 0.5,
                    -0.50, 0.50,
                    0.50, -0.50,
                    0, -0.5), # beta4 (first col) and beta5 (2nd col) set so that exposure groups have similar baseline weight values
                  nrow=4, byrow=T) 

# Create arrays to store data from for loops below.
store.1 = array(list(NULL), dim = c(length(unique(sigmas)), nrow(p.df), nrow(base.int)))
dim(store.1)

# START OF CODE TO CREATE DATA ---------------------------------------------
# NOTE: I cannot nest the who z score function into one of my functions.
#       Instead I will nest in for loops.
# -------------------------------------------------------------------------

# 2) Positive diff at baseline: 'high' group is 0.50 kg heavier than 'low' group at birth
# -------------------------------------------------------------------------------------
set.seed(1234)

# Run for loop to generate data and corresponding estimates from regression models
for(params in 1:nrow(p.df)){
  #params=5
  param.vals = p.df[params, 1:8] # extract out parameters for growth curves for each country
  study = p.df[params,9]
  for(sig in 1:length(unique(sigmas))){ # iterate over different levels of variance
    # sig=1
    sigval = unique(sigmas)[sig]
    for(basevals in 1:nrow(base.int)){ # extract out values for baseline and interaction for weight in model
      #basevals=3
      base.int.vals = base.int[basevals,]
      
      iter=1; num=100 # num is total number of people in simulation
        
        # 4c) Coefficients for Reed model. See http://stackoverflow.com/questions/16276667/using-apply-with-assign-in-r
        Vars <- c("a.1", "b.1", "c.1", "d.1", "a.2", "b.2", "c.2", "d.2")
        Dat <- as.vector(p.df[params,])
        for (i in 1:length(Vars)){
          assign(Vars[i], as.numeric(Dat[i]))
        }
        
        e = base.int.vals[1] # effect size for cov.1 (in this case it's in kg units)
        
        # create data following modified reed growth function with autocorrelated time
        # ..........................
        
        ### maximum number of possible observations
        m <- 13 # m time points
        n = num # n people (n/2 per gender group)
        
          # 2) Create error/residuals with autocorrelated time, rho=0.5
          # ..................................................
          
          ### simulate number of observations for each individual
          p <- round(runif(n,12,m)) # do this if you want missing observations
          
          ### simulate observation moments [time points] (assume everybody has 1st obs)
          time <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))-1 # for some reason this part doesn't work when I specify all 13 time points
          #length(time)
          ### set up data frame
          dat <- data.frame(id=rep(1:n, times=p), time=time)
          
          # see http://stats.stackexchange.com/questions/76999/simulating-longitudinal-lognormal-data-in-r
          ar.val=0.5 # this is the rho value
          sigma=sigval
          dat$eij <- unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2) * sigma))
          
          # 3) Make a dichotomous indicator variable to be in the growth model
          # .....................................................................
          count.1 = length(unique(dat$id))
          x.1 = ifelse(runif(count.1)<0.5,1,0)
          
          # 4) Make gender variable; 1=boy, 2=girl. According to who function required values for sex of child.
          gender = rbinom(count.1,1,0.5)+1 # make equally distributed boys and girls in data frame, boys=1, girls=2
          dat.2 = data.frame(id=1:count.1, x.1=x.1, gender=gender) 
          dat = merge(dat, dat.2, by="id") # Add the x.1 covariate onto the data frame dat
          
          # 4) Generate outcome for boys and girls, weight in kg, with autocorrelated errors -- following reed model
          # ..................................................................................
          
          # 4a) Make boy and girl indicators for each point in dat data frame
          dat$boys = grepl(1, dat$gender)
          head(dat)
          dat$yij = NA
          dat$yij.2 = NA
          
          # 4b) Coefficients for interaction terms 
          f.1 = base.int.vals[2]
          f.2 = 0
          f.3 = 0
          
#          f.1.star = f.1 /( (1+9)/9 - 1) # Note: a one unit change in months is equivalent to a 0.11 change in this scale -- so 1/0.11*0.1 = 0.9 should convert to 0.1 change on time scale
    
          # 4d) Make weight outcomes according to model used in each paper -- Chirwa model different than Pizzi model.
          dat = within(dat, {
            # transformation of time because data includes weight at birth, Pizzi model only
            
            t.star = as.numeric((time+9)/9)
            
            yij = ifelse(boys==T, a.1 + b.1*t.star + c.1*log(t.star) + d.1/t.star + e*x.1 + f.1*t.star*x.1 + eij, # Pizzi
                         ifelse(boys==F, a.2 + b.2*t.star + c.2*log(t.star) + d.2/t.star + e*x.1 + f.1*t.star*x.1 + eij, # Pizzi
                                              NA))
            
            # no effect for x.1 group. Use this to get alpha
            yij.2 = ifelse(boys==T, a.1 + b.1*t.star + c.1*log(t.star) + d.1/t.star + eij, #Pizzi
                           ifelse(boys==F, a.2 + b.2*t.star + c.2*log(t.star) + d.2/t.star + eij, # Pizzi
                                  NA))
            
            
            time.2 = log(time+1)
            time.3 = 1/(time+1)
            
            t.star.2 = log(t.star)
            t.star.3 = 1/t.star
            
            x.1f = factor(x.1, labels=c("low", "high"))
          })
          
          # After creating data, convert to weight values to z-scores.
          # .............................................................
          
          # Note: this macro comes from http://www.who.int/childgrowth/software/en/
          # I ran all necessary steps in the program titled, "run-who.R", in this folder.
          zscore=igrowup.standard(mydf=dat, sex=gender, age=time, age.month=T, weight=yij)
          

          # zwei is the zscore for weight
          zscore = within(zscore, {pctile = round(pnorm(zwei)*100,0)})
          zscore$sig = sig
          zscore$params = params
          zscore$evals = basevals
  
        store.1[[sig, params, basevals]] = zscore
    }
  }
}

# END OF SECTION TO CREATE DATA ---------------------------------------------
# -------------------------------------------------------------------------

# Extract out data from the for loop
# ....................................................
dim(store.1)
class(store.1) # store.1 is an array.

# try http://stackoverflow.com/questions/10865337/r-how-to-get-a-value-of-a-multi-dimensional-array-by-a-vector-of-indices
# or
# http://r.789695.n4.nabble.com/How-to-quot-flatten-quot-a-multidimensional-array-into-a-dataframe-td4572108.html

vals.p = adply(store.1, 1:3) 
dim(vals.p)
class(vals.p)
head(vals.p)
with(vals.p, table(evals, params))

# labels for the population and gender groups
levels(factor(vals.p$gender))
levels(factor(vals.p$params))
levels(factor(vals.p$sig))
levels(factor(vals.p$evals))

vals.p = within(vals.p, {
  gender.f = factor(gender, labels=c("Male", "Female"))
  params.f = factor(params, labels=c("Portugal", "Italy", "Chile"))
  sigma.f = factor(sig, labels=c("0.75"))
  evals.f = factor(evals, labels=c("0, 0.5", 
                                   "-0.50 and 0.50", 
                                   "0.50 and -0.50", 
                                   "0 and -0.50"))
  })


summary(vals.p)
dim(vals.p)
names(vals.p)

# Check mean yij by group in females, sigma=0.25, time=0: across all countries and baseline differences
# ............................................................................................

ddply(vals.p[vals.p$gender==1 & vals.p$sig==1 & vals.p$time==0,], 
      .(evals.f, params, x.1f), 
      summarise, 
      mean=mean(yij, na.rm=T)) # Note that the values at baseline for the non-African models are not beta.4 because they use a diff model.

names(vals.p)
table(vals.p$evals)

getwd()
save(vals.p, file="vals.p.Rda")

# Check value of interaction term in model with Z-score response...
sub.1 = with(vals.p, vals.p[gender==1 & sig==1 & evals==3 & params==5,])
head(sub.1)
dim(sub.1)
      
m.1 = gls(zwei ~ time*x.1, 
          correlation=corAR1(form = ~ 1 | id),
          data=sub.1) # Z-scores
summary(m.1)$tTable
