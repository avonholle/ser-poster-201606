## @knitr plotg-baselines

# plot.growth-baselines-copy.R

# Plot the growth curves for each of the countries
# and accompanying interaction term change
# grid will be by country and difference at baseline between 2 groups

# data from plot.growth-data.R

# .....................................................

library(reshape2)
library(ggplot2)
library(car)
library(lattice)
library(nlme)
#library(sitar)
require(data.table)
require(plyr)

# packages for parallel processing
library(foreach)
library(doParallel)

library(gridExtra)
library(grid)

# Load data from plot.growth-data.R, vals.p object
# .................................


load("vals.p.Rda")

# re-label the values for the parameter combinations
levels(vals.p$evals.f)
levels(vals.p$x.1f)


vals.p = within(vals.p, {
  evals.f = factor(evals, labels=c("Unequal baseline \n + slope",
                                   "Equal baseline \n + slope",
                                   "Equal baseline \n - slope",
                                   "Unequal baseline \n - slope"))
  x.1f = factor(x.1f, labels=c("unexposed", "exposed"))
})

# Plot with re-transformed time values
# .........................................

# get subset with countries having parameters matching correct model.
sub.vals.p = vals.p[vals.p$gender==2 & 
                      vals.p$sig %in% c(1) &
                      vals.p$params.f %in% c("Chile"),]
head(sub.vals.p)

# pick out first five individuals from groups of exposure (x.1f) and model sim types (params.f)
test = sub.vals.p[,c("id", "x.1", "evals")]
head(test)

d = data.table(unique(test), key=c("id", "x.1", "evals"))
id.sub = d[, head(.SD,5), by = list(x.1, evals)]
nrow(id.sub)
id.sub$new.id = with(id.sub, paste(id, ".", x.1, ".", evals, sep=""))
id.sub$new.id

sub.vals.p$new.id = with(sub.vals.p, paste(id, ".", x.1, ".", evals, sep=""))
random.sub = sub.vals.p[sub.vals.p$new.id %in% id.sub$new.id,]
length(unique(random.sub$new.id))

random.sub$params.f = factor(random.sub$params.f)
levels(random.sub$params.f)

# Make plot for females (in Chile) here
# with random 5 people per exposure group and fitted line
p.f = ggplot(data = random.sub,
             aes(x=time, y=yij, colour=x.1f)) +
  geom_line(aes(group=id)) + 
  facet_grid(evals.f ~ params.f) +
  scale_colour_manual("Group",
                      values = c("unexposed" = "blue",
                                 "exposed" = "red")) +
  scale_y_continuous(breaks=seq(0,12,by=4)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  ylab("Weight (kg)") +
  xlab("Time (months)") +
  theme_bw(base_size=26) +
  theme(legend.position="bottom",
        strip.text.y = element_text(size = 20))
p.f


# NOTE: just plotting main effect so no variance .. the random 5 growth curves can give an idea of the variance
# transformation of time because data includes weight at birth, Pizzi model only

param.vals.girls=
  t(matrix(
    c(4.44, 5.13, 9.72,	6.76,	16.87,
      0.13,	0.14, 0.79,	-0.41,	3.20,
      1.16,	0.48, 1.33,	5.87,	-8.96,
      -1.36,	-2.43, -7.49,	-3.22,	-16.77),
    nrow=4,
    ncol=5, byrow=T))

p.df = param.vals.girls
p.df = as.data.frame(p.df)
colnames(p.df) = c("a.2", "b.2", "c.2", "d.2")
p.df

Vars <- c( "a.2", "b.2", "c.2", "d.2")

slope.1 = function(f) {
#  t.star = (t+9)/9 # adjusted time scale
  f.star = f /( (1+9)/9 - 1) 
  
  return(e=-f.star)
  
#  yij.exp = a + b*t.star + c*log(t.star) + d/t.star + e*exposure + f.star*t.star*exposure
#  yij.exp0 = a + b*t.star + c*log(t.star) + d/t.star 
#  yij.exp1 = a + b*t.star + c*log(t.star) + d/t.star + e + f.star*t.star
#  to have weight for unexposed = weight exposed at t=0 (and t.star=1)
#  need to have e = -f.star*t.star

  # would like baseline weight, at time, t, =0, to be equal around live birth; 
  # t=0 and t.star=1;
  # e = -f.star*1
  # if f = 0.15 -> f.star = 0.15/(10/9-1) = 1.35 -> e = 1.35*1 = 1.35
  # for f = -0.15 -> f.star = -0.15/(10/9-1) = -1.35 -> e = -1.35*1 = -1.35
}

effect.est = slope.1(0.15); effect.est # slope=0.15; e=-1.35

# make matrix with combination of different 
# baseline and slope values for 
# 1) higher baseline for exposure; positive slope
# 2) equal baseline for exposure; positive slope
# 3) equal baseline for exposure; negative slope
# 4) negative baseline for exposure; negative slope

base.int = matrix(c(-0.50, 0.15,
                    -1.35, 0.15,
                    1.35, -0.15,
                    0.50, -0.15), # beta4 is diff at baseline (first col) and 
                  # beta5 is diff in slopes (2nd col)
                  # set so that exposure groups have similar baseline weight values
                  nrow=4, byrow=T)

extra.vars = c("e", "f.1")
#see http://stackoverflow.com/questions/1376967/using-stat-function-and-facet-wrap-together-in-ggplot2-in-r
time <- with(sub.vals.p, seq(min(time), max(time), length = 100))
names(sub.vals.p) 

# get predicted growth curve by country (evals.f) and 
# growth scenario (params.f)

# make a separate data set with the values for predicted lines by the 
# facet_grid parameters.

cl<-makeCluster(2) # 2 workers
registerDoParallel(cl)

# want to iterate over all samples (params.f), simulation
# scenarios (evals.f) and have the appropriate model 
# parameters matching each of those combinations (in Dat and Dat.2)

# get coefficients for fitted lines for z-score by groups
coef.zwei =
  ddply(sub.vals.p, c("params", "evals", "x.1"), 
      function(x) {
        coef(glm(zwei~time, data=x))[1:2]
})

pred.lines = foreach(samps = unique(sub.vals.p$params), .combine="rbind") %do% {
  foreach(type = unique(sub.vals.p$evals), .combine="rbind") %do% {
    foreach(exposure = unique(sub.vals.p$x.1), .combine="rbind") %do% {
      
      # assign the intercept and slope for zwei regressed onto time in linear model
      a = with(coef.zwei, coef.zwei[params %in% samps & evals %in% type & x.1 %in% exposure, 4])
      b = with(coef.zwei, coef.zwei[params %in% samps & evals %in% type & x.1 %in% exposure, 5])
      
    # Assign model parameters according to sample they were drawn from
      #a.2, b.2, c.2, and d.2
    Dat <- as.vector(p.df[samps,])
    for (i in 1:length(Vars)){
      assign(Vars[i], as.numeric(Dat[i]))
      }
    
    # Assign different baseline and interaction types
    # e and f.1
    Dat.2 = as.vector(base.int[type,]) # correspond to rows different combination of baseline and interaction types
    for (i in 1:length(extra.vars)){
      assign(extra.vars[i], as.numeric(Dat.2[i]))
      }

    t.star = (time+9)/9 # adjusted time scale
    f.1.star = f.1 /( (1+9)/9 - 1) 
    
    return(
      data.frame(
               time=time,
               t.star = t.star,
               predicted = time,
               yij =      a.2 + 
                          b.2*t.star +
                          c.2*log(t.star) + 
                          d.2/t.star + 
                          e*exposure + 
                          f.1.star*(time+9)/9*exposure,
               zij = a + b*time,
               params = samps,
               evals = type,
               x.1 = exposure)
      )
    
    }
  }
}

dim(pred.lines)
head(pred.lines)

# reformat the x.1 and evals to get labels
pred.lines = within(pred.lines, {
  params.f = factor(params,
                    labels=levels(factor(sub.vals.p$params.f)))
  evals.f = factor(evals,
                   labels=levels(vals.p$evals.f))
  x.1f = factor(x.1,
                labels=levels(vals.p$x.1f))
})

# PLOT here...........................................

# test predicted lines here
ggplot(data=pred.lines, 
       aes(x=time, y=yij, colour=x.1f)) +
  geom_line(lwd=1) +
  facet_grid(evals.f ~ params.f) +
  scale_colour_manual("Group",
                      values = c("unexposed" = "blue","exposed" = "red")) +
  geom_hline(yintercept=10, lty=3, lwd=1) +
  theme_bw()


# Output to png file, females
png(file="sim-female-baselines.png", widt=2400, height=4800, res=300)

p.f +   
  geom_line(data=pred.lines, 
                  aes(x=time, y=yij, colour=x.1f),
                  lwd=2) +
  geom_hline(yintercept=10, lty=3, lwd=1.5)


dev.off()

# Z-scores ---------------------------------------------
# ---------------------------------------------------------

# Output to png file


# test predicted lines here
ggplot(data=pred.lines, 
       aes(x=time, y=zij, colour=x.1f)) +
  geom_line(lwd=1) +
  facet_grid(evals.f ~ params.f) +
  scale_colour_manual("Group",
                      values = c("unexposed" = "blue","exposed" = "red")) +
  theme_bw()

head(pred.lines)


png(file="sim-female-Z-baselines.png", res=300, height=1800, width=2400)

z.1 = ggplot(random.sub, 
       aes(x=time, y=zwei, colour=x.1f)) +
  geom_line(aes(group=id), lty=1) + 
  geom_hline(yintercept=0, lty=2) +
  ylab("Z-scores") +
  xlab("Time (months)") +
  facet_grid(evals.f ~ params.f) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_bw(base_size=26) +
  theme(legend.position="bottom",
        strip.text.y = element_text(size = 20))  +
  scale_colour_manual("Group",
                      values = c("unexposed" = "blue","exposed" = "red")) 


z.1 +
  geom_line(data=pred.lines, 
               aes(x=time, y=zij, colour=x.1f),
               lwd=2) 
dev.off()

# Combine z-scores and weight into one figure as in poster-results.Rmd
# ....................................................................

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(z.p)


weight.p = p.f +   
  geom_line(data=pred.lines, 
            aes(x=time, y=yij, colour=x.1f),
            lwd=2) +
  geom_hline(yintercept=10, lty=3, lwd=1.5)

z.p = z.1 +
  geom_line(data=pred.lines, 
            aes(x=time, y=zij, colour=x.1f),
            lwd=2) 


png(file="sim-female.png", width=4800, height=4300, res=300)

grid.arrange(weight.p + theme(legend.position="none"),
             z.p + theme(legend.position="none"),
             legend, 
             ncol=2, nrow=2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), 
             heights = c(2.5, 0.2))

dev.off()
