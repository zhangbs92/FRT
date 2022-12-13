### This program generates datasets and implements codes to run the models of the paper ###
### Fan RZ, Zhang YW, Albert PS, Liu AY, Wang YJ, and Xiong MM (2012) Longitudinal ###
### association analysis of quantitative traits. Genetic Epidemiology 36:856-869. ###
### The codes does needs library(epicalc) ###
### If you use the codes, please cite the above paper. Thanks. ###
###########################################################################################
### Plus, This program generates dataset Age_Sex_Geno_Trait_data.csv ###
########################################################################
library(MASS)
library(nlme)
### Functions
t <- seq(15, 70);
mut <- function(t)
{return(-34.2+81.7*log(0.3*(t+21.7)) ) }
func1 <- mut(t)
plot(t, func1, type="l")
mut1 <-function(t)
{return(29.5+5.2*log(0.25*(t+14.2)) ) }
plot(t, mut1(t), type="l")
sint <- function(t)
{return(100 + 50 * sin(0.2*t) )}
plot(t, sint(t), type="l")
expt <- function(t)
{return(110 * exp(0.0002* (t - 25)^2 ))}
plot(t, expt(t), type="l")
sigt <- function(t)
{return(5+0.5*t)}
ar1 <- function(theta,t,var)
     {
     n <- length(t)
     cov <- matrix(0,n,n)
     for (i in 1:n)
          for (j in 1:i)
     {
               cov[i,j]<-theta^abs(t[i]-t[j])*sqrt(var[i]*var[j])
               cov[j,i]<-cov[i,j]
                         }
               return(cov)
     }
arma <- function(p,q,t,var)
               {
               n <- length(t)
     cov <- matrix(0,n,n)
               for (i in 1:n)
                    for (j in 1:i)
{
               cov[i,j] <- p*q^abs(t[i]-t[j])*sqrt(var[i]*var[j])
               cov[j,i] <- cov[i,j]
                         }
               return(cov)
     }
############# Design parameters #####################
## number of individuals
n <- 200
famid <- c(1:n)
## age
##nobs <- sample(c(4:8), n, replace=T) ### Results for Table 4 ###
nobs <- sample(c(3:6), n, replace=T)
age <- vector( "numeric", sum(nobs) )
tmpid <- 1
a1 <- vector( "numeric", n )
visit <- matrix(0, sum(nobs), 1)
ID <- matrix(0, sum(nobs), 1)
for (i in 1:n)
     {
     a1[i] <- runif(1, 20, 65)
     age[tmpid] <- a1[i]
     visit[tmpid] <- 1
     ID[tmpid] <- i
     for (j in 1:(nobs[i]-1))
          {
          tmpid <- tmpid + 1
          r <- rbinom(1,1,0.5)
          age[tmpid] <- age[tmpid-1]+2*r+(1-r)*4
          visit[tmpid] <- j + 1
          ID[tmpid] <- i
          }
     tmpid <- tmpid + 1
     }
all <- cbind(ID, visit, age)
colnames(all) <- c("ID", "visit", "age")
all <- data.frame(all)
attach(all)
###simulate genotype
p <- 0.25
geno <- vector("numeric", n)
for (i in 1:n)
     {
     seed <- runif(1,0,1)
     if ( seed < p * p)
          geno[i] <- 0
          else
               {
               if ( seed < p * p + 2 * p * (1-p) )
                    geno[i] <- 1
                              else
                         eno[i] <- 2
               }
     }
genodata <- cbind(famid, geno)
colnames(genodata) <- c("famid", "geno")
genodata <- data.frame(genodata)
AgeGenodata <- merge(all, genodata, by.x="ID", by.y="famid")
#sort by famid, subid, visit #
nt <- aggregate(AgeGenodata[,3], by=list(ID), func=c("count"), FUN = mean)
sex <- vector( "numeric", n)
sex <- rbinom(n, 1, 0.5)
nt <- cbind(nt, sex)
AgeSexGenodata <- merge(AgeGenodata, nt, by.x=c("ID"), by.y=c("Group.1"))
######### simulate phenotype #######
sig_E <- 5
vare <- 10 ## variance of the random residuals
nrep <- 20#00 ### You may change it to 2000
beta_sex <- 5
alpha <- 0 ##genetic effect
ymat <- matrix(0, nrow(AgeSexGenodata), nrep)
for (h in 1:nrep)
     {
     E <- mvrnorm(n, rep(0,1), sig_E^2)
     tempid <- 1
     y <- vector("numeric", nrow(AgeSexGenodata))
     theta <- 0
     for (i in 1:n)
          {
          Ei <- E[i] #shared environmental
          temps <- sex[i]
          nobstemp <- nobs[i]
               time <- AgeSexGenodata[tempid:(tempid+nobstemp-1),3]
               #var<-sigt(time)^2
               #cov<-ar1(theta,time,var))
               #diag(cov)<-var
               #epsilon<-mvrnorm(1,rep(0,nobstemp),cov)
               epsilon <- mvrnorm(1, rep(0,nobstemp), diag(vare, nobstemp))
               genet <- alpha*AgeSexGenodata[tempid, 4]
          sext <- beta_sex*AgeSexGenodata[tempid, 6]
          ### the function mut is used to generate the traits below ###
               y[tempid:(tempid+nobstemp-1)] <- mut(time) + rep(Ei, nobstemp) + rep(genet,nobstemp) + rep(sext,nobstemp) + epsilon
          tempid <- tempid + nobstemp
          }
     ymat[,h] <- y
     print(h)
     }
###save the dataset of AgeSexGenodata ###
Age_Sex_Geno_Trait_data = cbind( AgeSexGenodata[, c(1:4, 6)], y)
Data.name <- sprintf("C:/NICHD/Research/software/Fan/Longitudinal_qtl_popu_web/Age_Sex_Geno_Trait_data.csv")
#Data.name <- sprintf("/data/fanr/NICHD/Research/paper/2012/Longitudinal_qtl_popu/Simulation/Age_Sex_Geno_Trait_data.csv")
write.csv(Age_Sex_Geno_Trait_data, file = Data.name, row.names = FALSE)
###
######################## Fit mixed effects model:
## set up basis functions #linear truncated polynomial base
x <- AgeSexGenodata[,3] - mean(AgeSexGenodata[,3]) #centered age
id <- AgeSexGenodata[,1]
nknots <- 15
temp <- seq(min(x), max(x), length=nknots+2) #check this
knots <- temp[2:(length(temp)-1)]
m <- length(x)
X1 <- cbind(rep(1,m), AgeSexGenodata[,6], x )
X2 <- cbind(rep(1,m), AgeSexGenodata[,6], x, x^2)
X3 <- cbind(rep(1,m), AgeSexGenodata[,6], x, x^2, x^3)
XT <- cbind(rep(1,m), AgeSexGenodata[,6], mut( AgeSexGenodata[,3] ) )
Z <- outer(x, knots, "-")
Z <- Z*(Z>0)
colnames(X1) <- c("intercept", "sex", "age_m")
colnames(X2) <- c("intercept", "sex", "age_m", "age_m^2")
colnames(X3) <- c("intercept", "sex", "age_m", "age_m^2", "age_m^3")
colnames(XT) <- c("intercept", "sex", "mut")
group <- rep(1, nrow(AgeSexGenodata))
genotype <- AgeSexGenodata[, 4]
###################################################################
################ start from here to fit model and compute power ###
##compute power for testing genetic effect
power1 <- a_genet1 <- b_sex1 <- vector("numeric", nrep)
power2 <- a_genet2 <- b_sex2 <- vector("numeric", nrep)
power3 <- a_genet3 <- b_sex3 <- vector("numeric", nrep)
power4 <- a_genet4 <- b_sex4 <- vector("numeric", nrep)
powerT <- a_genetT <- b_sexT <- vector("numeric", nrep)
for (h in 1:nrep)
     {
     y <- ymat[,h]
     fit11 <- lme(y ~ -1 + X1 + genotype,
                    random =list(group = pdIdent(~-1+Z),
                                   id = ~1),
                    correlation = corExp(form = ~ x | group/id),
                    method="ML")
     fit12 <- lme(y ~ -1 + X1,
                    random =list(group = pdIdent(~-1+Z),
                                   id = ~1),
                    correlation = corExp(form = ~ x | group/id),
                    method="ML")
     power1[h] <- anova(fit11,fit12)[2,9]
     b_sex1[h] <- summary(fit11)$tTable[3,1]
     a_genet1[h] <- summary(fit11)$tTable[4,1]
     #####
     fit21 <- lme(y ~ -1 + X1 + genotype,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     fit22 <- lme(y ~ -1 + X1,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     power2[h] <- anova(fit21,fit22)[2,9]
     b_sex2[h] <- summary(fit21)$tTable[3,1]
     a_genet2[h] <- summary(fit21)$tTable[4,1]
     #####
     fit31 <- lme(y ~ -1 + X2 + genotype,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     fit32 <- lme(y ~ -1 + X2,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     power3[h] <- anova(fit31,fit32)[2,9]
     b_sex3[h] <- summary(fit31)$tTable[3,1]
     a_genet3[h] <- summary(fit31)$tTable[4,1]
     #####
     fit41 <- lme(y ~ -1 + X3 + genotype,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     fit42 <- lme(y ~ -1 + X3,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     power4[h] <- anova(fit41,fit42)[2,9]
     b_sex4[h] <- summary(fit41)$tTable[3,1]
     a_genet4[h] <- summary(fit41)$tTable[4,1]
     #####
     fitT1 <- lme(y ~ -1 + XT + genotype,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     fitT2 <- lme(y ~ -1 + XT,
                    random =list(id = ~1),
                    correlation = corExp(form = ~ x | id),
                    method="ML")
     powerT[h] <- anova(fitT1,fitT2)[2,9]
     b_sexT[h] <- summary(fitT1)$tTable[3,1]
     a_genetT[h] <- summary(fitT1)$tTable[4,1]
     print(h)
     }
sum(power1 < 0.05)/nrep
sum(power1 < 0.01)/nrep
mean(b_sex1)
mean(a_genet1)
alpha - mean(a_genet1)
sum(power2 < 0.05)/nrep
sum(power2 < 0.01)/nrep
mean(b_sex2)
mean(a_genet2)
alpha - mean(a_genet2)
sum(power3 < 0.05)/nrep
sum(power3 < 0.01)/nrep
mean(b_sex3)
mean(a_genet3)
alpha - mean(a_genet3)
sum(power4 < 0.05)/nrep
sum(power4 < 0.01)/nrep
mean(b_sex4)
mean(a_genet4)
alpha - mean(a_genet4)
sum(powerT < 0.05)/nrep
sum(powerT < 0.01)/nrep
mean(b_sexT)
mean(a_genetT)
alpha - mean(a_genetT)
