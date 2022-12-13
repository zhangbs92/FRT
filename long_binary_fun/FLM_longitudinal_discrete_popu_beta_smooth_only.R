library(fda)
library(MASS)
library(Matrix)
library(nlme)
library(lattice)
library(car)
library(glmmTMB)

### This Function is written by Ruzong Fan, March, 2018 ###
FLM_longi_disc_popu_beta_smooth <- function(ped, cov, geno, pos, mu_nknots, correlation, order=4, beta_basis, base) 
   {
   #ped<-ped_sim;cov<-cov_sim;geno<-gen_sim;pos<-pos_sim;mu_nknots<-5;correlation<-correlation;order<-order;beta_basis<-5;base<-"bspline"
   #General checking
   if (ncol(geno)-2 != length(pos)) 
      {
      return(print("Length of pos is not equal to the number of geno"))
      }
  
   ### Read covarite data ###
   cov <- cov[, -c(1:2)]
   cov[is.na(cov)] <- 0 
   
   #### CREATE BASIS   
   if (base ==  "bspline")
      {
      betabasis <- create.bspline.basis(norder = order, nbasis = beta_basis)
      } else if (base == "fspline"){
         betabasis <- create.fourier.basis(c(0,1), nbasis = beta_basis)
         }

   ###  Read genotype data ###
   geno <- geno[, -c(1:2)]  
   
   #### CASE WTIH NO MISSING GENO
   if (sum(is.na(geno)) == 0)
      {
      dqr     <- qr(geno)
      qr_id   <- dqr$pivot[1:dqr$rank]
      geno    <- as.matrix(geno[, qr_id])
      pos     <- pos[qr_id]

      ### Normalize the region to [0,1] if needed  
      if (max(pos) > 1)
         {
         pos <- (pos - min(pos)) / (max(pos) - min(pos))
	       }
  
      B  <- eval.basis(pos, betabasis)	
      UJ <- as.matrix(geno) %*% B
      } else if (sum(is.na(geno)) > 0)
         {
         if (max(pos) > 1)
            {
            pos <- (pos - min(pos)) / (max(pos) - min(pos))
            }
        
            B  <- eval.basis(pos, betabasis)	
            UJ <- NULL
            for (i in 1:nrow(geno))
               {
               idx <- which(is.na(geno[i,])) 
               if (length(idx) == 0 || length(idx) == ncol(geno))
                  {
                  gi <- geno[i,]
                  Bi <- B
                  } else 
                     {
                     gi <- geno[i, -idx] 
                     Bi <- B[-idx,]
                     }
          
               gi_m <- matrix(gi, nrow = 1)
               tmp <-  unlist(gi_m) %*% Bi
               UJ  <- rbind(UJ, tmp)
               }
         }

   #### SET TRAIT TO NA IF WHOLE GENO IS MISSING, SO THAT "fitnull" BELOW WILL USE DATA WITHOUT SUBJECT WITH WHOLE GENO MISSING
   if (sum(is.na(UJ)) > 0)
      {
      tmp1 <- is.na(UJ)
      tmp2 <- apply(tmp1, 1, sum)
      ped[tmp2>0, "event"] <- NA
      }

   idx       <- is.na(ped[, "event"])
   trait <- ped[!idx, "event"]
   cov   <- data.frame(cov)[!idx,]
   ped   <- ped[!idx, ]
   UJ       <- UJ[!idx,]
   id    <- ped[, "ped"]
   time  <- ped[, "time"]
   n_t      <- length(time)
   group    <- rep(1,n_t)
   time2   <- numFactor(time)
      
   ### Make sure UJ has full rank of bbasis or fbasis ###
   UJdqr   <- qr(UJ)
   UJindex <- UJdqr$pivot[1:UJdqr$rank]
   UJ      <- UJ[, UJindex]
         
   ##Defining matrix terms for the mean ###
   t_knot    <- seq(min(time, na.rm = T), max(time, na.rm = T),length = mu_nknots + 2)  #check this
   knot_mean <- t_knot[-c(1, length(t_knot))]
   z1        <- outer(time, knot_mean, "-")
   z2        <- z1*(z1>0)
   z_mean    <- cbind(time, z2)
   
   ###
   pval <- list()
   
   if (correlation == "none")
      {
      fit_cov        <- glmmTMB(trait ~ as.matrix(cov)      + (time + 0 | id), family = binomial)
      fit_cov_geno   <- glmmTMB(trait ~ as.matrix(cov) + UJ + (time + 0 | id), family = binomial)
      
      #treat mu(t) as fixed
      fit_cov_fix        <- glmmTMB(trait ~ as.matrix(cov) + z_mean      + (time + 0 | id), family = binomial)
      fit_cov_fix_geno   <- glmmTMB(trait ~ as.matrix(cov) + z_mean + UJ + (time + 0 | id), family = binomial)
      
      #treat mu(t) as random
      fit_cov_mix        <- glmmTMB(trait ~ as.matrix(cov) + (z_mean + 0 | group)      + (time + 0 | id), family = binomial)
      fit_cov_mix_geno   <- glmmTMB(trait ~ as.matrix(cov) + (z_mean + 0 | group) + UJ + (time + 0 | id), family = binomial)
      } else if (correlation == "ou")
         {
         fit_cov        <- glmmTMB(trait ~ as.matrix(cov)      + ou(time2 + 0 | id), family = binomial)
         fit_cov_geno   <- glmmTMB(trait ~ as.matrix(cov) + UJ + ou(time2 + 0 | id), family = binomial)
      
         #treat mu(t) as fixed
         fit_cov_fix        <- glmmTMB(trait ~ as.matrix(cov) + z_mean      + ou(time2 + 0 | id), family = binomial)
         fit_cov_fix_geno   <- glmmTMB(trait ~ as.matrix(cov) + z_mean + UJ + ou(time2 + 0 | id), family = binomial)
      
         #treat mu(t) as random
         fit_cov_mix        <- glmmTMB(trait ~ as.matrix(cov) + (z_mean + 0 | group)      + ou(time2 + 0 | id), family = binomial)
         fit_cov_mix_geno   <- glmmTMB(trait ~ as.matrix(cov) + (z_mean + 0 | group) + UJ + ou(time2 + 0 | id), family = binomial)
         } else if (correlation == "gauss")
            {
            fit_cov        <- glmmTMB(trait ~ as.matrix(cov)      + gau(time2 + 0 | id), family = binomial)
            fit_cov_geno   <- glmmTMB(trait ~ as.matrix(cov) + UJ + gau(time2 + 0 | id), family = binomial)
      
            #treat mu(t) as fixed 
            fit_cov_fix        <- glmmTMB(trait ~ as.matrix(cov) + z_mean      + gau(time2 + 0 | id), family = binomial)
            fit_cov_fix_geno   <- glmmTMB(trait ~ as.matrix(cov) + z_mean + UJ + gau(time2 + 0 | id), family = binomial)
      
            #treat mu(t) as random
            fit_cov_mix        <- glmmTMB(trait ~ as.matrix(cov) + (z_mean + 0 | group)      + gau(time2 + 0 | id), family = binomial)
            fit_cov_mix_geno   <- glmmTMB(trait ~ as.matrix(cov) + (z_mean + 0 | group) + UJ + gau(time2 + 0 | id), family = binomial)
            }
   
   ### 
   tryCatch(
      error = function(cnd)
         {
         pval$fit_cov_geno <- NA
         },
         pval$fit_cov_geno <- anova(fit_cov, fit_cov_geno)[2,8])
   
   tryCatch(
      error = function(cnd)
         {
         pval$fit_cov_mu_fixed_geno <- NA
         },
         pval$fit_cov_mu_fixed_geno <- anova(fit_cov_fix, fit_cov_fix_geno)[2,8])
   
   tryCatch(
      error = function(cnd)
         {
         pval$fit_cov_mu_mixed_geno <- NA
         },
         pval$fit_cov_mu_mixed_geno <- anova(fit_cov_mix, fit_cov_mix_geno)[2,8])
   
   pval     	
   }

