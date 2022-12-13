#Read in data
ped_dt <- read.csv("c:/GUMC/Research/software/Fan/FLM_longitudinal_discrete_popu/Example_data/ped.csv", header = TRUE)
cov_dt <- read.csv("c:/GUMC/Research/software/Fan/FLM_longitudinal_discrete_popu/Example_data/covariate.csv", header = TRUE)
gen_dt <- read.csv("c:/GUMC/Research/software/Fan/FLM_longitudinal_discrete_popu/Example_data/geno.csv", header = TRUE)
pos_dt <- read.csv("c:/GUMC/Research/software/Fan/FLM_longitudinal_discrete_popu/Example_data/pos.csv", header = TRUE)
pos_dt2 <- c(pos_dt$x)
#Read in functions
source("c:/GUMC/Research/software/Fan/FLM_longitudinal_discrete_popu/FLM_longitudinal_discrete_popu_beta_smooth_only.R")
source("c:/GUMC/Research/software/Fan/FLM_longitudinal_discrete_popu/FLM_longitudinal_discrete_popu_fixed_model.R")

pval <- list()

Beta_only_B <- FLM_longi_disc_popu_beta_smooth(ped = ped_dt, cov = cov_dt, geno = gen_dt, 
                  pos = pos_dt2, mu_nknots = 10, correlation = "ou", order = 4, beta_basis = 10, base = "bspline")

Beta_only_F <- FLM_longi_disc_popu_beta_smooth(ped = ped_dt, cov = cov_dt, geno = gen_dt, 
                  pos = pos_dt2, mu_nknots = 10, correlation = "ou", order = 4, beta_basis = 10, base = "fspline")

Fixed_B <- FLM_longi_disc_popu_fixed(ped = ped_dt, cov = cov_dt, geno = gen_dt, 
              pos = pos_dt2, mu_nknots = 10, correlation = "ou", order = 4, beta_basis = 10, geno_basis = 10, base = "bspline")

Fixed_F <- FLM_longi_disc_popu_fixed(ped = ped_dt, cov = cov_dt, geno = gen_dt, 
              pos = pos_dt2, mu_nknots = 10, correlation = "ou", order = 4, beta_basis = 10, geno_basis = 10, base = "fspline")

pval$Beta_only_b_mu_fixed <- Beta_only_B$fit_cov_geno
pval$Beta_only_f_mu_fixed <- Beta_only_F$fit_cov_geno
pval$Fixed_b_mu_fixed <- Fixed_B$fit_cov_geno
pval$Fixed_f_mu_fixed <- Fixed_F$fit_cov_geno

pval
