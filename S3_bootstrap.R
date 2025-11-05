#generate bootstrapped confidence intervals
#Daniel Suh

library(here)

source(here("base","src.R"))

#read data
m5 <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))

m5_len <- readRDS(here("processed_data", "mle", "m5_body_size_fit.rds"))


m5_coef <- coef(m5)

m5_f <- as.numeric(m5_coef[1])
m5_arr <- as.numeric(m5_coef[3])
m5_h <- as.numeric(m5_coef[5])
m5_w <- as.numeric(m5_coef[6])
m5_sd_est <- as.numeric(m5_coef[9])

m5_len_coef <- coef(m5_len)

m5_u <- as.numeric(m5_len_coef[1])/100000
m5_m <- as.numeric(m5_len_coef[2])
m5_arr_u <- as.numeric(m5_len_coef[3])
m5_rho <- as.numeric(m5_len_coef[4])
m5_phi <- as.numeric(m5_len_coef[5])


if(file.exists(here("processed_data", "supp_bootstraps.rds")) == TRUE &
   file.exists(here("processed_data", "supp_bootstrap_quantiles.rds")) == TRUE) {
  
  message("Bootstrap data already exist.")
  
} else {
  


source(here("03_infection.R"))

  
dataset %<>% filter(exposed==T)
life_data <- dataset
fora_data <- data



m5_len_sim <- function(R, time, f, u, m, length, gamma, Z, ref_t, arr_t_f, arr_t_u, h, temp, resource, w, rho, phi){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=u*exp(m*length),
              length = length,
              gamma = gamma,
              Z = Z,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              arr_t_u = arr_t_u,
              h = h,
              temp = temp,
              resource = resource,
              w = w,
              rho = rho,
              phi = phi)
  output <- as.data.frame(lsoda(y=xstart, times, m5_num_sol, params))
  
  if(Z == 0) {
    end_data <- slice_max(output, time)[,2]
  }
  else {
    end_data <- slice_max(output, time)[,4]  
  }
  
  return(end_data)
}

m5_len_ll <- function(f, u, m, arr_t_f, arr_t_u, h, w, rho, phi, sd_est){
  
    
    data <- boot_data
    
    R_end <- as.data.frame(mapply(m5_len_sim, 
                                  R=data$amt_init*fora_vol/1000, 
                                  time=data$time/60/24, 
                                  f=f/fora_vol,
                                  u = 0, 
                                  m = m, 
                                  length=data$mm, 
                                  gamma=gamma, 
                                  Z=0,
                                  ref_t=ref_t,
                                  arr_t_f=arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=h,
                                  temp = data$temp,
                                  resource = data$resource,
                                  w = w,
                                  rho = rho,
                                  phi = phi))
    colnames(R_end) <- "R_end"
    data$end <- R_end$R_end
    
    data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*fora_vol/1000))
    data %<>% mutate(sqrt_model_end = sqrt(end),
                     sqrt_data_end = sqrt(amt_rem*vol/1000))
    
    nll_sum <- 0
    for(i in 1:length(treatment_IDs)){
      treatment_data <- data %>% filter(treatment_ID == treatment_IDs[i])
      nll <- dnorm(treatment_data$sqrt_data_end,
                   mean = mean(treatment_data$sqrt_model_end),
                   sd = sd_est,
                   log = T)
      nll_sum <- -sum(nll) + nll_sum
    }
    
    
    exp_data <- boot_dataset %>% filter(exposed==TRUE)
    
    
    I_end <- as.data.frame(mapply(m5_len_sim, 
                                  R=exp_data$resource*life_vol/1000, 
                                  time=exp_data$time, 
                                  f=f/life_vol, 
                                  u=u/100000, 
                                  m = m,
                                  length = exp_data$life_mm, 
                                  gamma=gamma, 
                                  Z=spore_conc*life_vol,
                                  ref_t=ref_t,
                                  arr_t_f=arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=h,
                                  temp=exp_data$temp,
                                  resource=exp_data$resource,
                                  w = w,
                                  rho = rho,
                                  phi = phi))
    
    colnames(I_end) <- "I_end"
    
    exp_data %<>% cbind(I_end)
    
    
    exp_data %<>% 
      mutate(ll = mapply(inf_out, inf_status = inf_status, I_end=I_end))
    nll_sum <- nll_sum + -sum(exp_data$ll)
    return(nll_sum)  
    
    
  
}


# Run Bootstrap

iterations <- 1

boot_01 <- list()

start_time <- Sys.time()

for(i in 1:iterations){
  iter_time <- Sys.time()
  
  boot_data <- fora_data %>% 
    group_by(treatment_ID) %>%
    slice_sample(., prop=1, replace=T) %>%
    ungroup()
  
  boot_dataset <- life_data %>% 
    group_by(treatment) %>%
    slice_sample(., prop=1, replace=T) %>%
    ungroup()
  
  boot_01[[i]] <-mle2(m5_len_ll, start=list(u = m5_u*100000,
                                        f=m5_f, 
                                        m = m5_m,
                                        arr_t_f=m5_arr, 
                                        arr_t_u=m5_arr_u,
                                        h=m5_h,
                                        w=m5_w,
                                        rho=m5_rho,
                                        phi=m5_phi,
                                        sd_est=m5_sd_est), 
                      control=list(parscale = c(u = m5_u*100000,
                                                f=m5_f, 
                                                m = m5_m,
                                                arr_t_f=m5_arr, 
                                                arr_t_u=m5_arr_u,
                                                h=m5_h,
                                                w=m5_w,
                                                rho=m5_rho,
                                                phi=m5_phi,
                                                sd_est=m5_sd_est),
                                   maxit=5000,
                                   reltol=0.00001),
                      skip.hessian=F, 
                      method="Nelder-Mead")
  print(i)
  end_time <- Sys.time()
  print(end_time)
  print(end_time-iter_time)
  print(end_time-start_time)
}

#saveRDS(boot_01, file = here("processed_data", "supp_bootstraps.rds"))

}



#format bootstrap data


mod_boot <- readRDS(here("processed_data", "supp_bootstraps.rds"))

tib_length <- length(mod_boot)

mod_boot_coefs <- tibble(f = rep_len(0, tib_length),
                         u = rep_len(0, tib_length),
                         m = rep_len(0, tib_length),
                         arr = rep_len(0, tib_length),
                         arr_u = rep_len(0, tib_length),
                         h = rep_len(0, tib_length),
                         w = rep_len(0, tib_length),
                         rho = rep_len(0, tib_length),
                         phi = rep_len(0, tib_length),
                         sd_est = rep_len(0, tib_length),
                         warning = NA)


for(i in 1:length(mod_boot)){
  if(mod_boot[[i]]@details[4] == 52){
    print(paste(i, "convergence_failure", sep="_"))
    mod_boot_coefs$f[i] <- NA
    mod_boot_coefs$u[i] <- NA
    mod_boot_coefs$m[i] <- NA
    mod_boot_coefs$arr[i] <- NA
    mod_boot_coefs$arr_u[i] <- NA
    mod_boot_coefs$h[i] <- NA
    mod_boot_coefs$w[i] <- NA
    mod_boot_coefs$rho[i] <- NA
    mod_boot_coefs$phi[i] <- NA
    mod_boot_coefs$sd_est[i] <- NA
    mod_boot_coefs$warning[i] <- "52"
  }
  else if(mod_boot[[i]]@details[4] == 10){
    print(paste(i, "convergence_failure", sep="_"))
    mod_boot_coefs$f[i] <- NA
    mod_boot_coefs$u[i] <- NA
    mod_boot_coefs$m[i] <- NA
    mod_boot_coefs$arr[i] <- NA
    mod_boot_coefs$arr_u[i] <- NA
    mod_boot_coefs$h[i] <- NA
    mod_boot_coefs$w[i] <- NA
    mod_boot_coefs$rho[i] <- NA
    mod_boot_coefs$phi[i] <- NA
    mod_boot_coefs$sd_est[i] <- NA
    mod_boot_coefs$warning[i] <- "10"
  }
  else if(mod_boot[[i]]@min == 0){
    print(paste(i, "zero_ll", sep="_"))
    mod_boot_coefs$f[i] <- NA
    mod_boot_coefs$u[i] <- NA
    mod_boot_coefs$m[i] <- NA
    mod_boot_coefs$arr[i] <- NA
    mod_boot_coefs$arr_u[i] <- NA
    mod_boot_coefs$h[i] <- NA
    mod_boot_coefs$w[i] <- NA
    mod_boot_coefs$rho[i] <- NA
    mod_boot_coefs$phi[i] <- NA
    mod_boot_coefs$sd_est[i] <- NA
    mod_boot_coefs$warning[i] <- "zero_ll"
  }
  else if(mod_boot[[i]]@details[4] == 1){
    print(paste(i, "max_iter_reached", sep="_"))
    boot_coef <- coef(mod_boot[[i]])
    mod_boot_coefs$f[i] <- as.numeric(boot_coef[1])
    mod_boot_coefs$u[i] <- as.numeric(boot_coef[2])
    mod_boot_coefs$m[i] <- as.numeric(boot_coef[3])
    mod_boot_coefs$arr[i] <- as.numeric(boot_coef[4])
    mod_boot_coefs$arr_u[i] <- as.numeric(boot_coef[5])
    mod_boot_coefs$h[i] <- as.numeric(boot_coef[6])
    mod_boot_coefs$w[i] <- as.numeric(boot_coef[7])
    mod_boot_coefs$rho[i] <- as.numeric(boot_coef[8])
    mod_boot_coefs$phi[i] <- as.numeric(boot_coef[9])
    mod_boot_coefs$sd_est[i] <- as.numeric(boot_coef[10])
    mod_boot_coefs$warning[i] <- "1"
  }
  else{
    boot_coef <- coef(mod_boot[[i]])
    mod_boot_coefs$f[i] <- as.numeric(boot_coef[1])
    mod_boot_coefs$u[i] <- as.numeric(boot_coef[2])
    mod_boot_coefs$m[i] <- as.numeric(boot_coef[3])
    mod_boot_coefs$arr[i] <- as.numeric(boot_coef[4])
    mod_boot_coefs$arr_u[i] <- as.numeric(boot_coef[5])
    mod_boot_coefs$h[i] <- as.numeric(boot_coef[6])
    mod_boot_coefs$w[i] <- as.numeric(boot_coef[7])
    mod_boot_coefs$rho[i] <- as.numeric(boot_coef[8])
    mod_boot_coefs$phi[i] <- as.numeric(boot_coef[9])
    mod_boot_coefs$sd_est[i] <- as.numeric(boot_coef[10])
    mod_boot_coefs$warning[i] <- "0"
  }
}

mod_boot_coefs %<>% filter(warning == 0)

mod_quantiles <- tibble(id = c("est", "mean", "lower.025", "upper.975"),
                        f = NA,
                        u = NA,
                        m = NA,
                        arr = NA,
                        arr_u = NA,
                        h = NA,
                        w = NA,
                        rho = NA,
                        phi = NA,
                        sd_est = NA)


mod_quantiles$f[1] <- m5_f
mod_quantiles$f[2] <- mean(mod_boot_coefs$f)
mod_quantiles$f[3] <- quantile(mod_boot_coefs$f, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$f[4] <- quantile(mod_boot_coefs$f, probs=seq(0.025, 0.975, 0.95))[2]



mod_quantiles$u[1] <- m5_u
mod_quantiles$u[2] <- mean(mod_boot_coefs$u)/100000
mod_quantiles$u[3] <- quantile(mod_boot_coefs$u, probs=seq(0.025, 0.975, 0.95))[1]/100000
mod_quantiles$u[4] <- quantile(mod_boot_coefs$u, probs=seq(0.025, 0.975, 0.95))[2]/100000


mod_quantiles$m[1] <- m5_m
mod_quantiles$m[2] <- mean(mod_boot_coefs$m)
mod_quantiles$m[3] <- quantile(mod_boot_coefs$m, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$m[4] <- quantile(mod_boot_coefs$m, probs=seq(0.025, 0.975, 0.95))[2]



mod_quantiles$arr[1] <- m5_arr
mod_quantiles$arr[2] <- mean(mod_boot_coefs$arr)
mod_quantiles$arr[3] <- quantile(mod_boot_coefs$arr, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$arr[4] <- quantile(mod_boot_coefs$arr, probs=seq(0.025, 0.975, 0.95))[2]


mod_quantiles$arr_u[1] <- m5_arr_u
mod_quantiles$arr_u[2] <- mean(mod_boot_coefs$arr_u)
mod_quantiles$arr_u[3] <- quantile(mod_boot_coefs$arr_u, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$arr_u[4] <- quantile(mod_boot_coefs$arr_u, probs=seq(0.025, 0.975, 0.95))[2]


mod_quantiles$h[1] <- m5_h
mod_quantiles$h[2] <- mean(mod_boot_coefs$h)
mod_quantiles$h[3] <- quantile(mod_boot_coefs$h, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$h[4] <- quantile(mod_boot_coefs$h, probs=seq(0.025, 0.975, 0.95))[2]


mod_quantiles$w[1] <- m5_w
mod_quantiles$w[2] <- mean(mod_boot_coefs$w)
mod_quantiles$w[3] <- quantile(mod_boot_coefs$w, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$w[4] <- quantile(mod_boot_coefs$w, probs=seq(0.025, 0.975, 0.95))[2]



mod_quantiles$rho[1] <- m5_rho
mod_quantiles$rho[2] <- mean(mod_boot_coefs$rho)
mod_quantiles$rho[3] <- quantile(mod_boot_coefs$rho, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$rho[4] <- quantile(mod_boot_coefs$rho, probs=seq(0.025, 0.975, 0.95))[2]


mod_quantiles$phi[1] <- m5_phi
mod_quantiles$phi[2] <- mean(mod_boot_coefs$phi)
mod_quantiles$phi[3] <- quantile(mod_boot_coefs$phi, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$phi[4] <- quantile(mod_boot_coefs$phi, probs=seq(0.025, 0.975, 0.95))[2]


mod_quantiles$sd_est[1] <- m5_sd_est
mod_quantiles$sd_est[2] <- mean(mod_boot_coefs$sd_est)
mod_quantiles$sd_est[3] <- quantile(mod_boot_coefs$sd_est, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$sd_est[4] <- quantile(mod_boot_coefs$sd_est, probs=seq(0.025, 0.975, 0.95))[2]


theme_set(theme_bw(base_size = 8))

mod_quantiles %<>% pivot_longer(cols = f:sd_est) %>% pivot_wider(., names_from="id")

#saveRDS(mod_quantiles, file = here("processed_data", "supp_bootstrap_quantiles.rds"))


mod_quantiles_2F <- readRDS(file = here("processed_data", "supp_bootstrap_quantiles.rds"))

mod_quantiles_2E <- readRDS(file = here("processed_data", "m2E_bootstrap_quantiles.rds"))

mod_quantiles_combine <- 
  mod_quantiles %>% mutate(model = "2F") %>% 
  rbind(., mod_quantiles_2E %>% mutate(model = "2E"))

#Parameter estimates and CI for supplemental model



supp_ci_est_combine <- 
  mod_quantiles_combine %>% 
  ggplot(aes(x = model, 
             y = est,
             color = model)) + 
  geom_point() + 
  geom_linerange(aes(ymin = lower.025, 
                     ymax = upper.975)) + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = c("#E66100", "#4B0092")) + 
  labs(x = "",
       y = "",
       title = "Maximum Likelihood Estimates (95% CI)",
       color = "Model") + 
  facet_wrap(~name, scales = "free",
             nrow = 2) +
  theme_bw(base_size = 8) 



if(dir.exists(here("figures")) == FALSE) {
  message("Welcome! Let's make some room for figures.")
  dir.create(here("figures")) 
} else {
  message("/figures exists! Proceeeding to save.")
}




ggsave(here("figures", "supp_ci_est_combine.png"), 
       plot = supp_ci_est_combine,
       width = 7,
       height = 4,
       units = "in")

