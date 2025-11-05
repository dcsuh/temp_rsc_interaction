#Model fitting for supplemental model - size-dependent interactive u (model 2F)


library(here)

source(here("base","src.R"))
source(here("03_infection.R"))


if(file.exists(here("processed_data", "mle", "m5_body_size_fit.rds")) == TRUE) {
  
  message("file already exists")
  
} else {

m5_num_sol <- function(t, x, params){
  R <- x[1]
  S <- x[2]
  I <- x[3]
  Z <- x[4]
  
  with(as.list(params),{
    dR <- -R*(S+I)*(f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dS <- -Z*S*u*exp(m*length)*exp(arr_t_u*(1/ref_t - 1/temp))*exp(resource*rho)*exp(resource*temp*phi)*
      (f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dI <- Z*S*u*exp(m*length)*exp(arr_t_u*(1/ref_t - 1/temp))*exp(resource*rho)*exp(resource*temp*phi)*
      (f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dZ <- -Z*(S+I)*(f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    res <- c(dR, dS, dI, dZ)
    list(res)}
  )
}


m5_len_sim <- function(R, time, f, u, m, length, gamma, Z, ref_t, arr_t_f, arr_t_u, h, temp, resource, w, rho, phi){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=u,
              m = m,
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



m5_len_ll <- function(u, m, arr_t_u, h, rho, phi){
  
  if(#f<0 | f > f_max | h<0 | u<0 | arr_t_f<0 | arr_t_u<0 | rho < 0 | phi > 0 | w > 0 | 
     (u/100000*exp(m*min(length_summ$life_mm))*exp(arr_t_u*((1/15)-(1/15)))*exp(0.1/1000*rho)*exp(0.1/1000*15*phi)) <0 |
     (u/100000*exp(m*max(length_summ$life_mm))*exp(arr_t_u*((1/15)-(1/25)))*exp(1.0/1000*rho)*exp(1.0/1000*25*phi)) <0) {
    
    return(Inf)
    
  } 
  
  else {  
    cat((u/100000*exp(m*min(length_summ$life_mm))*exp(arr_t_u*((1/15)-(1/15)))*exp(0.1/1000*rho)*exp(0.1/1000*15*phi)),
        (u/100000*exp(m*max(length_summ$life_mm))*exp(arr_t_u*((1/15)-(1/25)))*exp(1.0/1000*rho)*exp(1.0/1000*25*phi)),
        "\n")
    
    
    R_end <- as.data.frame(mapply(m5_len_sim, 
                                  R=data$amt_init*fora_vol/1000, 
                                  time=data$time/60/24, 
                                  f=m5_f/fora_vol,
                                  u = 0, 
                                  m = m, 
                                  length=data$mm, 
                                  gamma=gamma, 
                                  Z=0,
                                  ref_t=ref_t,
                                  arr_t_f=m5_arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=m5_h,
                                  temp = data$temp,
                                  resource = data$resource,
                                  w = m5_w,
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
                   sd = m5_sd,
                   log = T)
      nll_sum <- -sum(nll) + nll_sum
    }
    
    
    exp_data <- dataset %>% filter(exposed==TRUE)
    
    
    I_end <- as.data.frame(mapply(m5_len_sim, 
                                  R=exp_data$resource*life_vol/1000, 
                                  time=exp_data$time, 
                                  f=m5_f/life_vol, 
                                  u=u/100000, 
                                  m = m,
                                  length = exp_data$life_mm, 
                                  gamma=gamma, 
                                  Z=spore_conc*life_vol,
                                  ref_t=ref_t,
                                  arr_t_f=m5_arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=m5_h,
                                  temp=exp_data$temp,
                                  resource=exp_data$resource,
                                  w = m5_w,
                                  rho = rho,
                                  phi = phi))
    
    colnames(I_end) <- "I_end"
    
    exp_data %<>% cbind(I_end)
    
    
    exp_data %<>% 
      mutate(ll = mapply(inf_out, inf_status = inf_status, I_end=I_end))
    nll_sum <- nll_sum + -sum(exp_data$ll)
    return(nll_sum)  
    
    
  }
}


#use fits from best infection model (2E)
m5 <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))

m5_coef <- coef(m5)
m5_f <- as.numeric(m5_coef[1])
m5_u <- as.numeric(m5_coef[2])
m5_arr_t_f <- as.numeric(m5_coef[3])
m5_arr_t_u <- as.numeric(m5_coef[4])
m5_h <- as.numeric(m5_coef[5])
m5_w <- as.numeric(m5_coef[6])
m5_rho <- as.numeric(m5_coef[7])
m5_phi <- as.numeric(m5_coef[8])
m5_sd  <- as.numeric(m5_coef[9])



m5_len_fit <-
  mle2(m5_len_ll, start=list(u = 10, 
                             m = 0.1,
                             arr_t_u=30,
                             rho=0.1,
                             phi=-0.1),
       control=list(parscale = c(u=10, 
                                 m = 0.1,
                                 arr_t_u=10,
                                 rho=0.1,
                                 phi=0.1),
                    maxit=5000),
       skip.hessian=F, 
       method="Nelder-Mead")

saveRDS(m5_len_fit, file = here("processed_data", "mle", "m5_body_size_fit.rds"))
}

