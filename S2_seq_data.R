#simulate results across length for size-dependent model

library(here)

source(here("base","src.R"))
source(here("03_infection.R"))


m5 <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))

m5_len <- readRDS(here("processed_data", "mle", "m5_body_size_fit.rds"))

if(file.exists(here("processed_data", "seq_data", "supp_body_size_fit_data.rds")) == TRUE) {
  
  message("Simulated data already exist.")
  
} else {

temp_range <- seq(15,25,by=0.1)
resource_range <- seq(0.1,1,by=0.01)


seq_data <- tibble(temp=c(),
                   resource=c())

for (i in 1:length(temp_range)){
  tmp <- tibble(temp = temp_range[i],
                resource = seq(0.1,1,by=0.01))
  seq_data %<>% bind_rows(., tmp)
}



lm(mm ~ resource + temp, data = data)

length_coef <- coef(lm(mm ~ resource + temp, data = data)) #lengths from foraging rate experiment
#length_coef <- coef(lm(mm ~ resource + as.numeric(temp_id), data = lengths)) #5 day lengths from life table


interpolate_length <- function(resource, temp){
  length <- length_coef[1] + resource*length_coef[2] + temp*length_coef[3]
}

seq_data %<>% mutate(length = mapply(interpolate_length, 
                                     resource = resource,
                                     temp = temp))


# add interpolated lengths

length_range <- seq(min(length_summ$life_mm), max(length_summ$life_mm), by = 0.01)

length_seq <- expand_grid(temp = c(15, 20, 25), 
                          resource = c(0.1, 0.5, 1.0), 
                          length = length_range)

seq_data %<>% rbind(length_seq)

# add in experimental data

prevalence <- readRDS(here("processed_data", "prevalence.rds"))
prevalence %<>%  
  mutate(resource = as.numeric(as.character(resource))) %>%
  mutate(ID = paste(temp, resource, sep="_"))

seq_data %<>% left_join(., prevalence, by=c("temp", "resource"))


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

seq_data %<>% mutate(m5_len_I_end = mapply(m5_len_sim, 
                                      R=resource*life_vol/1000, 
                                      time=1, 
                                      f=m5_f/life_vol, 
                                      u = m5_u,
                                      m = m5_m,
                                      length=length, 
                                      gamma=gamma, 
                                      Z=200*life_vol,
                                      ref_t = ref_t,
                                      arr_t_f = m5_arr,
                                      arr_t_u = m5_arr_u,
                                      h = m5_h,
                                      temp = temp,
                                      resource = resource,
                                      w = m5_w,
                                      rho = m5_rho,
                                      phi = m5_phi))

seq_data %<>% mutate(m5_len_rate = (m5_f*exp(m5_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m5_f*exp(m5_arr*(1/ref_t - 1/temp))*
                          (length^gamma)*
                          m5_h*exp(m5_w*temp)*
                          resource/1000),
                     m5_len_susc = m5_u*
                       exp(m5_m*length)*
                       exp(m5_arr_u*(1/ref_t - 1/temp))*
                       exp(m5_rho*resource)*
                       exp(resource*temp*m5_phi))



# spores_consumed ---------------------------------------------------------


m1_sim_Z <- function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, h, temp){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=0,
              length = length,
              gamma = gamma,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              h = h,
              temp = temp)
  output <- as.data.frame(lsoda(y=xstart, times, m1_num_sol, params))
  
  end_data <- slice_min(output, time)[,5] - slice_max(output, time)[,5]
  
  
  return(end_data)
}



seq_data %<>% mutate(spores_consumed_m5_len = mapply(m1_sim_Z, 
                                                 R=resource*life_vol/1000, 
                                                 time=1, 
                                                 f=m5_f/life_vol, 
                                                 u = 0, 
                                                 length=length, 
                                                 gamma=gamma, 
                                                 Z=200*life_vol,
                                                 ref_t = ref_t,
                                                 arr_t_f = m5_arr,
                                                 h = m5_h*exp(m5_w*temp),
                                                 temp = temp))


u_seq <- readRDS(file=here("processed_data", "seq_data", "infection_fit_data.rds"))

u_seq %<>% right_join(., seq_data)

m5_no_len <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))
m5_f_beta <- coef(m5_no_len)[1]
m5_arr_beta <- coef(m5_no_len)[3]
m5_h_beta <- coef(m5_no_len)[5]
m5_w_beta <- coef(m5_no_len)[6]

u_seq %<>% mutate(spores_consumed_m5 = 
                    case_when(is.na(spores_consumed_m5) ~
                                              mapply(m1_sim_Z, 
                                                     R=resource*life_vol/1000, 
                                                     time=1, 
                                                     f=m5_f_beta/life_vol, 
                                                     u = 0, 
                                                     length=length, 
                                                     gamma=gamma, 
                                                     Z=200*life_vol,
                                                     ref_t = ref_t,
                                                     arr_t_f = m5_arr_beta,
                                                     h = m5_h_beta*exp(m5_w_beta*temp),
                                                     temp = temp),
                              .default = spores_consumed_m5),
                  m5_rate = 
                    case_when(is.na(m5_rate) ~
                                (m5_f_beta*exp(m5_arr_beta*(1/ref_t - 1/temp))*(length^gamma))/
                                (1+m5_f_beta*exp(m5_arr_beta*(1/ref_t - 1/temp))*
                                   (length^gamma)*
                                   m5_h_beta*exp(m5_w_beta*temp)*
                                   resource/1000),
                              .default = m5_rate))

u_seq %>% arrange(resource, temp) %>% View()

saveRDS(u_seq, file=here("processed_data", "seq_data", "supp_body_size_fit_data.rds"))

}



# make figures ------------------------------------------------------------



u_seq <- readRDS(file=here("processed_data", "seq_data", "supp_body_size_fit_data.rds"))


margin <- theme(plot.margin = unit(c(0,0,0,0), "cm"))


#spores consumed
spores_comp <-
u_seq %>% filter(temp %in% const_temp) %>% drop_na(m1_rate) %>%
  ggplot(., aes(x=resource, y=spores_consumed_m5, color=as.factor(temp), group = as.factor(temp))) + 
  geom_line(size = 1, linetype = "dashed") + 
  labs(x = "", 
       y = "",
       color = "Temp") + 
  geom_line(size = 1, alpha = 0.5, aes(y=spores_consumed_m5_len)) + 
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  scale_y_continuous(breaks = c(0, 2000, 4000, 6000), limits = c(0, 6100)) + 
  guides(color = guide_legend(reverse=T)) + 
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  theme(legend.position = "none", axis.text.x = element_blank()) + margin



#per-parasite susceptibility
u_comp <- 
u_seq %>% 
  filter(temp %in% const_temp) %>% drop_na(m1_rate) %>%
  ggplot(.,aes(x=resource, y=m5_susc, color=as.factor(temp))) + 
  geom_line(size = 1, alpha = 1, linetype = "dashed") + 
  labs(x = "", y = "", color = "Temp", 
       title = "") +
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  geom_line(size = 1, aes(y=m5_len_susc)) + 
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  scale_y_continuous(breaks = c(0, 0.0025, 0.005), limits = c(0, 0.005)) + 
  guides(color = guide_legend(reverse=T)) +
  theme(legend.position = "none", axis.text.x = element_blank()) + margin

u_comp_rsc <- 
u_seq %>%
  filter(resource %in% c(0.1, 0.5, 1)) %>% drop_na(m1_rate) %>%
  ggplot(.,aes(x=temp, y=m5_susc, color=as.factor(resource))) + 
  geom_line(size = 1, alpha = 1, linetype = "dashed") + 
  labs(x = "", y = "", color = "Resource", 
       title = "") +
  scale_color_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
  geom_line(size = 1, aes(y=m5_len_susc)) + 
  scale_x_continuous(breaks = c(15, 20, 25)) + 
  scale_y_continuous(breaks = c(0, 0.0025, 0.005), limits = c(0, 0.005)) + 
  guides(color = guide_legend(reverse=T)) +
  theme(legend.position = "none", axis.text.x = element_blank()) + margin
  


#prevalence
prev_comp <-
u_seq %>% filter(temp %in% const_temp) %>% drop_na(m1_rate) %>%
  ggplot(.,aes(x=resource, 
               y=m5_I_end, 
               color=as.factor(temp))) + 
  geom_line(size=1, alpha = 1, linetype = "dashed") +
  geom_line(size = 1, aes(y=m5_len_I_end)) + 
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  new_scale_color() + 
  geom_point(aes(x=resource, y=prev, color=as.factor(temp)), size=3) + 
  geom_linerange(aes(ymin = conf$lower, ymax=conf$upper, color=as.factor(temp))) +
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) + 
  labs(x="Resource (mg/L)", y="", title = "") +
  theme(legend.position = "none") + margin



# examine length effects --------------------------------------------------

#by resource
length_spores_rsc <-
u_seq %>% 
  filter(temp %in% const_temp &
           resource %in% resource_IDs) %>% 
  ggplot(.,aes(x=length, y=spores_consumed_m5_len, color = as.factor(resource))) + 
  geom_line(size = 1) + 
  scale_color_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
  geom_line(aes(y=spores_consumed_m5, color = as.factor(resource)), linetype = "dashed") + 
  scale_x_continuous(breaks = c(0.8, 1.0, 1.2)) + 
  guides(color = guide_legend(reverse=T)) + 
  labs(x = "Length (mm)", y = "Spores Consumed", color = "Resource", 
       title = "") +
  facet_wrap(~temp, nrow = 1)

length_spores <- 
u_seq %>% 
  filter(temp %in% const_temp &
           resource %in% resource_IDs) %>% 
  ggplot(.,aes(x=length, y=spores_consumed_m5_len, color = as.factor(temp))) + 
  geom_line(size = 1, alpha = 0.5) + 
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  geom_line(aes(y=spores_consumed_m5, color = as.factor(temp)), linetype = "dashed") + 
  scale_x_continuous(breaks = c(0.8, 1.0, 1.2)) + 
  guides(color = guide_legend(reverse=T)) + 
  labs(x = "", y = "", color = "Temp",
       title = "") +
  facet_wrap(~resource, nrow = 1) +
  theme(axis.text.x = element_blank(), legend.position = "none") + margin


make_length_plot_rsc <- function(data, var_alpha, var_beta, y_label){
  data %>% 
    filter(temp %in% const_temp &
             resource %in% resource_IDs) %>% 
    ggplot(.,aes(x=length, y=!!sym(var_alpha), color = as.factor(resource))) + 
    geom_line(size = 1) + 
    scale_color_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
    geom_hline(aes(yintercept=!!sym(var_beta), color = as.factor(resource)), linetype = "dashed") + 
    scale_x_continuous(breaks = c(0.8, 1.0, 1.2)) + 
    guides(color = guide_legend(reverse=T)) + 
    labs(x = "Length (mm)", y = y_label, color = "Resource", 
         title = "") +
    facet_wrap(~temp, nrow = 1)
}

length_prev_rsc <- make_length_plot_rsc(u_seq, "m5_len_I_end", "m5_I_end", "")
length_susc_rsc <- make_length_plot_rsc(u_seq, "m5_len_susc", "m5_susc", "")


#by temp
make_length_plot_temp <- function(data, var_alpha, var_beta, x_label, y_label){
  data %>% 
    filter(temp %in% const_temp &
             resource %in% resource_IDs) %>% 
    ggplot(.,aes(x=length, y=(!!sym(var_alpha)), color = as.factor(temp))) + 
    geom_line(size = 1) + 
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    geom_hline(aes(yintercept=!!sym(var_beta), color = as.factor(temp)), linetype = "dashed") + 
    scale_x_continuous(breaks = c(0.8, 1.0, 1.2)) + 
    guides(color = guide_legend(reverse=T)) + 
    labs(x = x_label, y = y_label, color = "Temp",
         title = "") +
    facet_wrap(~resource, nrow = 1) +
    theme(legend.position = "none")
}

length_prev <- make_length_plot_temp(u_seq, "m5_len_I_end", "m5_I_end", x_label = "Length (mm)", y_label = "") + 
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) + 
  margin

length_susc <- make_length_plot_temp(u_seq, "m5_len_susc", "m5_susc", x_label = "", y_label = "") + 
  scale_y_continuous(breaks = c(0, 0.0025, 0.005)) + 
  theme(axis.text.x = element_blank()) + margin


(spores_comp / u_comp / prev_comp) | (length_spores / length_susc / length_prev)

ggsave(here("figures", "supp_model_results.png"), 
       plot = last_plot(),
       width = 12,
       height = 12,
       units = "in")



