#Summary of length data
library(here)

source(here("base","src.R"))

length_summ <- readRDS(here("processed_data", "life_length_summ.rds"))
data_summ <- readRDS(here("processed_data", "foraging.rds")) 



# length comparison -------------------------------------------------------


all_lengths <- 
  data_summ %>%
  mutate(n = length_n) %>%
  dplyr::select(c(temp, resource, mm_mean, mm_sd, mm_se, n)) %>%
  mutate(assay = "Foraging") %>%
  rbind(., 
        length_summ %>%
          rename(temp = temp_id,
                 mm_mean = life_mm,
                 mm_sd = sd,
                 mm_se = se) %>%
          dplyr::select(temp, resource, mm_mean, mm_sd, mm_se, n) %>%
          mutate(assay = "Infection")) %>%
  mutate(ID = paste(temp, resource, sep = "_")) %>%
  arrange(temp, resource)

all_lengths %>%
  ggplot(., aes(x = as.factor(resource), y = mm_mean, color = assay, group = assay)) + 
  geom_point(position = position_dodge(width = 1)) + 
  geom_errorbar(aes(ymin = mm_mean - 2*mm_sd, ymax = mm_mean + 2*mm_sd), position = position_dodge(width = 1)) +
  facet_wrap(~temp) +
  labs(x = "Resource (mg/L)", y = "Body Length (mm)", color = "Assay")

ggsave(here("figures", "supp_body_length.png"), 
       plot = last_plot(),
       width = 10,
       height = 7,
       units = "in")
