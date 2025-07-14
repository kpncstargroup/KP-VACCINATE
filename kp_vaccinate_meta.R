############################
# Project: KP VACCINATE
# Author: Rishi Parikh
# Purpose: Meta analysis of KPNC and KPMAS results
#   - Load and clean aggregate results data
#   - Perform fixed effects meta analysis
#   - Plot combined cumulative incidence curves and primary 
#     comparison differences.
#   - Perform meta analyses for subgroup effects 
############################

#Install and load packages
list.of.packages <- c("haven","dplyr", "tidyr","broom","ggplot2","survival","survminer","marginaleffects","meta")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(survival)
library(haven)
library(survminer)
library(marginaleffects)
library(meta)

#Set filepath string
path = "~/kp_vaccinate/"

#Meta analysis of primary effects
################
#Load site-specific results
nc_results = read.csv(paste0(path,"output/primary_comparisons_kpnc.csv")) %>% mutate(site="KPNC")
mas_results = read.csv(paste0(path,"output/primary_comparisons_kpmas.csv")) %>% mutate(site="KPNC")
combined= bind_rows(nc_results,mas_results)


comps = c('comp1','comp2','comp3','comp4','comp5','comp6')

results = NULL
for (comp in comps) {
  data = filter(combined,term==comp)
  
  m.gen <- metagen(TE = estimate,
                   seTE = std.error,
                   studlab = site,
                   data = data,
                   fixed = T,
                   level=0.992,
                   level.comb=0.992
  )
  
  
  result = data.frame(
    comparison = comp,
    total = paste0(round(m.gen$TE.common*100,3), " (",round(m.gen$lower.common*100,3),", ",round(m.gen$upper.common*100,3),")"),
    kpnc = paste0(round(m.gen$TE[[2]]*100,3), " (",round(m.gen$lower[[2]]*100,3),", ",round(m.gen$upper[[2]]*100,3),")"),
    kpmas = paste0(round(m.gen$TE[[1]]*100,3), " (",round(m.gen$lower[[1]]*100,3),", ",round(m.gen$upper[[1]]*100,3),")")
  )
  results = rbind(results,result)
}

write.csv(results,paste0(path,"output/primary_comparisons_combined.csv"))

#Combined cumulative incidence curves
#############################

nc_model <- readRDS(paste0(path,"output/survfit_KPNC.Rds"))
mas_model <- readRDS(paste0(path,"output/survfit_KPMAS.Rds"))

kpnc_data = data.frame(
  time = nc_model$time,
  risk = nc_model$n.risk,
  event = nc_model$n.event,
  censor = nc_model$n.censor
) %>%
  mutate(strata=cumsum(time==1),
         group = factor(case_when(
           strata==1 ~ 'CVCV',
           strata==2 ~ 'CVUC',
           strata==3 ~ 'UCCV',
           strata==4 ~ 'UCUC'
         ), levels = c('UCUC','CVUC','UCCV','CVCV'))
  )

kpmas_data = data.frame(
  time = mas_model$time,
  risk = mas_model$n.risk,
  event = mas_model$n.event,
  censor = mas_model$n.censor
) %>%
  mutate(strata=cumsum(time==1),
         group = factor(case_when(
           strata==1 ~ 'CVCV',
           strata==2 ~ 'CVUC',
           strata==3 ~ 'UCCV',
           strata==4 ~ 'UCUC'
         ), levels = c('UCUC','CVUC','UCCV','CVCV'))
  )

all_data = left_join(kpnc_data,kpmas_data,by=c('time','group')) %>%
  rowwise() %>%
  mutate(risk = sum(risk.x,risk.y,na.rm=T),
         event = sum(event.x,event.y,na.rm=T),
         censor = sum(censor.x,censor.y,na.rm=T)) %>%
  dplyr::select(time,group,risk,event,censor)


#Expand back into individual-level
expand_df <- all_data %>%
  group_by(group, time) %>%
  reframe(
    status = c(rep(1, event), rep(0, censor))
  )

expand_mas <- kpmas_data %>%
  group_by(group, time) %>%
  reframe(
    status = c(rep(1, event), rep(0, censor))
  )

expand_nc <- kpnc_data %>%
  group_by(group, time) %>%
  reframe(
    status = c(rep(1, event), rep(0, censor))
  )

#Full model
model = survfit(Surv(time,status) ~ group,data=expand_df)
mas_model = survfit(Surv(time,status) ~ group,data=expand_mas)
nc_model = survfit(Surv(time,status) ~ group,data=expand_nc)

saveRDS(model,paste0(path,"output/survfit_combined.Rds"))

#Difference in cumulative incidence curves by co-primary comparison
#################
model_df = tidy(model) %>%
  mutate(cuminc=(1-estimate),
         group = gsub("group=", "", strata))

cuminc_wide <- model_df %>%
  select(time, group, cuminc, std.error) %>%
  pivot_wider(
    names_from = group,
    values_from = c(cuminc, std.error),
    names_sep = "."
  ) %>%
  mutate(comp1_diff = cuminc.CVCV - cuminc.UCUC,
         comp1_se = sqrt(std.error.CVCV^2 + std.error.UCUC^2),
         comp2_diff = cuminc.CVUC - cuminc.UCUC,
         comp2_se = sqrt(std.error.CVUC^2 + std.error.UCUC^2),
         comp3_diff = cuminc.UCCV - cuminc.UCUC,
         comp3_se = sqrt(std.error.UCCV^2 + std.error.UCUC^2),
         comp4_diff = cuminc.CVCV - cuminc.CVUC,
         comp4_se = sqrt(std.error.CVCV^2 + std.error.CVUC^2),
         comp5_diff = cuminc.CVUC - cuminc.UCCV,
         comp5_se = sqrt(std.error.CVUC^2 + std.error.UCCV^2)
  )

#Comp6- re-model 
expand_df = mutate(expand_df,comp6=ifelse(group=='UCUC',0,1))
model6 = survfit(Surv(time,status) ~ comp6,data=expand_df)
model_df6 = tidy(model6) %>%
  mutate(cuminc=(1-estimate),
         group = gsub("comp6=", "", strata),
         group = ifelse(group=="0",'UCUC','anyCV'))

cuminc_wide2 <- model_df6 %>%
  select(time, group, cuminc, std.error) %>%
  pivot_wider(
    names_from = group,
    values_from = c(cuminc, std.error),
    names_sep = "."
  ) %>%
  select(time,cuminc.anyCV,std.error.anyCV) %>%
  left_join(cuminc_wide,by="time") %>%
  mutate(comp6_diff = cuminc.anyCV - cuminc.UCUC,
         comp6_se = sqrt(std.error.anyCV^2 + std.error.UCUC^2))

plot_df <- cuminc_wide2 %>%
  select(time, starts_with("comp")) %>%
  pivot_longer(
    cols = starts_with("comp"),
    names_to = c("comparison", ".value"),
    names_pattern = "comp(\\d+)_(diff|se)"
  ) %>%
  mutate(
    comparison = paste("Comparison", comparison),
    ci_lower = diff - 2.41 * se,
    ci_upper = diff + 2.41 * se
  )

comparison_labels <- c(
  "Comparison 1:\n(Cardiovascular/Cardiovascular - Usual Care/Usual Care)",
  "Comparison 2:\n(Cardiovascular/Usual Care - Usual Care/Usual Care)",
  "Comparison 3:\n(Usual Care/Cardiovascular - Usual Care/Usual Care)",
  "Comparison 4:\n(Cardiovascular/Cardiovascular - Cardiovascular/Usual Care)",
  "Comparison 5:\n(Cardiovascular/Usual Care - Usual Care/Cardiovascular)",
  "Comparison 6:\n(Any Cardiovascular Message - Usual Care/Usual Care)"
)

plot_df <- plot_df %>%
  mutate(
    comparison = factor(comparison, levels = paste("Comparison", 1:6), labels = comparison_labels),
    diff = diff * 100,
    ci_lower = ci_lower * 100,
    ci_upper = ci_upper * 100
  )

ggplot(plot_df, aes(x = time, y = diff)) +
  geom_line(color = "steelblue") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "steelblue", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits=c(-1.5,1.5)) +
  facet_wrap(~ comparison, ncol = 2) +
  labs(
    x = "Days from Randomization",
    y = "Difference in Cumulative Incidence of Vaccination, % (99.2% CI)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10)
  )

ggsave(paste0(path,"output/cuminc_diff_",site,".png"),device="png",height=12,width=7,units="in",dpi=300)


#Meta analyze subgroup effects
###################### 
kpnc_subgroups = read.csv(paste0(path,"output/subgroup_effects_KPNC.csv")) %>%
  mutate(site='kpnc',
         term2=paste(term,group,level,sep="-"))


kpmas_subgroups = read.csv(paste0(path,"output/subgroup_effects_KPMAS.csv")) %>%
  mutate(site='kpma',
         term2=paste(term,group,level,sep="-"))

all_subgroups = bind_rows(kpnc_subgroups,kpmas_subgroups)


comps = unique(all_subgroups$term2)

results = NULL
for (comp in comps) {
  data = filter(all_subgroups,term2==comp)
  
  #Stratified estimates
  m.gen <- metagen(TE = estimate,
                   seTE = std.error,
                   studlab = site,
                   data = data,
                   fixed = T,
                   level=0.992,
                   level.comb=0.992
  )
  
  #Interaction term estimates
  m.gen2 <- metagen(TE = estimate.1,
                    seTE = std.error.1,
                    studlab = site,
                    data = data,
                    fixed = T,
                    level=0.992,
                    level.comb=0.992
  )
  
  result = data.frame(
    comparison = strsplit(comp,"-")[[1]][[1]],
    group = strsplit(comp,"-")[[1]][[2]],
    level = strsplit(comp,"-")[[1]][[3]],
    total = paste0(round(m.gen$TE.common*100,3), " (",round(m.gen$lower.common*100,3),", ",round(m.gen$upper.common*100,3),")"),
    kpnc = paste0(round(m.gen$TE[[2]]*100,3), " (",round(m.gen$lower[[2]]*100,3),", ",round(m.gen$upper[[2]]*100,3),")"),
    kpmas = paste0(round(m.gen$TE[[1]]*100,3), " (",round(m.gen$lower[[1]]*100,3),", ",round(m.gen$upper[[1]]*100,3),")"),
    estimate = m.gen$TE.common,
    conf.low = m.gen$lower.common,
    conf.high = m.gen$upper.common,
    overallse = m.gen$seTE.common,
    overallp = m.gen$pval.common,
    interactionp = m.gen2$pval.common
  )
  results = rbind(results,result)
}

write.csv(results,paste0(path,"output/subgroup_effects_combined.csv"))


