############################
# Project: KP VACCINATE
# Author: Rishi Parikh
# Purpose: Primary/secondary outcome analysis: to be run for each site separately
#   - Load and clean data
#   - Calculate marginal risk differences for each co-primary comparison
#   - Save and plot aggregate cumulative incidence function
#   - Plot difference in cumulative incidence for each co-primary comparison
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
site = "KPNC"

#Load SAS dataset
df = read_sas(paste0(path,"data/baseline_analytic.sas7bdat")) %>% 
  dplyr::select(Age,SEX,MC,OUTREACH_DT,GROUPNM,index_date,prior_season,demo_race,demo_hispanic,census_ndi,
                nudge_cvd,nudge_ischemic_hd,nudge_mi,nudge_hf,nudge_af,nudge_cerebro,nudge_pvd,comorb_diabetes,nonenglish,interpreter,lab_outpat_egfr,
                comp1,comp2,comp3,comp4,comp5,comp6,vax_days,vax_2024, vax_history)


df = df %>% mutate(outreachdt = as.factor(OUTREACH_DT),
                   mc = as.factor(MC),
                   female =factor(ifelse(SEX=='F',1,0),levels = c(0,1), labels = c('Male','Female')),
                   agec = factor(case_when(
                     Age < 55 ~ 1,
                     Age >=55 & Age <65 ~ 2,
                     Age >= 65 ~ 3),
                     levels = c(1,2,3),
                     labels = c('<55','55-64','>=65')),
                   race = factor(demo_race,levels =c('F','A','B','C','D','E','G'),
                                 labels = c('White','AmInd','Asian','Black','NHPI','Multiracial','Other')),
                   vax_history = factor(vax_history,levels = c('Non-Believer','Inconsistent','Believer','Unknown')),
                   ndi_4tile_22 = factor(case_when(
                     !is.na(census_ndi) &  census_ndi < -0.918206413481892 ~ 1,
                     !is.na(census_ndi) & -0.918206413481892  <= census_ndi & census_ndi < -0.342035691796302 ~ 2,
                     !is.na(census_ndi) & -0.342035691796302  <= census_ndi & census_ndi < 0.3807222509529160 ~ 3,
                     !is.na(census_ndi) & 0.3807222509529160  <= census_ndi ~ 4,
                     .default = NA 
                   ), levels = c(1,2,3,4)),
                   ckd = factor(ifelse(lab_outpat_egfr < 60,1,0), levels = c(0,1), labels =c('No CKD','CKD')),
                   
)

#Marginal effects for co-primary outcome comparisons
#################
#Binary indicator variables for each comparison
comparisons = c("comp1","comp2","comp3","comp4","comp5","comp6")

all_comparisons = function(comparisons) {
  
  for (i in seq_along(comparisons)) {
    formula = as.formula(paste0("vax_2024 ~ ",comparisons[i]," + outreachdt + mc"))
    model = glm(formula,data=df,family ='binomial')
    if (i==1) { result = avg_comparisons(model, variables=comparisons[i], conf_level=0.992)}
    else {result = rbind(result,avg_comparisons(model, variables=comparisons[i], conf_level=0.992))}
  }
  return(result)
}

results = all_comparisons(comparisons)

write.csv(results,paste0(path,"output/primary_comparisons_",site,".csv"))


#KM Curves by Group
##############################

model <- survfit(Surv(vax_days,vax_2024) ~ GROUPNM,data=df)

#Save model for meta analysis
saveRDS(model,paste0(path,"output/survfit_",site,".Rds"))

survplot <- ggsurvplot(
  model,
  data = df,
  fun = function(y) (1-y)*100,
  censor = F,
  conf.int = F,            
  pval = F,                
  risk.table = TRUE,          
  risk.table.height = 0.3,    
  xlab = "Days from Randomization",
  ylab = "Cumulative Incidence of Vaccination, %",
  palette = c("darkorange1","royalblue4","green4","red"),            
  legend.title = "Group",
  legend.labs = c("Cardiovascular/Cardiovascular","Cardiovascular/Usual Care","Usual Care/Cardiovascular","Usual Care/Usual Care"),
  tables.y.text = FALSE,
  ggtheme = theme_bw() +      
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
)

survplot$plot <- survplot$plot + theme(plot.title = element_text(hjust = 0),
                                       axis.title.y= element_text(margin=margin(t=0,r=10,b=0,l=0)),
                                       axis.title.x= element_text(margin=margin(t=15,r=0,b=0,l=0)),
                                       plot.margin = unit(c(5.5,5.5,5.5,20),"points"),
                                       legend.text=element_text(size=12),
                                       legend.title=element_text(size=12)) +
  scale_y_continuous(limits = c(0,50),expand = c(0,0)) + scale_x_continuous(limits =c(115,120),breaks = seq(115,120,1))

survplot$table <- survplot$table + theme(plot.title = element_blank(),
                                         plot.margin = unit(c(-20,5.5,5.5,20),"points"),
                                         panel.grid = element_blank())


ggsave(paste0(path,"output/km_cuminc_",site,".png"),device="png",height=7,width=9,units="in",dpi=300)


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
model6 = survfit(Surv(time,status) ~ comp6,data=df)
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


#Subgroup effects
#################

comparisons = c("comp1","comp2","comp3","comp4","comp5","comp6")

modifiers = c("agec","female","race","demo_hispanic","interpreter","vax_history",
              "prior_season","ndi_4tile_22","nudge_cvd","nudge_ischemic_hd","nudge_mi",
              "nudge_hf","nudge_af","nudge_cerebro","comorb_diabetes","ckd")

all_comparisons_em = function(comparisons,modifiers) {
  
  for (i in seq_along(comparisons)) {
    for (x in seq_along(modifiers)) {
      formula = as.formula(paste0("vax_2024 ~ ",comparisons[i],"*",modifiers[x]," + outreachdt + mc"))
      model = glm(formula,data=df,family ='binomial')
      if (i==1 & x==1) { 
        result = avg_comparisons(model, variables=comparisons[i], by=modifiers[x], conf_level=0.992)%>%
          rename(level := !!modifiers[x]) %>%
          mutate(level = as.character(level),
                 group = modifiers[x])
        pval = avg_comparisons(model, variables=comparisons[i], by=modifiers[x], hypothesis = "b* - b1 = 0", conf_level=0.992)
        result = cbind(result,pval)
      }
      else {
        newr = avg_comparisons(model, variables=comparisons[i], by=modifiers[x], conf_level=0.992) %>%
          mutate(group = modifiers[x]) %>%
          rename(level := !!modifiers[x])
        pval = avg_comparisons(model, variables=comparisons[i], by=modifiers[x], hypothesis = "b* - b1 = 0", conf_level=0.992)
        newr = cbind(newr,pval)
        result = rbind(result,newr)
      }
    }
  }
  return(result)
}

results_em = all_comparisons_em(comparisons,modifiers)

write.csv(results_em,paste0(path,"output/subgroup_effects_",site,".csv"))
