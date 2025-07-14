############################
# Project: KP VACCINATE
# Author: Rishi Parikh
# Purpose: Forest plots for subgroups analysis
#   - Use meta-analyzed subgroup effects
#   - Plot an overall effect row
#   - Display stratified estimates and interaction p-values
############################

#Install and load packages
list.of.packages <- c("dplyr", "tidyr","ggplot2","forcats","grid")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(grid) 

results = read.csv(paste0(path,"output/subgroup_effects_combined.csv"))

#Create nice labels in results
results = mutate(results,
                 level = ifelse(group=="ndi_4tile_22" & level=="1","Q1",level))
results =  mutate(results,
                  groupf = factor(group, levels=unique(results$group),
                                 labels = c("Site","Age","Sex","Self-reported race",
                                            "Hispanic ethnicity","Needs interpreter","Vaccination history",
                                            "Vaccinated in prior season",
                                            "Neighborhood deprivation index",
                                            "Chronic cardiovascular disease","Ischemic heart disease",
                                            "Myocardial infarction","Heart failure","Atrial fibrillation",
                                            "Cerebrovascular disease","Diabetes","Chronic kidney disease")),
                 levelf = factor(level, levels=unique(results$level),
                                 labels = c("KPNC","KPMAS","<55","55-64","\u226565","Male","Female","White","American Indian/Alaska Native",
                                            "Asian","Black","Native Hawaiian/Pacific Islander","Multi-racial","Other",
                                            "No","Yes","Non-believer","Inconsistent","Believer","Unknown","Quartile 1","Quartile 2","Quartile 3","Quartile 4",
                                            "No","Yes")))

#Overall effects
overall = read.csv(paste0(path,"output/primary_comparisons_combined.csv"))

overall_comp = filter(results,term=="comp1")


#Choose comparison and filter
comps = c('comp1','comp2','comp3','comp4','comp5','comp6')

for (comp in comps) {
  results_em = filter(results,comparison==comp)
  overall_comp = filter(overall,term==comp)
  
  
  overall_effect <- data.frame(
    group = "Overall Effect",
    level = "Overall",
    estimate = overall_comp$estimate,
    conf.low = overall_comp$conf.low,
    conf.high = overall_comp$conf.high,
    estimate_pct = round(estimate * 100, 2),
    conf_low_pct = round(conf.low * 100, 2),
    conf_high_pct = round(conf.high * 100, 2),
    formatted_estimate = paste0(estimate_pct, " (", conf_low_pct, ", ", conf_high_pct, ")"),
    formatted_p = ""
  )
  
  # Format estimates and interaction p-values
  results_em <- results_em %>%
    mutate(
      estimate_pct = round(estimate * 100, 2),
      conf_low_pct = round(conf.low * 100, 2),
      conf_high_pct = round(conf.high * 100, 2),
      formatted_estimate = paste0(estimate_pct, " (", conf_low_pct, ", ", conf_high_pct, ")"),
      
      # Format interaction p-value
      formatted_p = case_when(
        is.na(interactionp) ~ "Ref",
        interactionp < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", interactionp)
      ),
      order = row_number()
    )
  
  # Get unique groups
  groups <- as.data.frame(results_em$groupf) %>% group_by(results_em$groupf) %>% slice(1)
  groups <- groups$`results_em$groupf`
  
  # Create headers for each group (without estimate and p-values)
  headers_df <- data.frame(
    groupf = groups,
    levelf = NA,
    estimate = NA,
    overallse = NA,
    conf.low = NA,
    conf.high = NA,
    formatted_estimate = "",
    formatted_p = ""
  )
  
  # Combine data and sort
  combined_df <- bind_rows(
    headers_df,
    results_em %>% select(group, level, groupf, levelf, estimate, overallse, conf.low, conf.high, formatted_estimate, formatted_p,interactionp, order)) %>%
    mutate(order = ifelse(is.na(level),0,order)) %>%
    arrange(groupf, order) %>%
    mutate(
      row_label = case_when(
        is.na(level) ~ as.character(groupf),
        !is.na(level) ~ paste0("   ", levelf)
      ),
      is_header = is.na(level)
    ) %>%
    mutate(row_id = n():1)
  
  combined_df <- bind_rows(overall_effect, combined_df) %>%
    mutate(
      row_label = case_when(
        group == "Overall Effect" ~ "Overall Effect",
        is.na(level) ~ as.character(groupf),
        !is.na(level) ~ paste0("   ", levelf)
      ),
      is_header = is.na(level)
    ) %>%
    mutate(row_id = n():1)
  
  # Plot data excluding headers
  plot_data <- combined_df %>% filter(!is_header)
  overall_data <- plot_data %>% filter(group == "Overall Effect")
  left_margin = 0.07
  x_min <- -0.02
  x_max <- 0.02
  
  diamond_data <- overall_data %>%
    mutate(
      y = row_id,
      x_center = estimate,
      x_low = conf.low,
      x_high = conf.high
    ) %>%
    rowwise() %>%
    do({
      data.frame(
        x = c(.$x_low, .$x_center, .$x_high, .$x_center),
        y = c(.$y, .$y + 0.7, .$y, .$y - 0.7)
      )
    })
  
  # Create clipped error bar data and flags for arrows
  plot_data <- plot_data %>%
    mutate(
      ci_low_clipped = pmax(conf.low, x_min),
      ci_high_clipped = pmin(conf.high, x_max),
      left_arrow = conf.low < x_min,
      right_arrow = conf.high > x_max
    )
  
  forest_plot <- ggplot() +
    # Reference line
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Vertical line for overall effect
    geom_vline(data = overall_data, aes(xintercept = estimate), linetype = "dotted", color = "red") +
    # Clipped error bars
    geom_errorbarh(data = plot_data, 
                   aes(y = row_id, xmin = ci_low_clipped, xmax = ci_high_clipped), 
                   height = 0.2) +
    
    # Left arrows for clipped lower limits
    geom_segment(data = plot_data %>% filter(left_arrow),
                 aes(y = row_id, yend = row_id, x = x_min + 0.0005, xend = x_min),
                 arrow = arrow(length = unit(0.08, "cm"), ends = "first"), 
                 color = "black") +
    
    # Right arrows for clipped upper limits
    geom_segment(data = plot_data %>% filter(right_arrow),
                 aes(y = row_id, yend = row_id, x = x_max - 0.0005, xend = x_max),
                 arrow = arrow(length = unit(0.08, "cm"), ends = "last"), 
                 color = "black") +
    
    # Points
    geom_point(data = plot_data %>% filter(group != "Overall Effect"), 
               aes(y = row_id, x = estimate, size = 0.5/overallse), shape = 15) +
    
    # Overall diamond
    geom_polygon(data = diamond_data, aes(x = x, y = y), fill = "darkred", color = "darkred") +
    
    # Left-side variable labels
    geom_text(data = combined_df, 
              aes(y = row_id, x = -0.045, label = row_label, 
                  fontface = ifelse(is_header, "bold", "plain")), 
              hjust = 0, size = 3.5) +
    
    # Estimate column
    geom_text(data = combined_df, 
              aes(y = row_id, x = 0.035, label = formatted_estimate), 
              hjust = 1, size = 3.2) +
    
    # P-value column
    geom_text(data = combined_df, 
              aes(y = row_id, x = 0.037, label = formatted_p), 
              hjust = 0, size = 3.2) +
    
    # Column headers
    annotate("text", x = -0.045, y = max(combined_df$row_id) + 1.5, 
             label = "Subgroup", hjust = 0, fontface = "bold", size = 3.5) +
    annotate("text", x = 0.035, y = max(combined_df$row_id) + 1.5, 
             label = "Estimate (99.2% CI)", hjust = 1, fontface = "bold", size = 3.5) +
    annotate("text", x = 0.037, y = max(combined_df$row_id) + 1.5, 
             label = "Interaction\np-value", hjust = 0, fontface = "bold", size = 3.5) +
    
    # X-axis settings
    scale_size_continuous(range = c(2, 5), guide = "none") +
    scale_x_continuous(limits = c(-0.15, 0.15), breaks = seq(-0.03, 0.03, 0.005),
                       labels = scales::percent_format(accuracy = 0.1)) +
    coord_cartesian(xlim = c(-0.018, 0.018), clip = "off") +
    
    labs(x = "Marginal Percentage Point Increase (99.2% CI) in Influenza Vaccination") +
    
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.title.x = element_text(face = "bold"),
      plot.margin = margin(t = 20, r = 160, b = 20, l = 160, unit = "pt")
    )
  
  print(forest_plot)
  ggsave(paste0(path,"output/subgroup_forest_",comp,".png"),device="png",height=12,width=7.9,units="in",dpi=600)
  
}




