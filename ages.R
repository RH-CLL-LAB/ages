# This script calculates age-distributions and compares with IPIs
# date: "2025-01-02"
# author: christian brieghel

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
load_dataset(c('patient', 't_dalycare_diagnoses'))
load_dataset(RKKP_DATASETS)

# IPIs were calculated in IPI_raw.R; source first 
IPI = read_csv2('**censored**/IPI.csv')  %>% 
  mutate(IPI = factor(IPI, levels = c('Low', 'Intermediate', 'High', 'Very high')))

load_dalycare_icd10()
patient %>% nrow_npatients()

disease_var = c("All",   "FL", "CLL", "cHL",  "MZL", "LPL", "DLBCL","MCL", "MM")
ALL = c(ICD10.FL, ICD10.CLL, ICD10.HL, ICD10.MZL, ICD10.LPL, ICD10.DLBCL, ICD10.MCL, ICD10.MM)
list_disease = list(ALL, ICD10.FL, ICD10.CLL, ICD10.HL, ICD10.MZL, ICD10.LPL, ICD10.DLBCL, ICD10.MCL, ICD10.MM)
palette(BLOOD)
list_plot = list()
list_df = list()
for(i in 1:length(list_disease)){
  list_df[[i]] = t_dalycare_diagnoses %>% 
    filter(datasource == 'RKKP') %>% # only RKKP diagnoses to compare Ages vs IPIs fairly
    filter_first_diagnosis(list_disease[[i]], str_contains = F) %>% 
    mutate(age = diff_years(date_birth, date_diagnosis),
           ages = cut(age, c(-Inf, 40, 50, 60, 70, 80, Inf), right = FALSE),
           Disease = disease_var[i]) %>% 
    filter(time_dx_death >=0) %>% 
    left_join(IPI %>% select(patientid, Disease, IPI))
  
  list_plot[[i]] = KM_plot(survfit(Surv(time_dx_death, status) ~ ages, list_df[[i]]),
                           title = disease_var[i],
                           labs = c('<40', '40-49', '50-59', '60-69', '70-79', '≥80'),
                           breaks = 3,
                           xlim = c(0, 15),
                           pval = T)
}

list_df[[1]] %>% nrow_npatients()
list_df[[2]] %>% nrow_npatients()
list_df[[3]] %>% nrow_npatients() 

ENTIRE_COHORT = bind_rows(list_df[[2]],  
                          list_df[[3]],  
                          list_df[[4]],  
                          list_df[[5]],  
                          list_df[[6]], 
                          list_df[[7]], 
                          list_df[[8]], 
                          list_df[[9]]) %>% 
    mutate(Disease = factor(Disease, levels = c('cHL', 'DLBCL', 'FL', 'MZL', 'MCL', 'CLL', 'LPL', 'MM')))
ENTIRE_COHORT %>%  nrow_npatients()

ENTIRE_COHORT_n = ENTIRE_COHORT %>% 
  group_by(patientid) %>% 
  arrange(date_diagnosis) %>% 
  mutate(n = n()) %>% 
  slice(1) %>% 
  ungroup()

# multiple LCs?
ENTIRE_COHORT_n$n %>% table 

#slice to 1st LC
ENTIRE_COHORT.1 = ENTIRE_COHORT %>% 
  group_by(patientid) %>% 
  arrange(date_diagnosis) %>% 
  slice(1) %>% 
  ungroup() 

ENTIRE_COHORT.1 %>% nrow_npatients()

ENTIRE_COHORT.1$IPI %>% table
ENTIRE_COHORT.1$ages %>% table

ENTIRE_COHORT2 = bind_rows(ENTIRE_COHORT.1,
                           ENTIRE_COHORT.1 %>% 
                             mutate(Disease = 'All')) %>% 
  group_by(Disease) %>% 
  mutate(Median = median(age),
         l.iqr = quantile(age, 0.25),
         u.iqr = quantile(age, 0.75)) %>% 
  ungroup() %>% 
  mutate(Disease = factor(Disease, levels = c('All', 'cHL', 'DLBCL', 'FL', 'MZL', 'MCL', 'CLL', 'LPL', 'MM'))) %>% 
  mutate(ages = recode_factor(ages,
                              `[-Inf,40)` = '<40',
                              `[40,50)` = '40-49',  
                              `[50,60)` = '50-59',  
                              `[60,70)` = '60-69',   
                              `[70,80)`= '70-79',
                              `[80, Inf)` = '≥80'))

ENTIRE_COHORT2 %>% nrow_npatients()
Table_1 = utable(Disease ~ Q(age) +ages +IPI, ENTIRE_COHORT2) %>% summary()
Table_1 = Table_1[,-12]
write_csv2(Table_1, '**censored**/table_1.csv')

# Plot Ages and IPIs
list_plot_ipi = list()
list_plot_age2 = list()
for(i in 1:length(list_disease)){
  print(disease_var[i])
  
  list_plot_age2[[i]] =  KM_plot(survfit(Surv(time_dx_death, status) ~ ages, 
                                         ENTIRE_COHORT2 %>% filter(!is.na(IPI),
                                                                   Disease == disease_var[i])),
                             title = disease_var[i],
                             labs = c('<40', '40-49', '50-59', '60-69', '70-79', '≥80'),
                             breaks = 3,
                             xlim = c(0, 15),
                             pval = T)
  
  list_plot_ipi[[i]] = KM_plot(survfit(Surv(time_dx_death, status) ~ IPI, 
                                       ENTIRE_COHORT2 %>% filter(!is.na(IPI),
                                                                 Disease == disease_var[i])),
                                 title = disease_var[i],
                               labs = ENTIRE_COHORT2 %>% 
                                 filter(!is.na(IPI),
                                        Disease == disease_var[i]) %>% 
                                 arrange(IPI) %>% 
                                 pull(IPI) %>% 
                                 unique(),
                               breaks = 3,
                                 xlim = c(0, 15),
                                 pval = T)
}

if(SAVE){
  
  ggsave('**censored**/figure_s2_OS_IPI.png',
         arrange_ggsurvplots(list_plot_ipi, nrow = 3, ncol = 3, print = FALSE),
         height = 16,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S3a-age-distribution.png',
         ggplot() +
           geom_histogram(data = ENTIRE_COHORT2 %>% filter(Disease != 'All'), aes(age, fill = BLOOD[1])) +
           geom_vline(data = ENTIRE_COHORT2, aes(xintercept = Median)) +
           geom_vline(data = ENTIRE_COHORT2, aes(xintercept = l.iqr), lty = 4) +
           geom_vline(data = ENTIRE_COHORT2, aes(xintercept = u.iqr), lty = 4) +
           # geom_vline(xintercept = ENTIRE_COHORT2$Median) +
           facet_wrap(~Disease) +
           theme_classic() +
           theme(legend.position = "none"),
         height = 10,
         width = 10,
         dpi = 300)
  ggsave('**censored**/figure_s3b-age-distribution_all.png',
         ggplot() +
           geom_histogram(data = ENTIRE_COHORT2 %>% filter(Disease == 'All'), aes(age, fill = BLOOD[1])) +
           geom_vline(data = ENTIRE_COHORT2, aes(xintercept = Median)) +
           geom_vline(data = ENTIRE_COHORT2, aes(xintercept = l.iqr), lty = 4) +
           geom_vline(data = ENTIRE_COHORT2, aes(xintercept = u.iqr), lty = 4) +
           # geom_vline(xintercept = ENTIRE_COHORT2$Median) +
           facet_wrap(~Disease) +
           theme_classic() +
           theme(legend.position = "none"),
         height = 10,
         width = 10,
         dpi = 300)
  
  
  ggsave('**censored**/figure_S4_OS_age2.png',
         arrange_ggsurvplots(list_plot_age2, nrow = 3, ncol = 3, print = FALSE),
         height = 16,
         width = 19,
         dpi = 300) # same as Figure 1 but only pts /w available IPI

  ggsave('**censored**/figure_1_OSage.png',
         arrange_ggsurvplots(list_plot, nrow = 3, ncol = 3, print = FALSE),
         height = 16,
         width = 19,
         dpi = 300)
}

# Double loop for subgroups analyses of LC (8; j) x ages (6; i); Figures S5-S12
age_levels = ENTIRE_COHORT2$ages %>% levels
age_labels = c('<40', '40-49', '50-59', '60-69', '70-79', '≥80')
disease_levels = ENTIRE_COHORT2$Disease %>% levels
list_plot_ages = list()
for (j in 1:length(disease_levels)) {
  df = ENTIRE_COHORT2 %>% 
    filter(Disease == disease_levels[j],
           !is.na(IPI)) 
  for(i in 1:length(age_levels)){
    k = (j*6)-6
    print(k+i)
    # tryCatch({c_index = cindex(Surv(time_dx_death, status) ~ IPI, df %>% filter(ages == age_levels[i]))[3]}, error=function(e){})
    list_plot_ages[[i+k]] = KM_plot(survfit(Surv(time_dx_death, status) ~ IPI, df %>% filter(ages == age_levels[i])),
                                    title = paste(disease_levels[j], age_labels[i], 'years'),
                                    labs = df %>% filter(ages == age_levels[i]) %>% 
                                      arrange(IPI) %>% 
                                      pull(IPI) %>% 
                                      unique(),
                                    xlim = c(0, 15),
                                    breaks = 3,
                                    pval = T
                                    )
  }
}

if(SAVE){
  ggsave('**censored**/figure_S5-cHL_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[7]], list_plot_ages[[10]], 
                                  list_plot_ages[[8]], list_plot_ages[[11]],
                                  list_plot_ages[[9]], list_plot_ages[[12]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S6-DLBCL_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[13]], list_plot_ages[[16]], 
                                  list_plot_ages[[14]], list_plot_ages[[17]],
                                  list_plot_ages[[15]], list_plot_ages[[18]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S7-FL_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[19]], list_plot_ages[[22]], 
                                  list_plot_ages[[20]], list_plot_ages[[23]],
                                  list_plot_ages[[21]], list_plot_ages[[24]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S8-MZL_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[25]], list_plot_ages[[28]], 
                                  list_plot_ages[[26]], list_plot_ages[[29]],
                                  list_plot_ages[[27]], list_plot_ages[[30]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S9-MCL_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[31]], list_plot_ages[[34]], 
                                  list_plot_ages[[32]], list_plot_ages[[35]],
                                  list_plot_ages[[33]], list_plot_ages[[36]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S10-CLL_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[37]], list_plot_ages[[40]], 
                                  list_plot_ages[[38]], list_plot_ages[[41]],
                                  list_plot_ages[[39]], list_plot_ages[[42]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S11-LPL_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[43]], list_plot_ages[[46]], 
                                  list_plot_ages[[44]], list_plot_ages[[47]],
                                  list_plot_ages[[45]], list_plot_ages[[48]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
  ggsave('**censored**/figure_S12-MM_ages.png',
         arrange_ggsurvplots(list(list_plot_ages[[49]], list_plot_ages[[52]], 
                                  list_plot_ages[[50]], list_plot_ages[[53]],
                                  list_plot_ages[[51]], list_plot_ages[[54]]), 
                             nrow = 2, ncol = 3, print = FALSE),
         height = 10,
         width = 19,
         dpi = 300)
  
}

# Calculates 5-y OS estimates by Ages and IPIs
# reference: ggplot code: https://r-charts.com/correlation/heat-map-ggplot2/?utm_content=cmp-true

data_surv = tibble()
list_surv = list()
for(i in 1:length(disease_var)){
  print(disease_var[i])
  list_surv = summary(survfit(Surv(time_dx_death, status) ~ ages+IPI, ENTIRE_COHORT2 %>% 
                                filter(!is.na(IPI),
                                       Disease == disease_var[i])),
                      times = 5, extend = TRUE)
  data_surv = tibble(label = str_trim(list_surv$table %>% rownames(), 'both'),
                     os5 = round(list_surv$surv*100, 0),
                     n = list_surv$n,
                     n.risk = list_surv$n.risk) %>% 
    mutate(disease = disease_var[i], 
           x = str_trim(str_split_fixed(label, '\\,', 2)[,2], 'both'),
           x = str_split_fixed(x, '=', 2)[,2],
           y = str_trim(str_split_fixed(label, '\\,', 2)[,1],  'both'),
           y = str_split_fixed(y, '=', 2)[,2]) %>% 
    bind_rows(data_surv)
}

#### Save Figure 2 ####

if(SAVE){
  ggsave('**censored**/figure_2_5y_risk_IPIxAges.png',
         ggplot(data_surv %>% 
                  filter(n >= 10,
                         disease != 'All' | disease == 'All' & x != 'Very high') %>% 
                  mutate(disease = recode_factor(disease, 
                                                 All = 'All', cHL = 'cHL: IPS', DLBCL = 'DLBCL: R-IPI', 
                                                 FL = 'FL: FLIPI2', MZL = 'MZL: MALT-IPI', MCL = 'MCL: MIPI', 
                                                 CLL = 'CLL: CLL-IPI', LPL = 'LPL: IPSSWM', MM = 'MM: R-ISS'),
                         x = factor(x, c('Low', 'Intermediate', 'High', 'Very high'))), 
                aes(x = x, y = y, fill = os5)) +
           geom_tile(color = "white",
                     lwd = 1.5,
                     linetype = 1) +
           geom_text(aes(label = os5), color = "black", size = 4) +
           scale_fill_gradientn('5-y OS', colors = hcl.colors(20, "RdYlGn"))+
           labs(x = 'Subtype-specific IPI', y = 'Ages',
                title = '5-year OS estimates') +
           facet_wrap(~disease, scales='free') +
           theme_classic() +
           theme(legend.position = 'none'),
         height = 8,
         width = 10,
         dpi = 300)
}
