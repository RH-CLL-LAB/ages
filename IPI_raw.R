## This script gathers lab values and calculates IPIs/IPS/ISS in 8 LCs: cHL, DLBCL, FL, MZL, MCL, CLL, LPL, and MM.
# created: "2024-04-26"
# author: christian brieghel

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')

#### RKKP IPI ####
load_dataset(c(RKKP_DATASETS, 'LAB_IGHVIMGT'))
ALL_RKKP = c(RKKP_CLL$patientid, RKKP_DaMyDa$patientid, RKKP_LYFO$patientid) %>% unique()
ALL_RKKP %>% n_distinct()

## load lab values using common NPU codes
load_npu_common() # loads NPU vectors
BIOCHEMISTRY = load_biochemistry(c(NPU.ALB, NPU.B2M, NPU.LDH, NPU.IGM, NPU.HGB, NPU.LEU, NPU.LYM, NPU.LDH)) %>% 
  filter(patientid %in% ALL_RKKP) %>% 
  clean_lab_values()

lab_helper = function(data, lab_name){
  name = lab_name
  npu_code = get(paste0('NPU.', {{lab_name}}))
  data %>% 
  filter(NPU %in% npu_code) %>% 
    transmute(patientid, 
              'date_{name}' := as_date(samplingdate), 
              '{name}_value' := value2)
}

# gather lab values
B2M = BIOCHEMISTRY %>% 
  lab_helper(lab_name = 'B2M') %>% 
  mutate(B2M_lab = factor(ifelse(B2M_value < 4.0, '<4.0 mg/L', '>4.0 mg/L')))

ALB = BIOCHEMISTRY %>% 
  lab_helper(lab_name = 'ALB')

IGM = BIOCHEMISTRY %>% 
  lab_helper(lab_name = 'IGM')

HGB = BIOCHEMISTRY %>% 
  lab_helper(lab_name = 'HGB')

WBC =  BIOCHEMISTRY %>% 
  lab_helper(lab_name = 'LEU')

ALC = BIOCHEMISTRY %>% 
  lab_helper(lab_name = 'LYM')

LDH =  BIOCHEMISTRY %>% 
  lab_helper(lab_name = 'LDH')

#### CLL ####

RKKP_CLL_clean = RKKP_CLL %>% 
  clean_RKKP_CLL() %>% 
  left_join(LAB_IGHVIMGT %>% transmute(patientid, IGHVIMGT = factor(ighv)), 'patientid') %>% 
  mutate(IGHV = if_else(is.na(IGHV), IGHVIMGT, IGHV)) %>% 
  left_join(B2M, 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_B2M, value = B2M_value, interval_days = c(-30, 30), name = 'B2M') %>% 
  mutate(B2M = if_else(is.na(B2M), B2M_lab, B2M)) %>% 
  CLL_IPI() %>%
  mutate(CLL_IPI = if_else(IPI_score_minus_B2M >= 7, 'Very high', CLL_IPI),
         CLL_IPI = if_else(IPI_score_minus_IGHV >= 7, 'Very high', CLL_IPI),
         CLL_IPI = if_else(IPI_score_minus_B2M == 4, 'High', CLL_IPI),
         CLL_IPI = if_else(IPI_score_minus_IGHV == 4, 'High', CLL_IPI),
         CLL_IPI = factor(CLL_IPI, levels = c('Low', 'Intermediate', 'High', 'Very high')))

RKKP_CLL_clean$CLL_IPI %>% table(exclude = NULL)
RKKP_CLL_clean %>%  nrow_npatients()

#### DAMYDA ####
RKKP_DaMyDa_clean = RKKP_DaMyDa %>% 
  clean_RKKP_DAMYDA() %>% 
  left_join(B2M, 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_B2M, value = B2M_value,  interval_days = c(-30, 30)) %>% 
  mutate(B2M = if_else(is.na(B2M), B2M_value, B2M))  %>% 
  
  left_join(ALB, 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_ALB, value = ALB_value, interval_days = c(-30, 30), name = 'ALB_val') %>% 
  mutate(ALB = if_else(is.na(ALB), ALB_value, ALB)) %>% 
  
  left_join(LDH, 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_LDH, value = LDH_value, interval_days = c(-30, 30), name = 'LDH_val') %>% 
  mutate(LDH = if_else(is.na(LDH), LDH_value, LDH)) %>% 
  
  mutate(ISS_ext = ifelse(ALB >= 35 & B2M <3.5, 1, 2),
         ISS_ext = ifelse(B2M >= 5.5, 3, ISS_ext),
         ISS = ifelse(is.na(ISS), ISS_ext, ISS)) %>% 
  mutate(LDH_high = ifelse(age < 70 & LDH > 205, 'Yes', NA),
         LDH_high = ifelse(age < 70 & LDH <= 205, 'No', LDH_high),
         LDH_high = ifelse(age >= 70 & LDH > 255, 'Yes', LDH_high),
         LDH_high = ifelse(age >= 70 & LDH <= 255, 'No', LDH_high),
         FISH.score = ifelse(FISH_t4_14 =='Yes' | FISH_t14_16 =='Yes' | FISH_del17p =='Yes', 'Yes', 'No')) %>%
  mutate(RISS_Addon = ifelse(LDH_high =='Yes' | FISH.score =='Yes', 'Yes', 'No')) %>%
  mutate(RISS = ifelse(ISS == 3 & RISS_Addon == 'Yes', 3, NA),
         RISS = ifelse(ISS == 1 & RISS_Addon == "No", 1, RISS),
         RISS = ifelse(ISS == 3 & RISS_Addon == "No", 2, RISS),
         RISS = ifelse(ISS == 1 & RISS_Addon == "Yes", 2, RISS),
         RISS = ifelse(ISS == 2, 2, RISS)) 

RKKP_DaMyDa_clean$RISS %>% table(exclude = NULL)
RKKP_DaMyDa_clean %>%  nrow_npatients()

#### LYFO ####
RKKP_LYFO_clean = RKKP_LYFO %>% 
  clean_RKKP_LYFO() %>% 
  left_join(HGB, by = 'patientid', relationship = "many-to-many") %>%
  slice_closest_value(date_diagnosis, date_HGB, value = HGB_value, interval_days = c(-30,30), name = 'HGB') %>% 
  left_join(WBC, by = 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_LEU, value = LEU_value, interval_days = c(-30,30), name = 'WBC') %>% 
  left_join(ALC, by = 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_LYM, value = LYM_value, interval_days = c(-30,30), name = 'ALC') %>% 
  left_join(LDH,  by = 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_LDH, value = LDH_value, interval_days = c(-30,30), name = 'LDH') %>% 
  left_join(IGM,  by = 'patientid', relationship = "many-to-many") %>% 
  slice_closest_value(date_diagnosis, date_IGM, value = IGM_value, interval_days = c(-30,30), name = 'IGM')  %>%
  mutate(HB_diagnosis = ifelse(is.na(HB_diagnosis), valueHGB, HB_diagnosis),
         WBC_diagnosis = ifelse(is.na(WBC_diagnosis), valueWBC, WBC_diagnosis),
         ALC_diagnosis = ifelse(is.na(ALC_diagnosis), valueALC, ALC_diagnosis),
         LDH_diagnosis = ifelse(is.na(LDH_diagnosis), valueLDH, LDH_diagnosis),
         IgM_gL_diagnosis = ifelse(is.na(IgM_gL_diagnosis), valueIGM, IgM_gL_diagnosis))

#### Gather all IPIs  ####
## DALY-CARE methods paper CLEP 2025: Table S8
ALL_IPI = bind_rows(RKKP_CLL_clean %>% 
                      transmute(patientid, 
                                sex,
                                PS,
                                IPI = CLL_IPI,
                                Disease = 'CLL'),
                    RKKP_DaMyDa_clean %>% 
                      filter(subtype == 'DC900') %>% 
                      transmute(patientid, 
                                sex,
                                PS,
                                IPI = recode_factor(RISS,
                                                    `1` = 'Low',
                                                    `2` = 'Intermediate',
                                                    `3` = 'High'),
                                Disease = 'MM'),
                    RKKP_LYFO_clean %>% 
                      filter(subtype  == 'DLBCL') %>% 
                      transmute(patientid, 
                                sex,
                                PS = PS_diagnosis,
                                IPI = RIPI_diagnosis,
                                Disease = 'DLBCL'),
                    RKKP_LYFO_clean %>% 
                      filter(subtype  == 'FL') %>% 
                      transmute(patientid, 
                                sex,
                                PS = PS_diagnosis,
                                IPI = FLIPI2_diagnosis,
                                Disease = 'FL'),
                    RKKP_LYFO_clean %>% 
                      filter(subtype  == 'MCL') %>% 
                      MIPI(subtype = subtype) %>% 
                      transmute(patientid, 
                                sex,
                                PS = PS_diagnosis,
                                IPI = MIPI,
                                Disease = 'MCL'),
                    RKKP_LYFO_clean %>% 
                      filter(subtype  == 'WM') %>% 
                      IPSSWM() %>%
                      transmute(patientid, 
                                sex,
                                PS = PS_diagnosis,
                                IPI = IPSSWM,
                                Disease = 'LPL'),
                    RKKP_LYFO_clean %>% 
                      filter(subtype == 'MZL') %>% 
                      MALT_IPI() %>% 
                      transmute(patientid, 
                                sex,
                                PS = PS_diagnosis,
                                IPI = MALT_IPI,
                                Disease = 'MZL'),
                    RKKP_LYFO_clean %>% 
                      filter(subtype =='cHL') %>% 
                      IPS() %>%
                      mutate(IPS = cut(IPS.score, c(-Inf, 2, Inf), labels = c('Low', 'High'))) %>% 
                      transmute(patientid, 
                                sex,
                                PS = PS_diagnosis,
                                IPI = cut(IPS.score, c(-Inf, 2, Inf), labels = c('Low', 'High')),
                                Disease = 'cHL'))  %>% 
  mutate(IPI = factor(IPI, levels =  c('Low', 'Intermediate', 'High', 'Very high')))

ALL_IPI$IPI %>% table(exclude = NULL)
ALL_IPI %>% nrow_npatients()

setwd('**censored**/tables/')
if(SAVE == TRUE){
  write_csv2(ALL_IPI, paste0(getwd(), '/IPI.csv'))
}

# End of script: Go to ages.R