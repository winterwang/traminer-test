
# run Traminer with the example data --------------------------------------
library(tidyverse)

path <- c("C:/Users/C.Wang/Google ドライブ/Chao/Work with Thomas Laurent/2019-09-13/MEP-APAC-Pain_20181019/MEP-APAC-Pain_20181019")


Sippei <- read.csv(paste(path, "/04_med_sippei.csv", sep = ""))
DPC <- read.csv(paste(path, "/08_dpc_sippei.csv", sep = ""))


MED <- Sippei %>% 
  select(c(4, 3, 5, 9, 7))

names(MED) <- c("ID", "receDate", "ICD10", "Flag_uta", "ICDcoded")
MED <- MED %>% 
  filter(Flag_uta == 0 & ICDcoded != 999 & ICD10 != "") %>% 
  mutate(ICD10 = as.character(ICD10))
head(MED)

ID_MED <- as.character(unique(MED$ID))



DPC0 <- DPC %>% 
  select(c(4, 3, 5, 9, 7))
names(DPC0) <- c("ID", "receDate", "ICD10", "Flag_uta", "ICDcoded")


DPC0 <- DPC0 %>% 
  # select(c(4, 3, 5, 9, 7)) %>% 
  filter(Flag_uta == 0 & ICDcoded != 999 & ICD10 != "") %>% 
  mutate(ICD10 = as.character(ICD10))

head(DPC0)

ID_DPC <- as.character(unique(DPC0$ID))

table(ID_DPC %in% ID_MED)
 
# FALSE  TRUE 
#     1  5202  

MED[MED$ID == ID_DPC[1], ] # so only 1 patient in DPC data is not recored in MED data


# lets combine the two data sets ------------------------------------------

icd10 <- rbind(DPC0, MED)

icd10 <- icd10 %>% 
  arrange(ID, receDate)
head(icd10)



# permute disease for the same month ------------------------------------------

set.seed(1234)



icd10_sample <- icd10 %>% 
  mutate(From = substr(receDate, 1, 6)) %>% 
  group_by(ID, From) %>% 
  sample_n(1) 

icd10_sample <- icd10_sample %>% 
  group_by(ID) %>% 
  distinct(ICD10, .keep_all = TRUE)


icd10_sample$receDate <-  as.Date(as.character(icd10_sample$receDate), "%Y%m%d")

icd10_sample$To <- icd10_sample$receDate - 1


icd10_sam_ana <- icd10_sample %>% 
    mutate(To = lead(To)) %>% 
  select(c(1, 2, 3, 7))

icd10_sam_ana <- icd10_sam_ana %>% 
   select(c(1, 2, 4, 3)) %>% 
     group_by(ID) %>% 
     mutate(Index = row_number())


# the maximum of receDate 

max(icd10_sam_ana$receDate)
# [1] "2017-12-01"

# recode the final "To" date for each patient into 2017-12-31

icd10_sam_ana$To[is.na(icd10_sam_ana$To)] <- as.Date("2017-12-31")




# do traminer -------------------------------------------------------------

library(TraMineR)


LA.labels <- seqstatl(icd10_sam_ana$ICD10)

LA.states <- 1:length(LA.labels)


icd10_sam_ana <- icd10_sam_ana %>% 
  mutate(from = as.numeric(receDate), 
         to = as.numeric(To))
icd10_sam_ana <- as.data.frame(icd10_sam_ana)

icd.seq <- seqdef(icd10_sam_ana, var = c("ID", "from", "to", "ICD10"), informat = "SPELL", 
       states = LA.states, labels = LA.labels, process = FALSE) #too heavy
# sample 100 patients to do the sequence 

id <- icd10_sam_ana$ID
id_test <- sample(id, 200)

icd10_test <- icd10_sam_ana %>% 
  filter(ID %in% id_test)


LA.labels <- seqstatl(icd10_test$ICD10)

LA.states <- 1:length(LA.labels)


icd.seq <- seqdef(icd10_test, var = c("ID", "from", "to", "ICD10"), informat = "SPELL", 
                  states = LA.states, labels = LA.labels, process = FALSE) 


seqiplot(icd.seq, title = "Index plot (first 10 sequences)",
         withlegend = FALSE)


seqtab(icd.seq, tlim = 1:4)

# most patients' sequence are unique and therefore maybe the list of ICD10 should be defined. 


# lets try to focus on ICD10 code start with I  ---------------------------

# ie. cardiovascular diseases 

ICD_i_med <- MED %>% 
  filter(str_detect(ICD10, 'I'))
ICD_i_dpc <- DPC0 %>% 
  filter(str_detect(ICD10, 'I'))


# lets combine the two data sets ------------------------------------------

ICD_i <- rbind(ICD_i_med, ICD_i_dpc)

ICD_i <- ICD_i %>% 
  arrange(ID, receDate)
head(ICD_i)





# permute disease for the same month ------------------------------------------

set.seed(1234)



icd_i_sample <- ICD_i %>% 
  mutate(From = substr(receDate, 1, 6)) %>% 
  group_by(ID, From) %>% 
  sample_n(1) 

icd_i_sam_dist <- icd_i_sample %>% 
  group_by(ID) %>% 
  distinct(ICD10, .keep_all = TRUE)

# icd_i_sample %>%
#   group_by(ID) %>%
#   summarise(first=head(ICD10, 1), count=n_distinct(ICD10))


icd_i_sam_dist$receDate <-  as.Date(as.character(icd_i_sam_dist$receDate), "%Y%m%d")

icd_i_sam_dist$To <- icd_i_sam_dist$receDate - 1


icdi_distana <- icd_i_sam_dist %>% 
  mutate(To = lead(To)) %>% 
  select(c(1, 2, 3, 7))

icdi_distana <- icdi_distana %>% 
  select(c(1, 2, 4, 3)) %>% 
  group_by(ID) %>% 
  mutate(Index = row_number())


# the maximum and of receDate 

max(icdi_distana$receDate)
# [1] "2017-12-01"

min(icdi_distana$receDate)

# recode the final "To" date for each patient into 2017-12-31

icdi_distana$To[is.na(icdi_distana$To)] <- as.Date("2017-12-31")



# recode ID  --------------------------------------------------------------

# ie create group id as a new variable

icdi_distanaPID <- icdi_distana %>% 
  group_by(ID) %>% 
  mutate(PID = group_indices()) %>% 
  ungroup() %>% 
  select(c(6, 2, 3, 4, 5))

icdi_distanaPID <- icdi_distanaPID %>% 
  mutate(timediff = as.numeric(To - receDate)/365.25) %>% 
  group_by(PID) %>% 
  mutate(from = row_number()) %>% 
  ungroup() %>% 
  # mutate(to = round(timediff)) %>% 
  # mutate(to = if_else(to <= from, from + 1, to)) 
  mutate(to = from)
  
  
# do traminer -------------------------------------------------------------

library(TraMineR)

icdi_distanaPID <- icdi_distanaPID %>% 
  mutate(ICD10_2 = substr(ICD10, 1, 3)) 

icdi_distanaPID <- icdi_distanaPID %>% 
  mutate(CVD = if_else(grepl("I0[0-2]", 
                             icdi_distanaPID$ICD10_2), "Acute rheumatic fever", 
                if_else(grepl("I0[5-9]", 
                              icdi_distanaPID$ICD10_2), "Chronic rheumatic heart diseases", 
                 if_else(grepl("I1[0-5]", 
                               icdi_distanaPID$ICD10_2), "Hypertensive diseases", 
                  if_else(grepl("I2[0-5]", 
                                icdi_distanaPID$ICD10_2), "Ischaemic heart diseases", 
                  if_else(grepl("I2[6-8]", 
                                icdi_distanaPID$ICD10_2), "Pulmonary heart disease", 
                  if_else(grepl("I[3-5][0-9]", 
                                icdi_distanaPID$ICD10_2), "Other forms of heart disease", 
                  if_else(grepl("I6[0-9]", 
                                icdi_distanaPID$ICD10_2), "Cerebrovascular diseases", 
                  if_else(grepl("I7[0-9]", 
                                icdi_distanaPID$ICD10_2), "Diseases of arteries", 
                  if_else(grepl("I8[0-9]", 
                                icdi_distanaPID$ICD10_2), "Diseases of veins", 
                  if_else(grepl("I9[5-9]", 
                                icdi_distanaPID$ICD10_2), "Other and unspecified", "NA")))))))))))
  
  
LA.labels <- seqstatl(icdi_distanaPID$CVD)



# LA.states <- 1:length(LA.labels)
LA.states <- c("ARF", "CD", "CRHD", "DA", "DV", "HD", "IHD", "OU", "OHD", "PHD")

# icdi_distanaPID <- icdi_distanaPID %>% 
#   mutate(from = from, 
#          to = to)

icdi_distanaPID <- as.data.frame(icdi_distanaPID)

icd.seq <- seqdef(icdi_distanaPID, var = c("PID", "from", "to", "CVD"), informat = "SPELL", 
                  states = LA.states, labels = LA.labels, process = FALSE)


print(icd.seq[1:5, ], format = "SPS")

seqiplot(icd.seq, with.legend = "none")

seqtab(icd.seq, idxs = 1:4)

icd.trate <- seqtrate(icd.seq)
round(icd.trate, 2)


icd.om <- seqdist(icd.seq, method = "OM", indel = 1, sm = "TRATE")
library(cluster)
clusterward <- agnes(icd.om, diss = TRUE, method = "ward")
