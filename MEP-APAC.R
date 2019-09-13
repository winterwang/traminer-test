
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
# 1      5307 

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
