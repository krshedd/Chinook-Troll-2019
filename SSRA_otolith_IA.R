library(tidyverse)
source("~/../R/Functions.GCL.R")

# read in genotypes for SSRA otolith project
load_sillys(path = "../Genotypes/strata_postQA/", sillyvec = c("SpringRet2_10145_2019",
                                                               "SpringRet2_10146_2019"))
load_objects(path = "../Objects/")

# do we have dna_specimen_no
names(SpringRet2_10145_2019.gcl$attributes)

SpringRet2_10145_2019.gcl$attributes$`Dna Specimen No`

# get the fish ids she wants
lauras_fish <- as.numeric(readClipboard())

# how many in there?
table(lauras_fish %in% SpringRet2_10145_2019.gcl$attributes$`Dna Specimen No`)
table(lauras_fish %in% SpringRet2_10146_2019.gcl$attributes$`Dna Specimen No`)

# how many in asl?
asl <- read_csv("../ASL Data/spring_troll_samples_ASL.csv")
table(lauras_fish %in% asl$`Dna Specimen No`)

# breakout asl by sw and area
asl %>% 
  filter(`Dna Specimen No` %in% lauras_fish) %>% 
  count(`Stat Area`, `Stat Week`)

# get fish_ids
lauras_10145 <- SpringRet2_10145_2019.gcl$attributes$FK_FISH_ID[SpringRet2_10145_2019.gcl$attributes$`Dna Specimen No` %in% lauras_fish]
lauras_10146 <- SpringRet2_10146_2019.gcl$attributes$FK_FISH_ID[SpringRet2_10146_2019.gcl$attributes$`Dna Specimen No` %in% lauras_fish]

# pool
PoolCollections.GCL(
  collections = c("SpringRet2_10145_2019", "SpringRet2_10146_2019"),
  loci = GAPSLoci_reordered,
  IDs = list(
    "SpringRet2_10145_2019" = lauras_10145,
    "SpringRet2_10146_2019" = lauras_10146
  ),
  newname = "SSRA_otos"
)

spring_atts <- bind_rows(SpringRet2_10145_2019.gcl$attributes, SpringRet2_10146_2019.gcl$attributes) %>% 
  select(SillySource, `Dna Specimen No`)

str(SSRA_otos.gcl)

# create ONCOR files
gcl2Genepop.GCL(sillyvec = "SSRA_otos", path = "../ONCOR/Mixture/SSRA_otos.gen", loci = GAPSLoci, VialNums = TRUE, usat = TRUE)

# read output
rawSSRA <- scan(file = "../ONCOR/Output/SSRA_otos_33RG_IA.txt", what = '', sep = '\n', blank.lines.skip = FALSE)
skip <- grep(pattern = "PROBABILITY OF EACH INDIVIDUAL IN THE MIXTURE BELONGING TO EACH REPORTING GROUP", x = rawSSRA) + 1
SSRA_dat <- read_delim(file = "../ONCOR/Output/SSRA_otos_33RG_IA.txt", delim = "\t", skip = skip, trim_ws = TRUE) %>% 
  dplyr::rename(FishID = X1)


# 1st probability RG
SSRA_dat_1 <- SSRA_dat %>% 
  gather(RG, prob, -FishID) %>% 
  group_by(FishID) %>% 
  slice(which.max(prob)) %>% 
  arrange(FishID) %>% 
  ungroup()

# 2nd probability RG
SSRA_dat_2 <- SSRA_dat %>% 
  gather(RG, prob, -FishID) %>% 
  group_by(FishID) %>% 
  filter(rank(desc(prob)) == 2) %>% 
  arrange(FishID) %>% 
  ungroup()


# Join 1st and 2nd probability RG with SSRA attributes
SSRA_dat_tdy <- left_join(x = SSRA_dat_1, y = SSRA_dat_2, by = "FishID") %>% 
  dplyr::rename(group_1 = RG.x, prob_1 = prob.x, group_2 = RG.y, prob_2 = prob.y) %>% 
  left_join(spring_atts, by = c("FishID" = "SillySource")) %>% 
  select(`Dna Specimen No`, everything(), -FishID)


SSRA_dat_tdy %>% 
  count(group_1)

SSRA_dat_tdy %>% 
  arrange(prob_1)

write_csv(x = SSRA_dat_tdy, path = "../SSRA_otolith_project_potential_wild_fish.csv")
