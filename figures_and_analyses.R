# Alberto Rovellini
# University of Washington - School of Aquatic and Fishery Sciences
# 01/21/2025
# This script takes the output of the Atlantis runs for the evaluation of the Optimum Yield in the Gulf of Alaska and creates the analyses and figures for the associated manuscript
# The code is organized by Figure number, following the order of the main text. Supplementary figures are created throughout with their main-text counterparts (if any), and are ordered more loosely

# Set up the environment and read data ------------------------------------------------

library(tidyverse)
library(here)
library(tidyr)
library(readxl)
library(ggh4x)
library(viridis)
library(tidync)
library(ncdf4)
library(rbgm)
library(RColorBrewer)
library(patchwork)

# general settings
t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S') # set the clock to date plots
select <- dplyr::select
burnin <- 30 # years of burn-in
maxmult <- 4 # this is the maximum explored multiplier on FMSY

# set paths to Atlantis results folders
batch_ss <- "results/ss/processed/" # there are the result from the single-species runs (Step 1)
batch_ms_flat <- "results/ms/flat_results/" # these are the biomage, catch, and mort files
batch_ms_nc <- "results/ms/nc_results/" # these are the full out.nc files

# read in Groups.csv file
grps <- read.csv('data/GOA_Groups.csv')

# read maturity and selectivity information
selex <- read.csv("data/age_at_selex_new.csv", header = T) # age at selectivity
fspb <- read.csv("data/fspb.csv", header = T) # proportion of spawning biomass per age class
# reshape fspb
fspb <- fspb %>%
  pivot_longer(-Code, names_to = "Age", values_to = "fspb") %>%
  mutate(Age = gsub("X","",Age))

# read in lookup tables and f35 proxy values
oy_key <- read.csv("data/oy_key.csv")
f_lookup <- read.csv("data/f_lookup_OY_SS.csv")

# list Tier 3 stocks 
t3_fg <- f_lookup %>% pull(species) %>% unique() %>% sort()
t3_names <- grps %>% filter(Code %in% t3_fg) %>% pull(Name) # names for nc files pulling



# Figure 1. GOA harvest specifications ------------------------------------
# this spreadsheet is available for download from AKFIN Answers
specs <- read_excel("data/GOA_harvest specs_1986-2024.xlsx", 
                    sheet = 1,
                    na = "n/a",
                    n_max = 131)

# need to clean the data set
# drop asterisks and commas and turn to numeric
for(col in names(specs)){
  specs[[col]] <- gsub("\\*","", specs[[col]])
  specs[[col]] <- gsub(",","", specs[[col]])
}

# now handle column names
# pad years
colnames(specs) <- c("", "", rep(2024:1986, each = 3))
# collapse column names with the first row for pivot later
new_row <- rep(NA, ncol(specs))
for(i in 1:ncol(specs)){
  new_row[i] <- paste(specs[1,i], names(specs)[i], sep = "_")
  new_row[1:2] <- gsub("_","",new_row[1:2])
}

# set new colnames
colnames(specs) <- new_row
# drop old row 1
specs <- specs[-1,]

# now pad the species column
for(i in 1:nrow(specs)){
  if(is.na(specs[i,1])){
    specs[i,1] <- specs[i-1,1]
  }
}

# pivot longer, split spec and year, add Tier and Atlantis functional group
specs_long <- specs %>%
  pivot_longer(-c(Species, Area), names_to = "Spec_Year", values_to = "mt") %>%
  separate(Spec_Year, into = c("Spec", "Year"), sep = "_") %>%
  filter(Year > 1990, Area %in% c("Total","Total (GW)", "GW")) %>%
  mutate(mt = as.numeric(mt)) %>%
  drop_na() %>%
  mutate(Area = "GOA")

# add in total groundfish catch reconstructions
# Data is from "Catch Data" tab in AKFIN Answers using the following tags:
# * Year: 1991-2024
# * FMP Area: GOA
# * FMP Subarea: --Select Value--
# * Gear: --Select Value--
# * Species Group: --Select Value--
# * Choose a Report: "Detail with Processor/Vessel Characteristics"

catch_data <- read.csv("data/Groundfish Total Catch.csv", fileEncoding = 'UTF-8-BOM')

# there are a lot of non TAC species reported here, as well as by-catch species, species that are in the FMP but are not groundfish, etc.
# For the purpose of comparing to ABC/TAC plots, we will filter only the species that have a TAC in the harvest specifications
# Map species in the catch to species in the harvest specification data set
tac_key <- read.csv("data/tac_catch_key.csv", header = T)

# process the catch data so that it can be mapped to the harvest specification data
catch_data_short <- catch_data %>%
  select(Year, Species.Group.Name, Catch..mt.) %>%
  left_join(tac_key, by = c("Species.Group.Name" = "Catch_sp")) %>%
  group_by(Year, TAC_sp) %>%
  summarise(mt = sum(Catch..mt., na.rm = T)) %>%
  mutate(Area = "GOA", Spec = "Catch") %>%
  select(TAC_sp, Area, Spec, Year, mt) %>%
  rename(Species = TAC_sp)

specs_long <- specs_long %>% rbind(catch_data_short) %>% drop_na()

# order factors
specs_long$Spec <- factor(specs_long$Spec, levels = c("OFL", "ABC", "TAC", "Catch"))

# format stock names for the plot
specs_long$Species <- gsub("/","\n",specs_long$Species)
specs_long$Species <- gsub(" and "," and\n",specs_long$Species)
specs_long$Species <- gsub(" \\(","\n\\(",specs_long$Species)

# remove mollusks for this plot - small catch 
specs_long <- specs_long %>%
  filter(!Species %in% c("Octopus","Squid"))

# set some colors - this plot has lots of species
# this combination works OK in that the key stocks are readable enough, but it is not greyscale- nor colorblind-friendly
colors <- c(viridis(11)[2:10], inferno(11)[2:10], cividis(10)[2:9])

# make a bar chart
harvest_specs_fig <- specs_long %>%
  #filter(Tier == 3) %>%
  group_by(Year, Spec, Species) %>%
  summarise(mt = sum(mt, na.rm = T)) %>%
  ggplot(aes(x = Year, y = mt/1000, fill = Species))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = 800, linetype = "dashed", color = "red")+
  theme_bw()+
  scale_x_discrete(breaks = seq(1992,2024,2))+
  labs(x = "", y = "1000 mt", fill = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  # theme(legend.position="bottom",
  #       legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(ncol = 1))+
  facet_wrap(~Spec, nrow = 2)
harvest_specs_fig

ggsave("results/figures/FIGURE_1.jpeg", harvest_specs_fig, width = 8, height = 6.5, dpi = 600)

# mean recent catch for text
tt <- specs_long %>%
  filter(Spec == "Catch") %>%
  group_by(Year) %>%
  summarise(Catch = sum(mt)) 

# Figure 2. Single-species biomass and catch  --------------------------------------------------------------

# list files from the single species runs (step 1)
# these are produced by the script ss_processing.R
f_files <- list.files(batch_ss, full.names = T)

# create empty list to fill with data frame for the yield curve
f_df_ls <- list()

for(i in 1:length(f_files)){
  
  this_f_files <- f_files[[i]]
  
  # read all csv files
  f_ls <- list()
  for(j in 1:length(this_f_files)){
    this_file <- this_f_files[j]
    f_ls[[j]] <- read.csv(this_file)
  }
  
  # bind into a data frame
  f_df <- f_ls %>% bind_rows() 
  
  # clean up and format
  f_df <- f_df %>%
    pivot_longer(
      cols = c(mean_biom, mean_catch, biom_cv, catch_cv),
      names_to = "temp",
      values_to = "value"
    ) %>%
    mutate(
      Var = case_when(
        temp == "mean_biom" ~ "Biomass",
        temp == "mean_catch" ~ "Catch",
        temp == "biom_cv" ~ "Biomass",
        temp == "catch_cv" ~ "Catch"
      ),
      type = ifelse(grepl("mean", temp), "Mean", "CV")
    ) %>%
    select(-temp) %>%
    pivot_wider(
      names_from = type,
      values_from = value
    ) %>%
    select(Code, f, fidx, Var, Mean, CV) %>%
    left_join(grps %>% select(Code, LongName), by = 'Code')
  
  f_df_ls[[i]] <- f_df
  
}

f_df <- bind_rows(f_df_ls)

# produce a dataset of 35% B0, to be used to plot horizontal lines that will intersect the yield curve
b35 <- f_df %>%
  group_by(LongName) %>%
  slice_min(f) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(b35 = ifelse(Var== 'Biomass', Mean * 0.35, NA)) %>%
  ungroup() %>%
  select(LongName, Code, Var, b35)

# read in MSY information (from FMP)
# this classification is current as of 2023, note that stocks may change Tier over time
tier3 <- read_xlsx('data/msy.xlsx', sheet = 1, range = 'A3:J19') %>%
  select(Stock, FOFL) %>%
  set_names(c('Stock', 'FMSY'))

tier4_5 <- read_xlsx('data/msy.xlsx', sheet = 2, range = 'A3:I10') %>%
  select(`Stock/Stock complex`, `M or FMSY`)%>%
  set_names(c('Stock', 'FMSY'))

tier_3_4_5 <- rbind(tier3, tier4_5)

# make key
tier_3_4_5 <- tier_3_4_5 %>%
  mutate(Code = c('POL','COD','SBF','FFS','FFS','FFS','FFS','FFD',
                  'REX','REX','ATF','FHS','POP','RFS','RFS','RFP',
                  'FFS','RFD','RFD','RFD','RFD','THO','DOG')) %>%
  group_by(Code) %>%
  summarise(FMSY = mean(FMSY))

all_f <- tier_3_4_5

# find groups to plot
grp_to_plot <- unique(f_df$Code)

# bind FMSY information
fmsy <- data.frame('Code' = grp_to_plot) %>%
  left_join(all_f) %>%
  left_join(grps %>% select(Code, LongName))

# add halibut (M from IPHC assessment)
fmsy[fmsy$Code=='HAL',]$FMSY <- 0.2 # this is M

# get f that returned the highest yield, and level of depletion for that F
sp <- unique(f_df$LongName)

# create a data frame with the reference point estimates from Atlantis
# FMSY here is the F that returned the highest catch at equilibrium
# BFMSY is the corresponding B
# depletion is the fraction of the unexploited spawning biomass
atlantis_fmsy_ls <- list()

for(i in 1:length(sp)){
  
  this_f_df <- f_df %>% filter(LongName == sp[i])
  
  atlantis_fmsy <- this_f_df %>% filter(Var == 'Catch') %>%
    slice_max(Mean) %>%
    pull(f)
  
  b0 <- this_f_df %>%
    filter(Var == 'Biomass') %>%
    slice_min(f) %>%
    pull(Mean)
  
  b_fmsy <- this_f_df %>%
    filter(f == atlantis_fmsy, Var == 'Biomass') %>%
    pull(Mean)
  
  depletion_fmsy <- b_fmsy / b0
  
  # track which run it is along the ramp for each stock
  fidx_fmsy <- this_f_df %>%
    filter(f == atlantis_fmsy) %>%
    pull(fidx) %>%
    unique()
  
  atlantis_fmsy_ls[[i]] <- data.frame('LongName' = sp[i], 
                                      'atlantis_fmsy' = atlantis_fmsy,
                                      'b_fmsy' = b_fmsy,
                                      'depletion' = depletion_fmsy,
                                      'fidx' = fidx_fmsy)
}

atlantis_fmsy <- bind_rows(atlantis_fmsy_ls)

# save this for future calculations
# write.csv(atlantis_fmsy, "data/f35_vector_PROXY_OY_SS.csv", row.names = F)

# annotations for the plots (atlantis depletion)
annotations <- atlantis_fmsy %>% 
  mutate(depletion=round(depletion,digits=2), atlantis_fmsy=round(atlantis_fmsy,digits = 2))

# plot
f_df_ms <- f_df

f_df_ms$LongNamePlot <- gsub(" ", "\n", f_df_ms$LongName)
fmsy$LongNamePlot <- gsub(" ", "\n", fmsy$LongName)
atlantis_fmsy$LongNamePlot <- gsub(" ", "\n", atlantis_fmsy$LongName)
b35$LongNamePlot <- gsub(" ", "\n", b35$LongName)
annotations$LongNamePlot <- gsub(" ", "\n", annotations$LongName)

# selected groups for main text
key_grps <- grps %>% filter(Code %in% c("POL", "COD", "ATF", "HAL", "SBF", "POP", "FFS")) %>% pull(LongName)
p_ms <- f_df_ms %>%
  filter(LongName %in% key_grps) %>%
  ggplot(aes(x = f, y = Mean/1000))+
  geom_line()+
  geom_point(size = 1.5)+
  geom_vline(data = fmsy %>% filter(LongName %in% key_grps), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongName %in% key_grps), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = atlantis_fmsy %>% filter(LongName %in% key_grps) %>% mutate(Var = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongName %in% key_grps) %>% mutate(Var = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1.5,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~Var, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
p_ms
ggsave('results/figures/FIGURE_3.jpeg', p_ms, width = 7, height = 7, dpi = 600)

# make figures with all stocks for supplement
grp1 <- unique(f_df_ms$LongNamePlot)[1:6]
f_plot1 <- f_df_ms %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = f, y = Mean/1000))+
  geom_line()+
  geom_point(size = 2)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp1) %>% mutate(Var = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongNamePlot %in% grp1) %>% mutate(Var = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1.5,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~Var, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot1

grp2 <- unique(f_df_ms$LongNamePlot)[7:12]
f_plot2 <- f_df_ms %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = f, y = Mean/1000))+
  geom_line()+
  geom_point(size = 2)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2) %>% mutate(Var = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongNamePlot %in% grp2) %>% mutate(Var = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1.5,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~Var, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot2

ggsave('results/figures/FIGURE_S4_1.jpeg', f_plot1, width = 7, height = 7)
ggsave('results/figures/FIGURE_S4_2.jpeg', f_plot2, width = 7, height = 7)

####################################
# Analysis of end-of-run variability
####################################
# This has the purpose of providing a general sense of the variability and noise over the last 5 years of a run
# CV are computed for metrics over the last 5 years of each run
# This is no true measure of uncertainty because runs are deterministic and no alternative starts were explored
# Future work should include more comprehensive analyses of uncertainty
cv_df_ss <- f_df_ms %>%
  dplyr::select(LongName, f, Var, CV)

cv_p_ss <- cv_df_ss %>%
  filter(Var == "Biomass") %>%
  ggplot()+
  geom_point(aes(x = f, y = CV), alpha = 0.7)+
  theme_bw()+
  labs(x = "F", y = "CV", fill = "") +
  facet_wrap(~LongName, nrow = 4)
cv_p_ss

ggsave("results/figures/FIGURE_S6.jpeg", cv_p_ss, width = 8, height = 5)

# Figure 4. Global yield ---------------------------------------------------------------
# # list the rds files
f35_results <- list.files(file.path("results", "ms", "flat_results"),
                          pattern = ".rds",
                          full.names = TRUE)
# order them correctly
# reorder these based on the number in the filename
num_idx <- as.numeric(gsub("([0-9]+)-result\\.rds", "\\1",
                           c(list.files(batch_ms_flat, pattern = ".rds", full.names = F))))
f35_results <- f35_results[order(num_idx)]

# create a data frame with catch output from each run within Step 2 (multispecies runs)
catch_list <- list()
for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("results/ms/flat_results/", "", f35_results[i])))
  
  # run information based on the index
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  
  # extract tables from results
  this_result <- readRDS(f35_results[i])
  # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
  if(length(this_result)==1) {
    this_result <- this_result[[1]]
  }
  
  this_catch <- this_result[[3]]
  this_catch <- this_catch %>%
    dplyr::select(c(Time, all_of(t3_fg))) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt") %>%
    group_by(Code) %>%
    slice_max(Time, n = 5) %>%
    summarise(mean_catch = mean(catch_mt),
              catch_cv = sd(catch_mt) / mean(catch_mt)) %>%
    mutate(mult = this_mult,
           run = this_run,
           idx = this_idx)
  
  catch_list[[i]] <- this_catch
  
}

catch_df <- bind_rows(catch_list)

# reshape and calculate total
catch_df_long <- catch_df %>%
  dplyr::select(-catch_cv) %>% # drop uncertainty for the bar plot
  #pivot_longer(-c(run, mult, idx), names_to = "Code", values_to = "mt") %>%
  filter(Code != "HAL") %>%
  group_by(run, mult) %>%
  mutate(total_yield = sum(mean_catch ),
         prop = mean_catch / total_yield) %>%
  ungroup() %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# apply scaling of catch by biomass in AK, so that we only have AK catch now
# source("nc_for_scaling.R")
# catch_scalars <- bind_rows(lapply(f35_nc, get_catch_ak_scalar)) # DO NOT RUN
# save this as it takes so long to produce
# write.csv(catch_scalars, "data/catch_scalars.csv", row.names = F)
catch_scalars <- read.csv("data/catch_scalars.csv")

catch_df_long_ak <- catch_df_long %>%
  left_join(catch_scalars %>% 
              left_join(grps %>% 
                          dplyr::select(Code, Name))) %>%
  mutate(mt_ak = mean_catch * ak_prop)

# add scenario information
catch_df_long_ak <- catch_df_long_ak %>%
  mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"),
                                     "Arrowtooth underexploitation",
                                     "MFMSY varies for all focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
catch_df_long_ak$`F on\narrowtooth` <- factor(catch_df_long_ak$`F on\narrowtooth`,
                                              levels = c("MFMSY varies for all focal groups",
                                                         "Arrowtooth underexploitation"))

# spaces
catch_df_long_ak$LongNamePlot <- gsub(" - ", "\n", catch_df_long_ak$LongName)

global_yield_ms <- catch_df_long_ak %>%
  ggplot(aes(x = mult, y = mt_ak / 1000, fill = LongNamePlot))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_viridis_d()+
  scale_y_continuous(limits = c(0,900))+
  geom_hline(yintercept = 800, color = "red", linetype = "dashed")+
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Catch (1000 mt)", fill = "") +
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(nrow = 4))+
  facet_grid(Climate~`F on\narrowtooth`)
global_yield_ms

ggsave(paste0("results/figures/FIGURE_4.jpeg"), global_yield_ms, width = 6, height = 5, dpi = 600)

# make a table with max catch per scenario for the report
max_catch <- catch_df_long_ak %>%
  select(run, idx, mt_ak) %>%
  group_by(run, idx) %>%
  summarize(total_yield_ak = sum(mt_ak, na.rm = T)) %>%
  select(run, total_yield_ak) %>%
  distinct() %>%
  group_by(run) %>%
  slice_max(total_yield_ak)

# Figure 5. Stocks under B35 ----------------------------------------------
# How many stocks are below 35% B0 for each scenario?
# extract biomass and catch from the multispecies runs, i.e. Step 2
ms_yield_list <- list()

for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("results/ms/flat_results/", "", f35_results[i])))
  
  # run information based on the index
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  
  # extract tables from results
  this_result <- readRDS(f35_results[i])
  # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
  if(length(this_result)==1) {
    this_result <- this_result[[1]]
  }
  
  biomage <- this_result[[2]]
  catch <- this_result[[3]]
  mort <- this_result[[4]]
  
  # now extract data
  # SSB to plot and report in tables
  spawning_biomass <- biomage %>%
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code %in% t3_fg) %>%
    left_join(fspb, by = c('Code','Age')) %>%
    mutate(biomass_mt = biomass_mt * fspb) %>%
    group_by(Time,Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    group_by(Code) %>%
    slice_max(Time, n = 5) %>%
    summarise(mean_biom = mean(biomass_mt),
              biom_cv = sd(biomass_mt) / mean(biomass_mt))
  
  # total catch
  # taking mean of the last 5 years
  catch_vals <- catch %>%
    dplyr::select(c(Time, all_of(t3_fg))) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt") %>%
    group_by(Code) %>%
    slice_max(Time, n = 5) %>%
    summarise(mean_catch = mean(catch_mt),
              catch_cv = sd(catch_mt) / mean(catch_mt))
  
  # # calculate realized F after 1 year of data
  # For runs with a burn-in, this has to be the biomass at the end of the burn-in, when we start fishing with the new scalar
  # # get initial biomass for the selected age classes
  biom_age_t1 <- biomage %>%
    filter(Time == 365 * burnin) %>%# this is the burn-in years
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass = sum(biomass)) %>%
    ungroup() %>%
    filter(Code %in% t3_fg)
  #
  # # catch (one time step after biomass: how much did we catch in this time?)
  catch_t1 <- catch %>%
    select(Time, all_of(t3_fg)) %>%
    filter(Time == 365 * (burnin + 1)) %>% # careful - there is a small transition phase
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
    select(-Time)
  #
  # # calc realized f
  f_t1 <- biom_age_t1 %>% left_join(catch_t1, by = "Code") %>%
    mutate(exp_rate = catch/biomass,
           f = -log(1-exp_rate)) %>%#,
    select(Code, f)#, fidx)
  
  # # bind all
  f_frame <- f_t1 %>%
    left_join(spawning_biomass) %>%
    left_join(catch_vals) %>%
    mutate(run = this_run,
           mult = this_mult)
  
  # reshape
  f_frame <- f_frame %>%
    pivot_longer(
      cols = c(mean_biom, mean_catch, biom_cv, catch_cv),
      names_to = "temp",
      values_to = "value"
    ) %>%
    mutate(
      Var = case_when(
        temp == "mean_biom" ~ "Biomass",
        temp == "mean_catch" ~ "Catch",
        temp == "biom_cv" ~ "Biomass",
        temp == "catch_cv" ~ "Catch"
      ),
      type = ifelse(grepl("mean", temp), "Mean", "CV")
    ) %>%
    select(-temp) %>%
    pivot_wider(
      names_from = type,
      values_from = value
    ) %>%
    select(Code, f, run, mult, Var, Mean, CV) %>%
    left_join(grps %>% select(Code, LongName), by = 'Code')
  
  # add to multispecies yield list
  ms_yield_list[[i]] <- f_frame
}

ms_yield_df <- bind_rows(ms_yield_list)

# handle the NaN's from FHS
ms_yield_df <- as.data.frame(ms_yield_df)
ms_yield_df$f[is.nan(ms_yield_df$f)] <- NA

# get static B0 from Base Scenario
b0 <- ms_yield_df %>% filter(mult == 0, Var == "Biomass", run == "base") %>% dplyr::select(LongName, Mean) %>% rename(b0 = Mean)

below_target <- ms_yield_df %>%
  filter(Var == "Biomass") %>%
  filter(!(run %in% c("atf", "atf_climate") & LongName == "Arrowtooth flounder")) %>%
  left_join(b0, by = "LongName") %>%
  mutate(depletion = Mean / b0) %>%
  mutate(below_target = ifelse(depletion < 0.35, 1, 0)) %>%
  group_by(run, mult) %>%
  mutate(n_below_target = sum(below_target)) %>%
  mutate(prop_below_target = n_below_target / length(unique(LongName))) %>%
  ungroup() %>%
  select(run, mult, n_below_target, prop_below_target) %>%
  distinct()

# add scenario information
below_target <- below_target %>%
  mutate(`F on\narrowtooth` = ifelse(run %in% c("atf","atf_climate"), 
                                     "Arrowtooth underexploitation",
                                     "MFMSY varies for all focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
below_target$`F on\narrowtooth` <- factor(below_target$`F on\narrowtooth`, 
                                          levels = c("MFMSY varies for all focal groups",
                                                     "Arrowtooth underexploitation"))

p_below_target <- below_target %>%
  ggplot(aes(x = mult, y = n_below_target))+
  geom_col()+
  theme_bw()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Stocks with SSB < B35%") +
  scale_y_continuous(breaks = 0:12)+
  theme(panel.grid.minor = element_blank())+
  facet_grid(Climate~`F on\narrowtooth`)
p_below_target

ggsave(paste0("results/figures/FIGURE_5.jpeg"), p_below_target, width = 6, height = 4, dpi = 600)

# Figure 6. Biomass and catch curves ------------------------------------------------
# get max yield (proxy of MSY)
ymax_ms <- ms_yield_df %>%
  filter(Var == "Catch") %>%
  group_by(LongName, run) %>%
  slice_max(Mean) %>%
  ungroup() %>%
  dplyr::select(LongName, run, Mean, f) %>%
  rename(ymax = Mean)

ymax <- ymax_ms

# plot catch and biomass curves
to_plot <- ms_yield_df

# spaces
to_plot$LongNamePlot <- gsub(" ", "\n", to_plot$LongName)
ymax$LongNamePlot <- gsub(" ", "\n", ymax$LongName)

# rename scenarios and order them
to_plot <- to_plot %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                          "Arrowtooth\nunderexploitation",
                          "MFMSY varies for\nall focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585\n(2075-2085)", "Historical\n(1999)"))

ymax <- ymax %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                          "Arrowtooth\nunderexploitation",
                          "MFMSY varies for\nall focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585\n(2075-2085)", "Historical\n(1999)"))

# reorder ATF F                             
to_plot$Fishing <- factor(to_plot$Fishing,
                          levels = c("MFMSY varies for\nall focal groups",
                                     "Arrowtooth\nunderexploitation"))
ymax$Fishing <- factor(ymax$Fishing,
                       levels = c("MFMSY varies for\nall focal groups",
                                  "Arrowtooth\nunderexploitation"))

# key focal groups only for main text
key_grps <- grps %>% filter(Code %in% c("POL", "COD", "ATF", "HAL", "SBF", "POP", "FFS")) %>% pull(LongName)
f_plot_ms <- to_plot %>%
  filter(LongName %in% key_grps, Var == "Catch") %>%
  ggplot(aes(x = f, y = Mean/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% 
               filter(LongName %in% key_grps) %>% 
               filter(!(LongNamePlot == "Arrowtooth\nflounder" & Fishing == "Arrowtooth\nunderexploitation")), 
             aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
  facet_grid2(LongNamePlot~Var, scales = 'free', independent = 'all')+
  #facet_wrap(~ LongName, scales = "free", ncol = 1)+
  theme(strip.text.y = element_text(angle=0))
f_plot_ms

ggsave('results/figures/FIGURE_6.jpeg', f_plot_ms, width = 4.8, height = 6.5, dpi = 600)

# make figures for supplement (break into two sets)
grp1 <- unique(to_plot$LongNamePlot)[1:6]
f_plot1 <- to_plot %>%
  filter(LongNamePlot %in% grp1) %>%
  filter(Var == "Catch") %>%
  ggplot(aes(x = f, y = Mean/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% filter(LongNamePlot %in% grp1), aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
  facet_grid2(LongNamePlot~Var, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot1

grp2 <- unique(to_plot$LongNamePlot)[7:12]
f_plot2 <- to_plot %>%
  filter(LongNamePlot %in% grp2) %>%
  filter(Var == "Catch") %>%
  ggplot(aes(x = f, y = Mean/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% filter(LongNamePlot %in% grp2), aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
  facet_grid2(LongNamePlot~Var, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot2

ggsave('results/figures/FIGURE_S9_1.jpeg', f_plot1, width = 7.5, height = 7)
ggsave('results/figures/FIGURE_S9_2.jpeg', f_plot2, width = 7.5, height = 7)

####################################
# Analysis of end-of-run variability
####################################
cv_df <- to_plot %>%
  dplyr::select(LongName, Fishing, Climate, mult, Var, CV)

# by species
cv_p1 <- cv_df %>%
  filter(Var == "Biomass") %>%
  ggplot(aes(x = LongName, y = CV, group = factor(mult), color =factor(mult)))+
  stat_summary(fun.data = function(x) {
    return(c(y = mean(x), ymin = min(x), ymax = max(x)))
  }, 
  geom = "pointrange",
  position = position_dodge(width = 0.75),
  size = 0.3) +
  scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9)+
  theme_bw()+
  guides(color = guide_legend(ncol = 2))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  labs(x = "", y = "CV", color = expression(atop(MF[MSY], "multiplier")))
cv_p1

# by variable
cv_p2 <- cv_df %>%
  ggplot()+
  geom_boxplot(aes(x = mult, y = CV, group = interaction(mult,Var), fill = Var, color = Var), alpha = 0.7)+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  scale_fill_viridis_d(option = "inferno", begin = 0.2, end = 0.8)+
  scale_color_viridis_d(option = "inferno", begin = 0.2, end = 0.8)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "CV", fill = "") +
  guides(color = "none")+
  facet_grid(Climate~Fishing)
cv_p2

cv_combo <- cv_p1 / cv_p2 + # Stack plots vertically
  plot_layout(heights = c(1, 1.5)) +  # 2:1 ratio between plots
  plot_annotation(tag_levels = 'A') # Add letters A, B

# Save the combined plot
ggsave("results/figures/FIGURE_S10.jpeg", cv_combo, height = 8.5, width = 8, dpi = 600)

# Figure 7: top predators and forage fish ------------------------------
top_preds <- c("SSL","PIN","DOL","BDF","BSF")
forage <- c("CAP","SAN","HER","EUL","FOS")
other_fg <- c(top_preds, forage)

ms_other_list <- list()

for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("results/ms/flat_results/", "", f35_results[i])))
  
  # run information based on the index
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  
  # extract tables from results
  this_result <- readRDS(f35_results[i])
  # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
  if(length(this_result)==1) {
    this_result <- this_result[[1]]
  }
  
  biomage <- this_result[[2]]
  
  # now extract data
  # SSB to plot and report in tables
  other_biomass <- biomage %>% 
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code %in% other_fg) %>%
    group_by(Time, Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    group_by(Code) %>%
    slice_max(Time, n = 5) %>%
    summarise(mean_biom = mean(biomass_mt),
              biom_cv = sd(biomass_mt) / mean(biomass_mt)) %>%
    mutate(run = this_run,
           mult = this_mult)
  
  # add to multispecies yield list
  ms_other_list[[i]] <- other_biomass
}

ms_other_df <- bind_rows(ms_other_list) %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# get b0
# static to "base"
b0_other <- ms_other_df %>% filter(mult == 0, run == "base") %>% dplyr::select(LongName, mean_biom) %>% rename(b0 = mean_biom)

ms_other_df <- ms_other_df %>%
  left_join(b0_other, by = c("LongName")) %>%
  mutate(biomchange = (mean_biom - b0)/b0 * 100)

# add factors for plot
ms_other_df <- ms_other_df %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                          "Arrowtooth\nunderexploitation",
                          "MFMSY varies for\nall focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
ms_other_df$Fishing <- factor(ms_other_df$Fishing,
                              levels = c("MFMSY varies for\nall focal groups",
                                         "Arrowtooth\nunderexploitation"))

# handle dash
ms_other_df$LongName <- gsub(" - "," ",ms_other_df$LongName)

# handle long names for the facet for predators, they are too wide
ms_other_df$LongNamePlot <- gsub(" ","\n",ms_other_df$LongName)

# order groups
ms_other_df$LongNamePlot <- factor(ms_other_df$LongNamePlot, levels = c(
  "Steller\nsea\nlion", 
  "Other\npinnipeds",
  "Dolphins",
  "Seabirds\nsurface\nfish",
  "Seabirds\ndiving\nfish",
  "Capelin",
  "Sandlance",
  "Pacific\nherring",
  "Forage\nfish\nslope",
  "Eulachon"
))

####################################
# Analysis of end-of-run variability
####################################
cv_df_other <- ms_other_df %>%
  dplyr::select(Code, LongName, Fishing, Climate, mult, biom_cv)

# split forage from top predators
cv_df_other <- cv_df_other %>%
  mutate(Guild = case_when(
    Code %in% top_preds ~ "Predator",
    TRUE ~ "Forage"
  ))

cv_p3 <- cv_df_other %>%
  #filter(mult > 0) %>%
  ggplot()+
  geom_boxplot(aes(x = mult, y = biom_cv, group = interaction(mult,Guild), fill = Guild, color = Guild), alpha = 0.7)+
  scale_x_continuous(breaks = seq(0,4,0.5))+
  scale_fill_viridis_d(option = "inferno", begin = 0.2, end = 0.8)+
  scale_color_viridis_d(option = "inferno", begin = 0.2, end = 0.8)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "CV", fill = "") +
  guides(color = "none")+
  facet_grid(Climate~Fishing)
cv_p3
ggsave("results/figures/FIGURE_S14.jpeg", cv_p3, width = 7, height = 4.5)

# for each predator, identify the main prey species from dietcheck (in baseline)
# sum up total prey biomass
# express changes from B0
# This is qualitative, but it demonstrate a likely trophic link and its effects

# read in diet comps for the base model and get the last 5 years
base_diet <- read.table("results/diets/output_0DietCheck.txt", sep = " ", header = T)

# put in long format
diet_long_other <- base_diet %>%
  mutate(Time = Time / 365) %>%
  filter(Time > 75 & Time <=80) %>% # 
  group_by(Predator) %>%
  summarise(across(KWT:DR, mean)) %>%
  ungroup() %>%
  pivot_longer(-Predator, names_to = 'Prey', values_to = 'Prop') %>%
  left_join((grps %>% select(Code, Name, LongName)), by = c('Predator'='Code')) %>%
  rename(Predator_Name = Name, Predator_LongName = LongName) %>%
  select(-Predator) %>%
  left_join((grps %>% select(Code, Name, LongName)), by = c('Prey'='Code')) %>%
  rename(Prey_Name = Name, Prey_LongName = LongName) %>%
  select(Prop, Predator_Name, Predator_LongName, Prey_Name, Prey_LongName)%>%
  filter(Prop > 0)

# for each top predator, which are the prey species?
top_pred_names <- grps %>% filter(Code %in% top_preds) %>% pull(Name) %>% as.character()

prey_per_predator <- list()
for(i in 1:length(top_pred_names)){
  
  this_top_pred <- top_pred_names[i]
  fav_prey <- diet_long_other %>%
    filter(Predator_Name == this_top_pred) %>%
    pull(Prey_Name)
  
  # df
  prey_per_predator[[i]] <- data.frame("Predator" = this_top_pred, "Prey" = fav_prey)
  
}
prey_per_predator <- bind_rows(prey_per_predator)

# for each prey, what are the predators?
forage_names <- grps %>% filter(Code %in% forage) %>% pull(Name) %>% as.character()

predator_per_prey <- list()
for(i in 1:length(forage_names)){
  
  this_forage <- forage_names[i]
  fav_predator <- diet_long_other %>%
    filter(Prey_Name == this_forage) %>%
    pull(Predator_Name)
  
  # df
  predator_per_prey[[i]] <- data.frame("Prey" = this_forage, "Predator" = fav_predator)
  
}
predator_per_prey <- bind_rows(predator_per_prey)

# now back to the biomasses
# for each predator or prey, loop over the results to extract, for each run, the total terminal biomass of the group of prey (or predators)
preds_and_prey <- c(top_pred_names, forage_names)

diet_biomass <- lapply(1:length(preds_and_prey), function(i) {
  # Create an inner list
  rep(list(NULL), length(f35_results))
})

for(i in 1:length(preds_and_prey)){
  this_sp <- preds_and_prey[i]
  
  # identify the groups that are the favorite prey or predator
  if(this_sp %in% top_pred_names){
    diet_grps <- prey_per_predator %>% filter(Predator == this_sp) %>% pull(Prey)
  } else {
    diet_grps <- predator_per_prey %>% filter(Prey == this_sp) %>% pull(Predator)
  }
  
  # bring in codes again as that's what the output works with
  diet_codes <- grps %>% filter(Name %in% diet_grps) %>% pull(Code)
  
  # now loop over results
  for(j in 1:length(f35_results)){
    
    print(paste("Doing", f35_results[j]))
    
    # grab the index from the file name
    this_idx <- as.numeric(gsub("-result.rds", "", gsub("results/ms/flat_results/", "", f35_results[j])))
    
    # run information based on the index
    this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
    this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
    
    # extract tables from results
    this_result <- readRDS(f35_results[j])
    # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
    if(length(this_result)==1) {
      this_result <- this_result[[1]]
    }
    
    biomage <- this_result[[2]]
    
    # now extract data
    # not bothering with the CV here
    # SSB to plot and report in tables
    this_diet_biomass <- biomage %>% 
      slice_tail(n = 5) %>% # use last 5 years
      summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
      ungroup() %>%
      pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
      separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
      filter(Code %in% diet_codes) %>%
      #group_by(Code) %>%
      summarise(biomass_of_prey_or_pred = sum(biomass_mt)) %>%
      #ungroup() %>%
      mutate(Name = this_sp, run = this_run, mult = this_mult)
    
    # add to list
    diet_biomass[[i]][[j]] <- this_diet_biomass
    
  }
  
  diet_biomass[[i]] <- bind_rows(diet_biomass[[i]])
  
}

diet_biomass <- bind_rows(diet_biomass)

# now need to rescale to b0 of the respective scenarios
b0_for_diets <- diet_biomass %>% filter(mult == 0) %>% dplyr::select(Name, run, biomass_of_prey_or_pred) %>% rename(b0 = biomass_of_prey_or_pred)

diet_biomass_scalars <- diet_biomass %>%
  left_join(b0_for_diets, by = c("Name", "run")) %>%
  mutate(scalar = (biomass_of_prey_or_pred - b0)/b0*100) %>%
  select(Name, run, mult, scalar) %>%
  left_join(grps %>% select(Name, LongName))

# fix dash
diet_biomass_scalars$LongName <- gsub(" - ", " ", diet_biomass_scalars$LongName)

# now join this to the ms_other_df frame
ms_other_df_diet <- ms_other_df %>%
  left_join(diet_biomass_scalars, by = c("LongName","run","mult"))

# now plot
other_plot_forage_diets <- ms_other_df_diet %>%
  filter(Code %in% c("CAP","SAN","HER","FOS","EUL")) %>%
  ggplot(aes(x = mult, y = mean_biom / 1000, fill = scalar, shape = Fishing))+
  geom_point(color = "black", size = 1.5)+
  scale_shape_manual(values = c(21,24))+
  colorspace::scale_fill_continuous_divergingx(palette = 'PRGn', mid = 0) + 
  geom_vline(xintercept = 1, color = 'black', linetype = "dashed", linewidth = 0.35)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Biomass (1000 mt)", fill = "Change in total predator\nbiomass from unfished (%)")+
  #guides(fill=guide_legend(order=1), shape=guide_legend(order=2))+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
other_plot_forage_diets

other_plot_top_diets <- ms_other_df_diet %>%
  filter(Code %in% c("DOL","SSL","PIN","BDF","BSF")) %>%
  ggplot(aes(x = mult, y = mean_biom / 1000, fill = scalar, shape = Fishing))+
  geom_point(color = "black", size = 1.5)+
  scale_shape_manual(values = c(21,24))+
  colorspace::scale_fill_continuous_divergingx(palette = 'PRGn', mid = 0) + 
  geom_vline(xintercept = 1, color = 'black', linetype = "dashed", linewidth = 0.35)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Biomass (1000 mt)", fill = "Change in total prey\nbiomass from unfished (%)")+
  #guides(fill=guide_legend(order=1), shape=guide_legend(order=2))+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
other_plot_top_diets

# combine into one figure
TL_combo <- other_plot_forage_diets / other_plot_top_diets + # Stack plots vertically
  plot_layout(heights = c(1, 1)) +  # 2:1 ratio between plots
  plot_annotation(tag_levels = 'a') # Add letters A, B

# Save the combined plot
ggsave("results/figures/FIGURE_7.jpeg", TL_combo, height = 8, width = 8, dpi = 600)

#########################
# SUPPLEMENTARY FIGURES #
#########################

# Diets -----------------------------------------
diet_runs <- data.frame("idx" = c(0,1,2,3,4),
                        "run" = c("Base model,\ncalibration fishing",
                                  "MFMSY varies for\nall focal groups,\nhistorical climate",
                                  "Arrowtooth\nunderexploitation,\nhistorical climate",
                                  "MFMSY varies for\nall focal groups,\nssp585",
                                  "Arrowtooth\nunderexploitation,\nssp585"))
diet_long_cohort <- list()

for(r in 1:nrow(diet_runs)){
  
  this_run_no <- diet_runs[r,1]
  this_run_lab <- diet_runs[r,2]
  
  diet <- read.table(paste0("results/diets/output_", this_run_no, "DietCheck.txt"), sep = " ", header = T)
  
  diet_long_cohort[[r]] <- diet %>%
    mutate(Time = Time / 365) %>%
    filter(Time > 75 & Time <=80) %>% # focus on last 5 years of the run
    group_by(Predator, Cohort) %>%
    summarise(across(KWT:DR, mean)) %>%
    ungroup() %>%
    pivot_longer(-c(Predator, Cohort), names_to = 'Prey', values_to = 'Prop') %>%
    left_join((grps %>% select(Code, Name, LongName)), by = c('Predator'='Code')) %>%
    rename(Predator_Name = Name, Predator_LongName = LongName) %>%
    select(-Predator) %>%
    left_join((grps %>% select(Code, Name, LongName)), by = c('Prey'='Code')) %>%
    rename(Prey_Name = Name, Prey_LongName = LongName) %>%
    select(Prop, Predator_Name, Predator_LongName, Cohort, Prey_Name, Prey_LongName)%>%
    filter(Prop > 0.005) %>%
    mutate(idx = this_run_no, run = this_run_lab)
  
}

diet_long_cohort <- bind_rows(diet_long_cohort)

pred_names <- grps %>% filter(Code %in% c("ATF", "HAL", top_preds)) %>% pull(Name)

diet_long_preds <- diet_long_cohort %>% filter(Predator_Name %in% pred_names)

# reorder factors for plotting
diet_long_preds$run <- factor(diet_long_preds$run, levels = c("Base model,\ncalibration fishing",
                                                              "MFMSY varies for\nall focal groups,\nhistorical climate",
                                                              "Arrowtooth\nunderexploitation,\nhistorical climate",
                                                              "MFMSY varies for\nall focal groups,\nssp585",
                                                              "Arrowtooth\nunderexploitation,\nssp585"))

# plot arrowtooth flounder alone for figure S1.3
atf_diet <- diet_long_preds %>%
  filter(Predator_Name == "Arrowtooth_flounder") %>%
  drop_na()

colourCount <- length(unique(atf_diet$Prey_Name))
colors <- c(viridis(11)[2:10], inferno(11)[2:10])#, cividis(6))

# plot ATF in the base model only (Fig. S1.3)
p_atf_diet_base <- atf_diet %>%
  filter(run == "Base model,\ncalibration fishing") %>%
  ggplot(aes(x = Cohort+1, y = Prop * 100, fill = Prey_LongName))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_x_continuous(breaks = 1:10)+
  scale_fill_manual(values = colors)+
  theme_bw()+
  labs(x = '', y = "Diet preference (%)", fill = "Prey")
p_atf_diet_base
ggsave("results/figures/diet_plots/FIGURE_S2.jpeg", p_atf_diet_base, width = 6, height = 5)

# plot across scenarios (except Base model)
p_atf_diet <- atf_diet %>%
  filter(run != "Base model,\ncalibration fishing") %>%
  ggplot(aes(x = Cohort+1, y = Prop * 100, fill = Prey_LongName))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_x_continuous(breaks = 1:10)+
  scale_fill_manual(values = colors)+
  theme_bw()+
  labs(x = '', y = "Diet preference (%)", fill = "Prey")+
  facet_wrap(~run)
p_atf_diet
ggsave("results/figures/diet_plots/FIGURE_S12.jpeg", p_atf_diet, width = 8, height = 6)

# Halibut
hal_diet <- diet_long_preds %>%
  filter(Predator_Name == "Halibut") %>%
  drop_na()

colourCount <- length(unique(hal_diet$Prey_Name))
colors <- c(viridis(9)[2:8], inferno(8)[2:7])#, cividis(6))

# plot Base diets for Halibut
p_hal_diet <- hal_diet %>%
  filter(run == "Base model,\ncalibration fishing") %>%
  ggplot(aes(x = Cohort+1, y = Prop * 100, fill = Prey_LongName))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_x_continuous(breaks = 1:10)+
  scale_fill_manual(values = colors)+
  theme_bw()+
  labs(x = '', y = "Diet preference (%)", fill = "Prey")
p_hal_diet
ggsave("results/figures/diet_plots/FIGURE_S3.jpeg", p_hal_diet, width = 6, height = 5)

# plot all together
# drop ATF here
diet_long_preds <- diet_long_preds %>%
  filter(!Predator_Name %in% c("Arrowtooth_flounder", "Halibut"))

# pred longname for facets
diet_long_preds$Predator_LongNamePlot <- gsub(" ", "\n", diet_long_preds$Predator_LongName)

colourCount <- length(unique(diet_long_preds$Prey_Name))
colors <- c(viridis(8)[2:8], inferno(8)[2:8])#, cividis(6))

p_all_diet <- diet_long_preds %>%
  ggplot(aes(x = Cohort+1, y = Prop * 100, fill = Prey_LongName))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_x_continuous(breaks = 1:10)+
  scale_fill_manual(values = colors)+
  theme_bw()+
  labs(x = 'Age class', y = "Diet composition (%)", fill = "Prey")+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(nrow = 4))+
  facet_grid2(Predator_LongNamePlot~run)+
  theme(strip.text.y = element_text(angle=0))
p_all_diet

ggsave("results/figures/diet_plots/FIGURE_S13.jpeg", p_all_diet, width = 8.5, height = 6)

# Production curves --------------------------------------------------
# These refer to the single-species runs in Step 1
ss_yield_long <- f_df %>%
  select(-CV)

# get b0 from SS curves
b0 <- ss_yield_long %>% 
  group_by(Code) %>%
  slice_min(f) %>% 
  ungroup() %>%
  filter(Var == "Biomass") %>% 
  dplyr::select(LongName, Mean) %>% 
  rename(b0 = Mean)

# get max yield
ymax <- ss_yield_long %>% 
  filter(Var == "Catch") %>% 
  group_by(LongName) %>%
  slice_max(Mean) %>%
  ungroup() %>%
  dplyr::select(LongName, Mean, f) %>% 
  rename(ymax = Mean) 

# we are plotting yield fraction against depletion
yield_func <- ss_yield_long %>%
  drop_na() %>%
  dplyr::select(LongName, Var, f, Mean) %>%
  pivot_wider(id_cols = c(LongName, f), names_from = Var, values_from = Mean) %>%
  left_join(b0, by = c("LongName")) %>% # if you are keep static reference point
  mutate(depletion = Biomass / b0) %>%
  left_join(ymax %>% select(-f), by = c("LongName")) %>%
  mutate(yfrac = Catch / ymax) %>%
  dplyr::select(LongName, yfrac, depletion, f)

# prepare data frames to write the following quantities on the plot:
# final depletion, final yield fraction, biomass corresponding to final depletion, biomass corresponding to final yield fraction
yfun_terminal <- yield_func %>%
  group_by(LongName) %>%
  slice_max(f) %>% # for each group, get highest f
  ungroup() %>%
  mutate(depletion = depletion * 100,
         yfrac = yfrac * 100) # turn depletion and yield to percentages for easier interpretation

annotations <- b0 %>%
  left_join(ymax, by = "LongName") %>%
  left_join(yfun_terminal, by = "LongName") %>%
  mutate(catch = ymax / 100 * yfrac / 1000,
         ssb = b0 / 100 * depletion / 1000) 

# plot
yield_func_plot <- yield_func %>%
  #filter(LongName %in% c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", "Pacific halibut")) %>%
  ggplot(aes(x = depletion, y = yfrac))+
  geom_point()+
  geom_line()+
  scale_x_reverse()+
  geom_text(data = annotations,
            aes(x = 0.5, y = 0.5, hjust=0.5, vjust=1,
                label=paste0('D(%)=', round(depletion,2),
                             '\n',
                             'YF(%)=', round(yfrac, 2),
                             '\n',
                             'SSB(1000mt)=', round(ssb, 2),
                             '\n',
                             'Catch(1000mt)=', round(catch, 2))),
            size = 3)+
  theme_bw()+
  labs(x = "Depletion", y = "Yield fraction")+
  facet_wrap(~LongName)
yield_func_plot

ggsave("results/figures/FIGURE_S5.jpeg", yield_func_plot, width = 8, height = 6.5)

# Biomass and catch ratio of focal grps vs total groundfish ---------------
# start from the base model (modified from Rovellini et al. 2024)
biom_base <- read.table("data/base_model_results/outputBiomIndx.txt", sep = " ", header = T)
biom_base <- biom_base %>%
  select(Time:DR) %>%
  pivot_longer(-Time, names_to = "Code", values_to = "mt")

# all fmp groups (i.e., no halibut)
all_fmp <- c("SHD", "SHP", "DOG", "POL", "COD", "ATF", "FHS", "REX", "FFS", "FFD", "SKL", "SKB", "SKO", "SBF", "POP", "RFS", "RFP", "RFD", "THO", "DFS", "DFD", "SCU")

biom_fmp <- biom_base %>%
  filter(Code %in% all_fmp) %>%
  mutate(is_focal = ifelse(Code %in% t3_fg, 1, 0))

# now proportion of is focal biomass over time
biom_fmp %>%
  group_by(Code) %>%
  slice_max(Time, n = 5) %>% # last 5 years
  group_by(Time, is_focal) %>%
  summarise(by_focal = sum(mt)) %>%
  group_by(Time) %>%
  mutate(tot = sum(by_focal)) %>%
  ungroup() %>%
  mutate(prop_focal = by_focal / tot) %>%
  filter(is_focal == 1)
# If you consider cephalopods, ~60% of the biomass in base run is from the focal groups, slowly declining
# if you only count fish, >75%

# now catch
catch_base <- read.table("data/base_model_results/outputCatch.txt", sep = " ", header = T)
catch_base <- catch_base %>%
  select(Time:BIV) %>%
  pivot_longer(-Time, names_to = "Code", values_to = "mt")

catch_fmp <- catch_base %>%
  filter(Code %in% all_fmp) %>%
  mutate(is_focal = ifelse(Code %in% t3_fg, 1, 0))

catch_fmp %>%
  group_by(Code) %>%
  slice_max(Time, n = 5) %>%
  group_by(Time, is_focal) %>%
  summarise(by_focal = sum(mt)) %>%
  group_by(Time) %>%
  mutate(tot = sum(by_focal)) %>%
  ungroup() %>%
  mutate(prop_focal = by_focal / tot) %>%
  filter(is_focal == 1)
# >90% of the catch is from the focal groups, fairly stable

# now do the same for the MS runs from step 2
# this has the purpose of showing the relation between the focal groups and total biomass
# source("nc_for_scaling.R")
# catch_scalars <- bind_rows(lapply(f35_nc, get_catch_ak_scalar)) # this is very slow, it can be optimized in many ways
# save this as it takes so long to produce
# write.csv(catch_scalars, "data/catch_scalars_fmp.csv", row.names = F)
catch_scalars <- read.csv("data/catch_scalars_fmp.csv")

all_yield_list <- list()

for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("results/ms/flat_results/", "", f35_results[i])))
  
  # run information based on the index
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  
  # extract tables from results
  this_result <- readRDS(f35_results[i])
  # the packaging of the RDS object was different between the eScience runs and the batch (doAzureParallel) runs
  if(length(this_result)==1) {
    this_result <- this_result[[1]]
  }
  
  biomage <- this_result[[2]]
  catch <- this_result[[3]]
  
  # bring in AK vs BC scalars
  this_scalar = catch_scalars %>% filter(idx == this_idx)
  
  # now extract data
  # get the end-of-run proportion of focal grop biomass / total gf biomass
  total_biomass <- biomage %>%
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code %in% all_fmp) %>%
    filter(Code != "HAL") %>% # drop halibut for this, it is not in OY
    group_by(Time,Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>% # sum across cohorts
    ungroup() %>%
    left_join(grps %>% select(Code, Name), by = "Code") %>%
    left_join(this_scalar, by = "Name") %>%
    mutate(biomass_mt = biomass_mt * ak_prop) %>%
    mutate(is_focal = ifelse(Code %in% t3_fg, 1, 0)) %>%
    group_by(Time, is_focal) %>%
    summarise(step1 = sum(biomass_mt, na.rm = T)) %>%
    group_by(Time) %>%
    mutate(step2 = sum(step1)) %>%
    ungroup() %>%
    mutate(prop = step1 / step2) %>%
    filter(is_focal == 1) %>%
    slice_max(Time, n = 5) %>%
    summarise(mean_prop = mean(prop),
              mean_focal = mean(step1),
              mean_total = mean(step2),
              cv_prop = sd(prop) / mean(prop),
              cv_focal = sd(step1) / mean(step1),
              cv_total = sd(step2) / mean(step2)) %>%
    mutate(Var = "Biomass")
  
  # total catch
  total_catch <- catch %>%
    dplyr::select(c(Time, all_of(all_fmp))) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt") %>%
    filter(Code %in% all_fmp) %>%
    filter(Code != "HAL") %>% # drop halibut for this, it is not in OY
    left_join(grps %>% select(Code, Name), by = "Code") %>%
    left_join(this_scalar, by = "Name") %>%
    mutate(catch_mt = catch_mt * ak_prop) %>%
    mutate(is_focal = ifelse(Code %in% t3_fg, 1, 0)) %>%
    group_by(Time, is_focal) %>%
    summarise(step1 = sum(catch_mt, na.rm = T)) %>%
    group_by(Time) %>%
    mutate(step2 = sum(step1)) %>%
    ungroup() %>%
    mutate(prop = step1 / step2) %>%
    filter(is_focal == 1) %>%
    slice_max(Time, n = 5) %>%
    summarise(mean_prop = mean(prop),
              mean_focal = mean(step1),
              mean_total = mean(step2),
              cv_prop = sd(prop) / mean(prop),
              cv_focal = sd(step1) / mean(step1),
              cv_total = sd(step2) / mean(step2)) %>%
    mutate(Var = "Catch")
  
  # # bind all
  prop_df <- rbind(total_biomass, total_catch) %>%
    mutate(run = this_run,
           mult = this_mult,
           idx = this_idx)
  
  # add to multispecies yield list
  all_yield_list[[i]] <- prop_df
}

all_yield_df <- bind_rows(all_yield_list)

fmp_key <- data.frame("run" = c("base","atf","climate","atf_climate"),
                      "run_lab" = c("MFMSY varies for\nall focal groups,\nhistorical climate",
                                    "Arrowtooth\nunderexploitation,\nhistorical climate",
                                    "MFMSY varies for\nall focal groups,\nssp585",
                                    "Arrowtooth\nunderexploitation,\nssp585"))

all_yield_df <- all_yield_df %>% left_join(fmp_key)

# reorder factors for plotting
all_yield_df$run_lab <- factor(all_yield_df$run_lab, levels = c("Base model,\ncalibration fishing",
                                                                "MFMSY varies for\nall focal groups,\nhistorical climate",
                                                                "Arrowtooth\nunderexploitation,\nhistorical climate",
                                                                "MFMSY varies for\nall focal groups,\nssp585",
                                                                "Arrowtooth\nunderexploitation,\nssp585"))

# make a plot
all_yield_df_2 <- all_yield_df %>%
  select(run_lab, mult, Var, mean_focal, mean_total, cv_focal, cv_total) %>%
  pivot_longer(-c(run_lab,mult,Var), names_to = "type_grp", values_to = "mt") %>%
  separate(type_grp, into = c("type", "grp"), sep = "_") %>%
  pivot_wider(names_from = type, values_from = mt) %>%
  mutate(grp = gsub("focal","Focal groups", grp),
         grp = gsub("total", "All GOA FMP species", grp))

p_fmp_2 <- all_yield_df_2 %>%
  filter(mult > 0) %>%
  ggplot(aes(x = mult, y = mean/1000, shape = grp))+
  geom_point(size = 1.5)+
  scale_shape_manual(values = c(1,2)) +
  geom_hline(yintercept = 800, color = "red", linetype = "dashed")+
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35)+
  labs(x = expression(MF[MSY] ~ "multiplier"), 
       y = "1000 mt", 
       color = "",
       shape = "") +
  scale_x_continuous(limits = c(0,4))+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  facet_grid2(Var~run_lab, scales = "free_y")
p_fmp_2
ggsave("results/figures/FIGURE_S7.jpeg", p_fmp_2, width = 8, height = 4.5)

# Walters plot ------------------------------------------------------
# create a plot akin to Walters et al. (2005) Fig. 3, except we do not organize it by TL
# treat this as model-wide MSY
ss_msy <- f_df %>%
  mutate(experiment = "ss") %>%
  select(Code, LongName, f, fidx, experiment, Var, Mean) %>%
  filter(Var == "Catch") %>%
  group_by(Code, LongName) %>%
  slice_max(Mean) %>%
  ungroup() %>%
  select(LongName, Mean) %>%
  rename(mt_ss = Mean)

# this is the scenario where ATF is kept at low fishing pressure, all other groups are varied
ms_msy <- catch_df %>%
  select(-catch_cv) %>%
  group_by(run, mult) %>%
  mutate(total_yield = sum(mean_catch),
         prop = mean_catch / total_yield) %>%
  ungroup() %>%
  left_join(grps %>% select(Code, LongName), by = "Code") %>%
  filter(mult == 1, run %in% c("base","atf")) %>%
  select(LongName, run, mean_catch) %>%
  rename(mt_ms = mean_catch)

ms_vs_ss_walters <- ms_msy %>%
  left_join(ss_msy) %>%
  mutate(ratio = mt_ms / mt_ss)

# reorder levels
levs <- levels(factor(reorder(ms_vs_ss_walters[ms_vs_ss_walters$run=="base",]$LongName, 
                              -ms_vs_ss_walters[ms_vs_ss_walters$run=="base",]$ratio)))

ms_vs_ss_walters$LongName <- factor(ms_vs_ss_walters$LongName, levels = levs)

# rename runs
ms_vs_ss_walters <- ms_vs_ss_walters %>%
  mutate(Fishing = ifelse(run == "atf", 
                          "Arrowtooth underexploitation",
                          "MFMSY varies for all focal groups"))

# 
ms_vs_ss_walters$Fishing <- factor(ms_vs_ss_walters$Fishing,
                                   levels = c("MFMSY varies for all focal groups",
                                              "Arrowtooth underexploitation"))

# plot
walters_plot <- ms_vs_ss_walters %>%
  mutate(ratio = ifelse(run == "atf" & LongName == "Arrowtooth flounder", NA, ratio)) %>%
  #filter(LongName != "Arrowtooth flounder") %>%
  ggplot(aes(x=LongName, y=ratio)) + 
  geom_hline(yintercept = 1, color = 'grey', linetype = 'dashed') +
  geom_point(stat='identity', fill="black", size=3)  +
  geom_segment(aes(y = 1,
                   x = LongName,
                   yend = ratio,
                   xend = LongName),
               linewidth = 1) +
  theme_bw() +
  labs(x = '', y = 'Step 2 MSY / Step 1 MSY') + 
  guides(color="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Fishing, nrow = 2)
walters_plot
ggsave("results/figures/FIGURE_S8.jpeg", walters_plot, width = 5, height = 7)

# Numbers at age from nc files --------------------------------------------
# expected to decline and be fairly close to 0 for older age classes when SSB is near 0
# # get the nc files
f35_nc <- c(list.files(batch_ms_nc, pattern = ".nc", full.names = T))
# reorder these based on the number in the filename
num_idx <- as.numeric(gsub("output_([0-9]+)\\.nc", "\\1",
                           c(list.files(batch_ms_nc, pattern = ".nc", full.names = F))))
f35_nc <- f35_nc[order(num_idx)]

# function sum over depth layers in each array slice
collapse_array <- function(mat){
  mat2 <- apply(mat, 3, colSums)
  mat3 <- data.frame(t(mat2))
  colnames(mat3) <- 0:108
  mat3
}

fl <- 'data/GOA_WGS84_V4_final.bgm'
bgm <- rbgm::read_bgm(fl)
goa_sf <- rbgm::box_sf(bgm)
boundary_boxes <- goa_sf %>% sf::st_set_geometry(NULL) %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes
# function to set values in the boundary boxes to NA
setNA <- function(mat) {
  mat2 <- mat
  if(length(dim(mat2))==3) mat2[,(boundary_boxes+1),]<-NA
  if(length(dim(mat2))==2) mat2[(boundary_boxes+1),] <- NA
  mat2
}

# get t3 names, make sure you maintain the same order as t3_fg
t3_names <- grps %>% 
  filter(Code %in% t3_fg) %>% 
  mutate(Code = factor(Code, levels = t3_fg)) %>%
  arrange(Code) %>%
  pull(Name) 

extract_naa <- function(ncfile){
  
  # get run number and corresponding multiplier for F35
  this_idx <- as.numeric(gsub(".*output_", "", gsub(".nc","", ncfile)))
  this_mult <- oy_key %>% filter(idx == this_idx) %>% pull(mult)
  this_run <- oy_key %>% filter(idx == this_idx) %>% pull(run)
  
  this_ncfile <- tidync(ncfile)
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata,varid = "t") %>% as.numeric
  tyrs <- ts/(60*60*24*365)
  
  # do one fg at a time, then bring them back together
  naa_frame <- data.frame()
  for (i in 1:length(t3_names)){
    
    fg <- t3_names[i]
    
    # Get numbers by box
    abun_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    abun1 <- purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_ncdata) %>% 
      lapply(setNA) %>%
      purrr::map(apply,MARGIN=3,FUN=sum,na.rm=T) %>% 
      bind_cols() %>% 
      suppressMessages() %>% 
      set_names(abun_vars$name) %>% 
      mutate(t=tyrs)
    
    abun2 <- abun1 %>%
      pivot_longer(cols = -t,names_to = 'age_group',values_to = 'abun') %>%
      mutate(age=parse_number(age_group)) %>%
      mutate(year = ceiling(t)) %>%
      group_by(year, age_group, age) %>%
      summarise(abun = mean(abun)) %>%
      ungroup() %>%
      mutate(Name = t3_names[i]) %>%
      dplyr::select(year, Name, age, abun)
    
    # get end of the time series (last 5 years average)
    abun3 <- abun2 %>%
      slice_max(year, n = 5) %>%
      group_by(Name, age) %>%
      summarize(abun = mean(abun)) %>%
      ungroup()
    
    naa_frame <- rbind(naa_frame, abun3)
    
  }
  
  # add multiplier for the run
  naa_frame <- naa_frame %>%
    mutate(mult = this_mult,
           run = this_run)
  
  return(naa_frame)
  
}

# apply function to the nc files
naa <- bind_rows(lapply(f35_nc, extract_naa))

# bring in long names
naa <- naa %>%
  left_join(grps %>% dplyr::select(Name, LongName), by = "Name")

# spaces
naa$LongNamePlot <- gsub(" ", "\n", naa$LongName)

# add scenario information
naa <- naa %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"), 
                          "Arrowtooth\nunderexploitation",
                          "MFMSY varies for\nall focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

# reorder ATF F
naa$Fishing <- factor(naa$Fishing, 
                      levels = c("MFMSY varies for\nall focal groups",
                                 "Arrowtooth\nunderexploitation"))


# add column for age as factor
naa$`Age class` <- factor(naa$age)

# plot
grp1 <- unique(naa$LongNamePlot)[1:6]
naa_plot1 <- naa %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = mult, y = abun/1000000, color = `Age class`, linetype = Fishing))+
  geom_line()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = 'Individuals (millions)')+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
naa_plot1

grp2 <- unique(naa$LongNamePlot)[7:12]
naa_plot2 <- naa %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = mult, y = abun/1000000, color = `Age class`, linetype = Fishing))+
  geom_line()+
  geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
  scale_color_viridis_d()+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = 'Individuals (millions)')+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))
naa_plot2

# make a figure
ggsave('results/figures/FIGURE_S11_1.jpeg', naa_plot1, width = 7, height = 7)
ggsave('results/figures/FIGURE_S11_2.jpeg', naa_plot2, width = 7, height = 7)
