# Alberto Rovellini
# 01/16/2024
# This code takes output of the MFMSY permutation runs from 4 scenarios and plots:
# Production functions
# Catch and biomass curves
# Numbers at age
# biomass of forage fish and predators

library(tidyverse)
library(here)
library(tidyr)
library(readxl)
library(ggh4x)
library(viridis)
library(tidync)
library(ncdf4)

# Set up env and read data ------------------------------------------------

burnin <- 30 # years of burn-in

# identify which data we want to work on
batch_res <- "results/ms/flat_results/" # these are the biomage, catch, and mort files
batch_nc <- "results/ms/nc_results/" # these are the full out.nc files

ss_job <- "results/ss/" # this is SS runs
maxmult <- 4 # this is the full range of explored F

# set the clock to date plots
t <- format(Sys.time(),'%Y-%m-%d %H-%M-%S')

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

# list the rds files
f35_results <- list.files(file.path("results", "ms", "flat_results"), 
                                         pattern = ".rds", 
                                         full.names = TRUE)
# order them correctly
# reorder these based on the number in the filename
num_idx <- as.numeric(gsub("([0-9]+)-result\\.rds", "\\1", 
                           c(list.files(batch_res, pattern = ".rds", full.names = F))))
f35_results <- f35_results[order(num_idx)]

# get the nc files
f35_nc <- c(list.files(batch_nc, pattern = ".nc", full.names = T))
# reorder these based on the number in the filename
num_idx <- as.numeric(gsub("output_([0-9]+)\\.nc", "\\1", 
                           c(list.files(batch_nc, pattern = ".nc", full.names = F))))
f35_nc <- f35_nc[order(num_idx)]

# extract biomass and catch from the MS runs
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
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code %in% t3_fg) %>%
    left_join(fspb, by = c('Code','Age')) %>%
    mutate(biomass_mt = biomass_mt * fspb) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup()
  
  # total catch
  # taking mean of the last 5 years
  catch_vals <- catch %>% 
    slice_tail(n = 5) %>%
    summarise(across(all_of(t3_fg), ~mean(.x, na.rm = T))) %>%
    pivot_longer(cols = everything(), names_to = "Code", values_to = "catch_mt")
  
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
    #fidx = fidx) %>% # need this for joining later on
    select(Code, f)#, fidx) 
  
  # # bind all
  f_frame <- f_t1 %>%
    left_join(spawning_biomass) %>%
    left_join(catch_vals) %>%
    mutate(run = this_run,
           mult = this_mult)
  
  # add to multispecies yield list
  ms_yield_list[[i]] <- f_frame
}

ms_yield_df <- bind_rows(ms_yield_list)

# add species long name and reshape catch and biomass and add a label that this is the ms approach
ms_yield_long <- df_mult <- ms_yield_df %>% # one of these is for later plots that need mult
  left_join(grps %>% select(Code, LongName), by = "Code") %>%
  rename(Biomass = biomass_mt, Catch = catch_mt) %>%
  pivot_longer(cols = -c(Code, LongName, f, mult, run), names_to = "type", values_to = "mt") %>%
  select(Code, LongName, run, f, mult, type, mt)

# Reference points
# get two data frames: one for b0 and one for maximum yield
# b0
# For climate scenarios leave it fixed to base conditions
b0 <- ms_yield_long %>% filter(mult == 0, type == "Biomass", run == "base") %>% dplyr::select(LongName, mt) %>% rename(b0 = mt)

# max yield
ymax_ms <- ms_yield_long %>% 
  filter(type == "Catch") %>% 
  group_by(LongName, run) %>%
  slice_max(mt) %>%
  ungroup() %>%
  dplyr::select(LongName, run, mt, f) %>% 
  rename(ymax = mt) 

ymax <- ymax_ms

# handle the NaN's from FHS
ms_yield_long <- as.data.frame(ms_yield_long)
ms_yield_long$f[is.nan(ms_yield_long$f)] <- NA

# Figure 2. Single-species biomass and catch  --------------------------------------------------------------

f_files <- list.files(ss_job, full.names = T)

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
  f_df <- f_ls %>% bind_rows() %>% rename(Biomass = biomass, Catch = catch)
  
  # clean up and format
  f_df <- f_df %>%
    pivot_longer(-c(Code, f, fidx), values_to = 'mt', names_to = 'type') %>%
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
  mutate(b35 = ifelse(type== 'Biomass', mt * 0.35, NA)) %>%
  ungroup() %>%
  select(LongName, Code, type, b35)

# read in MSY information (from FMP)
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
to_plot <- unique(f_df$Code)

# bind FMSY information
fmsy <- data.frame('Code' = to_plot) %>%
  left_join(all_f) %>%
  left_join(grps %>% select(Code, LongName))

# add halibut (M from IPHC assessment)
fmsy[fmsy$Code=='HAL',]$FMSY <- 0.2 # this is M

# get f that returned the highest yield, and level of depletion for that F
sp <- unique(f_df$LongName)

atlantis_fmsy_ls <- list()

for(i in 1:length(sp)){
  
  this_f_df <- f_df %>% filter(LongName == sp[i])
  
  atlantis_fmsy <- this_f_df %>% filter(type == 'Catch') %>%
    slice_max(mt) %>%
    pull(f)
  
  b0 <- this_f_df %>%
    filter(type == 'Biomass') %>%
    slice_min(f) %>%
    pull(mt)
  
  b_fmsy <- this_f_df %>%
    filter(f == atlantis_fmsy, type == 'Biomass') %>%
    pull(mt)
  
  depletion_fmsy <- b_fmsy / b0
  
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
# write.csv(atlantis_fmsy, "NOAA_Azure/data/f35_vector_PROXY_OY_SS.csv", row.names = F)

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
  ggplot(aes(x = f, y = mt/1000))+
  geom_line()+
  geom_point(size = 1.5)+
  geom_vline(data = fmsy %>% filter(LongName %in% key_grps), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongName %in% key_grps), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = atlantis_fmsy %>% filter(LongName %in% key_grps) %>% mutate(type = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongName %in% key_grps) %>% mutate(type = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1.5,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
ggsave(paste0('results/figures/biom_catch_key_stocks.png'), p_ms, width = 7, height = 7)

# make figures with all stocks for supplement
grp1 <- unique(f_df_ms$LongNamePlot)[1:6]
f_plot1 <- f_df_ms %>%
  filter(LongNamePlot %in% grp1) %>%
  ggplot(aes(x = f, y = mt/1000))+
  geom_line()+
  geom_point(size = 2)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp1), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp1) %>% mutate(type = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongNamePlot %in% grp1) %>% mutate(type = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1.5,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

grp2 <- unique(f_df_ms$LongNamePlot)[7:12]
f_plot2 <- f_df_ms %>%
  filter(LongNamePlot %in% grp2) %>%
  ggplot(aes(x = f, y = mt/1000))+
  geom_line()+
  geom_point(size = 2)+
  geom_vline(data = fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = FMSY, group = LongNamePlot), linetype = 'dashed', color = 'orange')+
  geom_vline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2), aes(xintercept = atlantis_fmsy, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_hline(data = atlantis_fmsy %>% filter(LongNamePlot %in% grp2) %>% mutate(type = 'Biomass'),
             aes(yintercept = b_fmsy/1000, group = LongNamePlot), linetype = 'dashed', color = 'blue')+
  geom_text(data = annotations %>% filter(LongNamePlot %in% grp2) %>% mutate(type = 'Biomass'),
            aes(x=Inf,y=Inf,hjust=1,vjust=1.5,label=paste0('Depletion=',depletion)), color = 'blue')+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = '1000\'s of tons')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))

ggsave(paste0('results/figures/yield_curves',t,'_OY_1.png'), f_plot1, width = 7, height = 7)
ggsave(paste0('results/figures/yield_curves',t,'_OY_2.png'), f_plot2, width = 7, height = 7)

# Figure 3. Global yield ---------------------------------------------------------------
# list
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
    slice_tail(n = 5) %>%
    summarise(across(all_of(t3_fg), ~mean(.x, na.rm = T))) %>%
    mutate(mult = this_mult,
           run = this_run,
           idx = this_idx)
  
  catch_list[[i]] <- this_catch
  
}

catch_df <- bind_rows(catch_list)

# reshape and calculate total
catch_df_long <- catch_df %>%
  pivot_longer(-c(run, mult, idx), names_to = "Code", values_to = "mt") %>%
  filter(Code != "HAL") %>%
  group_by(run, mult) %>%
  mutate(total_yield = sum(mt),
         prop = mt / total_yield) %>%
  ungroup() %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# apply scaling of catch by biomass in AK, so that we only have AK catch now
catch_scalars <- read.csv("data/catch_scalars.csv")

catch_df_long_ak <- catch_df_long %>%
  left_join(catch_scalars %>% 
              left_join(grps %>% 
                          dplyr::select(Code, Name))) %>%
  mutate(mt_ak = mt * ak_prop)

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
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Catch (1000 mt)", fill = "") +
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(nrow = 4))+
  facet_grid(Climate~`F on\narrowtooth`)
global_yield_ms

ggsave(paste0("results/figures/global_yield_ms_AK.png"), global_yield_ms, width = 6, height = 5)

# make a table with max catch per scenario for the report
max_catch <- catch_df_long_ak %>%
  select(run, idx, mt_ak) %>%
  group_by(run, idx) %>%
  summarize(total_yield_ak = sum(mt_ak, na.rm = T)) %>%
  select(run, total_yield_ak) %>%
  distinct() %>%
  group_by(run) %>%
  slice_max(total_yield_ak)

# Figure 4. Biomass and catch curves ------------------------------------------------
# plot catch and biomass curves
to_plot <- ms_yield_long

# spaces
to_plot$LongNamePlot <- gsub(" ", "\n", to_plot$LongName)
ymax$LongNamePlot <- gsub(" ", "\n", ymax$LongName)
ymax$type <- "Catch"

# rename scenarios and order them
to_plot <- to_plot %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                          "Arrowtooth\nunderexploitation",
                          "MFMSY varies for\nall focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

ymax <- ymax %>%
  mutate(Fishing = ifelse(run %in% c("atf","atf_climate"),
                          "Arrowtooth\nunderexploitation",
                          "MFMSY varies for\nall focal groups"),
         Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))

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
  filter(LongName %in% key_grps, type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  # geom_point(size = 1.6)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% 
               filter(LongName %in% key_grps) %>% 
               filter(!(LongNamePlot == "Arrowtooth\nflounder" & Fishing == "Arrowtooth\nunderexploitation")), 
             aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  #facet_wrap(~ LongName, scales = "free", ncol = 1)+
  theme(strip.text.y = element_text(angle=0))
f_plot_ms

ggsave(paste0('results/figures/catch',t,'_MS_ms.png'), f_plot_ms, width = 6, height = 6)

# POL and COD for AMSS
# amss_plot <- to_plot %>%
#   filter(Code %in% c("POL","COD"), type == "Catch") %>%
#   ggplot(aes(x = f, y = mt/1000, color = Climate, linetype = Fishing))+
#   geom_line(linewidth = 1)+
#   # geom_point(size = 1.6)+
#   scale_color_viridis_d(begin = 0.2, end = 0.8)+
#   geom_vline(data = ymax %>% 
#                filter(LongName %in% c("Walleye pollock","Pacific cod")) %>% 
#                filter(!(LongNamePlot == "Arrowtooth\nflounder" & Fishing == "Arrowtooth\nunderexploitation")), 
#              aes(xintercept = f, color = Climate, linetype = Fishing))+
#   theme_bw()+
#   scale_y_continuous(limits = c(0, NA))+
#   labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
#   facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
#   #facet_wrap(~ LongName, scales = "free", ncol = 1)+
#   theme(strip.text.y = element_text(angle=0))
# amss_plot
# ggsave("amss_catch.png",amss_plot,width=5,height = 5.5)

# make figures for supplement (break into two sets)
grp1 <- unique(to_plot$LongNamePlot)[1:6]
f_plot1 <- to_plot %>%
  filter(LongNamePlot %in% grp1) %>%
  filter(type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% filter(LongNamePlot %in% grp1), aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot1

grp2 <- unique(to_plot$LongNamePlot)[7:12]
f_plot2 <- to_plot %>%
  filter(LongNamePlot %in% grp2) %>%
  filter(type == "Catch") %>%
  ggplot(aes(x = f, y = mt/1000, color = Climate, linetype = Fishing))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(begin = 0.2, end = 0.8)+
  geom_vline(data = ymax %>% filter(LongNamePlot %in% grp2), aes(xintercept = f, color = Climate, linetype = Fishing))+
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))+
  labs(x = 'Fishing mortality (F)', y = 'Catch (1000 mt)')+
  facet_grid2(LongNamePlot~type, scales = 'free', independent = 'all')+
  theme(strip.text.y = element_text(angle=0))
f_plot2

ggsave(paste0('results/figures/biomass_catch',t,'_MS_1.png'), f_plot1, width = 7.5, height = 7)
ggsave(paste0('results/figures/biomass_catch',t,'_MS_2.png'), f_plot2, width = 7.5, height = 7)

# Figures 5 and 6: top predators and forage fish ------------------------------

top_preds <- c("SSL","PIN","DOL","BDF","BSF")
forage <- c("CAP","SAN","HER","EUL","FOS")
other_fg <- c(top_preds, forage)

ms_other_list <- list()

for(i in 1:length(f35_results)){
  
  print(paste("Doing", f35_results[i]))
  
  # grab the index from the file name
  this_idx <- as.numeric(gsub("-result.rds", "", gsub("results/ms/flat_results//", "", f35_results[i])))
  
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
    slice_tail(n = 5) %>% # use last xxx years
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    # separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code %in% other_fg) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup() %>%
    mutate(run = this_run,
           mult = this_mult)
  
  # add to multispecies yield list
  ms_other_list[[i]] <- other_biomass
}

ms_other_df <- bind_rows(ms_other_list) %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# get b0
# static to "base"
b0_other <- ms_other_df %>% filter(mult == 0, run == "base") %>% dplyr::select(LongName, biomass_mt) %>% rename(b0 = biomass_mt)

ms_other_df <- ms_other_df %>%
  left_join(b0_other, by = c("LongName")) %>%
  mutate(biomchange = (biomass_mt - b0)/b0 * 100)

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


# for each predator, identify the main prey species from dietcheck (in baseline)
# sum up total prey biomass
# express changes from B0
# could do the same for prey (but then pred species would be a lot, everyone eats CAP, but then changes should concern the most abundant ones)
# This is qualitative, but it demonstrate a likely trophic link and its effects

# for now, "baseline" is run 3, but replace with real base run when it's done (close enough)
base_diet <- read.table("data/output_1556DietCheck.txt", sep = " ", header = T)

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
  # Create an inner list of length X
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
  
  # bring in codes again as that's what the output works with...
  diet_codes <- grps %>% filter(Name %in% diet_grps) %>% pull(Code)
  
  # now loop over reuslts
  for(j in 1:length(f35_results)){
    
    print(paste("Doing", f35_results[j]))
    
    # grab the index from the file name
    this_idx <- as.numeric(gsub("-result.rds", "", gsub("results/ms/flat_results//", "", f35_results[j])))
    
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
    # SSB to plot and report in tables
    this_diet_biomass <- biomage %>% 
      slice_tail(n = 5) %>% # use last xxx years
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

# now need to rescale to b0, where b0 is for each scenario
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
other_plot_top_diets <- ms_other_df_diet %>%
  filter(Code %in% c("DOL","SSL","PIN","BDF","BSF")) %>%
  ggplot(aes(x = mult, y = biomass_mt / 1000, fill = scalar, shape = Fishing))+
  geom_point(color = "black", size = 1.5)+
  scale_shape_manual(values = c(21,24))+
  colorspace::scale_fill_continuous_divergingx(palette = 'PRGn', mid = 0) + 
  geom_vline(xintercept = 1, color = 'black', linetype = "dashed", linewidth = 0.35)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Biomass (1000 mt)", fill = "Change in total prey\nbiomass from unfished (%)")+
  #guides(fill=guide_legend(order=1), shape=guide_legend(order=2))+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))

other_plot_forage_diets <- ms_other_df_diet %>%
  filter(Code %in% c("CAP","SAN","HER","FOS","EUL")) %>%
  ggplot(aes(x = mult, y = biomass_mt / 1000, fill = scalar, shape = Fishing))+
  geom_point(color = "black", size = 1.5)+
  scale_shape_manual(values = c(21,24))+
  colorspace::scale_fill_continuous_divergingx(palette = 'PRGn', mid = 0) + 
  geom_vline(xintercept = 1, color = 'black', linetype = "dashed", linewidth = 0.35)+
  theme_bw()+
  labs(x = expression(MF[MSY] ~ "multiplier"), y = "Biomass (1000 mt)", fill = "Change in total predator\nbiomass from unfished (%)")+
  #guides(fill=guide_legend(order=1), shape=guide_legend(order=2))+
  facet_grid2(LongNamePlot~Climate, scales = 'free')+
  theme(strip.text.y = element_text(angle=0))

ggsave(paste0("results/figures/other_top_diets.png"), other_plot_top_diets, width = 7, height = 4.05)
ggsave(paste0("results/figures/other_forage_diets.png"), other_plot_forage_diets, width = 7, height = 4.05)

#########################
# SUPPLEMENTARY FIGURES #
#########################

# S1.1. Harvest sepcifications --------------------------------------------
grps <- read.csv("data/GOA_Groups.csv")

specs <- read_excel("data/GOA_harvest specs_1986-2024.xlsx", 
                    sheet = 1,
                    na = "n/a",
                    n_max = 131)

# lots of cleaning to do
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
  filter(Year > 1990, Area == "Total") %>%
  mutate(mt = as.numeric(mt))

species <- sort(unique(specs_long$Species))
key <- data.frame("Species" = species,
                  "Tier" = c(3,5,3,3,3,5,3,4,4,3,3,3,3,3,3,3,3,3,3),
                  "Code" = c("ATF","SKB","FFD","RFP","FHS","SKL","RFS","FFS","RFD","COD","POP","RFP","POL","REX","RFS","SBF","FFS","RFS","RFS"))

specs_long <- specs_long %>%
  left_join(key, by = "Species") %>%
  left_join(grps %>% select(Code, LongName), by = "Code")

# order factors
specs_long$Spec <- factor(specs_long$Spec, levels = c("OFL", "ABC", "TAC"))


# make a bar chart
harvest_specs_fig <- specs_long %>%
  filter(Tier == 3) %>%
  group_by(Year, Spec, LongName) %>%
  summarise(mt = sum(mt, na.rm = T)) %>%
  ggplot(aes(x = Year, y = mt/1000, fill = LongName))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_viridis_d()+
  geom_hline(yintercept = 800, linetype = "dashed", color = "red")+
  theme_bw()+
  scale_x_discrete(breaks = seq(1992,2024,2))+
  labs(x = "", y = "1000 mt", fill = "Stock")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(nrow = 4))+
  facet_grid(~Spec)

ggsave("results/figures/harvest_specs_S1.png", harvest_specs_fig, width = 8, height = 4)


# S1.3. Arrowtooth flounder diets -----------------------------------------

library(RColorBrewer)
diet_runs <- data.frame("idx" = c(1556,13,25,39,51),
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
    filter(Time > 75 & Time <=80) %>% # <= (ceiling(max(Time))-5)) %>% # focus on last 5 years of the run
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

pred_names <- grps %>% filter(Code %in% c("ATF", top_preds)) %>% pull(Name)

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
  filter(run == "Base model,\ncalibration fishing") %>%
  drop_na()

colourCount <- length(unique(atf_diet$Prey_Name))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

p_atf_diet <- atf_diet %>%
  ggplot(aes(x = Cohort+1, y = Prop * 100, fill = Prey_LongName))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_x_continuous(breaks = 1:10)+
  scale_fill_manual(values = getPalette(colourCount))+
  theme_bw()+
  labs(x = '', y = "Diet preference (%)", fill = "Prey")
p_atf_diet

ggsave("results/figures/diet_plots/ATF_diet_S3.png", p_atf_diet, width = 6, height = 6)


# plot all together for Figure S1.10
# drop ATF here
diet_long_preds <- diet_long_preds %>%
  filter(Predator_Name != "Arrowtooth_flounder")

# pred longname for facets
diet_long_preds$Predator_LongNamePlot <- gsub(" ", "\n", diet_long_preds$Predator_LongName)

colourCount <- length(unique(diet_long_preds$Prey_Name))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

p_all_diet <- diet_long_preds %>%
  ggplot(aes(x = Cohort+1, y = Prop * 100, fill = Prey_LongName))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_x_continuous(breaks = 1:10)+
  scale_fill_manual(values = getPalette(colourCount))+
  theme_bw()+
  labs(x = 'Age class', y = "Diet composition (%)", fill = "Prey")+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'))+
  guides(fill = guide_legend(nrow = 4))+
  facet_grid2(Predator_LongNamePlot~run)+
  theme(strip.text.y = element_text(angle=0))
p_all_diet

ggsave("results/figures/diet_plots/all_S10.png", p_all_diet, width = 8.5, height = 6)


# Production curves S1.5 --------------------------------------------------
ss_yield_long <- f_df %>%
  select(Code, LongName, f, fidx, type, mt)

# get b0 from SS curves
b0 <- ss_yield_long %>% 
  group_by(Code) %>%
  slice_min(f) %>% 
  ungroup() %>%
  filter(type == "Biomass") %>% 
  dplyr::select(LongName, mt) %>% 
  rename(b0 = mt)

# get max yield
ymax <- ss_yield_long %>% 
  filter(type == "Catch") %>% 
  group_by(LongName) %>%
  slice_max(mt) %>%
  ungroup() %>%
  dplyr::select(LongName, mt, f) %>% 
  rename(ymax = mt) 

# we are plotting yield fraction against depletion
yield_func <- ss_yield_long %>%
  drop_na() %>%
  dplyr::select(LongName, type, f, mt) %>%
  pivot_wider(id_cols = c(LongName, f), names_from = type, values_from = mt) %>%
  left_join(b0, by = c("LongName")) %>% # if you are keep static reference point
  mutate(depletion = Biomass / b0) %>%
  left_join(ymax %>% select(-f), by = c("LongName")) %>%
  mutate(yfrac = Catch / ymax) %>%
  #mutate(experiment = ifelse(experiment == "ms", "Multispecies", "Single-species")) %>%
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

ggsave(paste0("results/figures/yield_functions_S1.5.png"), yield_func_plot, width = 8, height = 6.5)


# S1.6, stocks below B35% -------------------------------------------------
# How many stocks are below 35% B0 for each scenario?
# Use static B0 from Base Scenario for this
b0 <- ms_yield_long %>% filter(mult == 0, type == "Biomass", run == "base") %>% dplyr::select(LongName, mt) %>% rename(b0 = mt)

below_target <- ms_yield_long %>%
  filter(type == "Biomass") %>%
  filter(!(run %in% c("atf", "atf_climate") & LongName == "Arrowtooth flounder")) %>%
  left_join(b0, by = "LongName") %>%
  mutate(depletion = mt / b0) %>%
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
  facet_grid(Climate~`F on\narrowtooth`)
p_below_target

ggsave(paste0("results/figures/below_target_S1.6.png"), p_below_target, width = 6, height = 4)

# S1.7. Walters plot ------------------------------------------------------
# create a plot akin to Walters et al. (2005) Fig. 3, except we do not organize it by TL for now
# treat this as model-wide MSY
ss_msy <- f_df %>%
  mutate(experiment = "ss") %>%
  select(Code, LongName, f, fidx, experiment, type, mt) %>%
  filter(type == "Catch") %>%
  group_by(Code, LongName) %>%
  slice_max(mt) %>%
  ungroup() %>%
  select(LongName, mt) %>%
  rename(mt_ss = mt)

# this is the scenario where ATF is kept at low fishing pressure, all other groups are varied
ms_msy <- catch_df %>%
  pivot_longer(-c(run, mult, idx), names_to = "Code", values_to = "mt") %>%
  group_by(run, mult) %>%
  mutate(total_yield = sum(mt),
         prop = mt / total_yield) %>%
  ungroup() %>%
  left_join(grps %>% select(Code, LongName), by = "Code") %>%
  filter(mult == 1, run %in% c("base","atf")) %>%
  select(LongName, run, mt) %>%
  rename(mt_ms = mt)

ms_vs_ss_walters <- ms_msy %>%
  left_join(ss_msy) %>%
  mutate(ratio = mt_ms / mt_ss)

# reorder levels
#ms_vs_ss_walters$LongNamePlot <- gsub(" - ", "\n", ms_vs_ss_walters$LongName)
# ms_vs_ss_walters$LongNamePlot <- reorder(ms_vs_ss_walters$LongNamePlot, -ms_vs_ss_walters$ratio)
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

# Comparing SS to (either) MS scenario:
# Across all groups, MS MSY is higher than SS MSY.
# Biggest differences for groups that are more top-down controlled (pollock, FHS, Cod, FFS)
# Groups that have smallest difference are higher trophic levels or groups that are parameterized to be less predated upon
# It makes sense though that there will always be less predators in the MS runs
# This shows that it varies by species, but there is a strong top-down control in the system

ggsave("results/figures/walters_plot_S1.7.png", walters_plot, width = 5, height = 7)

# SXX Numbers at age from nc files --------------------------------------------

# expected to decline and be fairly close to 0 for older age classes when SSB is near 0
# should that not be the case, there is an iddues with how we count biomass

# function sum over depth layers in each array slice
collapse_array <- function(mat){
  mat2 <- apply(mat, 3, colSums)
  mat3 <- data.frame(t(mat2))
  colnames(mat3) <- 0:108
  mat3
}

fl <- 'NOAA_Azure/data/GOA_WGS84_V4_final.bgm'
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
  pull(Name) #%>%
#sort()

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
naa <- bind_rows(lapply(f35_nc, extract_naa)) # this is slow with 40

# bring in long names
naa <- naa %>%
  left_join(grps %>% select(Name, LongName), by = "Name")

# spaces
naa$LongNamePlot <- gsub(" ", "\n", naa$LongName)

# add scenario information
# pretty verbose for WFC now
# naa <- naa %>%
#   mutate(Fishing = ifelse(run %in% c("atf","atf_climate"), 
#                                      "1/4 FOFL on arrowtooth flounder,\nMFMSY varying for all other stocks", 
#                                      "MFMSY varying for all stocks"),
#          Climate = ifelse(run %in% c("climate","atf_climate"), "ssp585 (2075-2085)", "Historical (1999)"))
# 
# # reorder ATF F
# naa$Fishing <- factor(naa$Fishing, 
#                                           levels = c("MFMSY varying for all stocks",
#                                                      "1/4 FOFL on arrowtooth flounder,\nMFMSY varying for all other stocks"))

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
# ggsave(paste0('NOAA_Azure/results/figures/oy_paper/naa',t,'_1.png'), naa_plot1, width = 7, height = 7)
# ggsave(paste0('NOAA_Azure/results/figures/oy_paper/naa',t,'_2.png'), naa_plot2, width = 7, height = 7)

# key groups
# naa_plot_ms <- naa %>%
#   filter(LongName %in% key_grps) %>%
#   ggplot(aes(x = mult, y = abun/1000000, color = `Age class`, linetype = Fishing))+
#   geom_line()+
#   geom_vline(xintercept = 1, color = "black", linetype = "dotted")+
#   scale_color_viridis_d()+
#   theme_bw()+
#   labs(x = expression(F[MSY] ~ "multiplier"), y = 'Individuals (millions)')+
#   facet_grid2(LongNamePlot~Climate, scales = 'free')+
#   theme(strip.text.y = element_text(angle=0))
# naa_plot_ms

# ggsave(paste0("NOAA_Azure/results/figures/oy_paper/NAA_ms.png"), naa_plot_ms, width = 8, height = 6.5)
