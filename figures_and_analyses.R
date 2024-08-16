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

# read in lookup tables
oy_key <- read.csv("data/oy_key.csv")
f_lookup <- read.csv("data/f_lookup_OY_SS.csv")

# list Tier 3 stocks 
t3_fg <- f_lookup %>% pull(species) %>% unique() %>% sort()
t3_names <- grps %>% filter(Code %in% t3_fg) %>% pull(Name) # names for nc files pulling

# list the rds files
f35_results <- c(list.files(batch_res, pattern = ".rds", full.names = T))

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

# make figures for supplement
# plot
f_df_ms <- f_df

f_df_ms$LongNamePlot <- gsub(" ", "\n", f_df_ms$LongName)
fmsy$LongNamePlot <- gsub(" ", "\n", fmsy$LongName)
atlantis_fmsy$LongNamePlot <- gsub(" ", "\n", atlantis_fmsy$LongName)
b35$LongNamePlot <- gsub(" ", "\n", b35$LongName)
annotations$LongNamePlot <- gsub(" ", "\n", annotations$LongName)

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

######################################################################
