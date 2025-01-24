# Alberto Rovellini
# 01/22/2025
# this code pulls terminal biomass and catch from each run, and it calculates realized F based on catch in the last 5 years
# these variables are stored in a table with one row and written to a csv file, to be processed later

library(dplyr)
library(tidyr)
library(here)

outdir <- "results/ss/processed/"
dir.create(outdir)

burnin <- 30 # years of burn-in

folder_path <- "results/ss/raw/"
# list the results files
results_list <- list.files(folder_path, full.names = T)

# list the data
grp_path <- here('data/GOA_Groups.csv') # functional groups
fspb_path <- here('data/fspb.csv') # proportion mature at age
selex_path <- here('data/age_at_selex_new.csv') # selectivity pattern, after adjusting startage
lookup_path <- here('data/f_lookup_OY_SS.csv')

atlantis_fg <- read.csv(grp_path)
fspb <- read.csv(fspb_path, header = T) # proportion mature
# reshape fspb
fspb <- fspb %>%
  pivot_longer(-Code, names_to = "Age", values_to = "fspb") %>%
  mutate(Age = gsub("X","",Age))

selex <- read.csv(selex_path, header = T) # age at selectivity
f_lookup <- read.csv(lookup_path) # lookup for species and F per run

for(i in 1:length(results_list)){
  
  print(paste("Doing",i))
  
  # load results
  this_result <- readRDS(results_list[i])
  
  # this idx - CAREFUL! This is not the same as i because of alphabetical sorting of the reuslts files
  this_idx <- as.numeric(names(this_result))
  
  # first identify the species and the level of fishing for this run. These are unrelated from runname
  sp <- f_lookup %>% filter(idx == this_idx) %>% pull(species)
  fidx <- f_lookup %>% filter(idx == this_idx) %>% pull(mult_idx)
  
  # rename the result object in the list to avoid problems with indexing
  names(this_result) <- "res"
  
  # extract tables from results
  biomage <- this_result$res[[paste0("biomage_",this_idx)]]
  catch <- this_result$res[[paste0("catch_",this_idx)]]
  
  # now extract data
  # SSB to plot and report in tables
  spawning_biomass <- biomage %>% 
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate(Code.Age, into = c('Code', 'Age'), sep = '\\.') %>%
    filter(Code == sp) %>%
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
    dplyr::select(c(Time, all_of(sp))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'catch_mt') %>%
    slice_max(Time, n = 5) %>%
    group_by(Code) %>%
    summarise(mean_catch = mean(catch_mt),
              catch_cv = sd(catch_mt) / mean(catch_mt))
  
  # # calculate realized F after 1 year of data
  # For runs with a burn-in, this has to be the biomass at the end of the burn-in, when we start fishing with the new scalar
  # # get initial biomass for the selected age classes
  biom_age_t1 <- biomage %>% 
    filter(Time == 365 * burnin) %>%# this is the burn-in years
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class_selex)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass = sum(biomass)) %>%
    ungroup() %>% 
    filter(Code == sp)
  # 
  # # catch (one time step after biomass: how much did we catch in this time?)
  catch_t1 <- catch %>% 
    select(Time, all_of(sp)) %>% 
    filter(Time == 365 * (burnin + 1)) %>% # careful - there is a small transition phase
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
    select(-Time)
  # 
  # # calc realized f
  f_t1 <- biom_age_t1 %>% left_join(catch_t1) %>%
    mutate(exp_rate = catch/biomass,
           f = -log(1-exp_rate),
           fidx = fidx) %>% # need this for joining later on
    select(Code, f, fidx) 
  
  # # bind all
  f_frame <- f_t1 %>%
    left_join(spawning_biomass) %>%
    left_join(catch_vals)
  
  # write out to be then brought together with all other runs
  write.csv(f_frame, paste(outdir,paste(sp,fidx,'f.csv',sep='_'), sep = "/"), row.names = F)
  
}
