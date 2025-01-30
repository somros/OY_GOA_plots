# Alberto Rovellini
# 3/15/2024
# Function to pull biomass proportion in AK over the model domain from out.nc files
# This is needed because we want to plot global T3 groundfish yield against the OY cap, but the OY cap is AK only
# We do not have spatial output for catch (because we save only some output from the large batch runs)

# Why does this work? Because mFC extracts catch proportionally to biomass and there are no further limitations to fishing
# See scripts comparing_biom_catch_spatial.R and check_catch_outputs.R for details on why this works
# Whn we diverge from this simple setup of no fleets, no spatial closures, and mFC, then this will no longer apply

# This function will:
# open an out.nc file for each of the MS runs
# pull WAA and NAA from nc file
# get biomass per cell
# collapse over depth
# filter to only selected age classes
# get proportions per box over the total for the selected biomass
# sum up over AK boxes
# average for the last 5 years (for consistency with how catch is being handled)
# return one scalar per group per run

get_catch_ak_scalar <- function(nc_file){
  
  print(paste("Doing", nc_file))
  
  this_idx <- as.numeric(gsub(".nc","", gsub("results.*output_", "", nc_file)))
  
  out <- tidync(nc_file)
  this.nc <- ncdf4::nc_open(nc_file)
  
  # make an empty list to populate with biomass for each group
  catch_prop_ls <- list()
  
  for(n in 1:length(fmp_names)){
    
    print(paste("Doing", fmp_names[n]))
    
    # args for the function below
    fg <- fmp_names[n] # this needs to use the "Name" to pull from the NC file
    
    #Extract from the output .nc file the appropriate reserve N time series variables
    resN_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
      filter(grepl("_ResN",name)) %>% # filter for reserve N
      filter(grepl(fg,name)) # filter for specific functional group
    
    #Extract from the output .nc file the appropriate structural N time series variables
    strucN_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
      filter(grepl("_StructN",name)) %>% # filter for structural N
      filter(grepl(fg,name)) # filter for specific functional group
    
    # Get numbers by box
    abun_vars <- hyper_vars(out) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    if(nrow(resN_vars)==0) {return("no data.")}
    else {
      # # Actually pull the data from the .nc
      # here we can collapse the depth layers but need to keep the boxes
      resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this.nc) 
      strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this.nc)
      nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this.nc) #numbers by age group,box,layer,time
      
      # Formula for biomass in each cell will be (resN+structN)*nums
      # this step maintains the spatial and temporal structure of the nc arrays
      biom <- list()
      for(i in 1:length(resN)){
        biom[[i]] <- (resN[[i]] + strucN[[i]]) * nums[[i]] * 20 * 5.7 / 1000000000 # go from mgN to tons
      }
      
      # now collapse over depth
      biom_box <- lapply(biom, function(x) apply(x, c(2, 3), sum))
      
      # now turn to data frame
      # create empty list
      biom_box_ls <- list()
      
      # Loop over each matrix in the collapsed list to fill the data frame (i.e. over each age class)
      for(i in 1:length(biom_box)) {
        # Get the current matrix
        mat <- biom_box[[i]]
        
        # Convert matrix to a data frame
        mat_df <- as.data.frame((mat))
        
        # add box_id
        mat_df <- mat_df %>% mutate(box_id = 0:108)
        
        # reshape
        mat_df_long <- mat_df %>%
          pivot_longer(-box_id, names_to = "ts", values_to = "mt")
        
        # turn ts column to integer
        mat_df_long <- mat_df_long %>%
          mutate(ts = gsub("V","",ts)) %>%
          mutate(ts = as.numeric(ts)) %>%
          mutate(ts = ts - 1) # start numbering ts from 0
        
        # add age
        mat_df_long <- mat_df_long %>%
          mutate(age = i-1) # number from 0 for consistency with age at selex and age mat
        
        biom_box_ls[[i]] <- mat_df_long
        
      }
      
      # turn list to data frame
      biom_box_df <- bind_rows(biom_box_ls)
      
      # to avoid immense tables when we do it for all species (and once per run), work out the proportions here
      # do it by age class first
      # this will be used to scale the biomass of each age class
      # I am not even convinced we need this - if we are keeping everything to total population level except for the global yield plot
      # 92 is the first box in BC
      # biom_ak_prop <- biom_box_df %>%
      #   mutate(area = ifelse(box_id < 92, "ak", "bc")) %>%
      #   group_by(ts, age) %>% # get total biomass by time step by age across the model domain
      #   mutate(mt_tot = sum(mt)) %>%
      #   group_by(ts, age, area) %>%
      #   mutate(mt_area = sum(mt)) %>%
      #   ungroup() %>%
      #   mutate(prop = mt_area / mt_tot) %>%
      #   select(ts, age, area, prop) %>%
      #   distinct()
      
      # view:
      # biom_ak_prop %>%
      #   filter(area == "ak") %>%
      #   ggplot(aes(x = ts, y = prop, color = factor(age))) +
      #   geom_line()
      
      # now bring in selex, filter at or above it, sum, and get prop for total
      catch_prop_from_ak <- biom_box_df %>%
        mutate(area = ifelse(box_id < 92, "ak", "bc")) %>% # define areas based on boxes
        mutate(Name = fg) %>% # add name
        left_join(grps %>% dplyr::select(Code, Name), by = "Name") %>% 
        left_join(selex, by = "Code") %>% # bring in age at selectivity (which counts from 0)
        mutate(idx = age - age_class_selex) %>% # is age >= age at selectivity?
        filter(idx >= 0) %>% # keep only age classes >= age at selex
        group_by(ts) %>% # get total biomass by time step across the model domain (aggregate age classes)
        mutate(mt_tot = sum(mt)) %>%
        group_by(ts, area) %>%
        mutate(mt_area = sum(mt)) %>%
        ungroup() %>%
        mutate(prop = mt_area / mt_tot) %>%
        select(ts, Name, area, prop) %>%
        distinct()
      
      # handle time:
      # drop t0 for easier filtering
      # subsample at every 5 time steps to have annual values
      # keep the last 5 of the series only
      # average
      # produce one value per run per species
      
      # NOTE: the MS batches have annual output in the biology file too, so no need to re-index
      catch_prop_from_ak <- catch_prop_from_ak %>%
        filter(area == "ak") %>% # keep proportion in AK only
        filter(ts > 0) %>% # drop t0
        # filter(ts %in% seq(5,250,5)) %>% # resample to have annual time steps instead of 73 days
        # mutate(ts = ts / 5) %>% # reindex the time step accordingly 
        slice_tail(n = 5) %>% # keep end of the run for consistency with the catch sampling
        group_by(Name) %>%
        summarize(ak_prop = mean(prop))
      
    }
    
    catch_prop_ls[[n]] <- catch_prop_from_ak
    
  }
  
  # now bind all species
  catch_prop_df <- bind_rows(catch_prop_ls)
  
  # add run
  catch_prop_df <- catch_prop_df %>% mutate(idx = this_idx) 
  
  # output this data frame from each run
  return(catch_prop_df)
  
}
