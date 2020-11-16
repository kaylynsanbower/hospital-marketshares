# This script fits the community detection model 
# Author: Kaylyn Sanbower; Code based on this repository: 
# https://github.com/graveja0/health-care-markets

# Date: September 14, 2020

# Import packages and run other support files  ------------------------

suppressWarnings(suppressMessages(source(here::here("/Data-Code/support/load_packages.R"))))

# source(here("Data-code/support/map_theme.R"))
# source(here("Data-code/support/move_ak_hi.R"))
# source(here("Data-code/support/get_geographic_info.R"))
# source(here("Data-code/support/get_contiguous_areas.R"))


# Import and restructure necessary data -----------------------------------

# Crosswalk of state fips to state abreviation
fips_to_state <- read_rds(here("Data/Output/county-fips-cw.rds")) %>% 
  mutate(statefp = str_sub(fips_code,1,2)) %>% 
  select(statefp,state) %>% unique()

# bring in county-fips cw with state and county fips separated
fips_xw <- read_rds(here("/Data/Output/county-fips-cw.rds")) %>%
  mutate(
    state_fips = stringr::str_sub(fips_code, 1,2),
    county_fips = stringr::str_sub(fips_code, 3, 5))

# Load the hospital-county patient sharing file
df_hosp_fips <- read_rds(here("Data/output/hospital-county-patient-data.rds"))

# Create Functions --------------------------------------------------------

# Create function to convert dataframe to bipartite matrix
convert_to_bipartite <- function(df,id) {
  id <- enquo(id)
  nn <- df %>% pull(!!id)
  foo <- df %>% select(-!!id) %>%
    as.matrix()
  rownames(foo) <- nn
  foo
}

# Create bipartite (bp) matrices -----------------------------------------------

# Contiguous county bp data frame to restrcit CMS data to only the fips_codes in the map
bp_contig_1 <- read_rds(here("Data/Output/tidy-mapping-files/df_county_info.rds")) %>% 
  tibble::as_tibble() %>% 
  mutate(fips_code = str_pad(paste0(geoid),width=5,pad="0")) %>% 
  select(fips_code, starts_with("contig_")) %>% 
  #filter(str_sub(fips_code,1,2)=="06") %>% 
  #filter(!(fips_code %in% setdiff(.$fips_code,rownames(up_fips)))) %>% 
  gather(key,fips_contig,-fips_code) %>% 
  #filter(!(fips_contig %in% setdiff(.$fips_contig,rownames(up_fips)))) %>% 
  filter(!is.na(fips_contig) & !is.na(fips_code)) %>% 
  select(fips_code,fips_contig) %>% 
  mutate(contig = 1) %>% 
  spread(fips_contig,contig)


# Set minimum number and share
# This is relevant for mapping
minimum_share = .1 # he had 0.10
minimum_number = 10 # he had 10

# create bp matrix (fips x hospitals), based on share of patients
bp_fips_hosp <-
  df_hosp_fips %>%
  group_by(fips_code) %>%
  mutate(share_of_patients = total_cases / sum(total_cases, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(connected = as.integer(share_of_patients >= minimum_share))  %>%
  mutate(share = ifelse(connected==1,share_of_patients,0)) %>% 
  select(fips_code, prvnumgrp, connected) %>%
  inner_join(bp_contig_1 %>% select(fips_code),"fips_code") %>% 
  spread(prvnumgrp, connected) %>%
  convert_to_bipartite(id = fips_code)


# Fill NA with zeros 
 bp_fips_hosp[is.na(bp_fips_hosp)] <- 0
 
# quick view of df 
 bp_fips_hosp[1:10,1:10]
 

# Create unipartite (up) matrices -----------------------------------------

# Create unipartite matrix out of fips x hosp matrix
up_fips <-bp_fips_hosp %*% t(bp_fips_hosp)
 
# Make sure that the matrix is symmetric (it should be)
isSymmetric(up_fips)

# ---- create the contiguous county unipartite matrix

# convert bp df to matrix form
bp_contig <- bp_contig_1 %>% 
  convert_to_bipartite(id = fips_code)

# replace NAs with 0 
bp_contig[is.na(bp_contig)] <- 0

up_contig <- bp_contig %*% t(bp_contig)
up_contig <- up_contig[rownames(up_fips),colnames(up_fips)]
up_contig[up_contig>0] <- 1


# Community Detection -----------------------------------------------------

# Restart
tn_contig_graph_adj_noloops <- 
  graph_from_adjacency_matrix(up_fips, weighted = TRUE) %>%
  simplify(., remove.loops = TRUE)

# Run cluster_walktrap on this network. 
initial_communities <-
  walktrap.community(tn_contig_graph_adj_noloops,
                     steps = 2,
                     merges = TRUE,
                     modularity = TRUE,
                     membership = TRUE) 

market <- membership(initial_communities)
df_walktrap <- bind_cols(fips_code = names(market), mkt = market) %>% 
  mutate(statefp = str_sub(fips_code,1,2))
df_walktrap %>% select(mkt) %>% unique() %>% 
  dim()

setwd("/Users/kaylynsanbower/Dropbox/CMS-mrktshr/Data/Output/mkt_dfs")
write.csv(df_walktrap, file = "combined-hospital-variables.csv")




# Create the igraph object from the adjacency matrix -- we'll use this graph for every CD method 
comm_graph <- 
  graph_from_adjacency_matrix(up_fips, weighted = TRUE, mode = "undirected") %>%
  simplify(., remove.loops = TRUE)


# Functions to create desired objects from CD models ----------------------

# function to run multiple CD algorithms on our graph, saves membership list results
mult_cd_methods <- function(grph) {
  mkt_walktrap <<-membership(walktrap.community(grph,steps = 2,merges = TRUE, membership = TRUE))
  mkt_fastgreedy <<- membership(fastgreedy.community(grph, merges = TRUE, membership = TRUE))
  mkt_louvain <<- membership(cluster_louvain(grph))
  mkt_infomap <<- membership(infomap.community(grph))
  mkt_labelprop <<- membership(label.propagation.community(grph))
}

# function to create a dataframe out of the list of markets 
df_from_mkts <- function(mkt_name){
  samplemarket <<- bind_cols(fips_code = names(mkt_name), mkt = mkt_name) %>% 
    mutate(statefp = str_sub(fips_code,1,2))
}


# Run functions on CD algorithm inputs ------------------------------------

# run each of the CD methods 
mult_cd_methods(comm_graph)

# create list of memberships and list of CD algorithm names
list_of_markets <- list(mkt_walktrap, mkt_fastgreedy, mkt_louvain, mkt_infomap, mkt_labelprop)
market_names <- c('walktrap', 'fastgreedy', 'louvain', 'infomap', 'labelprop')


# loop through list of markets (i.e. membership objects) and create dataframes 
n=1
for (i in list_of_markets){
  df_from_mkts(i)
  assign(paste('df_',market_names[[n]],sep=""), samplemarket)
  n = n + 1
}

# Summarize modularity and number of markets by CD algorithm.
mod_meas <- c()
for (i in list_of_markets){
  new_mod <- modularity(comm_graph, i)
  mod_meas <- append(mod_meas,new_mod)
}

number_mkts <- c()
for (i in list_of_markets){
  new_mkts <- length(unique(i))
  number_mkts <- append(number_mkts,new_mkts)
}

# Create df 
summ_of_cd_methods <- data.frame(market_names, mod_meas, number_mkts)
View(summ_of_cd_methods)


write_rds(df_louvain,path = here("Data/Output/markets-louvain-fips.rds"))


# playing around ----------------------------------------------------------


is.hierarchical(comms_fastgreedy)

plot_dendrogram(comms_fastgreedy, mode = "hclust")
is.hierarchical(comms_fastgreedy)















# Map Related Code --------------------------------------------------------

# 
# sf_walktrap <- 
#   sf::read_sf(here("Data/Input/CB/cb_2017_us_county_5m.shp")) %>% 
#   sf::st_transform(crs ="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96") %>% 
#   janitor::clean_names() %>% 
#   left_join(fips_to_state,"statefp") %>% 
#   filter(state %in% states) %>% 
#   move_ak_hi(state = state) %>% 
#   mutate(fips_code = geoid) %>% 
#   inner_join(df_walktrap, "fips_code") %>% 
#   group_by(mkt) %>% 
#   summarise() %>% 
#   ungroup() %>% 
#   st_simplify(dTolerance = 100)  %>% 
#   # left_join(get_contiguous(shp = ., id = commuting_zone_id_2010) %>% 
#   #mutate(commuting_zone_id_2010 = as.integer(commuting_zone_id_2010)), "commuting_zone_id_2010") %>% 
#   left_join(
#     df_walktrap %>% 
#       select(mkt,statefp) %>% 
#       left_join(fips_to_state,"statefp") %>% 
#       select(-statefp) %>% 
#       unique() %>% 
#       group_by(mkt) %>% 
#       mutate(n = paste0("state_",str_pad(row_number(),width = 2, pad="0"))) %>% 
#       spread(n,state)
#     ,"mkt") %>% 
#   rename(walktrap_id = mkt) %>% 
#   st_transform(crs = 4326)
# 
# sf_state <- read_sf(here("output/tidy-mapping-files/state/01_state-shape-file.shp")) %>% 
#   st_transform(crs = 4326)
# 
# states_to_map = states[-grep("HI|AK",states)] #c("TN","AL","GA","NC","VA","AR","MS","KY")
# sf_walktrap %>% 
#   filter(state_01 %in% states_to_map  | state_02 %in% states_to_map | state_03 %in% states_to_map | state_04 %in% states_to_map |
#            state_05 %in% states_to_map | state_06 %in% states_to_map | state_07 %in% states_to_map) %>% 
#   mutate(test = as.factor(sample(1:100,replace=TRUE,nrow(.)))) %>% 
#   ggplot() + geom_sf(aes(fill=test)) + theme_bw() + coord_sf(datum=NA) +
#   remove_all_axes + 
#   theme(legend.position = "none") + 
#   geom_sf(data = sf_state %>% filter(stusps %in% states_to_map), alpha = 0,lwd=.7,colour = "black") 
# 
