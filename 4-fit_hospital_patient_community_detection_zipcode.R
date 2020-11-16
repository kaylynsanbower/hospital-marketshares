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

# Load the hospital-county patient sharing file -- zip code level
df_hosp_zip <- read_rds(here("Data/output/hospital-county-patient-data-zip.rds"))

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

# Contiguous county bp data frame to restrcit CMS data to only the zip-codes in the map
bp_contig_1 <- read_rds(here("Data/Output/tidy-mapping-files/df_zip_info.rds")) %>% 
  tibble::as_tibble() %>% 
  mutate(zip_code = str_pad(paste0(geoid10),width=5,pad="0")) %>% 
  select(zip_code, starts_with("contig_")) %>% 
  gather(key,zip_contig,-zip_code) %>% 
  filter(!is.na(zip_contig) & !is.na(zip_code)) %>% 
  select(zip_code,zip_contig) %>% 
  mutate(contig = 1) %>% 
  spread(zip_contig,contig)


bp_contig <- bp_contig_1 %>% 
  convert_to_bipartite(id = zip_code)

# replace NAs with 0 
bp_contig[is.na(bp_contig)] <- 0

# Set minimum number and share
# This is relevant for mapping
minimum_share = 0 # he had 0.10
minimum_number = 0 # he had 10

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

# create the contiguous county unipartite matrix
up_contig <- bp_contig %*% t(bp_contig)
up_contig <- up_contig[rownames(up_fips),colnames(up_fips)]
up_contig[up_contig>0] <- 1



# Community Detection -----------------------------------------------------

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

mkt_infomap <- membership(infomap.community(comm_graph))


# function to create a dataframe out of the list of markets 
df_from_mkts <- function(mkt_name){
  samplemarket <<- bind_cols(fips_code = names(mkt_name), mkt = mkt_name) %>% 
    mutate(statefp = str_sub(fips_code,1,2))
}

df_from_mkts(mkt_infomap)
assign(df_infomap, samplemarket)



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




# playing around ----------------------------------------------------------


is.hierarchical(comms_fastgreedy)

plot_dendrogram(comms_fastgreedy, mode = "hclust")


is.hierarchical(comms_fastgreedy)





##### Might need these thigns? tbd tomorrow 
new_names <- list(df_walktrap, df_fastgreedy, df_louvain, df_infomap, df_labelprop)

df_list <- list(df_1, df_2, df_3, df_4, df_5)













# Map Related Code --------------------------------------------------------
# df_walktrap %>% filter(grepl("^06",fips_code)) %>% select(mkt) %>% unique() %>% 
dim()
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
