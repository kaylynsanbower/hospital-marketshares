# This script constructs county map data
# Author: Kaylyn Sanbower; Code based on this repository: 
# https://github.com/graveja0/health-care-markets

# Date: September 9, 2020


# Import packages and run other support files  ------------------------

suppressWarnings(suppressMessages(source(here::here("/Data-Code/support/load_packages.R"))))
#source(here("Data-Code/support/move_ak_hi.R"))
source(here("Data-Code/support/get_geographic_info.R"))
#source(here("Data-Code/support/map_theme.R"))

# Import Data -------------------------------------------------------------
# Import data from this site: 
# https://www.census.gov/geographies/mapping-files/2017/geo/carto-boundary-file.html

county_map <- readOGR(dsn=here("Data/Input/CB/cb_2017_us_county_5m.shp"),
                      layer = "cb_2017_us_county_5m",verbose = FALSE) 

zip_map <- readOGR(dsn=here("Data/Input/CB-zip/cb_2017_us_zcta510_500k.shp"),
                      layer = "cb_2017_us_zcta510_500k",verbose = FALSE) 


# Get Gegographic Information (e.g., centroid, contiguous geographies, etc.)
df_county_info <- 
  county_map %>% 
  subset(GEOID != "99") %>% 
  get_geograhic_info()

write_rds(df_county_info,here("Data/Output/tidy-mapping-files","df_county_info.rds"))


# need to modify the get_geographic_info function to match this. 
# change GEOID to GEOID10
# df_zip_info <- 
#   zip_map %>% 
#   subset(GEOID10 != "99") %>% 
#   get_geograhic_info()
# 
# write_rds(df_zip_info,here("Data/Output/tidy-mapping-files","df_zip_info.rds"))



# 
# # Create a ggplot-friendly map data frame.
# # The AK/HI command isn't working.
# #county_map <- move_ak_hi(county_map,type="county")
# county_map$fips_code <- paste0(county_map$GEOID)
# county_map$lng <- unlist(lapply(county_map@polygons, function(dt) dt@labpt[1]))
# county_map$lat <- unlist(lapply(county_map@polygons, function(dt) dt@labpt[2]))
# # Project to albers
# county_map <- spTransform(county_map,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))
# 
# # simplify the polgons a tad (tweak 0.00001 to your liking)
# #simplify_polygon = FALSE
# #if (simplify_polygon) df_map <- gSimplify(df_map, tol = 0.00001)
# county_map <- gBuffer(county_map, byid=TRUE, width=0)
# 
# df_county_map = fortify(county_map,region = "fips_code") %>%
#   rename(fips_code = id) %>%
#   dplyr::select(fips_code,everything()) %>% 
#   tbl_df()
# 
# write_rds(df_county_map,here("Data/Output/tidy-mapping-files","df_county.rds"))
# 
# df_county_map %>%
#   filter(grepl("^48",fips_code)) %>%
#   tbl_df() %>%
#   mutate(test = factor(sample(1:10,nrow(.),replace=TRUE))) %>%
#   ggplot() +
#   aes(long,lat,group=group) +
#   geom_polygon(aes(fill = test)) +
#   geom_path(color="black") +
#   coord_equal() +
#   ggthemes::theme_tufte() +
#   theme(legend.position = "none") +
#   remove_all_axes


