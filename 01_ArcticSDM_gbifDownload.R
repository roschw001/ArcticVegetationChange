library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

wd <- "ArcticSDM/"

## Spatial extent
ecoreg <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

bbox_grid <- ecoreg %>% st_bbox(crs = 4326) %>% st_as_sfc() %>%
  st_make_grid(n = c(60, 30)) %>% st_sf() %>% filter(apply(st_intersects(., ecoreg, sparse = FALSE), 1, any)) %>%
  mutate(grid_id = 1:nrow(.)) %>% dplyr::select(grid_id)

# ggplot() +
#   geom_sf(data = ecoreg %>% dplyr::select('ECO_ID_U')) +
#   geom_sf(data = bbox_grid, fill = NA)

save_dir   <- glue::glue("{wd}data/GBIF/")

for(i in 1:nrow(bbox_grid)) {
  
  if(!file.exists(glue::glue("{save_dir}all/{bbox_grid[i,] %>% pull('grid_id')}_gbif_all.csv"))) {
    
    poly <- ecoreg %>% 
      st_intersection(bbox_grid[i,]) %>% st_geometry() %>% st_union() %>% suppressMessages() %>% suppressWarnings()
    
    bbox <- poly %>% st_bbox(crs = 4326) %>% st_as_sfc() %>% sf::st_as_text()
    
    file_list <- occ_download(pred_within(bbox), format = "SIMPLE_CSV")
    repeat{
      wait <- tryCatch(occ_download_wait(file_list[1]), error = function(e) NULL)
      if(!is.null(wait)) break
    }
    
    dwnl <- tryCatch(occ_download_get(file_list[1], overwrite = T) %>%
      occ_download_import() %>% suppressMessages(), error = function(e) NULL)
    
    if(!is.null(dwnl)) {
      
      inPoly <- dwnl[c(dwnl %>% st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
                            st_geometry() %>% st_intersects(., poly, sparse = FALSE)),] %>% suppressMessages()
  
      write_csv(inPoly, 
                file = glue::glue("{save_dir}/all/{bbox_grid[i,] %>% pull('grid_id')}_gbif_all.csv"))
      
      write_csv(inPoly %>% filter(kingdom == "Plantae"), 
                file = glue::glue("{save_dir}/Plantea/{bbox_grid[i,] %>% pull('grid_id')}_gbif_plantea.csv"))
      
      unlink(list.files(pattern = ".zip", full.names = T))
      
    } else {
      
      inPoly <- tibble(error = TRUE)
      
      write_csv(inPoly, 
                file = glue::glue("{save_dir}all/{bbox_grid[i,] %>% pull('grid_id')}_gbif_all.csv"))
      
      write_csv(inPoly, 
                file = glue::glue("{save_dir}Plantea/{bbox_grid[i,] %>% pull('grid_id')}_gbif_plantea.csv"))
    }
  }
  
}
