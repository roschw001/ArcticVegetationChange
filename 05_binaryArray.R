### Distribution extraction
library(sf)
sf_use_s2(FALSE)
library(stars)
library(abind)
library(tidyverse)
library(terra)

{
  sdm_wd   <- "sdm_dispersal/"
  out_wd   <- "sdm_arrays/"
  map_data <- "ArtcticSDM_data/"
}

spTable <- tibble(species = list.files(sdm_wd)) %>%
  left_join(tibble(fls = list.files(sdm_wd, pattern = "disp_stars", recursive = T)) %>%
              mutate(species = sapply(strsplit(fls, "/"), function(x) x[[1]]), predict = TRUE) %>%
              dplyr::select(-fls), by = "species") %>% filter(!is.na(predict))
save(spTable, file = glue::glue("{out_wd}/spTable.rda"))

### Grid
proj    <- "+proj=laea +lon_0=-170 +lat_0=90"
ecoreg  <- st_read(glue::glue("{map_data}Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0) %>%
  dplyr::select(ECO_NAME) %>% st_transform(proj)

map_wrld <- st_read(glue::glue("{map_data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
  st_transform(proj) %>% st_intersection(ecoreg %>% st_union()) %>% dplyr::select('admin')

grid <- st_make_grid(ecoreg, cellsize = 5000) %>% st_centroid() %>% st_sf() %>%
  mutate(Ecoregion = ecoreg$ECO_NAME[apply(st_intersects(., ecoreg, sparse = FALSE), 1, function(x) ifelse(any(x), which(x), NA))]) %>%
  filter(Ecoregion!="<NA>") %>%
  mutate(Admin = map_wrld$admin[apply(st_intersects(., map_wrld, sparse = FALSE), 1, function(x) ifelse(any(x), which(x), NA))]) %>%
  filter(Admin!="<NA>") %>% dplyr::select("Admin", "Ecoregion") %>% mutate(Id = 1:nrow(.), .before = "Admin")

# save(grid, file = glue::glue("{out_wd}/grid_5km.rda"))
plot(grid %>% dplyr::select(Ecoregion), pch = 16, cex = 0.6)
# load(glue::glue("{out_wd}/grid_5km.rda"))

binarArray <- array(dim = c(nrow(spTable), nrow(grid), 4, 3))

for(sp in spTable$species) {
  
  current   <- rast(glue::glue("{sdm_wd}{sp}/current.tif")) %>% st_as_stars()
  curr_extr <- st_extract(current, grid$geometry %>% st_transform(st_crs(current))) %>% st_as_sf() %>%
    st_drop_geometry()
  
  load(glue::glue("{sdm_wd}/{sp}/disp_stars.rda"))
  extr_pred <- st_extract(disp_stars, grid$geometry %>% st_transform(st_crs(disp_stars))) %>% st_as_sf() %>%
    st_drop_geometry()
  
  listOut   <- lapply(1:3, function(x) {
    out <- as.numeric(unlist(cbind(curr_extr[,1], extr_pred[,rbind(c(1:3), c(4:6), c(7:9))[x,]])))
    ifelse(is.na(out), 0, out)
  })
  
  binarArray[which(spTable$species==sp), , , 1] <- listOut[[1]]
  binarArray[which(spTable$species==sp), , , 2] <- listOut[[2]]
  binarArray[which(spTable$species==sp), , , 3] <- listOut[[3]]
  
  rm(disp_stars)
}

save(binarArray, file = glue::glue("{out_wd}/binaryArray.rda"))
