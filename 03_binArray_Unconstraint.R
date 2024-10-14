suppressMessages({
  library(tidyverse)
  library(sf)
  sf_use_s2(FALSE)
  library(stars)
  library(terra)
})

run <- as.numeric(commandArgs(trailingOnly = T))

{
  sdm_wd <- "ArcticSDM/SDM_Results/"
  out_wd <- "Arrays_Unconstraint/"
}

spResults <- tibble(species = list.files(sdm_wd)) %>%
  left_join(tibble(fls = list.files(sdm_wd, pattern = "distRestriction.tif", recursive = T)) %>%
              mutate(species = sapply(strsplit(fls, "/"), function(x) x[[1]]), predict = TRUE) %>%
              dplyr::select(-fls), by = "species") %>% filter(!is.na(predict))

sp <- spResults$species[run]

#### Array Map
load("data/grid_25km.rda")
##############

if(!file.exists(glue::glue("{out_wd}/{sp}"))) {
  dir.create(glue::glue("{out_wd}/{sp}"))
}

### distance
proj <- "+proj=laea +lon_0=-170 +lat_0=90"
load(glue::glue("{sdm_wd}/{sp}/tmp/modelTab.rda"))

current <- read_stars(glue::glue("{sdm_wd}/{sp}/{sp}_MaxEnt_calibration.tif")) %>% setNames("present")

model_dist <- modelTab %>% filter(type == "occurance") %>% st_as_sf(coords = c("lon", "lat"), crs = proj) %>%
  filter(apply(st_distance(.), 1, function(x) min(x[x>0])) < 1500*1000)

ext_dist <- model_dist %>% mutate(p = 1) %>% dplyr::select(p) %>% st_rasterize(
  current %>% setNames('dist') %>% mutate(dist = NA)) %>% rast() %>% distance() %>% 
  st_as_stars() %>% st_extract(grid$geometry %>% st_transform(st_crs(current))) %>%
  st_as_sf() %>% st_drop_geometry() %>% suppressWarnings()

mean <- 1500*1000
sd   <- 250*1000

if(!file.exists(glue::glue("{out_wd}/{sp}/binArray_distConstraint.rda"))) {
  
  ### Maxent threshold
  maxentTrheshold <- as.numeric(strsplit(
    strsplit(readLines(glue::glue("{sdm_wd}/{sp}/MaxentModelOutput/species.html"))[12], "</th><th>")[[1]][5], "</td><td>")[[1]][20])
  
  load(glue::glue("{sdm_wd}/{sp}/Predictions/pred_stars.rda"))  
  
  curr_extr <- st_extract(current, grid$geometry %>% st_transform(st_crs(current))) %>% st_as_sf() %>%
    st_drop_geometry()
  
  extr_pred <- st_extract(pred_stars, grid$geometry %>% st_transform(st_crs(pred_stars))) %>% st_as_sf() %>%
    st_drop_geometry()
  
  listOut   <- lapply(1:3, function(x) {
    out <- as.numeric(unlist(cbind(curr_extr[,1], extr_pred[,rbind(c(1:3), c(4:6), c(7:9))[x,]]))) * (1 - pnorm(rep(ext_dist[,1], 4), mean, sd))
    ifelse(is.na(out), 0, ifelse(out>=maxentTrheshold, 1, 0))
  })
  
  binarArray <- array(dim = c(1, nrow(grid), 4, 3))
  binarArray[1, , , 1] <- listOut[[1]]
  binarArray[1, , , 2] <- listOut[[2]]
  binarArray[1, , , 3] <- listOut[[3]]
  
  save(binarArray, file = glue::glue("{out_wd}/{sp}/binArray_distConstraint.rda"))
}