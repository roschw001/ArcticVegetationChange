suppressMessages({
  library(tidyverse)
  library(sf)
  sf_use_s2(FALSE)
  library(stars)
  library(terra)
})

run <- as.numeric(commandArgs(trailingOnly = T))

{
  sdm_wd <- "/bioing/user/slisovsk/ArcticSDM/SDM_Results/"
  out_wd <- "/bioing/user/slisovsk/ArcticSDM/Arrays_dispersalVector_Sigmoid/"
}

speciesDisp <- suppressMessages(read_delim("/bioing/user/slisovsk/ArcticSDM/data/species_meters.csv", delim = ";"))

spResults <- (tibble(species = list.files(sdm_wd)) %>%
                left_join(tibble(fls = list.files(sdm_wd, pattern = "distRestriction.tif", recursive = T)) %>%
                            mutate(species = sapply(strsplit(fls, "/"), function(x) x[[1]]), predict = TRUE) %>%
                            dplyr::select(-fls), by = "species") %>% filter(!is.na(predict)) %>%
                left_join(speciesDisp, by = "species")) %>% filter(!is.na(meters_year_50), !is.na(meters_year_99))

sp <- spResults$species[run]


if(any(!file.exists(glue::glue("{out_wd}/{sp}/binArray_{c(0.5, 1, 2)}.rda")))) {
  
  
  #### Array Map
  load("/bioing/user/slisovsk/ArcticSDM/data/grid_25km.rda")
  ##############
  
  
  if(!file.exists(glue::glue("{out_wd}/{sp}"))) dir.create(glue::glue("{out_wd}/{sp}"))
  
    load(glue::glue("{sdm_wd}/{sp}/Predictions/pred_stars.rda"))
    
    ### Maxent threshold
    maxentTrheshold <- as.numeric(strsplit(
      strsplit(readLines(glue::glue("{sdm_wd}/{sp}/MaxentModelOutput/species.html"))[12], "</th><th>")[[1]][5], "</td><td>")[[1]][20])
    
    current <- read_stars(glue::glue("{sdm_wd}/{sp}/{sp}_MaxEnt_calibration_distRestriction.tif")) %>% setNames("present") %>%
      mutate(present = ifelse(present>=maxentTrheshold, 1, 0))
    
    distance <- current %>%
      mutate(present = ifelse(present==1, 1, NA)) %>% rast() %>% terra::distance() %>% st_as_stars() %>% setNames("distance") %>%
      suppressWarnings()
    
    ### Sigmoidal dispersal
    sigm_params <- lapply(c(0.5, 1, 2), function(d) {
      lapply(1:3, function(x) {
        p50 <- as.numeric(spResults %>% filter(species==sp) %>% dplyr::select(meters_year_50)) * 30 * x * d
        p99 <- as.numeric(spResults %>% filter(species==sp) %>% dplyr::select(meters_year_99)) * 30 * x * d
        s <- seq(0.0001, 0.1, length = 1000)
        c2  <- s[which.min(abs(0.99 - sapply(s, function(x) (1 / (1 + exp(-x*(p99 - p50)))))))]
        tibble(p50, c2)
      }) %>% Reduce("rbind",.)})
    
    
    for(d in 1:length(sigm_params)) {
      if(!file.exists(glue::glue("{out_wd}/{sp}/binArray_{c('0.5', '1', '2')[d]}.rda"))) {
        
        sppList <- lapply(1:length(pred_stars), function(spp) {
          tmp <- pred_stars[spp] %>% split(drop = TRUE)
          lapply(1:length(tmp), function(x) {
            
            p_dist <- distance %>% 
              mutate(runif = runif(prod(dim(.))),
                     p     = ifelse(distance<25, 1, 
                                    runif >= (1 / (1 + exp(-as.numeric(sigm_params[[d]][x,2])*(distance - as.numeric(sigm_params[[d]][x,1]))))))) %>%
              dplyr::select(p) %>% st_set_crs(st_crs(pred_stars))

            (rast(tmp[x]) * rast(p_dist)) %>% st_as_stars() %>% st_set_crs(st_crs(pred_stars)) %>%
              setNames("tmp") %>% mutate(tmp = ifelse(tmp>=maxentTrheshold, 1, 0)) %>%
              setNames(names(tmp)[x]) 
            }) %>% do.call("c", .) %>% merge(name = 'years') %>% setNames(names(pred_stars)[spp])
          })
        
        disp_stars <- do.call('c', sppList)
        
        curr_extr <- st_extract(current, grid$geometry %>% st_transform(st_crs(current))) %>% st_as_sf() %>%
          st_drop_geometry()
        
        extr_pred <- st_extract(disp_stars, grid$geometry %>% st_transform(st_crs(disp_stars))) %>% st_as_sf() %>%
          st_drop_geometry()
        
        listOut   <- lapply(1:3, function(x) {
          out <- as.numeric(unlist(cbind(curr_extr[,1], extr_pred[,rbind(c(1:3), c(4:6), c(7:9))[x,]])))
          ifelse(is.na(out), 0, out)
        })
        
        binarArray <- array(dim = c(1, nrow(grid), 4, 3))
        binarArray[1, , , 1] <- listOut[[1]]
        binarArray[1, , , 2] <- listOut[[2]]
        binarArray[1, , , 3] <- listOut[[3]]
        
        save(binarArray, file = glue::glue("{out_wd}/{sp}/binArray_{c('0.5', '1', '2')[d]}.rda"))
        gc()
      }
    }
  
}
