library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)
library(flexsdm)
library(e1071)
library(pROC)
library(patchwork)

wd        <- "ArcticSDM/"
gbif_data <- "ArcticSDM_data/"
map_data  <- "ArcticSDM/data/"
env_data  <- "ArcticSDM/environment/CHELSA data 5 km/Modern/"

out_wd    <- "ArcticSDM/SDM_Results"

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

### Maps
ecoreg   <- st_read(glue::glue("{map_data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

map_wrld <- st_read(glue::glue("{map_data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
  st_intersection(ecoreg %>% st_union()) %>% dplyr::select(c('name', 'admin'))

map  <- map_wrld %>% st_transform(proj)
# plot(map %>% dplyr::select(name))

### Resolution
res <- 5000

### Remove specific environmental variables
env_delete <- c("bio3", "swe", "scd")

##############
## GBIF Occ ##
##############

{
  ### Read in thinned occurrences
  thin_tab <- tibble(fls = list.files(glue::glue("{gbif_data}GBIF/GBIF_dataComp/ClusterThin"), pattern = ".rda")) %>%
    mutate(cell = as.numeric(sapply(strsplit(fls, "_"), function(x) x[[1]]))) %>%
    mutate(error = grepl("error", fls, fixed = T)) %>%
    filter(cell %in% ids)

  occTab_all <- parallel::mclapply(thin_tab$fls, function(x) {
      load(glue::glue("{gbif_data}GBIF/GBIF_dataComp/ClusterThin/{x}"))
      gbif_thin
    }, mc.cores = 1) %>% do.call("rbind", .) %>% st_as_sf(coords = c("lon", "lat"), crs = st_crs(map))

  # save(occTab_all, file = glue::glue("{map_data}GBIF_dataComp/occTab_all.rda"))
}
# load(glue::glue("{gbif_data}GBIF/GBIF_dataComp/occTab_all.rda"))

## reduce occTab
minSpecies <- 500
occTab <- occTab_all %>% filter(species %in% (occTab_all %>% group_by(species) %>% st_drop_geometry() %>% summarise(count = n()) %>%
                                                arrange(desc(count)) %>% filter(count>minSpecies) %>% pull(species))) %>%
  mutate(ecoreg = (ecoreg %>% pull(ECO_NAME))[apply(st_intersects(., ecoreg %>% dplyr::select(ECO_NAME) %>% st_transform(proj), sparse = FALSE), 1, function(x) which(x)[1])],
         .after = 'species')

spList <- occTab %>% st_drop_geometry() %>% group_by(species) %>% summarise(count = n()) %>%
  arrange(desc(count))


#######################
## SDM PreProcessing ##
#######################

### Prediction dataset
### Current
{
  env_list <- list.files(env_data)
  env_curr <- parallel::mclapply(env_list[-sapply(env_delete, function(x) which(grepl(x, env_list, fixed = TRUE)))], function(x) {
    read_stars(glue::glue("{env_data}{x}")) %>% setNames(strsplit(x, "_")[[1]][2])
  }, mc.cores = 3) %>% Reduce("c", .) %>% merge()
  
  predTab     <- env_curr %>% st_as_sf() %>% filter(apply(., 1, function(x) all(!is.na(x))))
}
### Future
{
  fut_wd   <- glue::glue("{gsub('Modern/', '', env_data)}Future_CMIP6")
  fls      <- list.files(fut_wd, recursive = T)
  
  futureList   <- tibble(fls  = fls,
                         year = sapply(strsplit(fls, "/"), function(x) round(mean(as.numeric(unlist(strsplit(x[1], "-")))),0)),
                         spp  = sapply(strsplit(fls, "/"), function(x) x[2]),
                         var  = sapply(strsplit(fls, "_"), function(x) x[2])) %>% arrange(year, spp, var) %>% filter(!(var%in%env_delete)) %>%
    filter(!is.na(var)) %>% group_split(spp)
  
  env_future_spp <- lapply(futureList, function(x) {
    year_stars <- x %>% group_split(year) %>% lapply(function(y) {
      read_stars(glue::glue("{fut_wd}/{y$fls}")) %>% setNames(y$var) %>% st_as_sf() %>% st_drop_geometry() %>%
        as.matrix()
    }) %>% abind::abind(., along = 3)
  }) %>% abind::abind(., along = 4)

  # save(env_future_spp, file = glue::glue("{fut_wd}/predictionArray.rds"))
  # load(glue::glue("{fut_wd}/predictionArray.rds"))
}


### Species Loop
for(sp in spList$species) {
  
  if(!file.exists(glue::glue("{out_wd}/{sp}"))) {
    dir.create(glue::glue("{out_wd}/{sp}"))
  }
  
  ### Occurance
  occ_sp    <- occTab %>% filter(species == sp) %>%
    mutate(id = 1:nrow(.), .before = "phylum")
  
  ### Spatial thinning with distance
  {
    distance <- 150
    
    ssPts <- sapply(3:4, function(x) (st_bbox(map)[x]*2)/1000/distance)
    grd   <- st_sample(map %>% st_union(), ssPts[1]*ssPts[2], type = "random") %>%
      st_intersection(map %>% st_buffer(distance*1000)) %>% suppressWarnings()
    
    occ_sp_filter <- parallel::mclapply(1:length(grd), function(x) {
      occ_sp %>% mutate(dist = as.numeric(st_distance(., grd[x]))) %>%
        arrange(dist) %>% filter(dist<(distance*1000)*2) %>%
        slice(1) %>% dplyr::select(-dist)
    }, mc.cores = 15) %>% do.call("rbind", .) %>% filter(!duplicated(id))
    
    # ggplot() +
    # geom_sf(data = map) +
    # geom_sf(data = occ_sp, shape = 16, size = 1, alpha = 0.4) +
    # geom_sf(data = occ_sp_filter, shape = 16, size = 1, color = "orange")
    }
  
  ### Co-linearity
  {
    env  <- (occ_sp_filter %>% st_drop_geometry())[,-c(1:7)]
    vif  <- correct_colinvar(env, method = c("vif", th = "10"))
    
    collSpecies <- tibble(species = sp) %>%
      bind_cols(
        tibble(var = names(occ_sp_filter %>% st_drop_geometry())[-c(1:7)]) %>%
          mutate(incl = ifelse(var%in%vif$removed_variables, FALSE, TRUE)) %>%
          pull(incl) %>% as.matrix(nrow = 1) %>% t() %>%
          as_tibble() %>% setNames(names((occ_sp_filter %>% st_drop_geometry())[,-c(1:7)]))
      )
    
    bioClim_filter <- collSpecies[,-1] %>% pivot_longer(cols = names(collSpecies)[-1]) %>%
      filter(value) %>% pull(name)
  }
  
  ### Partitioning
  {
    maxPart_perc <- 25
    
    
    occ_sp_eval <- occ_sp %>% mutate(selected = id%in%occ_sp_filter$id) %>%
      group_split(ecoreg) %>% lapply(function(tmp) {
        if(nrow(tmp)>6) {
          if(((sum(!tmp$selected)/nrow(tmp))*100) > maxPart_perc) {
            tmp %>% filter(!selected) %>% dplyr::select(-selected) %>% sample_n(floor((nrow(tmp %>% filter(selected))*(maxPart_perc/100))))
          } else {
            tmp_unselected <- tmp %>% filter(!selected)
            tmp_selected   <- tmp %>% filter(selected) %>% sample_n(floor((maxPart_perc*nrow(tmp))/100) - nrow(tmp_unselected))
            tmp_unselected %>% bind_rows(tmp_selected) %>% dplyr::select(-selected) %>% arrange(id)
          }
        }
      }) %>% Reduce("rbind", .)
    
    occ_sp_filter_select <- occ_sp_filter %>% filter(!(id%in%occ_sp_eval$id))
    
    
    occ_sp_filtered_partition <- occ_sp_filter %>% mutate(type = "occurance", .before = 'bio1') %>%
      bind_rows(occ_sp_eval %>% mutate(type = "occ_eval", .before = 'bio1')) %>% mutate(p = 1, .before = "bio1")
    
    pl <- ggplot() +
      geom_sf(data = ecoreg %>% st_transform(st_crs(occ_sp)) %>% st_simplify(dTolerance = 5000), 
              mapping = aes(geometry = geometry, fill = ECO_NAME), show.legend = FALSE, alpha = 0.5, color = "transparent") +
      geom_sf(data = occ_sp_filtered_partition %>% filter(type=="occurance"), shape = 21, size = 2, alpha = 1, fill = "white") +
      geom_sf(data = occ_sp_eval, size = 2) +
      theme_light()
    
    
    ggsave(pl, filename = glue::glue("{out_wd}/{sp}/occuranceMap.png"), width = 20, height = 20)
    save(occ_sp_filtered_partition, file = glue::glue("{out_wd}/{sp}/Occurence_thinned.rda"))
  }
  
  ### Pseudo absence
  {
    bfr <- occ_sp %>% st_buffer(150000) %>% st_union() %>%
      st_sym_difference(occ_sp %>% st_buffer(50000) %>% st_union() %>% st_simplify(dTolerance = 10000))
    
    occ_other <- occTab %>% filter(species != sp) %>%
      filter(st_intersects(., bfr, sparse = F)[,1])
    
    ## svm
    svm_tab <- occ_sp_filtered_partition %>% 
      dplyr::select(ecoreg, p, unique(unlist(bioClim_filter))) %>%
      bind_rows(
        occ_other %>% mutate(p = 0, .before = 'bio1') %>% 
          dplyr::select(ecoreg, p, bioClim_filter)
      ) %>%
      filter(apply(.[,bioClim_filter], 1, function(x) all(!is.na(x))))
    
    svm_mod <- svm(svm_tab %>% st_drop_geometry() %>% dplyr::select(bioClim_filter), 
                   svm_tab %>% st_drop_geometry() %>% pull(p), type = "one-classification")
    abs_svm <- svm_tab %>% filter(p == 0 & !predict(svm_mod)) 
    
    ## k-means
    if(nrow(abs_svm)>nrow(occ_sp_filtered_partition %>% filter(type=="occurance"))) {
      for(r in seq(1, 0.1, length = 100)) {
        km <- tryCatch(kmeans(abs_svm %>% dplyr::select(bioClim_filter) %>% st_drop_geometry(),
                              centers = round(nrow(occ_sp_filtered_partition %>% filter(type=="occurance"))*r, 0)), error = function(e) NULL)
        if(!is.null(km)) break
      }
      abs_tab <- km$centers %>% as_tibble() %>% mutate(p = 0, .before = names(.)[1]) 
    } else {
      abs_tab <- abs_svm %>% dplyr::select(bioClim_filter) %>% st_drop_geometry() %>% as_tibble() %>% mutate(p = 0, .before = names(.)[1]) 
    }
    
    modelTab <- occ_sp_filtered_partition %>% mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2], .before = 'type') %>%
      st_drop_geometry() %>% dplyr::select(species, ecoreg, lon, lat, type, p, bioClim_filter) %>%
      bind_rows(
        tibble(species = rep(sp, nrow(abs_tab)), ecoreg = NA, lon = NA, lat = NA, type = "absence") %>%
          bind_cols(abs_tab)
      ) %>%
      bind_rows(
        abs_svm %>% mutate(species = sp, .before = "ecoreg") %>%
          mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2], type = "abs_eval", .before = "p") %>%
          st_drop_geometry() %>% sample_n(min(nrow(abs_svm), sum(occ_sp_filtered_partition$type=="occ_eval")))
      )
    
    save(modelTab, file = glue::glue("{out_wd}/{sp}/modelTab_filtered_partition_absence.rda"))
  }
  
  ### MaxEnt Model
  {
    args_list <- c('responsecurves=FALSE',
                   'jackknife=FALSE',
                   'pictures=FALSE',
                   'autofeature=FALSE',
                   'linear=TRUE',
                   'quadratic=TRUE',
                   'product=TRUE',
                   'threshold=TRUE',
                   'hinge=TRUE',
                   'betamultiplier=1')
    
    maxTab <- modelTab %>% filter(type%in%c('occurance', 'absence')) %>%
      dplyr::select(p, bioClim_filter) %>% na.omit()
    
    maxent_mod <- dismo::maxent(maxTab[,-1], p = maxTab$p, nbg = 0, args = args_list,
                                path = glue::glue("{out_wd}/{sp}/MaxentModelOutput"))
    
    pred_maxent <- predTab %>% mutate(pred = dismo::predict(maxent_mod, predTab %>% st_drop_geometry() %>% dplyr::select(bioClim_filter), type = "response"))
    
    
    st_dimensions(env_curr)$x$delta
    
    rastOut     <- st_rasterize(pred_maxent %>% dplyr::select(pred), st_as_stars(st_bbox(env_curr), 
                          nx = abs(st_dimensions(env_curr)$x$delta), 
                          ny = abs(st_dimensions(env_curr)$y$delta), values = NA_real_)) %>% setNames(glue::glue("{sp}_MaxEnt_current"))
  }
  
  ### Maps
  {
    outMap <- ggplot() +
      geom_stars(data = rastOut, downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      labs(title = glue::glue("{sp}: Current projection (MaxEnt)"), x = "", y = "") +
      geom_sf(data = modelTab %>% filter(type == "occurance") %>% st_as_sf(coords = c("lon", "lat"), crs = proj), pch = 16, cex = 0.25) +
      theme_void() +
      theme(plot.title    = element_text(size = 9, face = "bold.italic", hjust = 0.5),
            plot.subtitle = element_text(size = 8, face = "bold.italic", hjust = 0.5))
    
    ggsave(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration.png"), outMap, units = "cm", width = 18, height = 15, bg = "white")
  }
  
  ### Predictions
  {
    dir.create(glue::glue("{out_wd}/{sp}/Predictions"))
    save(rastOut, file = glue::glue("{out_wd}/{sp}/Predictions/pred_stars_current.rda"))
    
    ## future
    pred_stars_list <- lapply(futureList %>% Reduce("rbind",.) %>% pull(spp) %>% unique(), function(spp) {
      spp_run <- which((futureList %>% Reduce("rbind",.) %>% pull(spp) %>% unique())==spp)
      envs    <- which((futureList %>% Reduce("rbind",.) %>% pull(var) %>% unique())%in%bioClim_filter)
      years   <- futureList %>% Reduce("rbind",.) %>% pull(year) %>% unique()
      
      tmp <- matrix(nrow = dim(env_future_spp)[1], ncol = 3)
      for(i in 1:3) {
        invisible(gc())
        tmp[,i] <- dismo::predict(maxent_mod, env_future_spp[, envs, i, spp_run], type = "response")
      }
      
      out <- st_rasterize(predTab %>% dplyr::select(geometry) %>% mutate(p1 = tmp[,1], p2 = tmp[,2], p3 = tmp[,3]), 
                   st_as_stars(st_bbox(env_curr), nx = abs(st_dimensions(env_curr)$x$delta), 
                   ny = abs(st_dimensions(env_curr)$y$delta), values = NA_real_))
      
      
      st_set_dimensions(out %>% setNames(years) %>% merge(), names = c("x", "y", "years")) %>% setNames(glue::glue("{spp}"))
                   
    })
  }
  
  ### Maps
  {
    mp1 <- ggplot() +
      geom_stars(data = pred_stars_list[[1]][,,,1], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('2026')
    
    mp2 <- ggplot() +
      geom_stars(data = pred_stars_list[[1]][,,,2], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('2056')
    
    mp3 <- ggplot() +
      geom_stars(data = pred_stars_list[[1]][,,,3], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('2086')
    
    spp1 <- (mp1 / mp2 / mp3)
    
    mp1 <- ggplot() +
      geom_stars(data = pred_stars_list[[2]][,,,1], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('')
    
    mp2 <- ggplot() +
      geom_stars(data = pred_stars_list[[2]][,,,2], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('')
    
    mp3 <- ggplot() +
      geom_stars(data = pred_stars_list[[2]][,,,3], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('')
    
    spp2 <- (mp1 / mp2 / mp3)
    
    
    mp1 <- ggplot() +
      geom_stars(data = pred_stars_list[[3]][,,,1], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('')
    
    mp2 <- ggplot() +
      geom_stars(data = pred_stars_list[[3]][,,,2], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('')
    
    mp3 <- ggplot() +
      geom_stars(data = pred_stars_list[[3]][,,,3], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      theme_void() +
      coord_equal() +
      ggtitle('')
    
    spp3 <- (mp1 / mp2 / mp3)
    
    outPredMaps <- spp1 | spp2 | spp3
    
  }
  
  ### Save
  {
    ggsave(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_predictions.png"), outPredMaps, units = "cm", width = 20, height = 20, bg = "grey70")
    save(pred_stars_list, file = glue::glue("{out_wd}/{sp}/Predictions/pred_stars_predictions.rda"))
  }
  
}

