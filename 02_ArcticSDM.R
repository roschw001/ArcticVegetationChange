library(sf)
sf_use_s2(FALSE)
library(stars)
library(terra)
library(tidyverse)
library(flexsdm)
library(e1071)
library(pROC)
library(patchwork)
source("functions/functions.R")

{
  wd        <- "ArcticSDM/"
  gbif_data <- "ArcticSDM_data/"
  map_data  <- "ArtcticSDM_data/"
  env_data  <- "ArcticSDM/data/Modern/"
  out_wd    <- "ArcticSDM/Results/"
}


### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

### Maps
ecoreg   <- st_read(glue::glue("{map_data}Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

map_wrld <- st_read(glue::glue("{map_data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
  st_intersection(ecoreg %>% st_union()) %>% dplyr::select(c('name', 'admin'))

map  <- map_wrld %>% st_transform(proj)
# plot(map %>% dplyr::select(name))

### Resolution
res <- 5000

### Remove environmental variables
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
  # load(glue::glue("{map_data}GBIF_dataComp/occTab_all.rda"))
}

## reduce occTab
minSpecies <- 150
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
  
  predTab_curr <- env_curr %>% st_as_sf() %>% filter(apply(., 1, function(x) all(!is.na(x))))
}

### Future
{
  fut_wd   <- glue::glue("{gsub('Modern/', '', env_data)}/Future_CMIP6")
  fls      <- list.files(fut_wd, recursive = T)[grep(".tif", list.files(fut_wd, recursive = T))]
  
  futureList   <- tibble(fls  = fls,
                         year = sapply(strsplit(fls, "/"), function(x) round(mean(as.numeric(unlist(strsplit(x[[1]], "-")))),0)),
                         spp  = sapply(strsplit(fls, "/"), function(x) x[[2]]),
                         var  = sapply(strsplit(fls, "_"), function(x) x[2])) %>% arrange(year, spp, var) %>% filter(!(var%in%env_delete)) %>%
    filter(!is.na(var)) %>% group_split(spp)
  
  # env_future <- lapply(futureList, function(x) {
  #   year_stars <- x %>% group_split(year) %>% lapply(function(y) {
  #     read_stars(glue::glue("{fut_wd}/{y$fls}")) %>% setNames(y$var) %>% st_as_sf() %>% st_drop_geometry() %>%
  #       as.matrix()
  #   }) %>% abind::abind(., along = 3)
  # }) %>% abind::abind(., along = 4)
  
  # save(env_future, file = glue::glue("{fut_wd}/env_future.rds"))
  load(glue::glue("{fut_wd}/env_future.rds"))
  
  predTab_future0 <- read_stars(glue::glue("{fut_wd}/{futureList[[1]][1,1]}")) %>% setNames("pred_future") %>% st_as_sf()
}


### Species Loop
for(sp in spList$species) {
  
  tryCatch({
    
    if(!file.exists(glue::glue("{out_wd}/{sp}"))) {
      dir.create(glue::glue("{out_wd}/{sp}"))
      dir.create(glue::glue("{out_wd}/{sp}/tmp"))
      dir.create(glue::glue("{out_wd}/{sp}/Predictions"))
    }
    
    ### Occurance
    occ_sp    <- occTab %>% filter(species == sp) %>%
      mutate(id = 1:nrow(.), .before = "phylum")
    
    ### Spatial thinning with distance
    {
      distance <- 150
      
      ssPts <- sapply(3:4, function(x) (st_bbox(map)[x]*2)/1000/distance)
      
      if(!file.exists(glue::glue("{out_wd}/{sp}/tmp/occ_sp_filter.rda"))) {
        grd   <- st_sample(map %>% st_union(), ssPts[1]*ssPts[2], type = "random") %>%
          st_intersection(map %>% st_buffer(distance*1000)) %>% suppressWarnings()
        
        occ_sp_filter <- parallel::mclapply(1:length(grd), function(x) {
          occ_sp %>% mutate(dist = as.numeric(st_distance(., grd[x]))) %>%
            arrange(dist) %>% filter(dist<(distance*1000)*2) %>%
            slice(1) %>% dplyr::select(-dist)
        }, mc.cores = 5) %>% do.call("rbind", .) %>% filter(!duplicated(id))
        
        save(grd, file = glue::glue("{out_wd}/{sp}/tmp/grd.rda"))
        save(occ_sp_filter, file = glue::glue("{out_wd}/{sp}/tmp/occ_sp_filter.rda"))
        
      } else {
        load(glue::glue("{out_wd}/{sp}/tmp/grd.rda"))
        load(glue::glue("{out_wd}/{sp}/tmp/occ_sp_filter.rda"))
      }
      
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
      
      
      occurence_thinned <- occ_sp_filter %>% mutate(type = "occurance", .before = 'bio1') %>%
        bind_rows(occ_sp_eval %>% mutate(type = "occ_eval", .before = 'bio1')) %>% mutate(p = 1, .before = "bio1")
      
      pl <- ggplot() +
        geom_sf(data = ecoreg %>% st_transform(st_crs(occ_sp)) %>% st_simplify(dTolerance = 5000), 
                mapping = aes(geometry = geometry, fill = ECO_NAME), show.legend = FALSE, alpha = 0.5, color = "transparent") +
        geom_sf(data = occurence_thinned %>% filter(type=="occurance"), shape = 21, size = 2, alpha = 1, fill = "white") +
        geom_sf(data = occ_sp_eval, size = 2) +
        theme_light()
      
      
      ggsave(pl, filename = glue::glue("{out_wd}/{sp}/occuranceMap.png"), width = 20, height = 20)
      save(occurence_thinned, file = glue::glue("{out_wd}/{sp}/tmp/occurence_thinned.rda"))
    }
    
    ### Pseudo absence
    {
      if(!file.exists(glue::glue("{out_wd}/{sp}/tmp/modelTab.rda"))) {
        bfr <- occ_sp %>% st_buffer(150000) %>% st_union() %>%
          st_sym_difference(occ_sp %>% st_buffer(50000) %>% st_union() %>% st_simplify(dTolerance = 10000))
        
        occ_other <- occTab %>% filter(species != sp) %>%
          filter(st_intersects(., bfr, sparse = F)[,1])
        
        # ggplot() +
        #   geom_sf(data = map) +
        #   geom_sf(data = bfr, fill = "orange", alpha = 0.5) +
        #   geom_sf(data = occ_other, size = 0.5)
        
        ## svm
        svm_tab <- occurence_thinned %>% 
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
        if(nrow(abs_svm)>nrow(occurence_thinned %>% filter(type=="occurance"))) {
          for(r in seq(1, 0.1, length = 100)) {
            km <- tryCatch(kmeans(abs_svm %>% dplyr::select(bioClim_filter) %>% st_drop_geometry(),
                                  centers = round(nrow(occurence_thinned %>% filter(type=="occurance"))*r, 0)), error = function(e) NULL)
            if(!is.null(km)) break
          }
          abs_tab <- km$centers %>% as_tibble() %>% mutate(p = 0, .before = names(.)[1]) 
        } else {
          abs_tab <- abs_svm %>% dplyr::select(bioClim_filter) %>% st_drop_geometry() %>% as_tibble() %>% mutate(p = 0, .before = names(.)[1]) 
        }
        
        modelTab <- occurence_thinned %>% mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2], .before = 'type') %>%
          st_drop_geometry() %>% dplyr::select(species, ecoreg, lon, lat, type, p, bioClim_filter) %>%
          bind_rows(
            tibble(species = rep(sp, nrow(abs_tab)), ecoreg = NA, lon = NA, lat = NA, type = "absence") %>%
              bind_cols(abs_tab)
          ) %>%
          bind_rows(
            abs_svm %>% mutate(species = sp, .before = "ecoreg") %>%
              mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2], type = "abs_eval", .before = "p") %>%
              st_drop_geometry() %>% sample_n(min(nrow(abs_svm), sum(occurence_thinned$type=="occ_eval")))
          )
        
        save(modelTab, file = glue::glue("{out_wd}/{sp}/tmp/modelTab.rda"))
      } else load(glue::glue("{out_wd}/{sp}/tmp/modelTab.rda"))
    }
    
    ### MaxEnt Model
    {
      if(!file.exists(glue::glue("{out_wd}/{sp}/MaxentModelOutput/maxent_mod.rda"))) {
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
        
        save(maxent_mod, file = glue::glue("{out_wd}/{sp}/MaxentModelOutput/maxent_mod.rda"))
      } else load(glue::glue("{out_wd}/{sp}/MaxentModelOutput/maxent_mod.rda"))
      
      if(!file.exists(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration.tif"))) {
        
        pred_maxent <- predTab_curr %>% 
          mutate(pred = dismo::predict(maxent_mod, predTab_curr %>% st_drop_geometry() %>% dplyr::select(bioClim_filter), type = "response"))
        
        rastOut     <- st_rasterize(pred_maxent %>% dplyr::select(pred), st_as_stars(st_bbox(env_curr), nx = st_dimensions(env_curr)$x$delta, ny = st_dimensions(env_curr)$x$delta, values = NA_real_)) %>%
          setNames(glue::glue("{sp}_MaxEnt_current"))
        
        write_stars(rastOut, glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration.tif"))
        
      } else rastOut <- read_stars(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration.tif"))
      
    }
    
    ### Restricted predictions
    {
      if(!file.exists(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration_restricted.tif"))) {
        
        env_var <- parallel::mclapply(env_list[sapply(names(modelTab)[-c(1:6)], function(x) which(grepl(x, env_list, fixed = TRUE)))], function(x) {
          read_stars(glue::glue("{env_data}{x}")) %>% setNames(strsplit(x, "_")[[1]][2])
        }, mc.cores = 3) %>% Reduce("c", .) %>% rast()
        names(env_var) <- names(modelTab)[-c(1:6)]
        
        eval <- extra_eval(training_data = modelTab,
                           pr_ab = "p",
                           projection_data = env_var,
                           metric = "mahalanobis",
                           univar_comb = FALSE,
                           n_cores = 5, 
                           aggreg_factor = 6)
        
        truncated <- extra_truncate(
          suit = rastOut %>% rast(),
          extra = terra::resample(x = eval, y = rastOut %>% rast()),
          threshold = c(75, 100),
          trunc_value = 0)
        
        ### geographical distance
        distRast <- occ_sp %>% mutate(p = 1) %>% dplyr::select(p) %>% st_rasterize(
          rastOut %>% setNames('dist') %>% mutate(dist = NA)) %>% rast() %>% distance()
        
        mean <- 1000*1000
        sd   <- 250*1000
        
        distRast[] <- 1 - pnorm(distRast[], mean, sd)
        
        trunc_stars <- c((truncated[[1]] * distRast) %>% st_as_stars(),
                         (truncated[[2]] * distRast) %>% st_as_stars()) %>% st_set_crs(st_crs(rastOut))

        
        write_stars(trunc_stars[1] %>% setNames(glue::glue("{sp}_MaxEnt_current")), glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration_restricted_50.tif"))
        write_stars(trunc_stars[2] %>% setNames(glue::glue("{sp}_MaxEnt_current")), glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration_restricted_100.tif"))
        
        calibMap <- ggplot() +
          geom_stars(data = rastOut, downsample = 5, show.legend = F) +
          scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
          labs(title = glue::glue("{sp}: Current projection (MaxEnt)"), x = "", y = "") +
          geom_sf(data = modelTab %>% filter(type == "occurance") %>% st_as_sf(coords = c("lon", "lat"), crs = proj), pch = 16, cex = 0.75) +
          theme_void() +
          theme(plot.title    = element_text(size = 9, face = "bold.italic", hjust = 0.5),
                plot.subtitle = element_text(size = 8, face = "bold.italic", hjust = 0.5)) |
          ggplot() +
          geom_stars(data = trunc_stars[1], downsample = 5, show.legend = F) +
          scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
          labs(title = glue::glue("{sp}: Restricted (50)"), x = "", y = "") +
          theme_void() +
          theme(plot.title    = element_text(size = 9, face = "bold.italic", hjust = 0.5),
                plot.subtitle = element_text(size = 8, face = "bold.italic", hjust = 0.5)) |
          ggplot() +
          geom_stars(data = trunc_stars[2], downsample = 5, show.legend = F) +
          scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
          labs(title = glue::glue("{sp}: Restricted (100)"), x = "", y = "") +
          theme_void() +
          theme(plot.title    = element_text(size = 9, face = "bold.italic", hjust = 0.5),
                plot.subtitle = element_text(size = 8, face = "bold.italic", hjust = 0.5))
        
        ggsave(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration_restriction.png"), calibMap, units = "cm", width = 18*3, height = 15, bg = "white")
        
      }
    }
    
    ### AUC
    {
      # aucTab <- modelTab %>% filter(type %in% c("occ_eval", "abs_eval")) %>% st_as_sf(coords = c("lon", "lat"), crs = proj) %>%
      #   dplyr::select(ecoreg, p) %>% st_transform(st_crs(rastOut)) %>% 
      #   mutate(pred = st_extract(rastOut, .) %>% st_drop_geometry() %>% setNames("pred") %>% pull(pred)) %>% st_drop_geometry()
      # 
      # aucSpecies <- collSpecies %>% mutate(auc = (as.numeric(auc(aucTab$p, aucTab$pred)) %>% suppressMessages()),
      #                                      occCount = spList %>% filter(species==sp) %>% pull(count), .before = "bio1")
      # write_csv(collSpecies, glue::glue("{wd}Results/collinearityTabel.csv"), append = T)
    }
    
    ### Map Current
    {
      if(!file.exists(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_calibration.png"))) {
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
    }
    
    ### Predictions
    {
      ## future
      if(!file.exists(glue::glue("{out_wd}/{sp}/Predictions/pred_stars.rda"))) {
        envs    <- names(maxent_mod@presence)
        
        predList <- list()
        for(spp in 1:length(futureList %>% Reduce("rbind",.) %>% pull(spp) %>% unique())) {
          predList[[spp]] <- matrix(nrow = dim(env_future)[1], ncol = 3)
          predList[[spp]][,1] <- dismo::predict(maxent_mod, env_future[,envs,1,spp], type = "response"); gc()
          predList[[spp]][,2] <- dismo::predict(maxent_mod, env_future[,envs,2,spp], type = "response"); gc()
          predList[[spp]][,3] <- dismo::predict(maxent_mod, env_future[,envs,3,spp], type = "response"); gc()
        }
        predTab_future  <- predTab_future0 %>% bind_cols(do.call("cbind", predList) %>% as_tibble()) %>% dplyr::select(-pred_future)
        invisible(gc())
        
        starsList <- parallel::mclapply(1:ncol(predTab_future %>% st_drop_geometry()), function(x) {
          st_rasterize(predTab_future[,x], st_as_stars(st_bbox(env_curr), nx = st_dimensions(env_curr)$x$delta, ny = st_dimensions(env_curr)$x$delta, values = NA_real_))
        }, mc.cores = ncol(predTab_future %>% st_drop_geometry()))
        
        
        pred_stars <- lapply(tibble(ind   = 1:length(starsList), 
                                    spp   = rep(futureList %>% Reduce('rbind',.) %>% pull(spp) %>% unique(), each = 3),
                                    years = rep(futureList %>% Reduce('rbind',.) %>% pull(year) %>% unique(), 3)) %>% group_split(spp),
                             function(spp) {
                               do.call("c", starsList[spp$ind]) %>% setNames(as.factor(spp$years)) %>% merge(name = "years") %>% setNames(spp$spp[1])
                             }) %>% do.call("c", .)
        
        save(pred_stars, file = glue::glue("{out_wd}/{sp}/Predictions/pred_stars.rda"))
        rm(list = c('starsList', 'predTab_future', 'predList', 'maxent_mod'))
      } else load(glue::glue("{out_wd}/{sp}/Predictions/pred_stars.rda"))
    }
    
    ### Maps Future
    {
      spp_list <- lapply(tibble(spp = rep(1:3, each = 3), years = rep(1:3, 3)) %>% group_split(spp, years), function(x) {
        out <- ggplot() + geom_stars(data = pred_stars[x$spp][,,,x$years], downsample = 5, show.legend = F) +
          scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
          coord_equal() +
          theme_void()
      })
      
      library(gridExtra)
      spp_out <- do.call('grid.arrange', c(spp_list, nrow = 3, ncol = 3, left = "spp: [5.85, 3.70, 1.26]", top = "years: [2026, 2056, 2086]"))
      ggsave(glue::glue("{out_wd}/{sp}/Predictions/{sp}_MaxEnt_predictions.png"), spp_out, units = "cm", width = 20, height = 20, bg = "grey90") 
      dev.off()
    }
    
  }, error = function(e) NULL)
  
  gc()
}

