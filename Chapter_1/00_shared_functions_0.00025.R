# ==============================================================================
# COMPREHENSIVE SHARED FUNCTIONS - SAGUENAY FJORD RISK ASSESSMENT
# ==============================================================================
# Supporting ALL 6 stressor analyses:
#   1. Forestry (watersheds)
#   2. Agriculture (watersheds)
#   3. Navigation/Shipping (grid corridors)
#   4. WWTP (point sources)
#   5. Ice Fishing (point sources)
#   6. Coastal Development (point sources from raster)
#
# USAGE: Place in same folder as .Rmd files, then: source("00_shared_functions.R")
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION A: UNIVERSAL FUNCTIONS (ALL 6 STRESSORS USE THESE)
# ------------------------------------------------------------------------------
# ==============================================================================
# REQUIRED PACKAGES - LOAD IN THIS ORDER
# ==============================================================================
library(raster)      # Load FIRST
library(sf)
library(dplyr)
library(gdistance)

# Prevent terra conflicts if loaded elsewhere (safer check)
if ("package:terra" %in% search()) {
  detach("package:terra", unload = TRUE)
  library(raster)  # Reload raster to ensure proper namespace
}
# ==============================================================================
# Standard CRS and extent
target_crs <- "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs"
xmin <- 345000
ymin <- 5325000
xmax <- 450000
ymax <- 5370000

# ------------------------------------------------------------------------------
# A1. LOAD BATHYMETRY AND WATER MASK
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# A1. LOAD BATHYMETRY AND WATER MASK
# ------------------------------------------------------------------------------
load_bathymetry <- function(bathy_path = "Data/bathy.tif") {
  
  # Check if file exists
  if (!file.exists(bathy_path)) {
    stop(paste("Cannot find bathy.tif at:", bathy_path,
               "\nCurrent directory:", getwd(),
               "\nPlease provide correct path."))
  }
  
  bathy <- raster::raster(bathy_path)
  bathy_utm <- raster::projectRaster(bathy, crs = target_crs)
  saguenay_extent <- raster::extent(xmin, xmax, ymin, ymax)
  bathy_proj <- raster::crop(bathy_utm, saguenay_extent)
  bathy_200 <- raster::projectRaster(bathy_proj, crs = target_crs, res = c(200, 200))
  bathy_200[bathy_200 > 0] <- 0
  
  # Create water mask
  water_mask <- bathy_200
  water_mask[water_mask < 0] <- -1  # Water
  water_mask[water_mask >= 0] <- 0  # Land
  
  # Convert to polygons
  water_poly <- raster::rasterToPolygons(water_mask, dissolve = TRUE,
                                         fun = function(x) x == -1)
  water_sf <- st_as_sf(water_poly)
  water_sf <- st_transform(water_sf, crs = target_crs)
  
  return(list(
    bathy_200 = bathy_200,
    water_mask = water_mask,
    water_sf = water_sf
  ))
}
# ------------------------------------------------------------------------------
# A2. CREATE GRID SYSTEM
# ------------------------------------------------------------------------------
create_grid_system <- function(xmin, ymin, xmax, ymax, grid_res = 200) {
  bbox_utm <- list(rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax),
                         c(xmin, ymax), c(xmin, ymin))) %>%
    st_polygon() %>%
    st_sfc(., crs = target_crs)
  
  grid_points <- st_make_grid(bbox_utm, cellsize = grid_res, what = "centers") %>%
    st_geometry()
  
  grid_polygons <- st_make_grid(bbox_utm, cellsize = grid_res, what = "polygons") %>%
    st_geometry()
  
  grid_points_sf <- st_sf(geometry = grid_points)
  grid_polygons_sf <- st_sf(geometry = grid_polygons)
  
  return(list(
    grid_points_sf = grid_points_sf,
    grid_polygons_sf = grid_polygons_sf
  ))
}

# ------------------------------------------------------------------------------
# A3. CREATE STANDARD TRANSITION LAYER
# ------------------------------------------------------------------------------
create_transition_layer <- function(bathy_200) {
  fun_bathy <- function(x) {
    ifelse(x[1] > x[2],
           pmax(0.001, 1 - abs((x[2] - x[1]) / x[2])),
           0.45)
  }
  
  tr_bathy <- transition(bathy_200,
                        transitionFunction = fun_bathy,
                        directions = 8,
                        symm = FALSE) %>%
    geoCorrection(., type = "c")
  
  return(tr_bathy)
}

# ------------------------------------------------------------------------------
# A4. RELATIVE EXPOSURE FUNCTION (UNIVERSAL - WORKS FOR ALL SOURCE TYPES)
# ------------------------------------------------------------------------------
# This works for:
#   - Point sources (WWTP, Ice Fishing, Coastal Development)
#   - Line sources (if needed)
#   - Polygon sources (converts to centroids)
#   - Grid cells (Navigation corridors)

RelativeExposure <- function(source, transition, decay_type, grid_points_sf,
                             empty_raster) {
  
  if (!identical(st_crs(source), st_crs(grid_points_sf))) {
    stop("CRS mismatch between source and grid points")
  }
  
  source_coords <- st_coordinates(source)
  source_sp <- SpatialPoints(coords = source_coords,
                            proj4string = CRS(st_crs(source)$proj4string))
  
  grid_points_sp <- as(grid_points_sf, "Spatial")
  
  distance <- costDistance(transition, grid_points_sp, source_sp)
  
  if (all(is.na(distance))) {
    stop("All distance values are NA")
  }
  
  if (ncol(distance) == 1) {
    min_dist <- distance[, 1]
  } else {
    min_dist <- apply(distance, 1, min, na.rm = TRUE)
  }
  
  # Decay coefficient
  a_param <- switch(decay_type,
                    "I"   = -1,
                    "II"  = -0.1,
                    "III" = -0.01,
                    "IV"  = -0.001,
                    stop("Invalid decay type"))
  
  # Gaussian kernel with scaling factor
  exposure <- exp((0.00025 * a_param) * (min_dist)^2)
  
  coords <- st_coordinates(grid_points_sf)
  sp_exposure <- SpatialPointsDataFrame(
    coords = coords,
    data = data.frame(exposure = exposure),
    proj4string = CRS(st_crs(empty_raster)$proj4string)
  )
  
  exposure_rast <- rasterize(sp_exposure, empty_raster,
                            field = "exposure", fun = mean, background = NA)
  
  max_value <- cellStats(exposure_rast, stat = "max", na.rm = TRUE)
  if (!is.na(max_value) && max_value > 0) {
    exposure_rast[] <- exposure_rast[] / max_value
  }
  
  return(exposure_rast)
}

# ------------------------------------------------------------------------------
# A5. STANDARD EXPOSURE PLOT
# ------------------------------------------------------------------------------
plot_exposure <- function(exposure_raster, title, xmin, xmax, ymin, ymax) {
  plot(exposure_raster,
       main = title,
       xlab = "Easting (UTM Zone 19N)",
       ylab = "Northing (UTM Zone 19N)",
       xlim = c(xmin, xmax),
       ylim = c(ymin, ymax),
       axes = TRUE)
  
  scalebar(d = 10000,
           xy = NULL,
           type = "bar",
           divs = 2,
           below = "km",
           label = c(0, 5, 10))
}

# ------------------------------------------------------------------------------
# A6. PROJECT POINT TO NEAREST WATER CELL
# ------------------------------------------------------------------------------
# Use for point sources that may be on land (WWTP, Ice Fishing)

nearest_valid_cell <- function(point, raster) {
  valid_cells <- which(!is.na(values(raster)), arr.ind = TRUE)
  valid_coords <- xyFromCell(raster, valid_cells)
  distances <- sqrt((valid_coords[,1] - point[1])^2 +
                   (valid_coords[,2] - point[2])^2)
  nearest_idx <- which.min(distances)
  return(valid_coords[nearest_idx, , drop=FALSE])
}

# ==============================================================================
# SECTION B: WATERSHED-SPECIFIC FUNCTIONS (ONLY FORESTRY & AGRICULTURE)
# ==============================================================================

# ------------------------------------------------------------------------------
# B1. LOAD WATERSHEDS
# ------------------------------------------------------------------------------
load_watersheds <- function(watershed_path = "Data/corrected_watersheds.gpkg") {  # ← Added parameter
  
  # Check if file exists
  if (!file.exists(watershed_path)) {
    stop(paste("Cannot find watersheds at:", watershed_path,
               "\nCurrent directory:", getwd(),
               "\nPlease provide correct path."))
  }
  
  watersheds <- st_read(watershed_path, quiet = TRUE)
  watersheds <- st_transform(watersheds, crs = target_crs)
  return(watersheds)
}

# ------------------------------------------------------------------------------
# B2. FIND RIVER MOUTHS
# ------------------------------------------------------------------------------
find_river_mouths <- function(watersheds, water_sf) {
  if (!identical(st_crs(watersheds), st_crs(water_sf))) {
    water_sf <- st_transform(water_sf, st_crs(watersheds))
  }
  
  river_mouths_list <- list()
  
  for (i in 1:nrow(watersheds)) {
    watershed <- watersheds[i, ]
    river_name <- watershed$RIVER_NAME
    
    watershed <- st_make_valid(watershed)
    water_sf <- st_make_valid(water_sf)
    
    tryCatch({
      watershed_boundary <- st_boundary(watershed)
      intersection <- st_intersection(watershed_boundary, water_sf)
      
      if (nrow(intersection) > 0) {
        if (nrow(intersection) > 1) {
          intersection$length <- st_length(intersection)
          intersection <- intersection[which.max(intersection$length), ]
        }
        river_mouth_point <- st_centroid(intersection)
      } else {
        watershed_buffered <- st_buffer(watershed_boundary, 200)
        buffered_intersection <- st_intersection(watershed_buffered, water_sf)
        
        if (nrow(buffered_intersection) > 0) {
          intersection_centroid <- st_centroid(buffered_intersection)
          closest_point <- st_nearest_points(watershed_boundary, intersection_centroid)
          river_mouth_point <- st_cast(closest_point, "POINT")[1]
        } else {
          water_boundary <- st_boundary(water_sf)
          closest_points <- st_nearest_points(watershed_boundary, water_boundary)
          river_mouth_point <- st_cast(closest_points, "POINT")[1]
        }
      }
      
      # Move point into water
      water_centroid <- st_centroid(water_sf)
      current_coords <- st_coordinates(river_mouth_point)
      water_coords <- st_coordinates(water_centroid)
      
      direction <- water_coords - current_coords
      direction_norm <- direction / sqrt(sum(direction^2)) * 75
      
      final_coords <- current_coords + direction_norm
      final_point <- st_sfc(st_point(final_coords), crs = st_crs(river_mouth_point))
      
      river_mouth <- st_sf(
        RIVER_NAME = river_name,
        geometry = final_point
      )
      
      river_mouths_list[[i]] <- river_mouth
      
    }, error = function(e) {
      watershed_centroid <- st_centroid(watershed)
      water_centroid <- st_centroid(water_sf)
      
      centroid_coords <- st_coordinates(watershed_centroid)
      water_coords <- st_coordinates(water_centroid)
      direction <- (water_coords - centroid_coords) * 0.8
      
      final_coords <- centroid_coords + direction
      final_point <- st_sfc(st_point(final_coords), crs = st_crs(watershed_centroid))
      
      river_mouth <- st_sf(
        RIVER_NAME = river_name,
        geometry = final_point
      )
      
      river_mouths_list[[i]] <- river_mouth
    })
  }
  
  river_mouths <- do.call(rbind, river_mouths_list)
  river_mouths <- st_transform(river_mouths, crs = target_crs)
  return(river_mouths)
}

# ------------------------------------------------------------------------------
# B3. PLOT RIVER MOUTHS DIAGNOSTIC
# ------------------------------------------------------------------------------
plot_river_mouths_diagnostic <- function(watersheds, river_mouths_sf, water_sf,
                                        xmin, xmax, ymin, ymax) {
  n_rivers <- nrow(river_mouths_sf)
  river_colors <- rainbow(n_rivers, alpha = 0.7)
  names(river_colors) <- river_mouths_sf$RIVER_NAME
  
  plot(st_geometry(water_sf),
       col = "lightblue",
       border = "blue",
       main = "River Mouth Locations and Watershed Boundaries",
       xlab = "Easting (UTM Zone 19N)",
       ylab = "Northing (UTM Zone 19N)",
       xlim = c(xmin, xmax),
       ylim = c(ymin, ymax))
  
  plot(st_geometry(watersheds),
       add = TRUE,
       col = adjustcolor(river_colors[watersheds$RIVER_NAME], alpha = 0.3),
       border = river_colors[watersheds$RIVER_NAME],
       lwd = 1.5)
  
  points(st_coordinates(river_mouths_sf)[,1],
         st_coordinates(river_mouths_sf)[,2],
         pch = 21, cex = 2,
         bg = river_colors[river_mouths_sf$RIVER_NAME],
         col = "black", lwd = 2)
  
  text(st_coordinates(river_mouths_sf)[,1],
       st_coordinates(river_mouths_sf)[,2],
       labels = river_mouths_sf$RIVER_NAME,
       pos = 3, cex = 0.7, col = "black", font = 2)
  
  legend("topright",
         legend = c("Water Areas", "Watersheds", "River Mouths"),
         fill = c("lightblue", "gray", NA),
         pch = c(NA, NA, 21),
         pt.bg = c(NA, NA, "red"),
         pt.cex = c(NA, NA, 1.5),
         cex = 0.8,
         bg = "white")
}

# ==============================================================================
# SECTION C: CUSTOM TRANSITION FUNCTIONS (SPECIFIC STRESSORS)
# ==============================================================================

# ------------------------------------------------------------------------------
# C1. SHIPPING TRANSITION FUNCTION
# ------------------------------------------------------------------------------
# Use for Navigation analysis
# High conductance in water, blocked by land

create_shipping_transition <- function(bathy_200) {
  fun_shipping <- function(x) {
    ifelse(x[1] < 0 & x[2] < 0,
           0.9,    # Both in water
           ifelse(x[1] >= 0 | x[2] >= 0,
                  0.001,  # One or both on land
                  0.5))   # Boundary
  }
  
  tr_shipping <- transition(bathy_200,
                           transitionFunction = fun_shipping,
                           directions = 8,
                           symm = FALSE) %>%
    geoCorrection(., type = "c")
  
  return(tr_shipping)
}

# ------------------------------------------------------------------------------
# C2. WATER-ONLY TRANSITION FUNCTION
# ------------------------------------------------------------------------------
# Use for Coastal Development analysis
# Simplified: full conductance in water, none on land

create_water_transition <- function(bathy_200) {
  # Create water-only layer
  water_bathy <- bathy_200
  water_bathy[water_bathy >= 0] <- NA
  water_bathy[is.na(water_bathy)] <- 0
  
  fun_water <- function(x) {
    ifelse(is.na(x[1]) || is.na(x[2]), 0, 1)
  }
  
  tr_water <- transition(water_bathy,
                        transitionFunction = fun_water,
                        directions = 8,
                        symm = FALSE) %>%
    geoCorrection(., type = "c")
  
  return(tr_water)
}

# ==============================================================================
# SECTION D: HELPER FUNCTIONS (OPTIONAL BUT USEFUL)
# ==============================================================================

# ------------------------------------------------------------------------------
# D1. COMBINE WEIGHTED EXPOSURES
# ------------------------------------------------------------------------------
# Use when you have multiple sources with intensity weights
# Example: Multiple WWTPs, Ice fishing villages, Shipping corridors

combine_weighted_exposures <- function(exposure_list, weights, water_mask = NULL) {
  if (length(exposure_list) != length(weights)) {
    stop("Number of exposures must match number of weights")
  }
  
  # Normalize weights to sum to 1
  weights <- weights / sum(weights)
  
  # Initialize combined exposure
  combined <- exposure_list[[1]] * 0
  
  for (i in 1:length(exposure_list)) {
    weighted_exp <- exposure_list[[i]] * weights[i]
    weighted_exp[is.na(weighted_exp)] <- 0
    combined <- combined + weighted_exp
  }
  
  # Apply water mask if provided
  if (!is.null(water_mask)) {
    combined <- mask(combined, water_mask < 0)
  }
  
  # Normalize to 0-1
  max_val <- cellStats(combined, stat = "max", na.rm = TRUE)
  if (!is.na(max_val) && max_val > 0) {
    combined <- combined / max_val
  }
  
  combined[combined == 0] <- NA
  
  return(combined)
}

# ------------------------------------------------------------------------------
# D2. CREATE TEMPLATE RASTER
# ------------------------------------------------------------------------------
create_template_raster <- function(bathy_200) {
  template_raster <- raster(extent(bathy_200),
                           resolution = 200,
                           crs = projection(bathy_200))
  return(template_raster)
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================
#
# FORESTRY & AGRICULTURE (Watersheds):
# -------------------------------------
# source("00_shared_functions.R")
# bathy_data <- load_bathymetry()
# watersheds <- load_watersheds()
# river_mouths <- find_river_mouths(watersheds, bathy_data$water_sf)
# grid_system <- create_grid_system(xmin, ymin, xmax, ymax)
# tr_bathy <- create_transition_layer(bathy_data$bathy_200)
# template <- create_template_raster(bathy_data$bathy_200)
# exposure <- RelativeExposure(river_mouths, tr_bathy, "III",
#                             grid_system$grid_points_sf, template)
#
# NAVIGATION (Grid corridors):
# ----------------------------
# source("00_shared_functions.R")
# bathy_data <- load_bathymetry()
# grid_system <- create_grid_system(xmin, ymin, xmax, ymax)
# tr_shipping <- create_shipping_transition(bathy_data$bathy_200)
# template <- create_template_raster(bathy_data$bathy_200)
# # Load your shipping corridor data
# corridors <- st_read("shipping_corridors.gpkg")
# exposure <- RelativeExposure(corridors, tr_shipping, "III",
#                             grid_system$grid_points_sf, template)
#
# WWTP (Point sources):
# ---------------------
# source("00_shared_functions.R")
# bathy_data <- load_bathymetry()
# grid_system <- create_grid_system(xmin, ymin, xmax, ymax)
# tr_bathy <- create_transition_layer(bathy_data$bathy_200)
# template <- create_template_raster(bathy_data$bathy_200)
# # Load WWTP locations
# wwtp <- st_read("wwtp_locations.shp")
# # Project to water
# for(i in 1:nrow(wwtp)) {
#   coords <- nearest_valid_cell(st_coordinates(wwtp[i,]), bathy_data$bathy_200)
#   st_geometry(wwtp[i,]) <- st_sfc(st_point(coords), crs = target_crs)
# }
# # Calculate intensity_index for each WWTP
# wwtp$intensity_index <- calculate_your_intensity(...)
# # Calculate exposures
# exposures <- list()
# for(i in 1:nrow(wwtp)) {
#   exposures[[i]] <- RelativeExposure(wwtp[i,], tr_bathy, "IV",
#                                     grid_system$grid_points_sf, template)
# }
# combined <- combine_weighted_exposures(exposures, wwtp$intensity_index,
#                                       bathy_data$water_mask)
#
# ICE FISHING (Point sources):
# ---------------------------
# Same as WWTP, use your cabin counts or fishing effort as intensity_index
#
# COASTAL DEVELOPMENT (Points from raster):
# ----------------------------------------
# source("00_shared_functions.R")
# bathy_data <- load_bathymetry()
# grid_system <- create_grid_system(xmin, ymin, xmax, ymax)
# tr_water <- create_water_transition(bathy_data$bathy_200)
# template <- create_template_raster(bathy_data$bathy_200)
# # Load nightlight raster, threshold, convert to points
# light_points <- your_processing_function(...)
# light_points$intensity_index <- light_points$light_value / max_light
# exposure <- RelativeExposure(light_points, tr_water, "III",
#                             grid_system$grid_points_sf, template)
#
# ==============================================================================

cat("✓ Shared functions loaded successfully!\n")
cat("  - Section A: Universal functions (ALL 6 stressors)\n")
cat("  - Section B: Watershed functions (Forestry & Agriculture)\n")
cat("  - Section C: Custom transitions (Navigation & Coastal)\n")
cat("  - Section D: Helper functions (Combining exposures)\n")
