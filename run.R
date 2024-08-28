#https://reginalexavier.com/post/geosampa-do-laz-ao-tif/

Sys.setenv(LANG = "en")

#Installing and loading libs
libs <- c(
  "tidyverse", "sf", "tidygeocoder",
  "terra", "rayshader", "lidR", "ggrepel",
  "httr2", "magrittr", "devtools", "av", "gifski",
  "rayimage"
)

installed_libs <- libs %in% rownames(
  installed.packages()
)

if(any(installed_libs == FALSE)) {
  install.packages(libs[!installed_libs])
}

lapply(
  libs,
  library,
  character.only = TRUE
)

# Create custom cache folder to save geosampa data
create_cache_folder <- function() {
  cache_dirname <- "cache"
  path <- file.path(getwd(), cache_dirname)
  print(paste0("Creating cache folder at: ",  path))
  if(dir.exists(path)) {
    print("Cache folder already exists.")
    return(path)
  }
  dir.create(path, recursive = TRUE)
  print("Cache folder created succesfully.")
  return(path)
}

# download shapefile containing the quadricle spacings
download_quadricle_shapefile <- function(cache_dir) {
  dtm_grid_shapefile_url <- "https://geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadArquivo.aspx?orig=DownloadCamadas&arq=21_Articulacao%20de%20Imagens%5C%5CArticula%E7%E3o%20MDT%5C%5CShapefile%5C%5CSIRGAS_SHP_quadriculamdt&arqTipo=Shapefile" # nolint: line_length_linter.
  
  req <- httr2::request(dtm_grid_shapefile_url) |>
    httr2::req_user_agent("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36") |> # nolint: line_length_linter.
    httr2::req_progress() |>
    httr2::req_cache(cache_dir, max_age = Inf, debug = TRUE)
  resp <- req |> httr2::req_perform(path = file.path(cache_dir, "dtm_quadricle.zip"))
  shapefiles_path <- unzip(resp$body, exdir = cache_dir)
  return(shapefiles_path)
}

#Get which cells does the latlong intersects with
get_intersecting_cells <- function(lat_longs, quadricles, buffer_size) {
  buffer <- lat_longs |>
    terra::vect(geom = c("longitude", "latitude"), crs = "EPSG:4326") |>
    terra::project(quadricles) |>
    terra::buffer(buffer_size)
  
  terra::plot(quadricles)
  terra::plot(buffer, col = "blue", add = TRUE)
  
  intersections <- terra::intersect(buffer, quadricles)
  return(intersections$qmdt_cod)
}


#Download lidar data from the specific cell codes
download_lidar <- function(cell_codes, cache_dir) {
  cells <- tibble::tibble(cell_code = intersecting_cells_code) %>%
    dplyr::select(cell_code) %>%
    dplyr::mutate(
      url = paste0(
        "https://geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadArquivo.aspx?orig=DownloadMapaArticulacao&arq=MDS_2020%5C",
        cell_code,
        ".zip&arqTipo=MAPA_ARTICULACAO"
      ),
      file_path = file.path(cache_dir, "lidar", paste0(cell_code, ".zip"))
    )


  cache_path <- file.path(getwd(), "cache", "lidar")
  dir.create(cache_path, recursive = TRUE)
   
  reqs <- purrr::map(cells$url, function(item) {
    httr2::request(item) %>%
      httr2::req_user_agent("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36") %>% # nolint: line_length_linter.
      httr2::req_progress() %>%
      httr2::req_cache(cache_dir, max_age = Inf, debug = TRUE) %>%
      return()
  })

  try(resps <- reqs |>
    httr2::req_perform_parallel(
      paths = cells$file_path,
      on_error = "stop",
      progress = TRUE
    )
  )
  unziped <- purrr::map(resps, function(x) unzip(x$body[1], exdir = file.path(cache_dir, "lidar", "data")))
  print(unziped %>% unlist())
  return(unziped %>% unlist())
}

load_lidar <- function(dirname, point, buffer_size) {
  geometry <- terra::geom(point) |> tibble::as_tibble()

  # The of read filter syntax, lidR uses the lastools
  # https://lastools.github.io/download/las2las_README.md
  my_filter <- paste0(
    #"-first_only ",
    "-inside_circle ",
    geometry$x,
    " ",
    geometry$y,
    " ",
    buffer_size
  )
  print(my_filter)

  #select <- "xyzuicRGB"
  select <- "xyzucRGB"
  catalog <- lidR::readLAScatalog(dirname, select = select, filter = my_filter)
  clipped <- lidR::clip_circle(catalog, xcenter = geometry$x, ycenter = geometry$y, buffer_size)
  return(clipped)
}


#Build raster from LASCatalog
build_lidar_raster <- function(catalog, terrain_only = FALSE) {
  canopy_raster <- lidR::rasterize_canopy(
    catalog,
    keep_class = c(3L, 4L, 5L, 6L),
    res = 0.3,
    algorithm = lidR::p2r(subcircle = 0.3),
    parallel = TRUE,
    pkg = "terra"
  )

  terra::plot(canopy_raster)

  terrain_raster <- lidR::rasterize_terrain(
    las = catalog,
    res = 0.3,
    algorithm = lidR::tin(),
    parallel = TRUE,
    pkg = "terra"
  )

  terra::plot(terrain_raster)

  resampled_terrain <- terra::resample(terrain_raster, canopy_raster)
  total_raster <- terra::cover(canopy_raster, resampled_terrain)

  terra::plot(total_raster)
  if(terrain_only == TRUE) {
    return(resampled_terrain)
  }
  return(total_raster)
}

# Download orthophoto for specific cell_codes
donwload_orthophoto <- function(cell_codes, cache_dir) {
  ortho_cells <- tibble::tibble(cell_code = intersecting_cells_code) %>%
    dplyr::select(cell_code) %>%
    dplyr::mutate(
      url = paste0(
        "https://geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadArquivo.aspx?orig=DownloadMapaArticulacao&arq=ORTOFOTOS_2020_RGB%5C",
        cell_code,
        ".zip&arqTipo=MAPA_ARTICULACAO"
      ),
      file_path = file.path(cache_dir, "ortho", paste0(cell_code, ".zip"))
    )

  ortho_cells$url

  ortho_cache_path <- file.path(getwd(), "cache", "ortho")
  dir.create(ortho_cache_path, recursive = TRUE)

  reqs <- purrr::map(ortho_cells$url, function(item) {
    httr2::request(item) %>%
      httr2::req_user_agent("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36") %>% # nolint: line_length_linter.
      httr2::req_progress() %>%
      httr2::req_cache(cache_dir, max_age = Inf, debug = TRUE) %>%
      return()
  })

  try(resps <- reqs |>
    httr2::req_perform_parallel(
      paths = ortho_cells$file_path,
      on_error = "stop",
      progress = TRUE
    )
  )

  unziped <- purrr::map(resps, function(x) unzip(x$body[1], exdir = file.path(cache_dir, "ortho", "data")))
  unziped |> unlist() |> print()
  return(unziped |> unlist())
}

#load orthophotos
load_orthophotos <- function(file_paths, buffer) {
  jp2_files <- grep("\\.jp2$", ortho_files_path, ignore.case = TRUE, value = TRUE)
  if(jp2_files |> length() > 1) {
    spat_raster_collection <- terra::sprc(jp2_files)
    ortho_crop_collection <- terra::crop(spat_raster_collection, buffer, snap = "in")
    ortho_mosaic <- terra::mosaic(ortho_crop_collection)
    return(ortho_mosaic)
  } else {
    spat_raster <- terra::rast(jp2_files[1])
    ortho_cropped <- terra::crop(spat_raster, buffer, snap = "in", mask = TRUE)
    return(ortho_cropped)
  }
}
#load orthophotos

#Create overlay based on orthophoto and sampled to the lidar raster
create_ortho_overlay <- function(ortho_raster, lidar_raster) {
  overlay_path <- file.path(getwd(), "cache", "data")
  
  dir.create(overlay_path, recursive = TRUE)
  
  #ortho_resampled <- terra::resample(ortho_raster, lidar_raster)
  #ortho_resampled <- terra::aggregate(ortho_raster, fact = 5)
  # terra::align()
  terra::writeRaster(
    #ortho_resampled,
    ortho_raster,
    filename = file.path(getwd(), "cache", "data", "ortho.png"),
    overwrite = TRUE
  )
  img <- png::readPNG(file.path(getwd(), "cache", "data", "ortho.png"))
  return(img)
}



######### RUN #########
######### RUN #########
######### RUN #########
######### RUN #########
######### RUN #########

cache_dir <- create_cache_folder()

#quadricle shapefile is in SIRGAS 2000 23S
sirgas_2000_23s_crs <- "epsg:31983"
quadricles <- download_quadricle_shapefile(cache_dir) %>%
  dirname() %>%
  dplyr::first() %>%
  terra::vect(crs = sirgas_2000_23s_crs)


#Por algum motivo não pode ter acento no endereço, senão o geocoding não funciona, tem que rever isso.
addresses <- tibble::tribble(
  ~name,              ~addr,
  #"Masp",             "Av. Paulista, 1578 - Bela Vista, Sao Paulo - SP, 01310-200",
  #"Ed. Banespa",      "Rua Joao Bricola, 24 Sao Paulo - SP, 01014-900",
  "Edifício Copan",    "Av. Ipiranga, 200, Sao Paulo - SP, 01046-010",
  #"Ed. Ouro Preto",    "R. Novo Horizonte, 64, São Paulo - SP, 01244-020"
)
buffer_size <- 400
#geocode address
lat_longs <- addresses %>%
  tidygeocoder::geocode(addr, method = "osm", lat = latitude, long = longitude, verbose = TRUE)
lat_longs

intersecting_cells_code <- get_intersecting_cells(lat_longs, quadricles, buffer_size)
intersecting_cells_code

#Prepere the lidar parts, get a lider raster as result
lidar_files <- download_lidar(intersecting_cells_code, cache_dir)
print(lidar_files)

point <- lat_longs %>%
  terra::vect(geom = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  terra::project(quadricles)

catalog <- load_lidar(file.path(cache_dir, "lidar", "data"), point, buffer_size)
catalog


# Plot intersecting cells
plot(catalog, color = "RGB")

total_lider_raster <- build_lidar_raster(catalog, terrain_only = T)

terra::plot(total_lider_raster)
#End preparing the lidar parts

#prepare orthophoto
ortho_files_path <- donwload_orthophoto(intersecting_cells_code, cache_dir)

buffer <- lat_longs %>%
  terra::vect(geom = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  terra::project(quadricles) %>%
  terra::buffer(buffer_size)

ortho_rast <- load_orthophotos(file.path(cache_dir, "ortho", "data"), buffer)

terra::plotRGB(ortho_rast)
img <- create_ortho_overlay(ortho_rast, total_lider_raster)

#render
elmat <- rayshader::raster_to_matrix(total_lider_raster)

elmat %>%
  rayshader::sphere_shade(texture = "imhof1", colorintensity = 1, zscale = 0.6) %>%
  #rayshader::height_shade() |>
  rayshader::add_overlay(
    img,
    alphalayer = 1
  ) |>
  rayshader::plot_3d(
    elmat,
    zscale = 0.6,
    fov = 0,
    theta = 45,
    zoom = 0.75,
    phi = 45,
    baseshape = "circle",
    windowsize = c(1000, 800),
    clear_previous = T,
    shadow = T
  )

rayshader::render_camera(
  phi = 30,
  zoom = 0.3,
  theta = 25,
)

pic <- jpeg::readJPEG("ed_copan.jpg")

pic_df <- data.frame(
  pic = c("ed_copan.jpg"),
  long = c(lat_longs$longitude),
  lat = c(lat_longs$latitude)
)

pic_vect <- pic_df %>%
  terra::vect(geom = c("long", "lat"), crs = "EPSG:4326") %>%
  terra::project(quadricles)

my_mesh <- rayvertex::xy_rect_mesh(
  material = rayvertex::material_list(texture_location = "./ed_copan.jpg"), scale = c(200, 200, 200)
)

rayshader::render_raymesh(
  extent = terra::ext(total_lider_raster), heightmap = elmat,  long = geom(pic_vect)[, "x"], lat =geom(pic_vect)[, "y"], altitude = 950, zscale = 0.6,
  raymesh = my_mesh, angle = c(180, 90, 0),clear_previous = T
)


############## Until here is fine ##########
# Now is if you want to render the point cloud on top of everything and create a movie
small_buffer <- lat_longs %>%
  terra::vect(geom = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  terra::project(quadricles) %>%
  terra::buffer(100)

#must be an sf and not a spatVector
cropped_catalog <- lidR::clip_roi(catalog, sf::st_as_sf(small_buffer))

catalog_points <- catalog
catalog_points <- cropped_catalog

pal <- colorRampPalette(height.colors(20))
col <- pal(30)[as.numeric(cut(catalog_points$Z, breaks = 30))]
rgb_colors <- rgb(catalog_points$R, catalog_points$G, catalog_points$B, maxColorValue = 99999)

#you can plot the point cloud on top of the rayshader using this.
# change the total lider raster for it to be just terrain raster and not canopy to see the point and not the polygons
rayshader::render_points(
  extent = terra::ext(total_lider_raster),
  long = catalog_points$X,
  lat = catalog_points$Y,
  altitude = catalog_points$Z,
  size = 2,
  color = rgb_colors,
  zscale = 0.6,
  clear_previous = F
)

render_scalebar(limits = c(0, 5, 10), label_unit = "km", position = "W", y = 0,
                scale_length = c(0.33, 1), clear_scalebar = TRUE)


phivechalf = 30 + 60 * 1/(1 + exp(seq(-7, 20, length.out = 180)/2))
phivecfull = c(phivechalf, rev(phivechalf))
thetavec = -90 + 45 * sin(seq(0,359,length.out = 360) * pi/180)
zoomvec = 0.1 + 0.5 * 1/(1 + exp(seq(-5, 20, length.out = 180)))
zoomvecfull = c(zoomvec, rev(zoomvec))


df <- data.frame(phivecfull, thetavec, zoomvecfull)


apply(df, 1, function(item) {
    
  rayshader::render_camera(
    phi = item[1],
    zoom = item[3],
    theta = item[2],
  )
  
})

#rayshader::render_movie("movie5", title_text = "Ed. Copan, São Paulo", zoom = 0.5)
rayshader::render_movie("movie13", title_text = "Ed. Copan, São Paulo", type = "custom", frames = 360, phi = phivecfull, zoom = zoomvecfull, theta = thetavec)


render_highquality(
  filename = "3d-copan.png",
  samples = 128,
  scale_text_size = 24,
  parallel = TRUE,
  preview = TRUE,
  light = TRUE,
  interactive = FALSE,
  width = 800,
  height = 800,
)
