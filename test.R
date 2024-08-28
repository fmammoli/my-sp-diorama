library(rayshader)
library(terra)
loadzip <- tempfile()

download.file("https://tylermw.com/data/dem_01.tif.zip", loadzip)
localtif <- terra::rast(unzip(loadzip, "dem_01.tif"))

unlink(loadzip)

elmat <- raster_to_matrix(localtif)

elmat %>%
  sphere_shade(texture = "desert") %>%
  add_water(detect_water(elmat), color = "lightblue") %>%
  add_shadow(cloud_shade(elmat, zscale = 10, start_altitude = 500, end_altitude = 1000,), 0) %>%
  plot_3d(elmat, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800),
          background="darkred")

render_camera(theta = 20, phi=40,zoom= 0.64, fov= 56 )
render_clouds(elmat, zscale = 10, start_altitude = 800, end_altitude = 1000, attenuation_coef = 2, clear_clouds = T)
render_snapshot(clear=TRUE)