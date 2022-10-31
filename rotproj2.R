source("R/workers.R")
source("R/noproj.R")
m <- do.call(cbind, maps::map(plot = F)[1:2])

centre <- geosphere::randomCoordinates(1)[1:2]
print(centre)

par(mar = rep(0.1, 4))
layout(matrix(c(1, 1, 1, 2, 3, 4), byrow = T, ncol = 3))
plot(rot(m,centre), pch = ".", asp = "")

plot(laea(m, centre = centre), pch = ".", asp = 1); title("laea", line = -2)
lines(reproj::reproj_xy(m, glue::glue("+proj=laea +lon_0={centre[1]} +lat_0={centre[2]}"), source = "OGC:CRS84"), col = "firebrick")
text(0, 0, lab = paste0(round(centre, digits = 1), collapse = ", "))

plot(aeqd(rot(m, centre)) * 6378137, pch = ".", asp = 1); title("aeqd", line = -2)
lines(reproj::reproj_xy(m, glue::glue("+proj=aeqd +lon_0={centre[1]} +lat_0={centre[2]}"), source = "OGC:CRS84"), col = "firebrick")

plot(cass(rot(m, centre)) * 6378137, pch = ".", asp = 1); title("cass", line = -2)
lines(reproj::reproj_xy(m, glue::glue("+proj=cass +R=6378137 +lon_0={centre[1]} +lat_0={centre[2]}"), source = "OGC:CRS84"), col = "firebrick")

#library(rgl)
#plot3d(cart(rot(m, centre)));spheres3d(0, 0, 0, radius = 6000000, col = "white"); lines3d(0, 0, c(-1, 1) * 6378137 * 1.2); rgl::rglwidget()
## this is different
#plot(reproj::reproj_xy(m, glue::glue("+proj=ob_tran +o_proj=longlat +o_lon_p={centre[1]} +o_lat_p={(90 - centre[2])}"), source = "OGC:CRS84"), asp = 1, pch = ".")

