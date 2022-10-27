
## from rgl
rotationMatrix <- function (angle, x, y, z, matrix)
{
  if (missing(matrix)) {
    if (angle == 0)
      return(identityMatrix())
    u <- c(x, y, z)/sqrt(x^2 + y^2 + z^2)
    cosa <- cos(angle)
    sina <- sin(angle)
    matrix <- (1 - cosa) * outer(u, u)
    matrix <- matrix + diag(3) * cosa
    matrix[1, 2] <- matrix[1, 2] - sina * u[3]
    matrix[1, 3] <- matrix[1, 3] + sina * u[2]
    matrix[2, 1] <- matrix[2, 1] + sina * u[3]
    matrix[2, 3] <- matrix[2, 3] - sina * u[1]
    matrix[3, 1] <- matrix[3, 1] - sina * u[2]
    matrix[3, 2] <- matrix[3, 2] + sina * u[1]
  }
  if (identical(all.equal(dim(matrix), c(3, 3)), TRUE))
    matrix <- cbind(rbind(matrix, c(0, 0, 0)), c(0, 0, 0,
                                                 1))
  return(matrix)
}
identityMatrix <- function () {diag(nrow = 4)}

## recentre angular coords lon,lat (degrees), centre is lon,lat
rot <- function(x, centre = c(0, 0)) {
  x <- cbind(cart(x), 1)
  if (abs(centre[1] ) > 1e-10) {
    x <- x %*% rotationMatrix(centre[1] * pi/180, 0, 0, 1)
  }
  if (abs(centre[2] ) > 1e-10) {
    x <- x %*% rotationMatrix(-centre[2] * pi/180, 0, 1, 0)
  }
  geod(x)
}
## xyz to lon lat (radians)
geod <- function(x) {
  if (ncol(x) == 2) x <- cbind(x, 0)
  if (ncol(x) > 3) x <- x[,1:3]
  lon <- atan2(x[, 2], x[, 1]) * 180/pi
  lat <- atan2(x[, 3], sqrt(x[, 1]^2 + x[, 2]^2)) * 180/pi
  cbind(lon, lat)
}

## lonlat (degrees) to xyz
cart <- function(x, radius = 6378137) {
  lon <- x[,1]
  lat <- x[,2]
  theta <- lon * pi/180
  phi <- lat * pi/180
  cbind(cos(phi) * cos(theta), cos(phi) * sin(theta), sin(phi)) * radius
}
