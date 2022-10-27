eccentricity <- function(a, b) {
  sqrt(1 - (b * b)/(a * a))
}

pj_qsfn <- function(sinphi, e, one_es) {
  EPSILON <- 1.0e-7

  #if (e >= EPSILON) {
    con = e * sinphi;
    div1 = 1.0 - con * con;
    div2 = 1.0 + con;

##
    res <- (one_es * (sinphi / div1 - (.5 / e) * log ((1. - con) / div2 )));
##      /* avoid zero division, fail gracefully */
    res[e< EPSILON] <- sinphi + sinphi;
    inf <- (div1 == 0.0 | div2 == 0.0)
    res[inf] <- Inf
    res
}
N_POLE = 0L
S_POLE = 1L
EQUIT  = 2L
OBLIQ  = 3L
EPS10   = 1.e-10
## not checking out of bounds
get_mode <- function(phi0) {
  t = abs(phi0);
  M_HALFPI = pi/2

  if (t > (M_HALFPI + EPS10 )) {
    stop("BAD lat")
  }
  if (abs(t - M_HALFPI) < EPS10) {
    if (phi0 < 0 ) return(S_POLE) else return(N_POLE)
  }

  if (phi0 < EPS10) return(EQUIT)
  return(OBLIQ)
}
WGS84 <- function() {
  list(a=6378137.0,      rf=298.257223563)
}
wgs <- WGS84()
a <- wgs[["a"]]; b <- wgs[["a"]]  - wgs[["a"]] * 1/wgs[["rf"]]

e <- eccentricity(a, b)
es <- e * e
one_es = 1.0 - es
if (one_es == 0.0) {
  #proj_log_error(P, _("Invalid eccentricity"));
  #proj_errno_set (P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
  #return PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE;
}
rone_es = 1.0/one_es
qp = pj_qsfn(1., e, one_es);
mmf = .5 / (1. - es);
#Q->apa = pj_authset(P->es);
#if (nullptr==Q->apa)
#  return destructor(P, PROJ_ERR_OTHER /*ENOMEM*/);
phi0 <- 0 #-42 * pi/180
lam0 <- 0 #147 * pi/180

mode <- EQUIT
if (mode == N_POLE) {

}
if (mode == S_POLE) {
  dd = 1.;
}
 if (mode == EQUIT) {
   rq = sqrt(.5 * qp)
      dd = 1. / (rq);
      xmf = 1.;
      ymf = .5 * qp;
}
if (mode == OBLIQ) {
        rq = sqrt(.5 * qp);
        sinphi = sin(phi0);
        sinb1 = pj_qsfn(sinphi, e, one_es) / qp;
        cosb1 = sqrt(1. - sinb1 * sinb1);
        dd = cos(phi0) / (sqrt(1. - es * sinphi * sinphi) *
                                  rq * cosb1);
        ymf = (xmf = rq) / dd;
        xmf = xmf * dd;
}


laea <- function(lon, lat = NULL) { #        /* Ellipsoidal, forward */
  xy <- xy.coords(lon, lat)
  lat <- xy$y
  lon <- xy$x
  bad <- is.na(lat)
  lp.phi <- lat[!bad] * pi/180
  lp.lam <- lon[!bad] * pi/180

  #double coslam, sinlam, sinphi, q, sinb=0.0, cosb=0.0, b=0.0;

    coslam = cos(lp.lam);
    sinlam = sin(lp.lam);
    sinphi = sin(lp.phi);
    q = pj_qsfn(sinphi, e, one_es);
    mode <- get_mode(0)
    if (mode == OBLIQ || mode == EQUIT) {
      sinb = q / qp;
      cosb2 = 1. - sinb * sinb;
      cosb = ifelse(cosb2 > 0, sqrt(cosb2), 0);
    }

if (mode == EQUIT) {

             b = 1. + cosb * coslam;
      }

    bad <- abs(b) < EPS10

          b = sqrt(2. / (1. + cosb * coslam));
          xy.y = b * sinb * ymf;
         # eqcon:
            xy.x = xmf * b * cosb * sinlam;
    out <- cbind(xy.x, xy.y)
    if (any(bad)) out[bad, ] <- NA
          out
}
aeqd <- function(lon, lat = NULL) {
  xy <- xy.coords(lon, lat)
  lat <- xy$y
  lon <- xy$x
  lp.phi <- lat * pi/180
  lp.lam <- lon * pi/180
  sinphi = sin(lp.phi);
  cosphi = cos(lp.phi);
  coslam = cos(lp.lam);
  # switch (Q->mode) {
  # case EQUIT:
  #   xy.y = cosphi * coslam;
  #   goto oblcon;
  #   case OBLIQ:
  Q.sinph0 <- 0
  Q.cosph0 <- 1
      xy.y = Q.sinph0 * sinphi + Q.cosph0 * cosphi * coslam;
   #plot(xy.y)
      # oblcon:
      #   if (fabs(fabs(xy.y) - 1.) < TOL)
      #     if (xy.y < 0.) {
      #       proj_errno_set(P, PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN);
      #       return xy;
      #     }
      # else
      #   return aeqd_e_forward(lp, P);
      # else {
        xy.y = acos(xy.y);
        xy.y = xy.y/ sin(xy.y);
        xy.x = xy.y * cosphi * sin(lp.lam);
        xy.y = xy.y * #(Q->mode == EQUIT) ? sinphi :
          Q.cosph0 * sinphi - Q.sinph0 * cosphi * coslam;

      cbind(xy.x, xy.y)
}


cass <- function(lon, lat = NULL) {
    xy <- xy.coords(lon, lat)
    lat <- xy$y
    lon <- xy$x
    lp.phi <- lat * pi/180
    lp.lam <- lon * pi/180

    xy.x  =  asin (cos (lp.phi) * sin (lp.lam));
    xy.y  =  atan2 (tan (lp.phi), cos (lp.lam)) - 0;
    cbind(xy.x, xy.y)
}





