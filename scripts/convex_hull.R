# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 20, 2022, rv15

# Purpose: Create convex hull of pts around given set of well locations,
# and then buffer this hull by a default value of 1% all around to create an 
# alternate, default boundary for the site of interest;
# returns a dataframe consisting of locations at vertices of expanded convex hull;
# can be used as alternate site boundary if user does not supply separate
# boundary file

# Revisions: rv15 allows alternate <locid> variable;
#  rv14 minor update;
#  rv13 does minor cleanup; update to library call;
#  rv12 makes minor changes;
#  rv11 fixes logical errors;
#  rv10 adjusts function name for consistency;
#  rv09 updates to latest <dplyr> along with minor fixes;
#  rv08 changes input to dataframe;
#  rv07 removes orientation() function as unnecessary;
#  rv06 fixes logic errors stemming from use of center of mass to compute expansion;
#  rv05 makes use of <dplyr>;
#  rv04 accounts for datasets with missing location coordinates;
#  rv03 adds 10% buffer to convex hull to ensure that no data pts
#  lie exactly on the alternate convex hull boundary


convex_hull <- function(f,locid=locid,east=east,north=north,bex=0.01) {
  library(tidyverse)
  library(rlang)
  library(purrr)
  
  locid <- ensym(locid)
  east <- ensym(east); north <- ensym(north)
  
# set up needed geometric point and vector functions  
  angle <- function(v1,v2) acos(sum(v1*v2)/(sqrt(sum(v1*v1))*sqrt(sum(v2*v2))))
  slope <- function(p1,p2) ifelse(p2[1]-p1[1]==0,Inf,(p2[2] - p1[2])/(p2[1] - p1[1]))
  yline <- function(x,p,m) m*(x - p[1]) + p[2]
  xline <- function(y,p,m) (1/m)*(y - p[2]) + p[1]
  linesect <- function(p1,p2,m1,m2) {
# handle case of vertical (infinite) slope m1
    if (is.infinite(m1) || is.nan(m1)) {
      x <- p1[1]; y <- p2[2] + m2*(x - p2[1])
    }
# handle case of vertical (infinite slope) m2
    if (is.infinite(m2) || is.nan(m2)) {
      x <- p2[1]; y <- p1[2]  + m1*(x - p1[1])
    }
# handles cases with non-vertical slopes
    x <- ((m1*p1[1] - m2*p2[1]) - (p1[2] - p2[2]))/(m1 - m2)
    y <- m1*(x - p1[1]) + p1[2]
    c(x,y)
  }
  dst <- function(p1,p2) sqrt(sum((p2 - p1)^2))
  
# identify unique well locations; eliminate missing locations
  loc <- f %>% select(!!locid,!!east,!!north) %>% filter(!is.na(!!east),!is.na(!!north)) %>% distinct()
# determine rows associated with vertices of strict convex hull;
# create subset of these pts; check for and remove duplicate planar pts in convex bound due to co-located wells at different depths
	pts <- with(loc,inject(chull(!!east,!!north)))
	convex.bnd <- loc[pts,] %>% select(!!east,!!north) %>% distinct()

# determine for each hull pt the line bisecting the angle formed with its two nearest hull nbrs
  convex.pts <- as.matrix(convex.bnd)
  n <- nrow(convex.bnd)
  E <- matrix(NA,n,2)
  for (i in 1:nrow(convex.bnd)) {
    if (i==1) {
     i1 <- n; i2 <- i; i3 <- i+1
    } else if (i==n) {
      i1 <- i-1; i2 <- i; i3 <- 1
    } else {
      i1 <- i-1; i2 <- i; i3 <- i+1
    }
    p1 <- convex.pts[i1,]; p2 <- convex.pts[i2,]; p3 <- convex.pts[i3,]
    v1 <- p1 - p2 ; v2 <- p3 - p2; v3 <- p1 - p3
    alph <- angle(v1,v2)/2; beta <- angle(v3,-v2)
    lx <- dst(p2,p3)
    
# use Law of Sines and Pythagorean theorem in two steps to deduce distance components from Q to R defined below
    dy <- lx*sin(alph)*sin(beta)/sin(alph + beta)
    dx <- dy/tan(alph)
# determine whether hull pt is above/below and to right/left of segment connecting nearest hull nbrs;
# also determine intersection point (R) between angle bisector and segment spanning p2 and p3 (v2)
    theta <- angle(v2,c(1,0))
    Rx <- abs(cos(theta)*dx); Ry <- abs(sin(theta)*dx)
    R <- p2 + sign(v2)*c(Rx,Ry)
    
# find intersection between line containing R perpendicular to v2 and segment spanning hull nbrs, p1 and p3 (v3);
# this intersection point (Q) will be along the angle bisector line on which the expansion point will be found;
# must first handle cases of possibly vertical (infinite) slopes for either point
    if (is.infinite(-1/slope(p2,p3))) {
      x <- R[1]; y <- p3[2] + slope(p1,p3)*(x - p3[1])
      Q <- c(x,y)
    } else if (is.infinite(slope(p1,p3))) {
      x <- p3[1]; y <- R[2] - (1/slope(p2,p3))*(x - R[1])
      Q <- c(x,y)
    } else {
      Q <- linesect(R,p3,-1/slope(p2,p3),slope(p1,p3))
    }
  # determine expansion length based on fixed percentage (bex) of diagonal of bounding box
    bx <- with(loc,inject(range(!!east))); by <- with(loc,inject(range(!!north)))
    dex <- bex*dst(c(bx[1],by[1]),c(bx[2],by[2]))
  # add expansion length along angle bisector line to find ith expansion point (Ei)
    es <- slope(p2,Q)

    if (es==Inf) {
      E[i,] <- ifelse((p2-Q)[2] >= 0,p2 + c(0,dex),p2 - c(0,dex))
    } else {
      sgn <- sign((p2-Q)[1])
      xe <- p2[1] + sgn*dex/sqrt(1 + es^2)
      E[i,] <- c(xe,yline(xe,p2,es))
    }
  }
# use expansion points as alternate convex boundary
  colnames(E) <- c('east','north')
  E <- as_tibble(E)
  E
}
