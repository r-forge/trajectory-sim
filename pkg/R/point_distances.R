## Several distance measures between points in spacetime,
## that can be used in the trajectories similarity measures.

## Each function takes two numeric vectors as input;
## The last element of each vector is regarded as time,
## the others as location coordinates.
## The input vectors should have identical lengths.

t1 <- matrix(c(0,1,5,0,2,3,1,2,4), 3)
t2 <- matrix(c(0,3,5,0,1,2,1,3,4), 3)

x1 <- c(0,0,0)
x2 <- c(2,0,0)

s <- Lp.norm(2)(x1,x2)
s
s <- Lp.norm2(2)(x1,x2)
s

"simple" <- function (p=2) {
  return (function(x1,x2){x1+x2})
  
}


"Lp.norm2" <- function (p=2) {
  if (p <= 0) {
    stop("The Lp norm is not defined for p <= 0")
  } else if (p < 1) {
    warning("The Lp norm is not convex for p < 1")
  }
  
  # The L_infinity norm is a special case
  if (is.infinite(p)) {
    return(function(x1, x2) {
      max(abs(x1[-length(x1)] - x2[-length(x2)]))
    })
  }
  return(function(x1, x2) {
    sum(abs(x1[-length(x1)] - x2[-length(x2)])^p)^(1/p)
  })
}


"Lp.norm" <- function (p=2) {
	if (p <= 0) {
		stop("The Lp norm is not defined for p <= 0")
	} else if (p < 1) {
		warning("The Lp norm is not convex for p < 1")
	}
	
	# The L_infinity norm is a special case
	if (is.infinite(p)) {
		return(function(x1, x2) {
			max(abs(x1[-length(x1)] - x2[-length(x2)]))
		})
	}
	return(function(x1, x2) {
		sum(abs(x1[-length(x1)] - x2[-length(x2)])^p)^(1/p)
	})
}
euclidian <- Lp.norm(2)
manhattan <- Lp.norm(1)



## Distance between lat/long coordinates (and time) on a sphere of radius R.
## Yields distance in same units as R. Default R is mean Earth radius in meters.
## Great circle distance is accurate to about 0.5% on Earth.
## Default lat/long units are degrees
"great.circle.dist" <- function(x1, x2, R=6371000, units=c("degrees","radians")) {
	if (length(x1) != 3 || length(x2) != 3) {
		stop("Points must have two coordinates and time")
	}
	units <- match.arg(units)
	if (units == "degrees") {
		x1 <- x1 * pi / 180
		x2 <- x2 * pi / 180
	}

	# The Haversine method for best accuracy and numerical stability
	a <- sin((x1[1]-x2[1])/2)^2 + cos(x1[1]) * cos(x2[1]) * sin((x1[2]-x2[2])/2)^2
	c <- 2 * atan2(sqrt(a), sqrt(1-a))
	R * c
}

