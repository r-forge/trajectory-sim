## Compute trajectory similarities according to the continuous Fréchet distance, under a convex distance measure.
## TODO: Check details of implementation against original paper

"frechet.decision.pairwise" <- function(T1, T2, epsilon=0, dp=euclidian) {
	# Set up the free space diagram
	fsd <- as.list(rep(NA, nrow(T1)*nrow(T2)))
	dim(fsd) <- c(nrow(T1), nrow(T2))	# T1 indexes rows, T2 columns.
	
	.frechet.free.edge(T1[11,], T2[10,], T2[11,], epsilon, dp)
	
}

## returns interval boundaries 0 <= a <= b <= 1 such that q1 + alpha (q2-q1) 
## has distance at most epsilon to p for a <= alpha <= b, or NULL if no 
## interval exists.
".frechet.free.edge" <- function(p, q1, q2, eps, dp) {
	a <- b <- NA
	if (dp(p,q1) <= eps) { a <- 0 }
	if (dp(p,q2) <= eps) { b <- 1 }
	if (any(is.na(c(a,b)))) {
		# Find where the distance reaches its minimum to set intervals for the
		# root finding step
		dalpha <- function(alpha) { dp(p, q1 + alpha*(q2-q1)) }
		min <- optimize(dalpha, interval=0:1)
		if (min$objective > eps) { return(NULL) }
		
		# Find roots of dalpha-eps to solve for a and/or b
		dae <- function(a) { dalpha(a)-eps }
		if (is.na(a)) { 
			a <- uniroot(dae, interval=c(0,min$minimum))$root
		}
		if (is.na(b)) {
			b <- uniroot(dae, interval=c(min$minimum,1))$root
		}
	}
	c(a,b)
}

## Find where the minimum distance from p to the segment (q1,q2) occurs
## under distance measure dp
".segment.min.dist" <- function(p, q1, q2, dp) {
	da <- function(alpha) { dp(p, q1 + alpha*(q2-q1)) }
	# perform golden section search
	ratio <- 1.5 - sqrt(1.25) # ~0.38166
	low <- 0
	high <- 1
	mid <- ratio
	dl <- da(low)
	dh <- da(high)
	dm <- da(mid)
	
	print(paste(dl, dm, dh))
	while (dm < min(dl, dh) && high-low > 10^-8) { # TODO: see what epsilon to use
		# Find where to evaluate the function again;
		# divide the bigger interval
		if (mid < (low+high)/2) {
			new <- mid + ratio * (high-mid) # golden section in (mid,high)
			dn <- da(new)
			if (dn < dm) { #minimum is in (dm, dh)
				low <- mid
				mid <- new
				dl <- dm
				dm <- dn
			} else {
				high <- new
				dh <- dn
			}
		} else {
			new <- mid - ratio * (mid-low)
			dn <- da(new)
			if (dn < dm) { #minimum is in (dm, dh)
				high <- mid
				mid <- new
				dh <- dm
				dm <- dn
			} else {
				low <- new
				dl <- dn
			}
		}
		print(paste(low, mid, high))
	}
	
}
