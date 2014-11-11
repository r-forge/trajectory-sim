## Compute trajectory similarities according to the continuous Fréchet distance, under a convex distance measure.

"frechet" <- function(trajectories, pd=euclidian) {
	trajectory.similarity(trajectories, implementation=frechet.pairwise, pd=pd, symmetric=TRUE, diagonal=0)
}

"frechet.decision" <- function(trajectories, epsilon=0, pd=euclidian) {
	trajectory.similarity(trajectories, implementation=frechet.decision,pairwise,
			epsilon=epsilon, pd=pd, symmetric=TRUE, diagonal=0)
}

"frechet.pairwise" <- function(T1, T2, pd=euclidian) {
	# Exponential search for upper bound on d_F
	epsilon <- 1.0
	while (!frechet.decision.pairwise(T1,T2,epsilon,pd)) {
		epsilon <- 2 * epsilon
	}

	# Binary search for exact value (up to machine precision)	
	low <- 0 # If d_F < 1 we may have d_F < epsilon/2
	high <- epsilon
	while (!(((low+high)/2) %in% c(low,high))) {
		epsilon <- (low+high)/2
		if (frechet.decision.pairwise(T1,T2,epsilon,pd)) {
			high <- epsilon
		} else {
			low  <- epsilon
		}
	}
	epsilon
}

"frechet.decision.pairwise" <- function(T1, T2, epsilon=0, pd=euclidian) {
	# Set up the free space diagram
	fsd <- as.list(rep(NA, (nrow(T1)-1)*(nrow(T2)-1)))
	dim(fsd) <- c(nrow(T1)-1, nrow(T2)-1)	# T1 indexes rows, T2 columns.
	
	## First check whether endpoints are free
	if (pd(T1[1,], T2[1,]) > epsilon
			|| pd(T1[nrow(T1),], T2[nrow(T2),]) > epsilon) { 
		return(FALSE)
	}
	for (i in 1:(nrow(T1)-1)) {
		for (j in 1:(nrow(T2)-1)) {
			free.L <- .frechet.free.edge(T2[j,], T1[i,], T1[i+1,], epsilon, pd)
			free.B <- .frechet.free.edge(T1[i,], T2[j,], T2[j+1,], epsilon, pd)
			## Compute reachable left and bottom boundary
			reach.L <- if(j==1) {
				# Left side of the diagram: only reachable if vertical segment
				# from (0,0) is completely free.
				if ((!is.null(free.L) && free.L[1]==0)
						&& (i==1 # && endpoint is known free
							|| (i>1 && 1.0 %in% fsd[[i-1,j]]$reach.L))) {
					free.L
				} else { NULL }
			} else {
				cl <- fsd[[i,j-1]] # left neighbour cell
				if (is.null(free.L)) { NULL }
				else if (!is.null(cl$reach.B)) { free.L }
				else if (!is.null(cl$reach.L) && cl$reach.L[1] <= free.L[2]) {
					c(max(cl$reach.L[1], free.L[1]), free.L[2])
				} else { NULL }
			}
			reach.B <- if(i==1) {
				# Bottom of the diagram: only reachable if horizontal segment
				# from (0,0) is completely free.
				if ((!is.null(free.B) && free.B[1]==0)
						&& (j==1  # && endpoint is known free
							|| (j>1 && 1.0 %in% fsd[[i,j-1]]$reach.B))) {
					free.B
				} else { NULL }
			} else {
				cb <- fsd[[i-1,j]] # bottom neighbour cell
				if (is.null(free.B)) { NULL }
				else if (!is.null(cb$reach.L)) { free.B }
				else if (!is.null(cb$reach.B) && cb$reach.B[1] <= free.B[2]) {
					c(max(cb$reach.B[1], free.B[1]), free.B[2])
				} else { NULL }
			}
			
			fsd[[i,j]] = list(free.L=free.L, free.B=free.B,
					reach.L=reach.L, reach.B=reach.B)
		} # for j
	} # for i
	# We know the endpoint is free, check whether the final cell is reachable
	!all(sapply(fsd[[nrow(T1)-1,nrow(T2)-1]][c("reach.L","reach.B")], is.null))
}

## returns interval boundaries 0 <= a <= b <= 1 such that q1 + alpha (q2-q1) 
## has distance at most epsilon to p for a <= alpha <= b, or NULL if no 
## interval exists.
".frechet.free.edge" <- function(p, q1, q2, eps, pd) {
	a <- b <- NA
	if (pd(p,q1) <= eps) { a <- 0 }
	if (pd(p,q2) <= eps) { b <- 1 }
	if (any(is.na(c(a,b)))) {
		# Find where the distance reaches its minimum to set intervals for the
		# root finding step
		dalpha <- function(alpha) { pd(p, q1 + alpha*(q2-q1)) }
		min <- optimize(dalpha, interval=0:1)
		# Check if the distance ever becomes <= epsilon
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

