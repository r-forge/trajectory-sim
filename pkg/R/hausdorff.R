## Compute trajectory similarities according to the Hausdorff distance, under a convex distance measure.

"hausdorff" <- function(trajectories, pd=euclidian) {
	dH <- hausdorff.directed(trajectories, pd)
	
	# Take the matrix of directed Hausdorff distances and its transpose
	# (directed Hausdorff distances of reversed pairs) and take maximum
	pmax(dH, t(dH))
}

"hausdorff.directed" <- function(trajectories, pd=euclidian) {
	trajectory.similarity(trajectories, implementation=hausdorff.directed.pairwise, pd=pd, symmetric=FALSE, diagonal=0)
}

"hausdorff.pairwise" <- function(T1, T2, pd=euclidian) {
	max(hausdorff.directed.pairwise(T1, T2, pd),
			hausdorff.directed.pairwise(T2, T1, pd))
}

# max_{p in T1} min_{q in T2} pd(p,q)
"hausdorff.directed.pairwise" <- function(T1, T2, pd=euclidian) {
	## Naive O(n^3) implementation.
	# TODO: Look for better algorithm for polylines in arbitrary dimension and 
	# arbitrary (possibly convex) distance metric.
	
	d <- 0
	# The Hausdorff distance is reached either at a vertex of T1...
	for (i in 1:nrow(T1)) {
		di <- Inf
		p <- T1[i,]
		for (j in 1:(nrow(T2)-1)) {
			# find the minimum distance from p along the segment (T2[j],T2[j+1])

			di <- min(di, .segment.distance(p, T2[j,], T2[j+1,], pd))
		}
		d <- max(d, di)
	}
		
	# ... or where an edge of T1 intersects the Voronoi diagram of T2's segments
	# Consider each segment of T1 with each pair of segments of T2
	## TODO: be careful about trajectory endpoints
	for (i in 1:(nrow(T1)-1)) {
		di <- Inf
		for (j1 in 1:(nrow(T2)-2)) {
			ds1 <- .segment.distance(T1[i,],   T2[j1,], T2[j1+1,], pd)
			de1 <- .segment.distance(T1[i+1,], T2[j1,], T2[j1+1,], pd)
			for (j2 in (j1+1):(nrow(T2)-1)) {
				# First figure out if the closest distance occurs at one of the endpoints of (T1[i], T1[i+1])
				ds2 <- .segment.distance(T1[i,],   T2[j2,], T2[j2+1,], pd)
				de2 <- .segment.distance(T1[i+1,], T2[j2,], T2[j2+1,], pd)
				
				# If the distance differences have opposite sign, they need checking
				if ((ds2-ds1)*(de2-de1) < 0) {
					# Find a root for the distance difference function
					interpolate <- function(alpha) {
						T1[i,] + alpha * (T1[i+1,] - T1[i,])
					}
					
					a <- uniroot(function(alpha) {
						p <- interpolate(alpha)
						.segment.distance(p, T2[j1,], T2[j1+1,], pd) -
								.segment.distance(p, T2[j2,], T2[j2+1,], pd)
					}, interval=0:1)$root
					di <- min(di, .segment.distance(interpolate(a),
							T2[j1,], T2[j1+1,], pd))
				}
			}
		}
		## Assumption: each trajectory is contained in a finite-size bounding box
		if (is.finite(di)) { d <- max(d, di) }
	}
	
	d
}

# Find the closest point on segment (q1, q2) from p under pointwise distance meaure pd
".segment.distance" <- function(p, q1, q2, pd) {
	md <- optimize(function(alpha) {
				pd(p, q1 + alpha * (q2 - q1))
			}, interval=0:1)
	md$objective
}
