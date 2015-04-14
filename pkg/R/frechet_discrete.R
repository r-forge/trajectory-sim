## Compute trajectory similarities according to the discrete Fréchet distance, under a convex distance measure.

"discrete.frechet" <- function(trajectories, pd=euclidian) {
	trajectory.similarity(trajectories, implementation=discrete.frechet.pairwise, pd=pd, symmetric=TRUE, diagonal=0)
}

"discrete.frechet.pairwise" <- function(T1, T2, pd=euclidian, ...) {
	traj.sim.dp(T1, T2, .discrete.frechet.step.fun, pd, ...)
}

## DP step function for discrete Fréchet
".discrete.frechet.step.fun" <- function(T1, T2, i, j, prev, pd) {
	if (length(prev) == 0) {
		## Lower left corner of table
		new("ts.dp.entry",
				value=pd(T1[i,], T2[j,]),
				pred=NULL)
	} else {
		# Select the predecessor with the smallest d_dF so far
		pm <- which.min(sapply(prev, slot, "value" ))
		p <- prev[[pm]]
		new("ts.dp.entry",
				value=max(p@value, pd(T1[i,], T2[j,])),
				pred=names(pm))
	}
}

