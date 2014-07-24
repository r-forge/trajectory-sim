## Compute trajectory similarities according to the discrete Fréchet distance, under a convex distance measure.

"discrete.frechet" <- function(trajectories, dp=euclidian) {
	trajectory.similarity(trajectories, implementation=discrete.frechet.pairwise, dp=dp, symmetric=TRUE, diagonal=0)
}

"discrete.frechet.pairwise" <- function(T1, T2, dp=euclidian, ...) {
	traj.sim.dp(T1, T2, .discrete.frechet.step.fun, dp, ...)
}

## DP step function for discrete Fréchet
".discrete.frechet.step.fun" <- function(i, j, prev=list()) {
	if (length(prev) == 0) {
		## Lower left corner of table
		new("ts.dp.entry",
				value=dp(T1[i,], T2[j,]),
				pred=NULL)
	} else {
		# Select the predecessor with the smallest d_dF so far
		pm <- which.min(sapply(prev, function(p) { p@value }))
		p <- prev[[pm]]
		new("ts.dp.entry",
				value=max(p@value, dp(T1[i,], T2[j,])),
				pred=names(pm))
	}
}

