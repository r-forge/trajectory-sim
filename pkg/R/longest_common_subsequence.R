## Compute trajectory similarities according to dynamic time warping
## TODO: Check details of implementation against original paper

"LCSS" <- function(trajectories, delta, epsilon, dp=Lp.norm(Inf)) {
	trajectory.similarity(trajectories, implementation=LCSS.pairwise, dp=dp,
			delta=delta, epsilon=epsilon, symmetric=TRUE, diagonal=0)
}

"LCSS.pairwise" <- function(T1, T2, delta, epsilon, dp=Lp.norm(Inf), ...) {
	traj.sim.dp(T1, T2, .LCSS.step.fun, dp, delta=delta, epsilon=epsilon, ...)
}

## DP step function for DTW
".LCSS.step.fun" <-  function(T1, T2, i, j, prev, dp, delta, epsilon) {
	## Take diagonal predecessor and add one to value
	if (dp(T1[i,], T2[j,]) < epsilon && abs(i-j) <= delta) {
		if ("D" %in% names(prev)) {
			p <- "D"
			pv <- prev[["D"]]@value
		} else {
			p <- NULL
			pv <- 0
		}
		pv <- pv+1
	} else {
		## Select only horiz. and vert. predecessors in a missing-friendly way
		prev <- prev[names(prev) %in% c("H","V")]
		if (length(prev) != 0) {
			pm <- which.max(sapply(prev, function(p) { p@value }))
			p  <- names(pm)
			pv <- prev[[pm]]@value
		} else {
			p <- NULL
			pv <- 0
		}
	}
	new("ts.dp.entry", value=pv, pred=p)
}

