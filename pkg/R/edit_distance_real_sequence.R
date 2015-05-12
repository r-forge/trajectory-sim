## Compute trajectory similarities according to edit distance on real sequence
## Follows definition in:
## Chen, Lei, M. Tamer Özsu, and Vincent Oria. 
## Robust and fast similarity search for moving object trajectories. 
## Proc. 2005 ACM SIGMOD intern. conf. Management of data.

"EDR" <- function(trajectories, epsilon, pd=Lp.norm(Inf)) {
	trajectory.similarity(trajectories, implementation=EDR.pairwise, pd=pd,
			epsilon=epsilon, symmetric=TRUE, diagonal=0)
}

"EDR.pairwise" <- function(T1, T2, epsilon, pd=Lp.norm(Inf), ...) {
	## Reverse T1 and T2, since the DP algorithm works back to front
	traj.sim.dp(T1[nrow(T1):1,], T2[nrow(T2):1,],
			.EDR.step.fun, pd, epsilon=epsilon, ...)
}

## DP step function for DTW
".EDR.step.fun" <-  function(T1, T2, i, j, prev, pd, epsilon) {
	if (length(prev) == 0) {
		return(new("ts.dp.entry", 
			value=sum(pd(T1[1,], T2[1,]) >= epsilon),
			pred=NULL))
	} else if (length(prev) == 1) {
		## First row or column entry
		p <- names(prev)
		pv <- prev[[1]]@value + 1
		
		## One may also make a diagonal move here.
		dc <- max(i,j)-1 + sum(pd(T1[i,], T2[j,]) >= epsilon)
		if (pv > dc) {
			p <- NULL
			pv <- dc
		}
	} else {
		## Generic entry
		pv <- sapply(prev, slot, "value")
		pv[c("H","V")] <- pv[c("H","V")]+1
		pv["D"] <- pv["D"] + sum(pd(
				T1[i,], 
				T2[j,]) >= epsilon)
		pv <- pv[which.min(pv)]
		p <- names(pv)
		names(pv) <- NULL
	}
	
	new("ts.dp.entry", value=pv, pred=p)
}

