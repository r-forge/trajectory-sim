## Generic code for trajectory distance measures based on dynamic programming.
## Each measure has to provide a "step function", that computes the value of
## the next cell in the DP table.

setClassUnion(".OptionalCharacter", c("character","NULL"))
setClass(Class = "ts.dp.entry",
	representation = representation(
			value = "numeric",
			pred  = ".OptionalCharacter",
			data  = "list"),
	prototype = prototype(
			value = as.double(NA),
			pred  = NULL),
	validity = function(object) {
		return(TRUE)
	}
)

setMethod(f = "as.numeric", 
		signature = c(x="ts.dp.entry"),
		function(x, ...) {
	x@value
})


# step.fun is a function that takes parameters
# 	T1, T2: trajectories
# 	i, j:   indices into T1 and T2
#	prev:   dyn prog table entries that precede this entry horizontally, vertically and/or diagonally
# 	pd:     a function that provides the distance between two location entries
# it returns the table entry for (i,j)
"traj.sim.dp" <- function(T1, T2, step.fun, pd=euclidian, ...,
		get.matching=FALSE, max.dt=Inf, 
		steps=list(H=c(-1,0), V=c(0,-1), D=c(-1,-1))) {
	dp.rows <- rep(list(NULL), nrow(T2))
	dp.emptyrow <- rep(list(NA), nrow(T1))
	
	## When max.dt is set, process only part of the column
	row.start <- 1 # At which row indices to start/stop processing the next column
	row.stop  <- 1 # First index that is NOT processed
	
	## Process the remaining columns
	for (j in 1:nrow(T2)) {
		dp.rows[[j]] <- dp.emptyrow
		
		## Search where to start processing this column; not before the previous
		while(row.start <= nrow(T1)
				&& T1[row.start,ncol(T1)] < T2[j,ncol(T2)] - max.dt) {
			row.start <- row.start + 1
		}
		while(row.stop < nrow(T1)
				&& T1[row.stop+1,ncol(T1)] <= T2[j,ncol(T2)] + max.dt) {
			row.stop <- row.stop + 1
		}
		
		for (i in row.start:row.stop) {
			preds <- lapply(steps, function(step) {
				i <- i + step[1]
				j <- j + step[2]
				if (i >= 1 && j >= 1) {
					dp.rows[[j]][[i]]
				} else {
					NA
				}
			})
			preds <- preds[!is.na(preds)] # Provide only predecessors that are set
			dp.rows[[j]][[i]] <- step.fun(T1, T2, i, j, preds, pd, ...)
		}

		if (!get.matching) {
			## Delete column that we no longer need to compute new entries,
			## in order to save space
			## Check step sizes to decide what can be removed
			keep.rows <- abs(min(0, sapply(steps, "[", 2)))
			if (j > keep.rows) {
				## Don't set to NULL: this shifts other entries in the list
				dp.rows[[j-keep.rows]] <- NA
			}
		}
	}
	
	res <- dp.rows[[nrow(T2)]][[nrow(T1)]]@value
	## Return matching, DP table values and the table itself, if requested.
	if (get.matching) {
		dp.table <- do.call(c, dp.rows)
		dim(dp.table) <- c(nrow(T1), nrow(T2))	# T1 indexes rows, T2 columns.
		dp.table <- t(dp.table) # Transpose to have T1 index horizontal axis
		attr(res, "matching") <- .matching(dp.table, steps)
		attr(res, "table.values") <- matrix(sapply(dp.table, function(entry) {
					as.numeric(entry) # Extracts value or passes through NA
				}), nrow(T2), nrow(T1))
		attr(res, "table") <- dp.table
	}
	res
}

".matching" <- function(dp.t, steps) {
	res <- matrix(NA, 2, sum(dim(dp.t))-1)
	pos <- c(ncol(dp.t), nrow(dp.t))
	l <- ncol(res)
	res[,l] <- pos

	entry <- dp.t[[pos[2],pos[1]]]
	while(!is.null(entry@pred)) {
		l <- l-1
		pos <- pos + steps[[entry@pred]]
		res[,l] <- pos
		entry <- dp.t[[pos[2],pos[1]]]
	}
	res[,l:ncol(res)]
}
