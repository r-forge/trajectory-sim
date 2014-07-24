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
	validity = function(object){
		if (!is.null(object@pred) && !(object@pred %in% c("H","V","D"))) {
			stop("'pred' must be NULL or set to 'H', 'V' or 'D', indicating direction of predecessor in table")
		}
		return(TRUE)
	}
)

setMethod(f = "as.numeric", 
		signature = c(x="ts.dp.entry"),
		function(x, ...) {
	x@value
})


# step.fun is a function that takes parameters
# 	i, j:   indices into T1 and T2
#	prev:   dyn prog table entries that precede this entry horizontally, vertically and/or diagonally
# It may assume the following variables to be defined in its environment:
# 	T1, T2: trajectories
# 	dp:     a function that provides the distance between two location entries
# it returns the table entry for (i,j)
"traj.sim.dp" <- function(T1, T2, step.fun, dp=euclidian, ...,
		get.matching=FALSE, max.dt=Inf, 
		steps=list(H=c(0,-1), V=c(-1,0), D=c(-1,-1))) {
	dp.cols <- rep(list(NULL), nrow(T2))
	dp.emptycol <- rep(list(NA), nrow(T1))
	col.start <- 1 # At which row indices to start/stop processing the next column
	col.stop  <- 1 # First index that is NOT processed
	
	# Set up an environment for the step function to avoid unnecessary argument passing
	env <- new.env(parent=environment())
	assign("T1", T1, envir=env)
	assign("T2", T2, envir=env)
	assign("dp", dp, envir=env)
	environment(step.fun) <- env
	
	## Process the remaining columns
	for (j in 1:nrow(T2)) {
		dp.cols[[j]] <- dp.emptycol
		
		## Search where to start processing this column; not before the previous
		while(col.start <= nrow(T1)
				&& T1[col.start,ncol(T1)] < T2[j,ncol(T2)] - max.dt) {
			col.start <- col.start + 1
		}
		while(col.stop <= nrow(T1)
				&& T1[col.stop,ncol(T1)] <= T2[j,ncol(T2)] + max.dt) {
			col.stop <- col.stop + 1
		}
		
		for (i in col.start:(col.stop-1)) {
			preds <- lapply(steps, function(step) {
				i <- i + step[1]
				j <- j + step[2]
				if (i >= 1 && j >= 1) {
					dp.cols[[j]][[i]]
				} else {
					NA
				}
			})
			preds <- preds[!is.na(preds)] # Provide only predecessors that are set
			dp.cols[[j]][[i]] <- step.fun(i, j, preds, ...)
		}
		if (!get.matching) {
			## Delete column that we no longer need to compute new entries
			dp.cols[[j-1]] <- NULL
		}
	}
	res <- dp.cols[[nrow(T2)]][[nrow(T1)]]@value
	## Return matching, DP table values and the table itself, if requested.
	if (get.matching) {
		dp.table <- do.call(c, dp.cols)
		dim(dp.table) <- dim(dp.table) <- c(nrow(T1), nrow(T2))	# T1 indexes rows, T2 columns.
		attr(res, "matching") <- .matching(dp.table, steps)
		attr(res, "table.values") <- matrix(sapply(dp.table, function(entry) {
					as.numeric(entry) # Extracts value or passes through NA
				}), nrow(T1), nrow(T2))
		attr(res, "table") <- dp.table
	}
	res
}

".matching" <- function(dp.t, steps) {
	res <- matrix(NA, 2, sum(dim(dp.t))-1)
	pos <- c(nrow(dp.t), ncol(dp.t))
	l <- ncol(res)
	res[,l] <- pos

	entry <- dp.t[[pos[1],pos[2]]]
	while(!is.null(entry@pred)) {
		l <- l-1
		pos <- pos + steps[[entry@pred]]
		res[,l] <- pos
		entry <- dp.t[[pos[1],pos[2]]]
	}
	res[,l:ncol(res)]
}
