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


# step.fun is a function that takes parameters
# 	i, j:   indices into T1 and T2
#	prev:   dyn prog table entries that precede this entry horizontally, vertically and/or diagonally
# It may assume the following variables to be defined in its environment:
# 	T1, T2: trajectories
# 	dp:     a function that provides the distance between two location entries
# it returns the table entry for (i,j)
"traj.sim.dp" <- function(T1, T2, step.fun, dp=euclidian) {
	dp.table <- as.list(rep(NA, nrow(T1)*nrow(T2)))
	dim(dp.table) <- c(nrow(T1), nrow(T2))	# T1 indexes rows, T2 columns.
	
	# Set up an environment for the step function to avoid unnecessary argument passing
	env <- new.env(parent=environment())
	assign("T1", T1, envir=env)
	assign("T2", T2, envir=env)
	assign("dp", dp, envir=env)
	environment(step.fun) <- env
	
	dp.table[[1,1]] <- step.fun(1, 1)
	for (i in 2:nrow(T1)) {
		dp.table[[i,1]] <- step.fun(i, 1, list(V=dp.table[[i-1,1]]))
	}
	for (j in 2:nrow(T2)) {
		dp.table[[1,j]] <- step.fun(1, j, list(H=dp.table[[1,j-1]]))
	}
	for (i in 2:nrow(T1)) {
		for (j in 2:nrow(T2)) {
			dp.table[[i,j]] <- step.fun(i, j, list(
					V=dp.table[[i-1,j]], H=dp.table[[i,j-1]], D=dp.table[[i-1,j-1]]))
		}
	}
	res <- dp.table[[nrow(T1),nrow(T2)]]@value
	## For now, always return matching, DP table values and the table itself for debugging.
	## TODO: Allow the user to choose what attributes to return and adapt the implementation to it.
	attr(res, "matching") <- .matching(dp.table)
	attr(res, "table.values") <- matrix(sapply(dp.table, function(entry) { entry@value }), nrow(T1), nrow(T2))
	attr(res, "table") <- dp.table
	res
}

".matching" <- function(dp.t) {
	res <- matrix(NA, 2, sum(dim(dp.t))-1)
	i <- nrow(dp.t)
	j <- ncol(dp.t)
	l <- ncol(res)
	res[,l] <- c(i,j)
	while(i > 1 || j > 1) {
		l <- l-1
		entry <- dp.t[[i,j]]
		p <- switch(entry@pred,
				V=c(i-1,j  ),
				H=c(i  ,j-1),
				D=c(i-1,j-1)
		)
		res[,l] <- p
		i <- p[1]
		j <- p[2]
	}
	res[,l:ncol(res)]
}
