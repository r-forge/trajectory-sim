## Generic functions for computing the similarity of a set of trajectories.
## These bring the input into a standard format and invoke the functions
## that compute pairwise similarity for specific measures.

setGeneric("trajectory.similarity", function(trajectories, implementation, ..., symmetric = FALSE, diagonal=NA) {
	standardGeneric("trajectory.similarity")
})

if ("move" %in% rownames(installed.packages()) && require(move)) {
	setMethod(f = "trajectory.similarity", 
			signature = c(trajectories = "MoveStack"),
			function(trajectories, ...) {
		## Convert the MoveStack into a list of individual trajectories,
		## represented as matrices with coordinate and time columns
		trs <- lapply(split(monkey.tr), function(tr) {
				cbind(tr@coords, datetime=as.double(tr@timestamps)) 
		})
		trajectory.similarity(trs, implementation, ..., symmetric=symmetric, diagonal=diagonal)
	})
}

setMethod(f = "trajectory.similarity", 
		signature = c(trajectories = "list", implementation="function"),
		function(trajectories, implementation, ..., symmetric = FALSE, diagonal=NA) {
	res <- matrix(NA, length(trajectories), length(trajectories),
			dimnames=list(names(trajectories), names(trajectories)))
	
	for (i in 1:(length(trajectories))) {
		## Figure out for which pairs (i,j) we actually need to compute pairwise distances
		## If the measure is symmetric or has a default for identical trajectories,
		## some can be skipped.
		j.start <- if (!symmetric) { 1 } else { i }
		j.vals <- j.start:length(trajectories)
		if (!is.na(diagonal)) { j.vals <- j.vals[j.vals != i] }
		
		for (j in j.vals) {
			res[i,j] <- implementation(trajectories[[i]], trajectories[[j]], ...)
		}
	}
	
	if (symmetric) {
		# Copy the upper triangle of the result into the lower triangle while transposing
		res[lower.tri(res)] <- t(res)[lower.tri(res)]
	}
	if (!is.na(diagonal)) {
		# Set the diagonal entries to a default value(s)
		diag(res) <- diagonal
	}
	res
})

