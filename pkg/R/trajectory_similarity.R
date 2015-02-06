## Generic functions for computing the similarity of a set of trajectories.
## These bring the input into a standard format and invoke the functions
## that compute pairwise similarity for specific measures.

setGeneric("trajectory.similarity", function(trajectories, implementation, ..., symmetric = FALSE, diagonal=NA) {
	standardGeneric("trajectory.similarity")
})

if ("move" %in% rownames(installed.packages()) && require(move)) {
	setMethod(f = "trajectory.similarity", 
			signature = c(trajectories = "MoveStack"),
			function(trajectories, ..., symmetric = FALSE, diagonal=NA) {
		## Convert the MoveStack into a list of individual trajectories,
		## represented as matrices with coordinate and time columns
		trs <- lapply(split(trajectories), function(tr) {
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
	
	## Figure out for which pairs of trajectories to compute similarity scores
	ij.compute <- lapply(1:length(trajectories), function (i) {
		j.start <- if (!symmetric) { 1 } else { i }
		j.vals <- j.start:length(trajectories)
		if (!is.na(diagonal)) { j.vals <- j.vals[j.vals != i] }

		sapply(j.vals, function(j) {
			c(i,j)
		})
	})
	ij.compute <- matrix(unlist(ij.compute), nrow=2)
	ij.compute <- split(ij.compute, col(ij.compute))
	
	similarities <- mclapply(ij.compute, function(ij, ...) {
		i <- ij[1]
		j <- ij[2]
		
		implementation(trajectories[[i]], trajectories[[j]], ...)
	}, mc.cores=detectCores(logical=T), ...)
	## Store results in matrix
	for (k in 1:length(ij.compute)) {
		ij <- ij.compute[[k]]
		res[ij[1], ij[2]] <- similarities[[k]]
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

