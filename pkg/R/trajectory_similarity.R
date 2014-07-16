## Generic functions for computing the similarity of a set of trajectories.
## These bring the input into a standard format and invoke the functions
## that compute pairwise similarity for specific measures.

setGeneric("trajectory.similarity", function(trajectories, implementation, point.distance, ...) {
	standardGeneric("trajectory.similarity")
})

if ("move" %in% rownames(installed.packages()) && require(move)) {
	setMethod(f = "trajectory.similarity", 
			signature = c(trajectories = "MoveStack"),
			function(trajectories, implementation, point.distance, ...) {
		## Convert the MoveStack into a list of individual trajectories,
		## represented as matrices with coordinate and time columns
		trs <- lapply(split(monkey.tr), function(tr) {
				cbind(tr@coords, datetime=as.double(tr@timestamps)) 
		})
		trajectory.similarity(trs, implementation, point.distance, ...)
	})
}

setMethod(f = "trajectory.similarity", 
		signature = c(trajectories = "list", implementation=closure),
		function(trajectories, implementation, point.distance, ...) {
	trajectory.similarity(split(trajectories), implementation, point.distance, ...)
})
