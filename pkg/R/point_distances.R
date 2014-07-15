## Several distance measures between points in spacetime,
## that can be used in the trajectories similarity measures.

## Each function takes two numeric vectors as input;
## The last element of each vector is regarded as time,
## the others as location coordinates.
## The input vectors should have identical lengths.

"Lp.norm" <- function (p=2) {
	if (p <= 0) {
		stop("The Lp norm is not defined for p <= 0")
	} else if (p < 1) {
		warning("The Lp norm is not convex for p < 1")
	}
	
	# The L_infinity norm is a special case
	if (is.infinite(p)) {
		return(function(x1, x2) {
			max(abs(x1[-length(x1)] - x2[-length(x2)]))
		})
	}
	function(x1, x2) {
		sum(abs(x1[-length(x1)] - x2[-length(x2)])^p)^(1/p)
	}
}
euclidian <- Lp.norm(2)
manhattan <- Lp.norm(1)
