## Implements the equal-time distance between a pair of trajectories
## Requires that each trajectory has the same temporal extent

"equal.time" <- function(trajectories, method=c("mean","max","euclidian"), pd=euclidian) {
	trajectory.similarity(trajectories, implementation=equal.time.pairwise, 
			pd=pd, method=match.arg(method), symmetric=TRUE, diagonal=0)
}

"equal.time.pairwise" <- function(T1, T2, method=c("mean","max","euclidian"), pd=euclidian) {
	method <- match.arg(method)

	t1 <- T1[,ncol(T1)]
	t2 <- T2[,ncol(T2)]
	ts <- c(t1, t2)
	ts <- sort(ts)
	
	dist <- mapply(pd, 
			data.frame(position(T1, ts)), 
			data.frame(position(T2, ts)))

	if (method=="euclidian" && !identical(t1, t2)) {
		stop("euclidian is only allowed for trajectories with identical sampling")
	}
	switch(method,
		mean=mean(dist),
		max=max(dist),
		euclidian=sqrt(0.5*sum(dist^2))) ## Every time step counted twice
}

## Linearly interpolates the position of Tr at the times in ts
"position" <- function(Tr, ts) {
	## Find where the times are in the trajectory
	i <- findInterval(ts, Tr[,ncol(Tr)], rightmost.closed=TRUE)
	
	mapply(function(t, i) {
		## Interpolate linearly between observation i and i+1 in Tr
		if (i <= 0 || i >= nrow(Tr)) { return(c(rep(NA, ncol(Tr)-1), t)) }
		
		a <- (t - Tr[i,ncol(Tr)]) / (Tr[i+1,ncol(Tr)] - Tr[i,ncol(Tr)])
		(1-a) * Tr[i,] + a * Tr[i+1,]
	}, ts, i)
}

