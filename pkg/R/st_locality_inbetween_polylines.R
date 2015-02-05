library(sp)

"simple.LIP" <- function(t1, t2) {
	## concatenate t1 with t2 reversed to construct the polygon
	p <- rbind(t1, t2[nrow(t2):1,])[,1:2] # drop timestamps
	
	# Return the polygon's area (abs accounts for possibly clockwise vertex order)
	0.5 * sum((p[,1] + p[c(2:nrow(p), 1)]) 
	              * (p[c(2:nrow(p), 1), 2] - p[, 2]))
}

"gen.LIP" <- function(t1,t2) {
	if (nrow(t1) == 1) { tmp <- t1; t1 <- t2; t2 <- tmp }
	if (nrow(t1) == 1) { return(0) } ## both trajectories are degenerated to points
	print(t1)
	print(t2)
	
	## TODO: split the trajectories in parts which generate simple polygons
	s1 <- 1
	s2 <- 1
	
	i1 <- 2
	i2 <- 1
	
	res <- 0
	
	## Establish in which direction the line will move initially
	dir <- .line.side(t1[s1,], t2[s2,], t1[i1,])
	if (dir == 0) { dir <- 1 } # Break degenerate positions arbitrarily
	
	## Walk along either of the trajectories as long as we can
	while (i1 < nrow(t1) || i2 < nrow(t2)) {
		print(cat(i1, "/", nrow(t1), ",", i2, "/", nrow(t2)))
		if (i1 < nrow(t1)) {
			d1 <- .line.side(t1[i1,], t2[i2,], t1[i1+1,])
			if (d1 * dir >= 0) { i1 <- i1 + 1 }
		} else { d1 <- -dir }
		
		if (i2 < nrow(t2)) {
			d2 <- .line.side(t1[i1,], t2[i2,], t2[i2+1,])
			if (d2 * dir >= 0) { i2 <- i2 + 1 }
		} else { d2 <- -dir }
		
		if (d1 == -dir && d2 == -dir) {
			## Break the polygon here
			res <- res + simple.LIP(t1[s1:i1,,drop=F], t2[s2:i2,,drop=F])
			
			s1 <- i1
			s2 <- i2
			dir <- -dir ## reverse direction
			
			print (cat("split at ", s1, s2))
		}
	}
	
	
	## Add the last polygon piece and return
	abs(res + simple.LIP(t1[s1:i1,,drop=F], t2[s2:i2,,drop=F]))
}

"LIP" <- function(t1, t2, weighted=TRUE) {
	## Split trajectories at pairwise intersections
	is <- red_blue.intersections(t1, t2)
	is <- is[order(is[,3]), , drop=F]
	## Coinciding endpoints are not proper intersections; drop them
# 	if (is[1,3] == t1[1,3] && is[1,4] == t2[1,3]) {
# 		is <- is[-1, ,drop=F]
# 	}
# 	if (is[nrow(is),3] == t1[nrow(t1),3] && is[nrow(is),4] == t2[nrow(t2),3]) {
# 		is <- is[-nrow(is), ,drop=F]
# 	}
	
	t1s <- split.at.intersections(t1, is[,-4, drop=F]) # drop timestamps for the other trajectory
	t2s <- split.at.intersections(t2, is[,-3, drop=F])
	
	# Compute gen.LIP for corresponding pairs of subtrajectories between intersections
	if (weighted) {
		weights <- (sapply(t1s, .tr.length) + sapply(t2s, .tr.length)) / 
				(.tr.length(t1) + .tr.length(t2))
	} else {
		weights <- 1
	}
	sum(weights*mapply(gen.LIP, t1s, t2s))
}

"red_blue.intersections" <- function(t1, t2) {
	## Simple O(n^2) implementation
	## TODO: how to deal with intersections of more than 2 segments?
	
	time1 <- t1[,3]
	t1 <- t1[,1:2] # drop timestamps
	time2 <- t2[,3]
	t2 <- t2[,1:2]
	is <- matrix(double(0), ncol=4)
	colnames(is) <- c('x','y','time1','time2')
	
	for (i in 1:(nrow(t1)-1)) {
		for (j in 1:(nrow(t2)-1)) {
			p <- segment.intersection.pairwise(t1[i:(i+1),], t2[j:(j+1),])
			if (!is.null(p)) {
				# Compute the times of the intersection
				p[3] <- time1[i] + p[3] * (time1[i+1]-time1[i])
				p[4] <- time2[j] + p[4] * (time2[j+1]-time2[j])
				is <- rbind(is, p)
			}
		}
	}
	is
}

"segment.intersection.pairwise" <- function(s1, s2) {
	## Segments are given as 2x2 matrices with one point per row
	## WARNING: Does not report intersections/overlap between collinear segments
	
	num <- det(s1) * (s2[1,]-s2[2,]) - 
			det(s2) * (s1[1,]-s1[2,])
	den <- (s1[1,1]-s1[2,1]) * (s2[1,2]-s2[2,2]) -
			(s1[1,2]-s1[2,2]) * (s2[1,1]-s2[2,1])
	p <- num/den # p is the intersection between both lines

	# If the segments are parallel, p has at least one coordinate at infinity
	
	## Check if p is within the segments
	s1r <- apply(s1, 2, range)
	s2r <- apply(s2, 2, range)
	if (all(is.finite(p)
			& p >= pmax(s1r[1,], s2r[1,]) 
			& p <= pmin(s1r[2,], s2r[2,]))) {
		## Find where on the segments p lies
		a1 <- (p-s1[1,])/(s1[2,]-s1[1,])
		a1 <- a1[is.finite(a1)][1] ## All finite parameters should be equal
		a2 <- (p-s2[1,])/(s2[2,]-s2[1,])
		a2 <- a2[is.finite(a2)][1]
		
		c(p, a1, a2)
	} else {
		NULL
	}
}

"split.at.intersections" <- function (tr, is) {
	is <- is[order(is[,3]), , drop=F] # Sort intersections chronologically
	
	i <- 1
	res <- list(matrix(double(0), 3))
	
	for (j in 1:nrow(tr)) {
		while (i <= nrow(is) && tr[j,3] >= is[i,3]) {
			res[[i]] <- cbind(res[[i]], is[i,])
			i <- i+1
			res[[i]] <- matrix(is[i-1,], 3)
		}
		res[[i]] <- cbind(res[[i]], tr[j,])
	}
	
	## the result contains transposed trajectories, revert
	## Make sure there are no duplicate points (in space and time)
	lapply(lapply(res, t), unique)
}

".tr.length" <- function(tr) {
	tr <- tr[,1:2, drop=F] # drop time information
# 	print(apply((tr[-1,] - tr[-nrow(tr),])^2, 1)
	sum(sqrt(apply((tr[-1, , drop=F] - tr[-nrow(tr), , drop=F])^2, 1, sum)))
}

## returns on which side of the line through l1 and l2 the point p lies
## -1 for right, 0 for on the line, 1 for left
".line.side" <- function(l1, l2, p) {
	sign( (l2[1] * p[2]) - (l1[1] * (p[2] - l2[2])) + (l1[2] * (p[1] - l2[1])) )
}

line.side(c(0,1), c(1,0), c(0,0))

t1 <- matrix(c(0,1,1,2, 0,0,2,2, 0,1,2,3), 4)
t2 <- matrix(c(0,0,2,2, 0,1,1,2, 0,1,2,3), 4)
# t1 <- matrix(c(0,2,2,0,   0,1,-1,0, 0,1,2,3), 4)
# t2 <- matrix(c(1,-1,-1,1, 0,-1,1,0, 0,1,2,3), 4)
t1 <- matrix(c(-1,4,0,0, 1,1,-3,2, 0,1,2,3), 4)
t2 <- matrix(c(-1,2,1,1, 0,0,-1,2, 0,1,2,3), 4)
# t1 <- matrix(c(-1,4,2, 1,1,-1, 0,1,2), 3)
# t2 <- matrix(c(-1,2,1.5, 0,0,-0.5, 0,1,2), 3)

LIP(t1,t2, weight=F)

