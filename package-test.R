library(trajectorySim)
## The following line unloads the package and reloads it again, to incorporate any changes
detach("package:trajectorySim", unload=TRUE); library(trajectorySim)

## Test some random pointwise distance measures

## 2d locations plus time
p <- c(1,1,0)
q <- c(2,2,0)
euclidian(p, q)

Lp.norm()(p,q)    # Lp.norm() with default parameter is Euclidian
Lp.norm(1)(p,q)   # Manhattan
Lp.norm(Inf)(p,q) # L_infinity

## 3d locations plus time
p <- c(1,1,1,0)
q <- c(2,2,2,0)
euclidian(p,q)

Lp.norm()(p,q)    # Lp.norm() with default parameter is Euclidian
Lp.norm(1)(p,q)   # Manhattan
Lp.norm(Inf)(p,q) # L_infinity


data(example_traj, package="trajectorySim")
T1 <- example.traj[[2]]#[9:13,]
T2 <- example.traj[[4]]#[9:12,]


## Frechet distance is known to be sqrt(2) for this example
T1 <- matrix(c(0,0,2,2,0,2,2,0,1,2,3,4), 4)
T2 <- matrix(c(0,1,1,3,0,1,0,0,1,2,3,4), 4)

dist <- equal.time.pairwise(T1, T2)
dist <- discrete.frechet.pairwise(T1,T2, get.matching=T)
dist <- LCSS.pairwise(T1, T2, 1, 2, get.matching=T)
dist <- DTW.pairwise(T1, T2, get.matching=T)
dist
tab <- attr(dist, "table")

vals <- attr(dist, "table.values")
image(1:nrow(vals), 1:ncol(vals), vals)
points(t(attr(dist, "matching")))
lines(t(attr(dist, "matching")))

plot(rbind(T1, T2))
lines(T1, type="l", col="red")
lines(T2, type="l", col="blue")
apply(attr(dist, "matching"), 2, function(m) {
	lines(rbind(T1[m[1],], T2[m[2],]))
})

## Compute pairwise DTW distances between all four trajectories
DTW(example.traj)

frechet.decision.pairwise(T1, T2, 1)
frechet.pairwise(T1, T2)


## Testing for LIP
t1 <- matrix(c(0,1,1,2, 0,0,2,2, 0,1,2,3), 4)
t2 <- matrix(c(0,0,2,2, 0,1,1,2, 0,1,2,3), 4)
# t1 <- matrix(c(0,2,2,0,   0,1,-1,0, 0,1,2,3), 4)
# t2 <- matrix(c(1,-1,-1,1, 0,-1,1,0, 0,1,2,3), 4)
#t1 <- matrix(c(-1,4,0,0, 1,1,-3,2, 0,1,2,3), 4)
#t2 <- matrix(c(-1,2,1,1, 0,0,-1,2, 0,1,2,3), 4)
# t1 <- matrix(c(-1,4,2, 1,1,-1, 0,1,2), 3)
# t2 <- matrix(c(-1,2,1.5, 0,0,-0.5, 0,1,2), 3)

LIP(t1,t2, weight=F)


