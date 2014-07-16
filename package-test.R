library(trajectorySim)
## The following line unloads the package and reloads it again, to incorporate any changes
detach("package:trajectorySim", unload=TRUE); library(trajectorySim)

## Test some random pointwise distance measures

## 2d locations plus time
p <- c(1,1,0)
q <- c(2,2,0)
euclidian(p, q)

Lp.norm()(p,q)
Lp.norm(1)(p,q)
Lp.norm(Inf)(p,q)

## 3d locations plus time
p <- c(1,1,1,0)
q <- c(2,2,2,0)
euclidian(p,q)

Lp.norm()(p,q)
Lp.norm(1)(p,q)
Lp.norm(Inf)(p,q)

data(example_traj, package="trajectorySim")
T1 <- example.traj[[2]]#[9:13,]
T2 <- example.traj[[4]]#[9:12,]

T1 <- matrix(c(0,0,2,2,0,2,2,0,1,2,3,4), 4)
T2 <- matrix(c(0,1,1,3,0,1,0,0,1,2,3,4), 4)

equal.time (T1, T2)

dtw <- DTW(list(T1, T2))[1,2]
dtw
tab <- attr(dtw, "table")

vals <- attr(dtw, "table.values")
image(1:nrow(vals), 1:ncol(vals), vals)
points(t(attr(dtw, "matching")))
lines(t(attr(dtw, "matching")))!

plot(rbind(T1, T2))
lines(T1, type="l", col="red")
lines(T2, type="l", col="blue")
apply(attr(dtw, "matching"), 2, function(m) {
	lines(rbind(T1[m[1],], T2[m[2],]))
})

## Compute pairwise DTW distances between all four trajectories
DTW(example.traj)

frechet.decision.pairwise(T1, T2, 50)
frechet.pairwise(T1, T2)