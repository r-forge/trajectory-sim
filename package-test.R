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
equal.time (example.traj[[1]], example.traj[[2]])