library(igraph)

# There are 3^6 possible combinations of three species interacting in one of six ways
## none 0/0
## competition -/-
## mutualisim +/+
## predation +/-
## amensalism -/0
## commensalism +/0

n.mot <- expand.grid(0, c(1,-1,0), c(1,-1,0), c(1,-1,0), 0, c(1,-1,0), c(1,-1,0), c(1,-1,0), 0)
motmat <- as.matrix(n.mot)
head(n.mot)

subg <- lapply(1:nrow(n.mot), function(x) matrix(unlist(n.mot[x,]), ncol =3))
iconn <- sapply(lapply(subg, function(x) graph.adjacency(abs(x))), is.connected)
subg2 <- apply(motmat, 1, paste, collapse = " ")

#######################################################################
### create network

build <- function(N, C){
  mat <- get.adjacency(erdos.renyi.game(N, p.or.m = C, type = "gnp", directed = FALSE), sparse = FALSE)
  return(mat)
}

# Define interaction types
typ_int <- function(x){
  if(x == "none"){
    return(c(0, 0))
  }else if(x == "pred"){
    return(sample(c(-1, 1)))
  }else if(x == "comp"){
    return(c(-1, -1))
  }else if(x == "mut"){
    return(c(1, 1))
  }else if(x == "amen"){
    return(c(-1, 0))
  }else if(x == "comm"){
    return(sample(c(0, 1)))
  }
}

# sample interaction types for the network
get_ints2 <- function(mat, ...){
  alli <- cbind(mat[lower.tri(mat)], t(mat)[lower.tri(mat)])
  alli2 <- alli
  
  isam <- sample(c("pred", "comp", "mut", "amen", "comm"), sum(alli[,1] == 1), replace = TRUE, ...)
  ilist <- rep("none", nrow(alli))
  ilist[which(alli[,1] == 1)] <- isam
  
  allints <- sapply(ilist, typ_int)
  
  m <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  
  m[lower.tri(m)] <- allints[1,]
  m2 <- t(m)
  m2[lower.tri(m)] <- allints[2,]
  return(t(m2))
  
}

# get distribution of interaction types from matrix
itypes <- function(x){
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  comp <- sum(i1 < 0 & i2 < 0)
  mut <- sum(i1 > 0 & i2 > 0)
  pred <- sum(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
  amens <- sum(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
  comm <- sum(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
  
  return(c(comp = comp, mut = mut, pred = pred, amens = amens, comm = comm))
}


imat <- get_ints2(build(20, .3))


dim(t(combn(20, 3)))

combin <- c()
for(i in 10:250){
  combin[i-9] <- nrow(t(combn(i, 3)))
}
plot(10:250, combin)
max(combin)


N <- seq(10,300,10)
Co <- seq(.1, .8, .1)

par <- expand.grid(N, Co)
timed <- c()
t.unit <- c()
timed.all <- c()
for(i in 1:nrow(par)){
  t0 <- Sys.time()
  imat <- get_ints2(build(par[i,1], par[i,2]))
  trips <- t(combn(nrow(imat), 3))
  t1 <- Sys.time()
  mot <- vector(length = nrow(trips))
  for(j in 1:nrow(trips)){
    subi <- paste(as.vector(imat[trips[j,], trips[j,]]), collapse = " ")
    mot[j] <- which(subg2 %in% subi)
  }
  t2 <- Sys.time()
  timed[i] <- t2 - t1
  t.unit[i] <- units(t2-t1)
  timed.all[i] <-t2 - t0
  (cat(i,";", timed[i], t.unit[i], ";", "N =", par[i,1],";", "C =", par[i,2], "\n"))
}

timed.2 <- timed
timed.2[t.unit == "mins"] <- timed[t.unit == "mins"]*60

# timing rises exponentially with size
plot(timed.2~par[,1])
# no effect of connectance
plot(timed.2~par[,2])
# does not seem to be a strong effect of connectance within size classes
plot(timed.2[par[,1] == 300]~par[par[,1] == 300,2])


df1 <- data.frame(N = par[,1], Co = par[,2], time = timed, unit = t.unit, secs = timed.2)
head(df1)

head(t(sapply(subg, itypes)))


tatoosh <- read.csv("./data/tatoosh.csv", header = F)
tatoosh <- as.matrix(tatoosh)
diag(tatoosh) <- 0
trips <- t(combn(nrow(tatoosh), 3))
t1 <- Sys.time()
mot <- vector(length = nrow(trips))
for(j in 1:nrow(trips)){
  subi <- paste(as.vector(tatoosh[trips[j,], trips[j,]]), collapse = " ")
  mot[j] <- which(subg2 %in% subi)
}
t2 <- Sys.time()