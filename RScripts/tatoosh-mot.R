library(igraph)
library(NetIndices)
library(deSolve)

tatoosh <- read.csv("./data/tatoosh.csv", header = F)
tatoosh <- as.matrix(tatoosh)

getmotlst <- function(adj){
  com <- combn(nrow(adj), 3)
  
  diag(adj) <- 0
  tbt <- lapply(1:ncol(com), function(x){adj[com[,x], com[,x]]})
  
  conn <- sapply(tbt, function(x) is.connected(graph.adjacency(abs(x))))
  
  m3 <- tbt[conn]
  
  return(m3)
}

idmotifs <- function(g){
  mo <- motifs(g)
  mo2 <- mo[c(5, 8, 12, 3, 7, 14, 9, 10, 6, 13, 16, 15, 11)]  
  names(mo2) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
  
  return(mo2)
}


t0 <- Sys.time()
mlst <- getmotlst(tatoosh)
t1 <- Sys.time()
t1-t0


## Create matrices for each of the possible 3 species configurations

s1<-matrix(c(0,1,0,-1,0,1,0,-1,0),nrow=3,ncol=3)
s2<-matrix(c(0,1,1,-1,0,1,-1,-1,0),nrow=3,ncol=3)

plot(graph.adjacency(abs(s1)))
plot(graph.adjacency(abs(s2)))

el2i <- get.edgelist(graph.adjacency(abs(s1)))
el3i <- get.edgelist(graph.adjacency(abs(s2)))

inames <- c("competition", "mutualism", "predation", "predation.alt", "ammensalism", "ammensalism.alt", "commensalism", "commensalism.alt", "none")
i3opts <- expand.grid(inames, inames, inames)
i2opts <- expand.grid(inames, inames)
i3opts <- i3opts[!apply(i3opts, 1, function(x) sum(x == "none")) >= 2,]

isign <- function(x){
  if(x == "competition"){
    return(c(-1, -1))
  }else if(x == "mutualism"){
    return(c(1, 1))
  }else if(x == "predation"){
    return(c(1, -1))
  }else if(x == "ammensalism"){
    return(c(0, -1))
  }else if(x == "commensalism"){
    return(c(0, 1))
  }else if(x == "predation.alt"){
    return(c(-1, 1))
  }else if(x == "ammensalism.alt"){
    return(c(-1, 0))
  }else if(x == "commensalism.alt"){
    return(c(1, 0))
  }else if(x == "none"){
    return(c(0, 0))
  }
}


glist <- list()
for(i in 1:nrow(i3opts)){
  el3s <- cbind(el3i, 0)
  el3s[c(1,3), 3] <- isign(i3opts[i,1])
  el3s[c(2,5), 3] <- isign(i3opts[i,2])
  el3s[c(4,6), 3] <- isign(i3opts[i,3])
  
  glist[[i]] <- el3s
}

i3list <- lapply(glist, function(x){
  m <- matrix(0, nrow = 3, ncol = 3)
  for(y in 1:nrow(x)){
    m[x[y,1], x[y,2]] <- x[y,3]
  }
  return(m)
})



sgl <- t(sapply(i3list, as.vector))


tatmotv <- t(sapply(mlst, as.vector))
tatmotv <- apply(tatmotv, 1, paste0, collapse = ",")

isoA <- apply(sgl, 1, paste0, collapse = ",")

id <- c()
for(i in 1:length(tatmotv)){
   wA <- isoA %in% tatmotv[i]

   if(any(wA)){
     id[i] <- which(wA)
   }else{
     id[i] <- NA
   }
}

sum(is.na(id))

#########################################################
#########################################################
#########################################################

tmot <- as.numeric(names(table(id)))
plot(as.vector(table(id)), sapply(lapply(alleq, "[[", 1), mean)[tmot])

#########################################################
#########################################################
#########################################################

# Basic Lotka-Volterra model
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
    
    list(dB)
  })
}

# Function to detect extinction (prevents negative abundances)
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
   
    return(c(states))
  })
}

ext2 <- function (times, states, parms){
  with(as.list(states), {
    states[states <= 10^-10] <- 0
    
    states <- states + states * (rbeta(1, .3, 6)*sample(c(-1,1), 1))
    
    return(c(states))
  })
}

sim_community <- function(times, state, parms, eq = lvmod, ex = ext1){
  out <- ode(state, times, parms = parms, func = eq, events = list(func = ex, time =  times))
  return(out[,-1])
}

ran.unif <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 1)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- -1 #runif(1, -1, 0)
  return(newmat)
}

eqS <- function(mat, iter, extf = ext1){
  eq <- c()
  plist <- list()
  for(i in 1:iter){
    cond <- FALSE
    while(!cond){
      istate <- runif(3, 0, 1)
      igr <- runif(3, .1, 1)
      par <- list(alpha = igr, m = ran.unif(mat))
      
      sc <- sim_community(1:1000, istate, par, ex = extf)
      
      cond <- nrow(sc) == 1000
    }
    plist[[i]] <- cbind(alpha = par$alpha, par$m) 
    eq[i] <- ifelse(nrow(sc) == 1000,sum(sc[1000,] != 0), NA)
  }
  return(list(eq, plist))
}

t0 <- Sys.time()
#alleq <- lapply(i3list, eqS, iter = 1000)
alleq <- list()
for(i in 1:length(i3list)){
  alleq[[i]] <- eqS(i3list[[i]], iter = 1000, extf = ext1)
  cat(rep(c(i, "\n"), 10))
}
t1 <- Sys.time()
t1-t0


tA <- Sys.time()
#alleq <- lapply(i3list, eqS, iter = 1000)
alleq1 <- list()
for(i in 1:length(i3list)){
  alleq1[[i]] <- eqS(i3list[[i]], iter = 1000, extf = ext2)
  cat(rep(c(i, "\n"), 10))
}
tB <- Sys.time()
tB-tA




qstab <- sapply(lapply(alleq, "[[", 1), function(x) sum(x == 3)/1000)
mab <- table(id)













# save.image("D://jjborrelli/motifINTS/neq-tri-int.Rdata")
# load("D://jjborrelli/motifINTS/neq-tri-int.Rdata")