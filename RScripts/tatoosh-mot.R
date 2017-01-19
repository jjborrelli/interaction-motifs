library(igraph)

tatoosh <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/rKeystone/tatoosh.csv", header = F)
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
test <- sapply(mlst, function(x) idmotifs(graph.adjacency(abs(x))))
dim(test)
rowSums(test)
t1 <- Sys.time()
t1-t0


## Create matrices for each of the possible 3 species configurations

s1<-matrix(c(0,1,0,-1,0,1,0,-1,0),nrow=3,ncol=3)
s2<-matrix(c(0,1,1,-1,0,1,-1,-1,0),nrow=3,ncol=3)

plot(graph.adjacency(abs(s1)))
plot(graph.adjacency(abs(s2)))

el2i <- get.edgelist(graph.adjacency(abs(s1)))
el3i <- get.edgelist(graph.adjacency(abs(s2)))

inames <- c("competition", "mutualism", "predation", "ammensalism", "commensalism")
i3opts <- expand.grid(inames, inames, inames)
i2opts <- expand.grid(inames, inames)

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
  }
}


glist <- list()
for(i in 1:nrow(i2opts)){
  el2s <- cbind(el2i, 0)
  el2s[1:2, 3] <- isign(i2opts[i,1])
  el2s[3:4, 3] <- isign(i2opts[i,2])
  
  glist[[i]] <- el2s
}

glist2 <- list()
for(i in 1:nrow(i3opts)){
  el3s <- cbind(el3i, 0)
  el3s[c(1,3), 3] <- isign(i3opts[i,1])
  el3s[c(2,5), 3] <- isign(i3opts[i,2])
  el3s[c(4,6), 3] <- isign(i3opts[i,3])
  
  glist2[[i]] <- el3s
}

glist

i2list <- lapply(glist, function(x){
  m <- matrix(0, nrow = 3, ncol = 3)
  for(y in 1:nrow(x)){
    m[x[y,1], x[y,2]] <- x[y,3]
  }
  return(m)
})

i2list.iso <- lapply(i2list, t)

i3list <- lapply(glist2, function(x){
  m <- matrix(0, nrow = 3, ncol = 3)
  for(y in 1:nrow(x)){
    m[x[y,1], x[y,2]] <- x[y,3]
  }
  return(m)
})

i3list.iso <- lapply(i3list, t)

sgl1 <- t(sapply(c(i2list, i3list), as.vector))
sgl2 <- t(sapply(c(i2list.iso, i3list.iso), as.vector))


tatmotv <- t(sapply(mlst, as.vector))
tatmotv <- apply(tatmotv, 1, paste0, collapse = ",")

isoA <- apply(sgl1, 1, paste0, collapse = ",")
isoB <- apply(sgl2, 1, paste0, collapse = ",")

id <- c()
for(i in 1:length(tatmotv)){
   wA <- isoA %in% tatmotv[15]
   wB <- isoB %in% tatmotv[i]
   
   if(any(wA)){
     id[i] <- which(wA)
   }else if(any(wB)){
     id[i] <- which(wB)
   }else{
     id[i] <- NA
   }
}

sum(is.na(id))
which(is.na(id[1:20]))
