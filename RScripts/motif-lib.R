n.mot <- expand.grid(0, c(1,-1,0), c(1,-1,0), c(1,-1,0), 0, c(1,-1,0), c(1,-1,0), c(1,-1,0), 0)
head(n.mot)

subg <- lapply(1:nrow(n.mot), function(x) matrix(n.mot[x,], ncol =3))
