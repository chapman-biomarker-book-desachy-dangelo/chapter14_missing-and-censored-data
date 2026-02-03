meanbygpf <- function(transxx,Ut,inveffectif){
  meanbygp2<-Rfast::mat.mult(transxx,Ut)
  meanbygp2<-Rfast::transpose(meanbygp2)
  meanbygp2<-meanbygp2*inveffectif
  return(meanbygp2)
}

mysvd <- function (don, q) {
  #from fast.svd (package corpcor)
  n <- dim(don)[1]
  p <- dim(don)[2]
  EDGE.RATIO <- 2
  if(EDGE.RATIO * n < p){
    don.tr <- transpose(don)
    don.tmp <- mat.mult(don,don.tr)
    svd.tmp <- La.svd(don.tmp,nu=q,nv=0)
    tol <- n * svd.tmp$d[1] * .Machine$double.eps
    Positive <- svd.tmp$d[seq.int(q)] > tol
    svd.tmp$d <- sqrt(svd.tmp$d[seq.int(q)[Positive]])
    svd.tmp$u <- submatrix(svd.tmp$u,rowStart = 1,rowEnd = n,colStart = 1,colEnd = sum(Positive))
    dontru <- mat.mult(don.tr,svd.tmp$u)
    svd.tmp$v <- mat.mult(dontru,
                          diag(1/svd.tmp$d,nrow = length(svd.tmp$d)))
    
  }else if  (n > EDGE.RATIO * p) {
    don.tmp <- Crossprod(don,don)
    svd.tmp <- svd(don.tmp,nu=0,nv=q)
    tol <- p * s$d[1] * .Machine$double.eps
    Positive <- svd.tmp$d[seq.int(q)] > tol
    svd.tmp$d <- sqrt(svd.tmp$d[seq.int(q)[Positive]])
    svd.tmp$v <- submatrix(svd.tmp$v,rowStart = 1,rowEnd = n,colStart = 1,colEnd = sum(Positive))
    svd.tmp$u <- mat.mult(don,
                          mat.mult(svd.tmp$v,
                                   diag(1/svd.tmp$d,
                                        nrow=length(svd.tmp$d))))
  }else {
    svd.tmp <- svd(don, nu=q,nv=q)
  }
  return(svd.tmp)
}

rkm.intern <- function(xx, k, q, U=NULL, thres = 10^(-3), maxiter = 100,printflag=FALSE){
  #initialisation
  xx.df <- as.data.frame(xx)
  transxx <- transpose(xx)
  Ut <- U;Utvec <- rowMaxs(Ut)
  effectif <- Table(Utvec)
  inveffectif <- 1/effectif
  sqrtinveffectif <- sqrt(inveffectif)
  critrkm <- rep(NA,maxiter)
  critrkm[1] <- Inf
  
  continue <- TRUE
  nbiter <- 1
  while(continue){
    if(printflag){cat(nbiter,"...")}
    #mise à jour
    
    ##A (pxq)
    meanbygp<-meanbygpf(transxx=transxx,Ut = Ut,inveffectif = inveffectif)
    res.svd <- mysvd(meanbygp*sqrt(effectif),q=q) #svd[(U'U)^(-1/2)U'X]
    
    prodsvd <- transpose(res.svd$u*sqrtinveffectif)*res.svd$d[seq.int(q)]#SigmaQ'(U'U)^(-1/2)
    prodsvd <- mat.mult(res.svd$v,prodsvd)#PSigmaQ'(U'U)^(-1/2)
    prodsvd<-transpose(prodsvd)#(U'U)^(-1/2)QSigmaP'
    
    At.tmp <- mysvd(prodsvd,
                    q = q)$v#svd[(U'U)^(-1/2)QSigmaP']
    
    ## expression de G dans A
    Attransxx.tmp <- mat.mult(xx,At.tmp)
    Gt.tmp<-mat.mult(Tcrossprod(diag(inveffectif),Ut),
                     Attransxx.tmp)
    ##U (nxk)
    Utvec.tmp <- colMins(dista(xnew = Gt.tmp,x = Attransxx.tmp))
    cond0 <- length(unique(Utvec.tmp))==k#pas de classes vides
    if(cond0){
      At<- At.tmp
      Attransxx <- Attransxx.tmp
      Utvec <- Utvec.tmp
      rm(list = c("At.tmp","Attransxx.tmp","Utvec.tmp","Gt.tmp"))
    }else if(!cond0 & (nbiter==1) ){
      At<- At.tmp
      Gt <- Gt.tmp
      UtGt<-mat.mult(Ut,Gt)
      recons.trans <- Tcrossprod(At,UtGt)
      break()
    }else{
      break()
    }
    
    Ut<-caret:::class2ind(as.factor(Utvec))
    effectif <- Table(Utvec)
    inveffectif <- 1/effectif
    sqrtinveffectif <- sqrt(inveffectif)
    ##F
    Gt <- mat.mult(Tcrossprod(diag(inveffectif),Ut),
                   Attransxx)
    #critere
    UtGt<-mat.mult(Ut,Gt)
    recons.trans <- Tcrossprod(At,UtGt)
    critrkm[nbiter+1] <-sum((transxx - recons.trans)^2)
    #convergence
    if(!cond0){warning(paste0(k-length(unique(Utvec)),"classe(s) vide(s)"))}
    cond1 <- nbiter < maxiter
    if(!cond1){warning(paste("pas de convergence en",maxiter,"iterations"))}
    ratio <- ifelse(nbiter>1,((critrkm[nbiter+1]-critrkm[nbiter])/critrkm[nbiter]),-Inf)
    cond2 <- (ratio <=0) & (abs(ratio)>=thres)#decroissance et pas suffisant
    continue <- cond0 & cond1 & cond2
    nbiter <- nbiter + 1
  }
  if(printflag){cat("done!\n")}
  res.out <- list(U=Ut, A=At, G=Gt,crit=critrkm, cluster=Utvec,recons.trans = recons.trans)
  return(res.out)
}


rkm <- function(xx, k, q, thres = 10^(-3), maxiter = 100,
                method.init=c("kmpp","sample","maximin","MacQueen"),
                maxiter.init=5,
                printflag = TRUE,
                nnodes=1,
                axes=NULL,
                part.init = NULL){
  if(is.data.frame(xx)){stop("xx must be a matrix")}
  if((q>(k-1))){stop("q is too large compared to k")}
  
  if (printflag){
    departtime<-Sys.time()
    cat("initialization...")}
  #init for U
  if(is.null(part.init)){
    part.init <- lapply(method.init,initrkm,k=k,xx=xx)
    names(part.init) <- paste(seq.int(length(method.init)),method.init)
    garde <- which(!sapply(part.init,is.null))
    
    part.init <- part.init[garde][!duplicated(part.init[garde])]
  }#else{
  U <- lapply(part.init,FUN = function(yy){
    res.out <- caret:::class2ind(as.factor(yy))
    return(res.out)
  }
  )
  # }
  
  
  if (printflag){
    cat(paste0("\n U done! [",round(difftime(Sys.time(),departtime,units="secs"),2)," sec]"))
    horloge<-Sys.time()
    cat("\n pre-run for RKM...")
  }
  
  #rkm
  if(nnodes==1){
    res.init <- lapply(U,FUN=function(U,xx,k,q,thres,maxiter.init){
      res.out <- rkm.intern(xx=xx, k=k, q=q, U=U, thres = thres, maxiter = maxiter.init)
      return(res.out)
    },xx=xx,k=k,q=q,thres=thres,maxiter.init=maxiter.init)
  }else{
    nnodes.intern<-min(nnodes, length(U),detectCores(logical = FALSE))
    cl <- makeCluster(nnodes.intern, type = "PSOCK")
    
    clusterEvalQ(cl, library("Rfast"))
    clusterExport(cl,
                  list("U","xx","k","q","thres","maxiter.init","rkm.intern","mysvd","meanbygpf"),
                  envir = environment())
    
    
    res.init <- parLapply(cl,U,fun=function(U,xx,k,q,thres,maxiter.init){
      res.out <- rkm.intern(xx, k, q, U=U, thres = thres, maxiter = maxiter.init)
      return(res.out)
    },xx=xx,k=k,q=q,thres=thres,maxiter.init=maxiter.init)
    
    
    stopCluster(cl)
  }
  
  #best init
  mincrit<-sapply(res.init,FUN = function(res){
    yy <- na.omit(res$crit)
    res.out <- yy[length(yy)]
    return(res.out)
  }, simplify =TRUE)
  
  winner <- res.init[[which.min(mincrit)]]
  if(printflag){
    cat(paste0("done! [",round(Sys.time()-horloge,2)," sec]"))
    cat(paste0("\n Elapsed time for initialization...[",round(difftime(Sys.time(),departtime,units="secs"),2)," sec]"))
    horloge <- Sys.time()
    cat("\nRKM is runing using the best starting values...")
  }
  #rkm
  res.rkm <- rkm.intern(xx, k=k, q=q, U=winner$U, thres = thres, maxiter = maxiter,printflag=FALSE)
  res.rkm$crit<- res.rkm$crit[!is.na(res.rkm$crit)]
  res.out <- list(res.rkm=res.rkm, res.init =res.init)
  if(printflag){
    cat(paste0("done! [",round(difftime(Sys.time(),horloge,units="secs"),2)," sec]\n"))
    cat(paste0("Total elpased time : ",round(difftime(Sys.time(),departtime,units="secs"),2)," sec\n"))
  }
  if(!is.null(axes)){
    proj <- mat.mult(xx,res.rkm$A[,axes])
    inertie<-round(100*colVars(proj)/sum(colVars(xx)),2)
    plot(proj,
         col=res.rkm$cluster,pch=16,
         xlab=paste("Dim",axes[1],"(",inertie[1],"%)"),
         ylab=paste("Dim",axes[2],"(",inertie[2],"%)")
         , main= "RKM map")
  }
  return(res.out)
}

initrkm <- function(xx,k,method){
  if(method=="kmpp"){
    res.init <- try(inaparc::kmpp(xx,k=k))
    cond1 <- cond2 <- !inherits(res.init, "try-error")
    if (cond1) {
      res.kmeans <- try(kmeans(xx, centers = res.init$v,iter.max = 1))
      cond2 <- !inherits(res.kmeans, "try-error")
    }
    if(cond2){
      res.out <- res.kmeans$cluster
    }
    
    if(!cond1|!cond2){
      res.out <- NULL
      return(res.out)
    }
  }else if(method=="maximin"){
    res.init <- try(inaparc::maximin(xx,k=k),silent = TRUE)
    cond1 <- cond2 <- !inherits(res.init, "try-error")
    if (cond1) {
      res.kmeans <- try(kmeans(xx, centers = res.init$v,iter.max = 1))
      cond2 <- !inherits(res.kmeans, "try-error")
    }
    if(!cond1|!cond2){
      res.out <- NULL
      return(res.out)
    }
  }else if (method == "sample") {
    #=forgy
    res.out <- sample(seq(k), size = dim(xx)[1], 
                      replace = TRUE)
  }else if (method == "MacQueen") {
    res.init <- sample(seq(dim(xx)[1]), size = k, 
                       replace = FALSE)
    res.kmeans <- try(kmeans(xx, centers = xx[res.init,],iter.max = 1))
    cond2 <- !inherits(res.kmeans, "try-error")
    if(cond2){
      res.out <- res.kmeans$cluster
    }else{
      res.out <- NULL
      return(res.out)
    }
  }else {
    stop("initialization method is unknown")
  }
  
  #verif nombre de classes
  if(length(unique(res.out))!=k){
    res.out <- NULL
    return(res.out)
  }
  
  return(res.out)
}

rkpod.intern <-function(xxna,
                        xximp,
                        k,
                        q,
                        nnodes = 1,
                        thres = 10^(-3),
                        maxiter = 100,
                        thres.rkm = 10^(-3),
                        maxiter.rkm = 20,
                        maxiter.init.rkm = 5,
                        method.init.rkm = rep("kmpp",20),
                        part.init=NULL,
                        printflag=TRUE){
  #initialisation
  ismissing<-which(is.na(as.matrix(xxna)))
  isnotmissing<-which(!is.na(as.matrix(xxna)))
  currentxx<- as.matrix(xximp)
  if(is.null(part.init)){
    part.init.intern <- lapply(X = method.init.rkm,FUN = initrkm,k=k,xx=xximp)
    names(part.init.intern) <- paste(seq.int(length(method.init.rkm)),method.init.rkm)
    garde <- which(!sapply(part.init.intern,is.null))
    part.init.intern <- part.init.intern[garde][!duplicated(part.init.intern[garde])]
  }else {
    part.init.intern <- part.init
  }
  continue <- TRUE
  nbiter <- 0
  critrkm <-rep(NA,maxiter)
  while (continue){
    nbiter <- nbiter+1
    if(printflag){cat(nbiter,"...")}
    #minimisation
    res.rkm <- rkm(xx=currentxx,
                   k=k,
                   q=q,
                   thres=thres.rkm,
                   maxiter.init=maxiter.init.rkm,
                   nnodes=nnodes,
                   part.init =part.init.intern,printflag = FALSE)
    #majoration
    
    UGAt <- transpose(res.rkm$res.rkm$recons.trans)
    currentxx[ismissing]<-UGAt[ismissing]
    critrkm[nbiter] <- sum(((xximp - UGAt)[isnotmissing])^2)
    
    #convergence
    cond1<- nbiter<maxiter
    cond2 <- ifelse(nbiter>1,yes = ((critrkm[nbiter-1]-critrkm[nbiter])/critrkm[nbiter-1])>=thres,no=TRUE)
    continue <- cond1 & cond2
    
    #init rkm
    part.init.intern <- lapply(method.init.rkm,initrkm,k=k,xx=currentxx)
    names(part.init.intern) <- paste(seq.int(length(method.init.rkm)),method.init.rkm)
    part.init.intern[[length(part.init.intern)+1]]<-res.rkm$res.rkm$cluster
    names(part.init.intern)[length(part.init.intern)] <- "Uprec"
    garde <- which(!sapply(part.init.intern,is.null))
    part.init.intern <- part.init.intern[garde][!duplicated(part.init.intern[garde])]
  }
  res.out <- list(U=res.rkm$res.rkm$U,G=res.rkm$res.rkm$G,A=res.rkm$res.rkm$A,
                  ximp=currentxx,
                  crit=critrkm[seq.int(nbiter)],
                  cluster=rowMaxs(res.rkm$res.rkm$U))
  return(res.out)
}

rkpod <- function(xxna,
                  method.init.rkpod = rep("MIPCA",10),
                  k,
                  q,
                  nnodes = 2,
                  thres = 10^(-4),
                  maxiter = 100,
                  thres.rkm = 10^(-3),
                  maxiter.rkm = 10,
                  maxiter.init.rkm = 5,
                  maxiter.init.rkpod = 5,
                  method.init.rkm = rep("kmpp",10),
                  printflag=TRUE,...){
  dots <- list(...)
  if("ncp" %in% names(dots)){
    ncp.mipca <- dots[["ncp"]]
  }else{
    ncp.mipca <- q
  }
  
  if(printflag){
    cat("Initialization...\n xmiss...")
  }
  departtime <- horloge <- Sys.time()
  res.imp <- vector(mode="list",length=length(method.init.rkpod ))
  names(res.imp) <- paste(seq.int(length(method.init.rkpod)),method.init.rkpod)
  for(method.tmp in unique (method.init.rkpod)){
    if(method.tmp%in%c("JM-DP","JM-GL","FCS-homo","FCS-hetero")){
      res.imp.tmp <- imputedata(data=xxna,
                                method=method.tmp,
                                nb.clust=k,
                                m = sum(method.init.rkpod==method.tmp),
                                Lstart = 100,
                                L = 20,
                                verbose = FALSE)$res.imp
    }else if (method.tmp=="MIPCA"){
      res.imp.tmp <-MIPCA(xxna,ncp = ncp.mipca,scale = FALSE,nboot = sum(method.init.rkpod==method.tmp),verbose = FALSE,method = "EM")$res.MI
      res.imp.tmp <- lapply(res.imp.tmp,as.matrix)
    }else if (method.tmp=="KNN"){
      res.imp.tmp <- replicate(n = sum(method.init.rkpod==method.tmp),expr = 
                                 VIM::kNN(xxna,
                                          k= max(1,round((dim(xxna)[1])/(10*k))),
                                          numFun = function(xx){sample(xx,size=1)},
                                          imp_var = FALSE),simplify = FALSE)
      res.imp.tmp <- lapply(res.imp.tmp,as.matrix)
    }
    
    res.imp[which(method.init.rkpod==method.tmp)] <- res.imp.tmp
    rm(list = "res.imp.tmp")
  }
  
  if (printflag){
    cat(paste0("done! [",round(difftime(Sys.time(),horloge,units="secs"),2)," sec]\n"))
    horloge <- Sys.time()
    cat(" pre-run for RKpod...")
  }
  if(nnodes> 1){
    nnodes.intern<-min(nnodes, length(res.imp),detectCores(logical = FALSE))
    cl <- makeCluster(nnodes.intern, type = "PSOCK")
    
    clusterExport(cl,
                  list("res.imp",
                       "xxna",
                       "k",
                       "q",
                       "thres",
                       "maxiter.init.rkpod","maxiter.init.rkm",
                       "thres.rkm",
                       "maxiter.rkm",
                       "method.init.rkm",
                       "rkm.intern",
                       "rkm",
                       "rkpod.intern",
                       "initrkm","mysvd","meanbygpf"),
                  envir = environment())
    
    clusterEvalQ(cl, library("Rfast"))
    clusterEvalQ(cl, library("inaparc"))
    
    res.init <- parLapply(cl,
                          X=res.imp,
                          fun = function(ximp,
                                         xxna,
                                         k,
                                         q,
                                         thres,
                                         maxiter,
                                         thres.rkm,
                                         maxiter.rkm,
                                         maxiter.init.rkm,
                                         method.init.rkm){
                            res.out <- rkpod.intern(xxna=xxna,
                                                    xximp=ximp,
                                                    k=k,
                                                    q=q,
                                                    nnodes = 1,
                                                    thres =  thres ,
                                                    maxiter = maxiter,
                                                    thres.rkm = thres.rkm,
                                                    maxiter.rkm = maxiter.rkm,
                                                    maxiter.init.rkm = maxiter.init.rkm,
                                                    method.init.rkm = method.init.rkm,
                                                    part.init=NULL,
                                                    printflag=TRUE)
                            return(res.out)
                          },
                          xxna=xxna,
                          k=k,
                          q=q,
                          thres =  thres,
                          maxiter = maxiter.init.rkpod,
                          thres.rkm = thres.rkm,
                          maxiter.rkm = maxiter.rkm,
                          maxiter.init.rkm = maxiter.init.rkm,
                          method.init.rkm = method.init.rkm)
    
    stopCluster(cl)
  }else if(nnodes==1){
    res.init <- lapply(res.imp,
                       FUN=function(ximp,
                                    xxna,
                                    k,
                                    q,
                                    nnodes,
                                    thres,
                                    maxiter,
                                    thres.rkm,
                                    maxiter.rkm,
                                    maxiter.init.rkm,
                                    method.init.rkm){
                         rkpod.intern(xxna=xxna,
                                      xximp=ximp,
                                      k=k,
                                      q=q,
                                      nnodes = nnodes,
                                      thres =  thres ,
                                      maxiter = maxiter,
                                      thres.rkm = thres.rkm,
                                      maxiter.rkm = maxiter.rkm,
                                      maxiter.init.rkm = maxiter.init.rkm,
                                      method.init.rkm = method.init.rkm,
                                      part.init=NULL,
                                      printflag=FALSE)
                       },xxna=xxna,
                       k=k,
                       q=q,
                       nnodes = nnodes,
                       thres =  thres,
                       maxiter = maxiter.init.rkpod,
                       thres.rkm = thres.rkm,
                       maxiter.rkm = maxiter.rkm,
                       maxiter.init.rkm = maxiter.init.rkm,
                       method.init.rkm = method.init.rkm)
  }
  names(res.init) <- names(res.imp)
  
  #winner
  
  mincrit<-sapply(res.init,FUN = function(res){
    yy <- na.omit(res$crit)
    res.out <- yy[length(yy)]
    return(res.out)
  }, simplify =TRUE)
  
  winner <- which.min(mincrit)
  # print(winner)
  if(printflag){
    cat(paste0("done! [",round(difftime(Sys.time(),horloge,units="secs"),2)," sec]"))
    cat(paste0("\n Elapsed time for initialization [",round(difftime(Sys.time(),departtime,units="secs"),2)," sec]"))
    horloge <- Sys.time()
    cat("\nRKpod is runing using the best starting values...")
  }
  #rkpod
  res.kpod <-  rkpod.intern(xxna=xxna,
                            xximp=res.init[[winner]]$ximp,
                            k=k,
                            q=q,
                            nnodes = nnodes,
                            thres =  thres ,
                            maxiter = maxiter,
                            thres.rkm = thres.rkm,
                            maxiter.rkm = maxiter.rkm,
                            maxiter.init.rkm = maxiter.init.rkm,
                            part.init=list(as.factor(res.init[[winner]]$cluster)),
                            printflag=TRUE)
  if(printflag){
    cat(paste0("done! [",round(difftime(Sys.time(),horloge,units="secs"),2)," sec]\n"))
    cat(paste0("Total elpased time : ",round(difftime(Sys.time(),departtime,units="secs"),2)," sec\n"))
  }
  
  res.out <- list(res.kpod = res.kpod, res.init=res.init)
  return(res.out)
}

plot.rkpod <- function(res.rkpod, axes = c(1,2)){
  G <- res.rkpod$res.kpod$G
  part.rkpod <- res.rkpod$res.kpod$cluster
  proj <- mat.mult(res.rkpod$res.kpod$ximp,res.rkpod$res.kpod$A[,axes])
  inertie<-round(100*colVars(proj)/sum(colVars(res.rkpod$res.kpod$ximp)),2)
  plot(proj,
       col=gray.colors(max(part.rkpod))[part.rkpod],
       pch=16,
       xlab=paste0("Dim",axes[1]," (",inertie[1],"%)"),
       ylab=paste0("Dim",axes[2]," (",inertie[2],"%)")
       , main= "RKM map")
  points(res.rkpod$res.kpod$G[,axes],col=gray.colors(nrow(res.rkpod$res.kpod$G)),pch=3,cex=3,lwd=2)
}