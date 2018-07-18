fwtnppermC <-
function (x, f, LocalPred = LinearPred, neighbours = 1, intercept = TRUE, 
    closest = FALSE, nkeep = 2, mod = sample(1:length(x), (length(x) - nkeep), FALSE)) 
{
    X <- x
    I <- intervals(X, "reflect")
    lengths <- lengthintervals(X, I, type = "midpoints", neighbours, closest)
    X <- as.row(X)
    f <- as.row(f)
    nkeep <- max(nkeep, 1)
    n <- length(X)
    removelist <- NULL
    lengthsremove <- NULL
    neighbrs <- list()
    gamlist <- list()
    alphalist <- list()
    schemehist <- NULL
    interhist <- NULL
    clolist <- NULL
    pointsin <- matrix(1:n, 1, n)
    pointsin <- pointsin[order(X)]

    Ialpha <- NULL

    # alterations for multifilter lifting.  This assumes nfilt = 2* neigbours
    # for a (2*neighbours+1)-tap filter length.  This is (for the moment) for 
    # closest = FALSE.  The routine should also be for a fixed trajectory so 
    # that the order of removal is the same for all filters.

    # EDIT: for the Hamilton et al. (2018), nfilt = 2, so changed to be fixed below.

#    nfilt <- 2*neighbours     
    nfilt <- 2
    ri<- neighbours + 1		# index of remove in filter

    coeff <- matrix(f,nrow = n, ncol = nfilt)	# should be f rep'd by col

    matno <- n - nkeep

    W <- diag(n)
    W <-array(rep(W,times=nfilt),c(n,n,nfilt))

    filtlist<-vector("list",matno)

    for (j in 1:matno) {
        remove <- mod[j]
        removelist[j] <- remove
        out <- getnbrs(X, remove, pointsin, neighbours, closest)
        nbrs <- out$n
        index <- out$index

	ds <- abs(X[nbrs] - X[remove])
	irregdeg <- max(ds)/min(ds)

        res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, 
            neighbours)
        if (length(res) == 2) {
            l <- res[[1]]
            clolist[j] <- res[[2]][[1]]
            nbrs <- res[[2]][[2]]
            index <- res[[2]][[3]]
        }
        else {
            l <- res
        }
        neighbrs[[j]] <- nbrs
        weights <- l[[1]]
        pred <- l[[2]]
        if (length(l) == 3) {
            scheme <- NULL
            int <- NULL
            details <- NULL
        }
        else {
            scheme <- l[[5]]
            int <- l[[4]]
            details <- l[[6]]
        }

	##  here come the modifications:
	##  note that the original weight vector may change a little
	
	# temporary fix to avoid edge prediction problems:
	weightsa<-append(weights,1,after=neighbours)	# adds unit weight into middle

	# just in case, do this the append way:
	nbrsa<-append(nbrs,remove,after=neighbours)	# adds remove into middle

	# update all coeff vectors by first updated c [replace] (see below) 
	# and apply filters to all coeff matrix
	
	m<-matrix(coeff[nbrs,],nrow=length(nbrs))
	tmpca<-insertRow(m,ri, v = coeff[remove,])

	filtmat<-orthpredfilters(weightsa)

	filtlist[[j]]<-filtmat

	# does all filters at once
       	coeff[remove,] <- colSums(filtmat * tmpca)

	# update accordingly: see above
	# new bit should be fine now too: use first filter (already removed invalid one above)

	if(length(nbrs)>1){
		w<- -filtmat[-ri,1]
	}
	else{
		w<-1	# single weight
	}
       	l1 <- PointsUpdate(X, coeff[,1], nbrs, index, remove, pointsin, w, lengths)
       	coeff[nbrs,] <- l1$coeff[nbrs]

        lengths <- l1$lengths
        r <- l1$r
        weights <- l1$weights
        N <- l1$N
        alpha <- l1$alpha

	# now update the lifting matrix  

	Wtmp<-array(sapply(1:nfilt, function(i) W[nbrsa,,i]*filtmat[,i]),dim=c(length(nbrsa),n,nfilt))
        W[remove,, ] <- apply(Wtmp, c(2,3), sum )	# should be n x nfilt
	
        W[nbrs,,]  <- W[nbrs,,1] + matrix(alpha) %*% W[remove,,1]

        lengthsremove[j] <- lengths[r]

	Ialpha[j]<-irregdeg

        gamlist[[j]] <- weights
        alphalist[[j]] <- alpha
        schemehist[j] <- scheme
        interhist[j] <- int
        lengths <- lengths[setdiff(1:length(pointsin), r)]
        pointsin <- setdiff(pointsin, remove)
    }

   coeffv<-NULL
   coeffv<-complex(real=coeff[,1],imaginary=coeff[,2])

   # set the scaling coeffs to be the ones from only the first transform (real)
   # and do the same for the lifting matrix

   W<-matrix(complex(real=W[,,1],imaginary=W[,,2]),nrow=n,ncol=n)

    N <- length(pointsin)
    return(list(coeff = coeff, lengthsremove = lengthsremove, pointsin = pointsin, removelist = removelist, 
        gamlist = filtlist, alphalist = alphalist, W = W, reo = match(1:n,c(pointsin,removelist)), coeffv=coeffv, Ialpha = Ialpha))
}
