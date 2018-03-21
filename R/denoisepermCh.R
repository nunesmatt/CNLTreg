denoisepermCh <-function (x, f, returnall = FALSE, verbose = FALSE, ...) 
{

sdvec<-rep(1,length(x))
    newcoeff <- NULL
    ndetlist <- list()
    tclist <- NULL
    out <- fwtnppermC(x, f, ...)

    W<-out$W			
    nfilt<- 2   		# symmetric neighbours for the moment

    lr <- out$lengthsremove
    rem <- out$removelist
    al <- artlev(lr, rem)
    levno <- length(al)

    Gpre = W %*% t(Conj(W))	

# do the real transform heterosced. variance estimation.  This will overrule any computation above

n<-length(x)

out1 <- fwtnpperm(x, f, do.W = FALSE, varonly = FALSE,mod=rem, ...)

    po1 <- out1$pointsin
    nonorcoeff1 <- out1$coeff
    lr1 <- out1$lengthsremove	# although these should be the same as the complex one above
    al1 <- artlev(lr1, rem)	#
    levno1 <- length(al1)
    detail <- y <- matrix(0, 1, n - length(po1))
    y <- x[setdiff(1:n, po1)]			# sort rem for speed?
    detail <- nonorcoeff1[setdiff(1:n, po1)]
    h <- heterovar(y, detail, al1)
    sdvech1 <- h$varvec1
    sdh1 <- NULL
    sdh1[setdiff(1:n, po1)] <- sdvech1
    sdh1[po1] <- NA

sdsq<-sdh1^2

# now proceed as before...


#    C = sdsq * W%*%t(W)		# should be no conjugation involved here
    C = sdsq * tcrossprod(W)		# should be no conjugation involved here
    G = sdsq * Gpre
    P = Conj(G) - t(Conj(C))%*%solve(G)%*%C

    Sigma<-array(0,dim=c(2,2,length(nonorcoeff1)))	# Vxx Vyy Vxy Vyx 
    Sigma[1,1,] <- diag(Re(G+C)/2)	# gets jk, jk value
    Sigma[2,2,] <- diag(Re(G-C)/2)
    Sigma[1,2,] <- -diag(Im(G-C)/2)	
    Sigma[2,1,] <- diag(Im(G+C)/2)      

	ll<-mthreshC(out$coeffv, Sigma, rem, out$pointsin, ali=al, verbose = verbose)
	newcoeff<-ll$coeffvt

fhat<-NULL
out$coeffv<-newcoeff
fhat<-solve(Re(W))%*%Re(newcoeff)

if (returnall) {
        return(list(fhat = fhat, w = out$W, al = al, sd = sqrt(sdsq)))
}
else {
        return(fhat)
}

}
