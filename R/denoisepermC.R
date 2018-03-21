denoisepermC <-function (x, f, returnall = FALSE, sdtype="adlift", verbose = FALSE, ...) 
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

#    Gpre = W %*% t(Conj(W))	
    Gpre = tcrossprod(W,Conj(W))	


if(sdtype!="adlift"){
    	indsd<-sqrt(diag(Gpre))	
}
else{
	indsd<-sqrt(diag(tcrossprod(Re(W))))
}

norcoeff <- out$coeffv/indsd
 
ndetlist <- norcoeff[al[[1]]]

if(sdtype!="adlift"){
  	sdsq <- mean(c(mad(Re(ndetlist)), mad(Im(ndetlist))))^2
}
else{
    	sdsq <- mad(Re(ndetlist))^2 
}


#    C = sdsq * W%*%t(W)		# should be no conjugation involved here

    C = sdsq * tcrossprod(W)		# should be no conjugation involved here

    G = sdsq * Gpre
    P = Conj(G) - t(Conj(C))%*%solve(G)%*%C

    Sigma<-array(0,dim=c(2,2,length(norcoeff)))	# Vxx Vyy Vxy Vyx 
    Sigma[1,1,] <- diag(Re(G+C)/2)	# gets jk, jk value
    Sigma[2,2,] <- diag(Re(G-C)/2)
    Sigma[1,2,] <- -diag(Im(G-C)/2)	
    Sigma[2,1,] <- diag(Im(G+C)/2)      

    ll<-mthreshC(out$coeffv,Sigma,rem,out$pointsin,ali=al, verbose = verbose)
    newcoeff<-ll$coeffvt

fhat<-NULL
out$coeffv<-newcoeff
fhat<-solve(Re(W))%*%Re(newcoeff)

if (returnall) {
        return(list(fhat = fhat, w = out$W, indsd = indsd, al = al, 
            sd = sd))
}
else {
        return(fhat)
}

}
