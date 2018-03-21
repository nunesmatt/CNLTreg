mthreshC <-function (coeffv, Sigma, rl ,po, ali, verbose = FALSE) 
{

n<-length(coeffv)

nsteps<-length(rl)	# number of lifting steps
nkeep <-length(po)	# number not lifted

nfilt<-2

coeffvt<-cbind(Re(coeffv),Im(coeffv))

chisq<-rep(NA,nsteps)

onenbrcheck<-function(m){abs(det(m)) < 1e-7}

nlev<-length(ali)
for (i in 1:nlev) {
	rows<-ali[[i]]
	nj<-length(rows)		
	chithresh <- 2 * log(nj) + (nfilt -2) * log(log(nj))
	cat("chi^2 threshold is:",chithresh,"\n")
	for(j in 1:nj){
		rem<-rows[j]
		dr<-matrix(coeffvt[rem,],ncol=1)

		# now check for singularity:
		m <- Sigma[,,rem]
		onenbr<-onenbrcheck(m)
		if(onenbr){	
			# "ignore correlation
			m[2,1]<-m[1,2]<-0
		}
		sigjinv<-solve(m)
		chisq[i]<-t(dr) %*% sigjinv %*% dr	# should be a 1x1 value to test
		if(verbose){
			cat("chisq stat is:",chisq[i],"\n")
		}
		# now threshold, a la Barber and Nason 04 P.9
		# shrink according to the modulus of the vector 
            	if (chisq[i] > 0){ 
                	shrink <- (max(chisq[i] - chithresh, 0))/chisq[i]
		}
            	else{
			shrink <- 0
		}
            	coeffvt[rem,] <- t(dr) * shrink
   	}
}

l<-list()
l$chi<-chisq
l$coeffvt<-complex(real=coeffvt[,1],imaginary=coeffvt[,2])

return(l)

}
