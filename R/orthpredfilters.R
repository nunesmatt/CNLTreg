orthpredfilters<-function(filter = c(0.5, 1, 0.5)){

n<-length(filter)

ri<-(n+1)/2	# only works for odd length filter

A<-matrix (NA,n,n)

if(n==2){	# nbr , remove
	A<-matrix(c(-1,1,1,1),2,2)
}

if (n==3){
	A[,1]<-c(1,-1,1)
	A[,2]<-filter
	A[,3]<-c(-(1+filter[3])/(filter[1]+1),-(filter[3]-filter[1])/(filter[1]+1),1)
}

if(n==2){
	QA<-A
}
else{
	#  orthogonal and equal norms (filter 1) 
	nn<-colSums(A^2)
	n1<-nn[2]		# should be the same as sum(filter^2)		
	QA<-A
	QA[,3:ncol(QA)]<-sqrt(n1)*QA[,3:ncol(QA)]/rep(sqrt(nn[3:ncol(QA)]),each=nrow(QA))		

# remove first column and minus
QA<-QA[,-1]
QA<- -QA
QA[ri,]<--QA[ri,]
}

return(QA)

}
