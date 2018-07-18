\name{mthreshC}
\alias{mthreshC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform 'multiwavelet style' level-dependent soft
thresholding for complex-valued wavelet coefficients
}
\description{
This function uses chi^2 statistics similar to Barber and Nason
     (2004) to threshold wavelet coefficients based on their magnitude
}
\usage{
mthreshC(coeffv, Sigma, rl, po, ali, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coeffv}{
A matrix of complex-valued wavelet coefficients (columns are
          real and imaginary parts of the coefficients respectively).
}
  \item{Sigma}{
An array of dimension \code{2 x 2 x n} describing the covariance
          between real and imaginary parts of the wavelet coefficients.
          In particular, \code{Sigma[,,i]} represents the covariance between
          real and imaginary parts of the \code{i}th lifted wavelet
          coefficient (see \env{rl} argument).
}
  \item{rl}{
The removelist (trajectory of lifted points) corresponding to
          a forward lifting transform.
}
  \item{po}{
A vector of indices describing the unlifted scaling
          coefficients in a forward lifting transform.
}
  \item{ali}{
A list of indices of observations, each entry corresponding
          to an 'artificial level' (finest to coarsest), see \code{\link[adlift]{artlev}}
          for more details.
}
  \item{verbose}{
 Indicates whether helpful messages should be printed to the
          console during the procedure.
}
}
\details{
The procedure in Downie and Silverman (1998) or Barber and Nason (2004) makes
     use of the magnitude of wavelet coefficients to threshold them.
     In particular, the covariance between the components of the
     wavelet coefficients (contained in \env{Sigma} is taken into account
     to compute a thresholding statistic, the distribution of which is
     chi-squared_2 distributed, see \code{\link[wavethresh]{cthresh}} for more details.  These
     statistics are then compared with level-dependent universal
     thresholds computed by counting the numnber of coefficients in
     specific artificial levels.
}
\value{
A list with the following components:
\item{chi}{the vector of chi-squared statistics used in the thresholding procedure.}
\item{coeffvt}{the matrix of thresholded coefficients, columns representing
          the real and imaginary components respectively.}
}
\references{
 Hamilton, J., Knight, M. I., Nunes, M. A. and Fryzlewicz (2018)
     Complex-valued wavelet lifting and applications. \emph{Technometrics}, \bold{60} (1), 48-60,DOI 10.1080/00401706.2017.1281846.
     Barber, S. and Nason, G. P. (2004) Real nonparametric regression
     using complex wavelets. \emph{J. Roy. Stat. Soc. B} \bold{66} (4), 927--939.
     Downie, T. R. and Silverman, B. W. (1998) The discrete multiple
     wavelet tranform and thresholding methods. \emph{IEEE Trans. Sig.
     Proc.} \bold{46} 2558--2561.

}
\author{
 Matt Nunes, Marina Knight
}

\seealso{
\code{\link[wavethresh]{cthresh}},
\code{\link{denoisepermC}},
\code{\link{denoisepermCh}}
}
\examples{

library(adlift)

set.seed(100)

# construct an (irregular) sampling structure:

x<-sort(runif(200))

g<-make.signal2("bumps",x=x)

# generate IID noise with a particular sd
noise<-rnorm(200,0,sd=0.5)

f<-g+noise

# perform forward complex lifting transform

out<-fwtnppermC(x,f,LocalPred=LinearPred,neigh=1)

# have a look at some of the coefficients

out$coeffv[1:10]

# extract lifting matrix and induced lifting variances
W <- out$W

Gpre<-tcrossprod(W,Conj(W))

indsd<-sqrt(diag(Gpre))

# now estimate noise sd using the first artificial level:

al<-artlev(out$lengthsremove,out$removelist)

fine<-(out$coeffv/indsd)[al[[1]]]

varest<-mad(Re(fine))^2

# now compute coefficient covariance structure, see
# Hamilton et al. (2018), Appendix B

C = varest * tcrossprod(W)
    G = varest * Gpre
    P = Conj(G) - t(Conj(C)) %*% solve(G) %*% C
    Sigma <- array(0, dim = c(2, 2, length(out$coeffv)))
    Sigma[1, 1, ] <- diag(Re(G + C)/2)
    Sigma[2, 2, ] <- diag(Re(G - C)/2)
    Sigma[1, 2, ] <- -diag(Im(G - C)/2)
    Sigma[2, 1, ] <- diag(Im(G + C)/2)

# now threshold complex coefficients according to this structure:

coeff.thresh<-mthreshC(out$coeffv,Sigma,out$removelist,out$pointsin,al)

# have a look at some of these coefficients

coeff.thresh$coeffv[1:10]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
