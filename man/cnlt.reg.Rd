\name{cnlt.reg}
\alias{cnlt.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs 'nondecimated' complex-valued wavelet lifting for signal denoising
}
\description{
The transform-threshold-invert procedure for signal denoising is dependent on the trajectory (lifting order) used in the forward lifting transform.  This procedure uses trajectory bootstrapping and averaging of estimates to gain denoising performance
}
\usage{
cnlt.reg(x, f, P, returnall = FALSE, nkeep = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of any length (not necessarily equally spaced) that gives the grid on which the signal is
observed.}
  \item{f}{Vector of the same length as \code{x} that gives the signal values corresponding to the
\code{x}-locations. }
  \item{P}{
Number of trajectories to be used by the nondecimated lifting algorithm.
}
\item{returnall}{Indicates whether the function returns useful variables or just the denoised datapoints.}
  \item{nkeep}{
Number of scaling points we want at the end of the transform. The usual choice is \code{nkeep}=2.
}
  \item{\dots}{Any other arguments to be passed to \code{\link{denoisepermC}}, see the function documentation for more details.}
}
\details{
	Essentially, this function applies the complex wavelet lifting
     denoising procedure \code{\link{denoisepermC}} \env{P} times, each with a different random lifting trajectory. 
     This results in \env{P} estimates of the (unknown) true
     signal. The average of these estimators is the proposed estimator.
}
\value{
 If \env{returnall=FALSE}, the estimate of the function after denoising.  If \env{returnall=TRUE},
a list with components:\cr

  \item{vec}{A matrix of dimension \code{P x (n - nkeep)}, each row corresponding to a different lifting trajectory.}
  \item{aveghat}{Estimated signal after removing the noise.}
}
\references{
Hamilton, J., Nunes, M. A., Knight, M. I. and Fryzlewicz, P. (2018) Complex-valued wavelet lifting and applicati$
\emph{Technometrics}, bold{60} (1), 48-60, DOI 10.1080/00401706.2017.1281846.\cr\cr
For the real-valued equivalent procedure, see also\cr
Knight, M. I. and Nason, G. P. (2009) A 'nondecimated' wavelet transform. \emph{Stat. Comput.}
\bold{19} (1), 1--16.
}
\author{
Matt Nunes
}
\section{Warning}{Using a large number of trajectories for long datasets could take a long time!}
\seealso{
\code{\link{denoisepermCh}},
\code{\link{fwtnppermC}},
\code{\link{mthreshC}},
\code{\link[nlt]{nlt}}
}
\examples{ 
    library(adlift)

    # construct an (irregular) observation grid
     x<-runif(256) 
     
     #construct the true, normally unknown, signal
     g<-make.signal2("blocks",x=x) 
     
     #generate noise with mean 0 and signal-to-noise ratio 5
     noise<-rnorm(256,mean=0,sd=sqrt(var(g))/5)
     
     #generate a noisy version of g 
     f<-g+noise 
     
     #decide on a number of random trajectories to be used (e.g. J=5 below), and apply
     # the nondecimated lifting transform to the noisy signal (x,f): 
\dontrun{
     est<-cnlt.reg(x,f,P=50,LocalPred=AdaptPred,neighbours=1,returnall=FALSE) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{regression}

