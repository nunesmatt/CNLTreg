\name{CNLTreg-package}
\alias{CNLTreg-package}
\alias{CNLTreg}
\docType{package}
\title{
\packageTitle{CNLTreg}
}
\description{
\packageDescription{CNLTreg}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{CNLTreg}
\packageIndices{CNLTreg}

The main routines of the package are \code{\link{denoisepermC}} and \code{\link{cnlt.reg}} which perform complex-valued lifting-based denoising, using a single or a multiple (chosen) number of lifting trajectories, respectively.  
}
\author{
\packageAuthor{CNLTreg}

Maintainer: \packageMaintainer{CNLTreg}
}
\references{
Hamilton, J., Nunes, M. A, Knight, M. I. and Fryzlewicz, P. (2018) Complex-valued wavelet lifting and applications. 
\emph{Technometrics}, \bold{60} (1), 48-60, DOI 10.1080/00401706.2017.1281846.\cr\cr
For related literature on the lifting methodology adopted in the technique, see\cr\cr
Nunes, M. A., Knight, M. I and Nason, G. P. (2006) Adaptive lifting for nonparametric regression. 
\emph{Stat. Comput.} \bold{16} (2), 143--159.\cr\cr
Knight, M. I. and Nason, G. P. (2009) A 'nondecimated' wavelet transform. \emph{Stat. Comput.} 
\bold{19} (1), 1--16. 
}

\keyword{ package }
\seealso{
	\code{\link[adlift]{denoise}}
	\code{\link[nlt]{denoiseperm}}
	\code{\link[nlt]{nlt}}
}
