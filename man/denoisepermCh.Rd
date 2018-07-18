\name{denoisepermCh}
\alias{denoisepermCh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Denoises a signal using the complex-valued lifting transform and multivariate soft thresholding and heteroscedastic variance computation
}
\description{
Denoises an input signal contaminated by noise. First the signal is decomposed using the complex-valued lifting
scheme (see \code{\link{fwtnppermC}}) using an order of point removal. A sliding window approach is then used
on these wavelet coefficients to estimate a local noise variance.
The resulting complex-valued wavelet
coefficients are then thresholded using a soft thresholding rule on the details' magnitude.
The transform is inverted and an estimate of the noisy signal is obtained.
}
\usage{
denoisepermCh(x, f, returnall = FALSE, verbose = FALSE, ...)
%denoisepermCh(x, f, pred = LinearPred, neigh = 1, int = TRUE, clo = FALSE, keep = 2, 
%per = sample(1:length(x), (length(x) - keep), FALSE), returnall = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of any length (not necessarily equally spaced) that gives the grid on which the signal is
observed.}
  \item{f}{Vector of the same length as \code{x} that gives the signal values corresponding to the
\code{x}-locations. }
%  \item{pred}{The type of regression to be used in the prediction step of the modified lifting algorithm.
%See \code{\link[adlift]{fwtnp}} for more details on the choices.}
%  \item{neigh}{Number of neighbours to be used in order to construct the neighbourhood of each point that
%has to be removed. If '\code{clo}=FALSE', this gives the number of neighbours on each side of the removed point.}
%  \item{int}{Specifies whether (\code{int}=TRUE) or not (\code{int}=FALSE) an intercept is to be used in the
%regression curve. For \code{pred}=\code{AdaptPred} or \code{AdaptNeigh}, the algorithm automatically makes
%this choice.}
%  \item{clo}{If (\code{clo}=TRUE) or (\code{clo}=FALSE), then at each step the neighbours are in closest,
%respectively symmetrical configuration.}
%  \item{keep}{Number of scaling points we want at the end of the transform. The usual choice is \code{keep}=2.}
%  \item{per}{Vector of length (length(\code{x})-\code{keep}) which gives the order of point removal in the
%lifting algorithm.}
  \item{returnall}{Indicates whether the function returns useful variables or just the denoised datapoints.}
  \item{verbose}{Indicates whether useful messages should be printed to the console during the procedure.}
\item{...}{Any other arguments to be passed to \code{\link[nlt]{fwtnpperm}} and \code{\link{fwtnppermC}}.}
}
\details{
After the complex lifting transform is applied, the wavelet coeficients are divided into artificial levels.
The details from the lifting scheme have different variances, and will therefore be normalized
to have the same variance as the noise, by using the lifting matrix.  A sliding window is used to compute a local 'heteroscedastic' noise variance by taking the MAD of those normalized details falling into the window, see Nunes et al. (2006) for more details.  
Given the noise estimates for each observation,  the normalized details can then be thresholded and un-normalized.  The transform is then inverted to
give an estimate of the signal.
}
\value{
 If \env{returnall=FALSE}, the estimate of the function after denoising.  If \env{returnall=TRUE},
a list with components:\cr

  \item{fhat}{Estimated signal after removing the noise.}
  \item{w}{This is the matrix associated to the modified lifting transform.}
  \item{al}{List giving the split of points between the artificial levels.}
  \item{sd}{Estimated heteroscedastic standard deviation of the noise.}
}
\references{
Hamilton, J., Nunes, M. A, Knight, M. I. and Fryzlewicz, P. (2018) Complex-valued wavelet lifting and applications. 
\emph{Technometrics}, \bold{60} (1), 48-60, DOI 10.1080/00401706.2017.1281846.\cr\cr
Nunes, M. A., Knight, M. I and Nason, G. P. (2006) Adaptive lifting for nonparametric regression.
\emph{Stat. Comput.} \bold{16} (2), 143--159.\cr\cr
Knight, M. I. and Nason, G. P. (2009) A 'nondecimated' wavelet transform. \emph{Stat. Comput.}
\bold{19} (1), 1--16.
}
\author{
Matt Nunes, Marina Knight
}

\seealso{
\code{\link{denoisepermC}},
\code{\link{fwtnppermC}},
\code{\link[nlt]{fwtnpperm}},
\code{\link[adlift]{heterovar}},
\code{\link{mthreshC}}
}
\examples{

library(MASS)	# where the motorcyle data lives

mcycleu<-mcycle[which(duplicated(mcycle$times)=='FALSE'),]
time<-mcycleu[,1]
accel<-mcycleu[,2]

set.seed(200)
est1<-denoisepermCh(time,accel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{regression}
