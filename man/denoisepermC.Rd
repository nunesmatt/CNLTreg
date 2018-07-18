\name{denoisepermC}
\alias{denoisepermC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Denoises a signal using the complex-valued lifting transform and multivariate soft thresholding
}
\description{
Denoises an input signal contaminated by noise. First the signal is decomposed using the complex-valued lifting 
scheme (see \code{\link{fwtnppermC}}) using an order of point removal. The resulting complex-valued wavelet
coefficients are then thresholded using a soft thresholding rule on the details' magnitude.  
The transform is inverted and an estimate of the noisy signal is obtained.
}
\usage{
denoisepermC(x, f, returnall = FALSE, sdtype = "adlift", verbose = FALSE, ...)
%denoisepermC(x, f, pred = LinearPred, neigh = 1, int = TRUE, clo = FALSE, keep = 2, 
%per = sample(1:length(x), (length(x) - keep), FALSE), returnall = FALSE, 
%sdtype = "adlift", verbose = FALSE)
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
  \item{sdtype}{Options are either \code{"adlift"} or \code{"complex"}, indicating whether the noise variance is estimated with the average of the mean absolute deviations of both real and imaginary components of the finest wavelet coefficients, or just the real component, as in \code{\link[adlift]{denoise}}.}
  \item{verbose}{Indicates whether useful messages should be printed to the console during the procedure.}
\item{...}{Any other arguments to be passed to \code{\link{fwtnppermC}}, see documentation for this function for more details.}
}
\details{
After the complex lifting transform is applied, the wavelet coeficients are divided into artificial levels. 
The details from the lifting scheme have different variances, and will therefore be normalized 
to have the same variance as the noise, by using the lifting matrix. Those normalized details falling into the finest artificial level will 
be used for estimating the standard deviation of the noise that contaminated the signal.  The variable \env{sdtype} is used for this estimate, see Appendix B of Hamilton et al. (2018) for more details.  Using this 
estimate,  the normalized details can then be thresholded and un-normalized.  The transform is then inverted to give an estimate of the signal. 
}
\value{
 If \env{returnall=FALSE}, the estimate of the function after denoising.  If \env{returnall=TRUE}, 
a list with components:\cr

  \item{fhat}{Estimated signal after removing the noise.}
  \item{w}{This is the matrix associated to the modified lifting transform.}
  \item{indsd}{Vector giving the standard deviations of the detail and scaling coefficients.}
  \item{al}{List giving the split of points between the artificial levels.}
  \item{sd}{Estimated standard deviation of the noise.}
}
\references{
Hamilton, J., Nunes, M. A., Knight, M. I. and Fryzlewicz, P. (2018) Complex-valued wavelet lifting and applications.
\emph{Technometrics}, \bold{60} (1), 48-60, DOI 10.1080/00401706.2017.1281846.
}
\author{
Matt Nunes, Marina Knight
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{denoisepermCh}},
\code{\link{fwtnppermC}},
\code{\link{mthreshC}},
\code{\link{cnlt.reg}}
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

     # perform the complex-valued denoising procedure to the noisy signal (x,f):
     est<-denoisepermC(x,f,LocalPred=AdaptPred,neigh=1,returnall=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{regression}
