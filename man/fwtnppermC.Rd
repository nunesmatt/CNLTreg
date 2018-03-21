\name{fwtnppermC}
\alias{fwtnppermC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward complex wavelet lifting transform
}
\description{
Performs the complex-valued lifting transform on a signal with grid \env{x} and corresponding function 
values \env{f}.
}
\usage{
fwtnppermC(x, f, LocalPred = LinearPred, neighbours = 1, 
intercept = TRUE, closest = FALSE, nkeep = 2, 
mod = sample(1:length(x), (length(x) - nkeep), FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of grid values. Can be of any length, not necessarily equally spaced.}
  \item{f}{A vector of function values corresponding to \env{x}. Must be of the same length as \env{x}.}
  \item{LocalPred}{The type of regression to be performed in the prediction lifting step. Possible options are \code{LinearPred}, \code{QuadPred}, \code{CubicPred}, \code{AdaptPred} and
\code{AdaptNeigh}.}
  \item{neighbours}{The number of neighbours over which the regression is performed at each step.
If \env{closest} is \code{FALSE}, then this in fact denotes the number of neighbours on each side of the removed point.}
  \item{intercept}{Indicates whether or not the regression prediction includes an intercept.}
  \item{closest}{Refers to the configuration of the chosen neighbours. If \env{closest} is \code{FALSE},
the neighbours will be chosen symmetrically around the removed point.
Otherwise, the closest neighbours (in distance) will be chosen.}
  \item{nkeep}{The number of scaling coefficients to be kept in the final representation of the initial signal.
 This must be at least two. }
  \item{mod}{Vector of length (length(\code{x})-\code{nkeep}). This gives
  the trajectory for the lifting algorithm to follow, i.e. it gives
  the order of point removal.}
}
\details{
Given \eqn{n} points on a line, \env{x}, each with a corresponding env{f} value this function computes
the complex-valued lifting transform of the (\env{x},\env{f}) data. This is similar in spirit to the 
real-valued lifting transform in \code{\link[nlt]{fwtnpperm}}, except that the algorithm constructs *two* 
orthogonally linked prediction filters, as in Section 2.2 of Hamilton et al. (2017).  A summary of the 
procedure is as follows:

Step One. Compute "integrals" associated to each point, representing the intervals that each gridpoint 
\code{x_i} spans. 

Then for each point index in the lifting trajectory \env{mod}, 

Step Two(a). The neighbours of the removed point are identified using the specified neighbour configuration.
The value of \env{f} at the removed point is predicted using the specified regression over the
neighbours, unless an adaptive procedure is chosen. In this case, the algorithm chooses the regression which
produces the minimal detail coefficient (in magnitude) from a range of regression types (see \code{\link[adlift]{AdaptPred}} or \code{\link[adlift]{AdaptNeigh}} for more information).  In either case, the regression
specifies a local filter of the function values over the neighbourhood, \code{L}.

Step Two(b). A second filter, \code{M}, is then constructed orthogonal to \code{L}, such that it has unit norm, see Hamilton et al. (2017) for more details. 

The differences between the removed point's \env{f} value and the predictions using the two filters are computed, which constitute the real and imaginary parts of the complex-valued wavelet coefficient.  This coefficient is then stored

Step Three. The integrals and the scaling function values (neighbouring \env{coeffv} values) are updated according to the filter \code{L}. 

The algorithm continues until all points in \env{mod} are removed.  

}
\value{
\item{coeff}{matrix of detail and scaling coefficients in the wavelet decomposition of the signal; first column: real component, second column: imaginary component.}
\item{lengthsremove}{vector of interval lengths corresponding to the points removed during the transform
(in \env{removelist}).}
\item{pointsin}{indices into \env{X} of the scaling coefficients in the wavelet decomposition. These are the
indices of the \env{X} values which remain after all points in \env{removelist} have been predicted and
removed. This has length \env{nkeep}.}
\item{removelist}{a vector of indices into \env{X} of the lifted coefficients during the transform
(in the order of removal).}
\item{gamlist}{a list of all the prediction weights used at each step of the transform; each list entry is a matrix of two rows, corresponding to the filters \code{L} and \code{M}.}
\item{alphalist}{a list of the update coefficients used in the update step of the decomposition.}
\item{W}{The complex-valued lifting matrix associated to the transform.}
\item{reo}{An index into the observations indicating a reordering to give \code{1:n}.  This is reported for convenience for other functions, and is not intended for use by the user.}
\item{coeffv}{vector of complex-valued detail and scaling coefficients in the wavelet decomposition of the signal; contains the same information as \env{coeff}.}
\item{Ialpha}{Vector of "irregularity degree" measures corresponding to each lifting step of the transform.  Note that this is returned for convenience in other functions, and is not intended for use by the user.}
}
\references{
Hamilton, J., Knight, Nunes, M. A. and Fryzlewicz (2017) Complex-valued wavelet lifting and applications.
\emph{Technometrics (to appear)} DOI 10.1080/00401706.2017.1281846.\cr\cr
For related literature on the lifting methodology adopted in the technique, see\cr\cr
Nunes, M. A., Knight, M. I and Nason, G. P. (2006) Adaptive lifting for nonparametric regression.
\emph{Stat. Comput.} \bold{16} (2), 143--159.\cr\cr
Knight, M. I. and Nason, G. P. (2009) A 'nondecimated' wavelet transform. \emph{Stat. Comput.}
\bold{19} (1), 1--16.
}
\author{
Matt Nunes, Marina Knight}

\seealso{
\code{\link[adlift]{AdaptNeigh}}, \code{\link[adlift]{AdaptPred}}, \code{\link[adlift]{CubicPred}}, 
\code{\link{denoisepermC}}, 
\code{\link{denoisepermCh}}, 
\code{\link[adlift]{LinearPred}}, 
\code{\link{orthpredfilters}}, 
\code{\link[adlift]{QuadPred}}     
}

\examples{
    library(adlift)

    # construct an (irregular) observation grid
     x<-runif(256)

     #construct a signal
     f<-make.signal2("blocks",x=x)

     fwd<-fwtnppermC(x,f,LocalPred=AdaptPred,neigh=1,closest=FALSE)

     # have a look at the complex-valued coefficients and the removal trajectory:

     fwd$coeffv

     fwd$removelist	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

