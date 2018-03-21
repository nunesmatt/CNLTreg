\name{orthpredfilters}
\alias{orthpredfilters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes orthogonal filters
}
\description{
Given a filter L, finds a second filter M, orthogonal to L and with unit norm
}
\usage{
orthpredfilters(filter = c(0.5, 1, 0.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{An initial filter L}
}
\details{
	See Hamilton et al. (2017), section 2.2.
}
\value{
A matrix with two rows, the first row corresponding to L, the second corresponding to the orthogonal filter M.
}
\references{
Hamilton, J., Nunes, M. A., Knight, M. I. and Fryzlewicz, P. (2017) Complex-valued wavelet lifting and applications.
\emph{Technometrics (to appear)} DOI 10.1080/00401706.2017.1281846.
}
\author{
Marina Knight, Matt Nunes
}

\section{Warning}{At present only works with odd length filters}

\seealso{
	\code{\link{fwtnppermC}}
}
\examples{

# create a vector representing a filter for one neighbour either side of a removed point
# (equally weighted):

L = c(0.5, 1, 0.5)

# now work out a unit-norm filter orthogonal to L

out <- orthpredfilters(L)

# M should be the second row:

out[2,]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
