\name{getDir}
\alias{getDir}
\title{Auxillary functions for the executive scripts}
\description{
  \code{getDir} creates a directory in case it does not exist and
  returns the directory name.
}
\usage{getDir(dir, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{Directory name}
  \item{\dots}{Other parameters passed to \code{\link{dir.create}}}
}
\details{
  Auxillary functions
}
\value{
  \code{getDir} returns the directory name
}
\author{Jitao David Zhang <jitao_david.zhang@roche.com>}
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
getDir(tempdir())
}
