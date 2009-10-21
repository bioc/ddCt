\name{replaceVectorByEquality}
\alias{replaceVectorByEquality}
\title{REPLACE ITEMS OF VECTOR BY EQUALITY}
\description{
  The function replces (or updates) the items of a given vector by
  checking the equality with the \code{target} parameter. If found, the
  item will be replaced by the \code{value} parameter. The length of
  both \code{target} and \code{value} must be the same and could be
  longer than 1, in which case the replace will be iterated.
}
\usage{
replaceVectorByEquality(vector, target, value)
}
\arguments{
  \item{vector}{A vector to be replaced. The items of the vector must be
  atom types, since the equality is checked by '=='.}
  \item{target}{targets to be replaced, could be either single or a
    vector}
  \item{value}{values to be replaced at the positions of targets, must
    be of the same length of \code{target}}
}
\details{
  A warning will be prompted if any item in the \code{target} cannot be found
}
\value{
  A vector of the same length as the parameter \code{vector}
}
\author{
  Jitao David Zhang
}
\seealso{
  \code{==} for checking equality.
}
\examples{
vector <- c("java", "perl", "python", "c#")
replaceVectorByEquality(vector, target="c#", value="c/c++")
replaceVectorByEquality(vector, target=c("c#","perl"), value=c("c/c++","R"))
}
