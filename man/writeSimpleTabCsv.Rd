\name{writeSimpleTabCsv}
\title{ Write a data frame into a tab delimited file }
\alias{writeSimpleTabCsv}
\description{
   Write a 'data.frame' into a tab delimited file (not quoted and no-row-name CSV file)
}
\usage{
writeSimpleTabCsv(x, file, ...)
}
\arguments{
  \item{x}{ 'data.frame' }
  \item{file}{ character. File name. }
  \item{...}{ Additional arguments passed onto the function }
}
\value{
  The function is called for its side effect: writing a file
}
\author{ Wolfgang Huber }
