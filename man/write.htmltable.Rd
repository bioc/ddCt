\name{write.htmltable}
\title{ Write a data frame into an html table within a html page }
\alias{write.htmltable}
\description{
   Write a 'data.frame' into an html table within a html page
}
\usage{
write.htmltable(x, filename, title = "", sortby = NULL, decreasing = TRUE, open = "wt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ 'data.frame' }
  \item{filename}{ character. File name. }
  \item{title}{ character. Title of html page }
  \item{sortby}{ character. Name of column by which to sort the table rows }
  \item{decreasing}{ logical. Should the sort order be increasing or decreasing? }
  \item{open}{ character. This argument is passed onto 'file' }
}
\value{
  The funciton is called for its side effect: writing a file
}
\author{ Wolfgang Huber }
