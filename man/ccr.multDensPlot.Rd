\name{ccr.multDensPlot}
\alias{ccr.multDensPlot}
\title{
    Mutiple shaded density plot
}
\description{
This functions plots multiple distribution densities with solid colors for the curves and shaded colors for underlying areas.
}
\usage{
ccr.multDensPlot(TOPLOT, COLS,
                 XLIMS, TITLE, LEGentries, XLAB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TOPLOT}{
A list of density object computed using the \code{density} function of the \code{stats} package.
}
  \item{COLS}{
  A vector of colors of the same length of \code{TOPLOT} that are used to plot the density curves. Alpha-reduced versions of these colors are used to fill the underlying areas.
}
  \item{XLIMS}{
    A vector of two numerical values optionally specifying x-axis limits (NULL by default).
}
  \item{TITLE}{
A string containing the plot title.
}
  \item{LEGentries}{
A vector of strings (one per each density in \code{TOPLOT}) specifying corresponding legend entries.
}
  \item{XLAB}{
  A string containing the x-axis label.
}
}
\author{
Francesco Iorio (francesco.iorio@fht.org)}
\examples{\dontshow{setwd(tempdir())}
## generating random data
x <- rnorm(1000, 0, 0.5)
y <- rnorm(1000, 2, 0.4)
z <- rnorm(1000, -1, 1.5)


## assembling kernel estimated distributions into a list
ToPlot<-list(x=density(x),y=density(y),z=density(z))

## density visualisation
ccr.multDensPlot(ToPlot,COLS = c('red','blue','gray'),
    TITLE = 'example',LEGentries = c('x','y','z'),
    XLIMS = c(-5,3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Assessment and Visualisation}% use one of  RShowDoc("KEYWORDS")

