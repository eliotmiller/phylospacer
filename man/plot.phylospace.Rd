\name{plot.phylospace}
\alias{plot.phylospace}
\title{Plot a phylospace object.}
\usage{
  \method{plot}{phylospace} (phylospace.object,
    climate.points, species.labels = FALSE,
    node.labels = FALSE, label.adjust = 0, x.label,
    y.label, x.limits, y.limits, ...)
}
\arguments{
  \item{phylospace.object}{Phylospace object. The results
  of a call to phylospace.}

  \item{climate.points}{Two column data frame (or matrix?)
  of x, y points. Columns in same order as species niche
  data frame. Here, each row refers to separate climate
  point.}

  \item{species.labels}{Logical indicating whether species
  names should be plotted.}

  \item{node.labels}{Logical indicating whether node names
  should be plotted.}

  \item{label.adjust}{Numeric value by which species labels
  should be offset from their points. Currently only works
  well if x & y axes are on similar scales.}

  \item{x.label}{Optional label for x axis. Defaults to
  "trait1".}

  \item{y.label}{Optional label for y axis. Defaults to
  "trait2".}

  \item{x.limits}{Optional x limits for the plot.}

  \item{y.limits}{Optional y limits for the plot.}

  \item{...}{Some arguments can be passed to graphical
  parameters, in particular "lwd". See: \code{\link{par}}}
}
\value{
  Returns a plotted phylospace object.
}
\description{
  Plots a phylospace object and, if desired, available
  climate space, species labels, and/or node labels.
}
\details{
  This function plots an object of class phylospace.
  Because currently the label adjust argument calculates
  how much to offset the labels from the species' points as
  a function of absolute slope (not slope as you see it
  plotted), this currently only works well if the x and y
  axes are on similar scales. This will be updated in
  future versions.
}
\examples{
#See examples under phylospace function
}
\references{
  Miller, E.T., A.E. Zanne, & R. E. Ricklefs. In press.
  Niche conservatism constrains Australian honeyeater
  assemblages in stressful environments. Ecology Letters.
}

