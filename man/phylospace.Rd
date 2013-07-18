\name{phylospace}
\alias{phylospace}
\title{Create a phylospace object.}
\usage{
  phylospace(ape.phylo, species.niches, node.niches,
    subset.node, subset.to.root = FALSE, jitter.level = 0,
    replacement.colors)
}
\arguments{
  \item{ape.phylo}{Phylo object}

  \item{species.niches}{Two column data frame (or matrix)
  with species names as row names and exactly matching
  phylogeny. Columns as traits of interest.}

  \item{node.niches}{Two column data frame (or matrix) with
  node names as row names and exactly matching phylogeny.
  Columns as traits of interest.}

  \item{subset.node}{Optional numeric value (possible
  character vector if nodes are so named?) specifying a
  node at which to subset the phylospace object.}

  \item{subset.to.root}{Logical argument used for when
  subset.node is specified. If TRUE, then all edges between
  the subset node and the root of the entire phylogeny will
  also be retained in the resulting phylospace object.}

  \item{jitter.level}{Optional numeric value to jitter both
  species and node values. Useful for when species' and/or
  nodes have very similar trait values.}

  \item{replacement.colors}{Optional vector of replacement
  colors by which to scale the phylospace. If not
  specified, defaults to blue to cyan to green to yellow to
  red. Otherwise, a vector of 15 numbers needs to be
  provided. This vector needs to be in the form of color1
  red, color1 green, color1 blue, color2 red, color2 green,
  color2 blue...color5 red, color5 green, color5 blue. The
  value of each of the 15 numbers needs to be between 0 and
  1, e.g. if you want color1 to be red, the first three
  numbers need to be 1,0,0, and if you want color2 to be
  yellow, the second three numbers need to be 1,1,0. Color1
  refers to color closest to root, while color5 is the
  color of the tips. See examples for futher details.}
}
\value{
  Returns a phylospace object, which is a list with five
  elements: the original phylo object; a data frame with 14
  columns describing the colors and bounds of all line
  segments; a character vector of all species retained; the
  identity of the subset node, if there is one; a vector of
  all descendant nodes from the subset node, and including
  also all ancestors of the subset node if subset.to.root =
  TRUE
}
\description{
  Given a phylo object, a data frame of species trait
  values, and another of reconstructed node values, returns
  a phylospace object.
}
\details{
  This function does not currently reconstruct ancestral
  states. These need to be calculated beforehand with any
  of a variety of functions, e.g. fastAnc() in phytools or
  ace() in ape. It is intended to take two data frames,
  with the species names of the species.niches frame
  exactly matching the phylogeny. Importantly, these need
  to be specified as the row names of the data frame, not
  as a separate column. Similarly, the reconstructed node
  names need to be the row names of the node.niches data
  frame. The subset.node must be a single number (or
  object) corresponding to a unique node in the tree. There
  are a variety of ways you could find this, including
  simply plotting the node labels onto a tree.
  Alternatively, consider the findMRCA function in
  phytools. Because species may have similar trait values,
  a quick and dirty jitter option is available, though this
  has the effect of shifting all points, internal nodes
  included. Thus, users may choose to slightly shift the
  points in question manually. The function is still in
  development version, please report bugs to me.
}
\examples{
library(ape)

#simulate tree with birth-death process
tree <- rbdtree(birth=0.1, death=0, Tmax=40)

#prune the phylogeny down to 50 species. if you get an error here, re-run.
tree <- drop.tip(tree, tip=51:length(tree$tip.label))

#simulate trait evolution up tree with Brownian motion process
trait1 <- rTraitCont(tree, model="BM")
trait2 <- rTraitCont(tree, model="BM")

#bind the traits together into a matrix
species.niches <- cbind(trait1, trait2)

#reconstruct ancestral states
nodes.trait1 <- ace(trait1, tree, type="continuous", method="REML")
nodes.trait2 <- ace(trait2, tree, type="continuous", method="REML")

node.niches <- cbind(nodes.trait1$ace, nodes.trait2$ace)

#simulate available climate space
climate.points <- cbind(rnorm(500, mean=0, sd=1), rnorm(100, mean=0, sd=1))

#calculate a phylospace object for the entire phylogeny
entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches)

#plot the entire phylospace, with species labels plotted and slightly offset from tips, and with available climate space in background
plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)

#plot the same phylospace as previous example, but with a green scale
entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, replacement.colors=c(0,0.1,0,0,0.4,0,0,0.6,0,0,0.8,0,0,1,0))
plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)

#plot the same phylospace as previous example, but with color scale from black to blue to purple to red to orange
black <- c(0,0,0)
blue <- c(0,0,1)
purple <- c(1,0,1)
red <- c(1,0,0)
orange <- c(1,0.5,0)
entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, replacement.colors=c(black, blue, purple, red, orange))
plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)

#plot the same phylospace as previous example, but in grayscale
entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, replacement.colors=c(0.95,0.95,0.95, 0.75,0.75,0.75, 0.5,0.5,0.5, 0.25,0.25,0.25, 0,0,0))
plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)

#Example of how one could create a series of panels for an animation. Here we are retaining all branches from each subset node to the root
cladeA <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, subset.node=63, subset.to.root=TRUE)
cladeB <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, subset.node=75, subset.to.root=TRUE)
cladeC <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, subset.node=90, subset.to.root=TRUE)

#These x & y limits can most likely be made more restrictive, depending on the simulated trait data. Keeping broad here, but modify if desired.
#Save each plot before calling next.
plot(phylospace.object=cladeA, species.labels=TRUE, label.adjust=0.1, x.label="Trait 1", y.label="Trait 2", x.limits=c(-2,2), y.limits=c(-2,2), lwd=2)
plot(phylospace.object=cladeB, species.labels=TRUE, label.adjust=0.1, x.label="Trait 1", y.label="Trait 2", x.limits=c(-2,2), y.limits=c(-2,2), lwd=2)
plot(phylospace.object=cladeC, species.labels=TRUE, label.adjust=0.1, x.label="Trait 1", y.label="Trait 2", x.limits=c(-2,2), y.limits=c(-2,2), lwd=2)

#Example of how to quickly plot phylospace for a small clade without the edges to the root.
plot(phylospace(tree, species.niches, node.niches, subset.node=85), lwd=2)
}
\references{
  Miller, E.T., A.E. Zanne, & R. E. Ricklefs. In press.
  Niche conservatism constrains Australian honeyeater
  assemblages in stressful environments. Ecology Letters.
}

