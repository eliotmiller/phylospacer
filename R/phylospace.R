#' Create a phylospace object.
#'
#' Given a phylo object, a data frame of species trait values, and another of
#' reconstructed node values, returns a phylospace object.
#'
#' @param ape.phylo Phylo object
#' @param species.niches Two column data frame (or matrix) with species names as row names  
#' and exactly matching phylogeny. Columns as traits of interest.
#' @param node.niches Two column data frame (or matrix) with node names as row names  
#' and exactly matching phylogeny. Columns as traits of interest.
#' @param subset.node Optional numeric value (possible character vector if nodes are so
#' named?) specifying a node at which to subset the phylospace object.
#' @param subset.to.root Logical argument used for when subset.node is specified. If TRUE,
#' then all edges between the subset node and the root of the entire phylogeny will also
#' be retained in the resulting phylospace object.
#' @param jitter.level Optional numeric value to jitter both species and node values. 
#' Useful for when species' and/or nodes have very similar trait values.
#' @param replacement.colors Optional vector of replacement colors by which to scale the
#' phylospace. If not specified, defaults to blue to cyan to green to yellow to red.
#' Otherwise, a vector of 15 numbers needs to be provided. This vector needs to be
#' in the form of color1 red, color1 green, color1 blue, color2 red, color2 green, color2
#' blue...color5 red, color5 green, color5 blue. The value of each of the 15 numbers needs
#' to be between 0 and 1, e.g. if you want color1 to be red, the first three numbers need
#' to be 1,0,0, and if you want color2 to be yellow, the second three numbers need to be
#' 1,1,0. Color1 refers to color closest to root, while color5 is the color of the tips. 
#' See examples for futher details. 
#' 
#' @details This function does not currently reconstruct ancestral states. These need to
#' be calculated beforehand with any of a variety of functions, e.g. fastAnc() in phytools
#' or ace() in ape. It is intended to take two data frames, with the species names of the
#' species.niches frame exactly matching the phylogeny. Importantly, these need to be
#' specified as the row names of the data frame, not as a separate column. Similarly, the
#' reconstructed node names need to be the row names of the node.niches data frame. The
#' subset.node must be a single number (or object) corresponding to a unique node in the  
#' tree. There are a variety of ways you could find this, including simply plotting the 
#' node labels onto a tree. Alternatively, consider the findMRCA function in phytools. 
#' Because species may have similar trait values, a quick and dirty jitter option is 
#' available, though this has the effect of shifting all points, internal nodes included. 
#' Thus, users may choose to slightly shift the points in question manually.
#' The function is still in development version, please report bugs to me. 
#'
#' @return Returns a phylospace object, which is a list with five elements: the original
#' phylo object; a data frame with 14 columns describing the colors and bounds of all line
#' segments; a character vector of all species retained; the identity of the subset
#' node, if there is one; a vector of all descendant nodes from the subset node, and
#' including also all ancestors of the subset node if subset.to.root = TRUE
#'
#' @export
#'
#' @import ape phylobase plotrix
#'
#' @references Miller, E.T., A.E. Zanne, & R. E. Ricklefs. In press. Niche conservatism 
#' constrains Australian honeyeater assemblages in stressful environments. Ecology Letters. 
#'
#' @examples
#' library(ape)
#'
#' #simulate tree with birth-death process
#' tree <- rbdtree(birth=0.1, death=0, Tmax=40)
#'
#' #prune the phylogeny down to 50 species. if you get an error here, re-run.
#' tree <- drop.tip(tree, tip=51:length(tree$tip.label))
#'
#' #simulate trait evolution up tree with Brownian motion process
#' trait1 <- rTraitCont(tree, model="BM")
#' trait2 <- rTraitCont(tree, model="BM")
#'
#' #bind the traits together into a matrix
#' species.niches <- cbind(trait1, trait2)
#'
#' #reconstruct ancestral states
#' nodes.trait1 <- ace(trait1, tree, type="continuous", method="REML")
#' nodes.trait2 <- ace(trait2, tree, type="continuous", method="REML")
#'
#' node.niches <- cbind(nodes.trait1$ace, nodes.trait2$ace)
#'
#' #simulate available climate space
#' climate.points <- cbind(rnorm(500, mean=0, sd=1), rnorm(100, mean=0, sd=1))
#'
#' #calculate a phylospace object for the entire phylogeny
#' entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches)
#'
#' #plot the entire phylospace, with species labels plotted and slightly offset from tips, and with available climate space in background
#' plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)
#'
#' #plot the same phylospace as previous example, but with a green scale
#' entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, replacement.colors=c(0,0.1,0,0,0.4,0,0,0.6,0,0,0.8,0,0,1,0))
#' plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)
#'
#' #plot the same phylospace as previous example, but with color scale from black to blue to purple to red to orange
#' black <- c(0,0,0)
#' blue <- c(0,0,1)
#' purple <- c(1,0,1)
#' red <- c(1,0,0)
#' orange <- c(1,0.5,0)
#' entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, replacement.colors=c(black, blue, purple, red, orange))
#' plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)
#'
#' #plot the same phylospace as previous example, but in grayscale
#' entire <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, replacement.colors=c(0.95,0.95,0.95, 0.75,0.75,0.75, 0.5,0.5,0.5, 0.25,0.25,0.25, 0,0,0))
#' plot(entire, species.labels=TRUE, climate.points=climate.points, label.adjust=0.05, lwd=2)
#'
#' #Example of how one could create a series of panels for an animation. Here we are retaining all branches from each subset node to the root
#' cladeA <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, subset.node=63, subset.to.root=TRUE)
#' cladeB <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, subset.node=75, subset.to.root=TRUE)
#' cladeC <- phylospace(tree, species.niches=species.niches, node.niches=node.niches, subset.node=90, subset.to.root=TRUE)
#'
#' #These x & y limits can most likely be made more restrictive, depending on the simulated trait data. Keeping broad here, but modify if desired.
#' #Save each plot before calling next.
#' plot(phylospace.object=cladeA, species.labels=TRUE, label.adjust=0.1, x.label="Trait 1", y.label="Trait 2", x.limits=c(-2,2), y.limits=c(-2,2), lwd=2)
#' plot(phylospace.object=cladeB, species.labels=TRUE, label.adjust=0.1, x.label="Trait 1", y.label="Trait 2", x.limits=c(-2,2), y.limits=c(-2,2), lwd=2)
#' plot(phylospace.object=cladeC, species.labels=TRUE, label.adjust=0.1, x.label="Trait 1", y.label="Trait 2", x.limits=c(-2,2), y.limits=c(-2,2), lwd=2)
#'
#' #Example of how to quickly plot phylospace for a small clade without the edges to the root.
#' plot(phylospace(tree, species.niches, node.niches, subset.node=85), lwd=2)

phylospace <- function(ape.phylo, species.niches, node.niches, subset.node, subset.to.root=FALSE, jitter.level=0, replacement.colors)
{	
	##convert to phylobase phylo. the suppress warnings command is because if one makes a tree ultrametric it can come with some unexpected parameters that phylobase doesn't know how to deal with. just ignores them and all is fine, but no reason to print warnings.	
	phylobase.phylo <- suppressWarnings(as(ape.phylo,"phylo4"))

	##derive a vector of species' names to subset later
	keep.species <- ape.phylo$tip.label

	##allow people to jitter the final results so if species have the exact same traits the tips can be distinguished. jitter will be normally set to zero
	species.niches <- cbind(jitter(species.niches[,1], factor=jitter.level), jitter(species.niches[,2], factor=jitter.level))
	node.niches <- cbind(jitter(node.niches[,1], factor=jitter.level), jitter(node.niches[,2], factor=jitter.level))
	
	##add the internal node values onto the end of a vector of the species values for each trait
	temp.trait1 <- c(species.niches[,1], node.niches[,1])
	temp.trait2 <- c(species.niches[,2], node.niches[,2])

	##combine those two vectors alonge with a column for the node names
	lookup.table <- data.frame(1:length(temp.trait1), temp.trait1, temp.trait2)
	row.names(lookup.table) <- NULL
	names(lookup.table) <- c("node","col1","col2")

	##make a new data frame that details the node to node connections. do not use the label, edge length or the node type columns
	segments.to.plot <- data.frame(phylobase.phylo@edge)
	
	##the root edge connects to node "0", which isn't in either species or node niches. however, we still want to plot the root on there
	##we would lose it on the merge command later because we merge on ancestor (i.e. 0) and lose the whole row
	##so, basically make a fake little branch that just goes from root to root (remember there will always be one branch less than number of nodes in a fully dichotomous tree)
	
	segments.to.plot$ancestor[segments.to.plot$ancestor==0] <- length(ape.phylo$tip.label)+1

	##merge segments.to.plot with the lookup table. note that the by.x & by.y arguments refer not to columns and row but to first dataframe (x) and second (y)
	##add in four new columns: x0, y0, x1, y1; to be used with drawing segments to connect the nodes. x = col1, y = col2. we will fill these new columns using the lookup table created above
	##first add in the x,y coordinates for the ancestor
	segments.to.plot <- merge(segments.to.plot, lookup.table, by.x="ancestor", by.y="node")
	names(segments.to.plot)[3:4] = c("x0","y0")

	##now add in the x,y coordinates for the descendant
	segments.to.plot <- merge(segments.to.plot, lookup.table, by.x="descendant", by.y="node")
	names(segments.to.plot)[5:6] = c("x1","y1")
	
	##color the branches by how deep they are in the phylogeny. for every node, will derive its distance from the root. will look at the distribution of these distances and break into 5 categories (blue, cyan, green, yellow, red). whether or not one uses an ultrametric tree here has a big influence on results

	all.dist <- dist.nodes(ape.phylo)

	##get the distances between the interior nodes and the root AND the tips and the root
	root.dist <- all.dist[length(ape.phylo$tip.label)+1, ]

	##use quantile to determine where breaks are. call tips "red" and divide the internal nodes four categories

	to.break <- root.dist[(length(ape.phylo$tip.label)+1):length(root.dist)]

	one <- quantile(to.break, 0.25)
	two <- quantile(to.break, 0.50)
	three <- quantile(to.break, 0.75)
	four <- quantile(to.break, 1)

	##make an empty character vector for use in coloring branches and fill based on distance of respective node to root. 
	node.color <- character(length(root.dist))
	names(node.color) <- names(root.dist)

	if(missing(replacement.colors))
	{
		##make an empty character vector for use in coloring branches and fill based on distance of respective node to root. 
		node.color[root.dist >= 0 & root.dist < one] <- "blue"
		node.color[root.dist >= one & root.dist < two] <- "cyan"
		node.color[root.dist >= two & root.dist < three] <- "green"
		node.color[root.dist >= three & root.dist <= four] <- "yellow"
		node.color[root.dist > four] <- "red"
	}
	else if(!missing(replacement.colors))
	{
		node.color[root.dist >= 0 & root.dist < one] <- "color1"
		node.color[root.dist >= one & root.dist < two] <- "color2"
		node.color[root.dist >= two & root.dist < three] <- "color3"
		node.color[root.dist >= three & root.dist <= four] <- "color4"
		node.color[root.dist > four] <- "color5"
	}

	##bind the character vector colors back into the segments to plot dataframe, first by ancestor, then by descendent
	temp <- as.data.frame(node.color)
	temp <- cbind(temp, 1:length(root.dist))
	names(temp)[2]="num"
	segments.to.plot <- merge(segments.to.plot, temp, by.x="ancestor", by.y="num")
	names(segments.to.plot)[7]="from"
	segments.to.plot <- merge(segments.to.plot, temp, by.x="descendant", by.y="num")
	names(segments.to.plot)[8]="to"

	##need to coerce these to characters for it to plot right
	segments.to.plot$from <- as.character(segments.to.plot$from)
	segments.to.plot$to <- as.character(segments.to.plot$to)

	##add six new columns here. The first three will specify in R,G,B space the color from the "from" column, the second are for the "to" column

	if(missing(replacement.colors))
	{
		segments.to.plot$from.r[segments.to.plot$from == "blue"] <- 0
		segments.to.plot$from.g[segments.to.plot$from == "blue"] <- 0
		segments.to.plot$from.b[segments.to.plot$from == "blue"] <- 1
		
		segments.to.plot$from.r[segments.to.plot$from == "cyan"] <- 0
		segments.to.plot$from.g[segments.to.plot$from == "cyan"] <- 1
		segments.to.plot$from.b[segments.to.plot$from == "cyan"] <- 1
		
		segments.to.plot$from.r[segments.to.plot$from == "green"] <- 0
		segments.to.plot$from.g[segments.to.plot$from == "green"] <- 1
		segments.to.plot$from.b[segments.to.plot$from == "green"] <- 0

		segments.to.plot$from.r[segments.to.plot$from == "yellow"] <- 1
		segments.to.plot$from.g[segments.to.plot$from == "yellow"] <- 1
		segments.to.plot$from.b[segments.to.plot$from == "yellow"] <- 0

		segments.to.plot$from.r[segments.to.plot$from == "red"] <- 1
		segments.to.plot$from.g[segments.to.plot$from == "red"] <- 0
		segments.to.plot$from.b[segments.to.plot$from == "red"] <- 0

		segments.to.plot$to.r[segments.to.plot$to == "blue"] <- 0
		segments.to.plot$to.g[segments.to.plot$to == "blue"] <- 0
		segments.to.plot$to.b[segments.to.plot$to == "blue"] <- 1
		
		segments.to.plot$to.r[segments.to.plot$to == "cyan"] <- 0
		segments.to.plot$to.g[segments.to.plot$to == "cyan"] <- 1
		segments.to.plot$to.b[segments.to.plot$to == "cyan"] <- 1
		
		segments.to.plot$to.r[segments.to.plot$to == "green"] <- 0
		segments.to.plot$to.g[segments.to.plot$to == "green"] <- 1
		segments.to.plot$to.b[segments.to.plot$to == "green"] <- 0

		segments.to.plot$to.r[segments.to.plot$to == "yellow"] <- 1
		segments.to.plot$to.g[segments.to.plot$to == "yellow"] <- 1
		segments.to.plot$to.b[segments.to.plot$to == "yellow"] <- 0

		segments.to.plot$to.r[segments.to.plot$to == "red"] <- 1
		segments.to.plot$to.g[segments.to.plot$to == "red"] <- 0
		segments.to.plot$to.b[segments.to.plot$to == "red"] <- 0
	}

	else
	{
		segments.to.plot$from.r[segments.to.plot$from == "color1"] <- replacement.colors[1]
		segments.to.plot$from.g[segments.to.plot$from == "color1"] <- replacement.colors[2]
		segments.to.plot$from.b[segments.to.plot$from == "color1"] <- replacement.colors[3]
		
		segments.to.plot$from.r[segments.to.plot$from == "color2"] <- replacement.colors[4]
		segments.to.plot$from.g[segments.to.plot$from == "color2"] <- replacement.colors[5]
		segments.to.plot$from.b[segments.to.plot$from == "color2"] <- replacement.colors[6]
		
		segments.to.plot$from.r[segments.to.plot$from == "color3"] <- replacement.colors[7]
		segments.to.plot$from.g[segments.to.plot$from == "color3"] <- replacement.colors[8]
		segments.to.plot$from.b[segments.to.plot$from == "color3"] <- replacement.colors[9]

		segments.to.plot$from.r[segments.to.plot$from == "color4"] <- replacement.colors[10]
		segments.to.plot$from.g[segments.to.plot$from == "color4"] <- replacement.colors[11]
		segments.to.plot$from.b[segments.to.plot$from == "color4"] <- replacement.colors[12]

		segments.to.plot$from.r[segments.to.plot$from == "color5"] <- replacement.colors[13]
		segments.to.plot$from.g[segments.to.plot$from == "color5"] <- replacement.colors[14]
		segments.to.plot$from.b[segments.to.plot$from == "color5"] <- replacement.colors[15]

		segments.to.plot$to.r[segments.to.plot$to == "color1"] <- replacement.colors[1]
		segments.to.plot$to.g[segments.to.plot$to == "color1"] <- replacement.colors[2]
		segments.to.plot$to.b[segments.to.plot$to == "color1"] <- replacement.colors[3]
		
		segments.to.plot$to.r[segments.to.plot$to == "color2"] <- replacement.colors[4]
		segments.to.plot$to.g[segments.to.plot$to == "color2"] <- replacement.colors[5]
		segments.to.plot$to.b[segments.to.plot$to == "color2"] <- replacement.colors[6]
		
		segments.to.plot$to.r[segments.to.plot$to == "color3"] <- replacement.colors[7]
		segments.to.plot$to.g[segments.to.plot$to == "color3"] <- replacement.colors[8]
		segments.to.plot$to.b[segments.to.plot$to == "color3"] <- replacement.colors[9]

		segments.to.plot$to.r[segments.to.plot$to == "color4"] <- replacement.colors[10]
		segments.to.plot$to.g[segments.to.plot$to == "color4"] <- replacement.colors[11]
		segments.to.plot$to.b[segments.to.plot$to == "color4"] <- replacement.colors[12]

		segments.to.plot$to.r[segments.to.plot$to == "color5"] <- replacement.colors[13]
		segments.to.plot$to.g[segments.to.plot$to == "color5"] <- replacement.colors[14]
		segments.to.plot$to.b[segments.to.plot$to == "color5"] <- replacement.colors[15]
	}

	##if only a specific clade within the whole phylogeny is to be plotted, subset segments to plot accordingly
	if(missing(subset.node))
	{
		subset.node <- "none"
		print("Retaining entire phylogeny")
	}
	else if(!missing(subset.node) & subset.to.root==FALSE)
	{
		keep.branches <- descendants(phylobase.phylo, subset.node, type="all")
		keep.species <- descendants(phylobase.phylo, subset.node, type="tips")
		segments.to.plot <- segments.to.plot[segments.to.plot$descendant %in% keep.branches, ]
	}
	else
	{
		keep.branches <- descendants(phylobase.phylo, subset.node, type="all")
		to.root <- ancestors(phylobase.phylo, subset.node, "ALL") ##calling it like this includes also the subset node
		keep.branches <- c(rev(to.root), keep.branches)
		keep.species <- descendants(phylobase.phylo, subset.node, type="tips")
		segments.to.plot <- segments.to.plot[segments.to.plot$descendant %in% keep.branches, ]
	}

	##switch the order ancestor and descendant appear in data frame for ease of reading
	
	segments.to.plot <- segments.to.plot[ , c("ancestor","descendant","x0","y0","x1","y1","from","to","from.r","from.g","from.b","to.r","to.g","to.b")]

	##derive a vector of species' names based on whatever taxa are left to be plotted
	name.vector <- phylobase.phylo@label[segments.to.plot$descendant][!is.na(phylobase.phylo@label[segments.to.plot$descendant])]

	##derive a vector of node names based on whatever nodes are left to be plotted
	node.vector <- segments.to.plot$descendant[segments.to.plot$descendant > length(ape.phylo$tip.label)]

	if(missing(replacement.colors))
	{
		replacement.colors <- "default"
	}
	else
	{
		replacement.colors <- replacement.colors
	} 

	output <- list(ape.phylo=ape.phylo, segments.to.plot=segments.to.plot, name.vector=name.vector, subset.node=subset.node, node.vector=node.vector, replacement.colors=replacement.colors)

	class(output) <- "phylospace"

	return(output)
}
