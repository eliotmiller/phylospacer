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
#' #prune the phylogeny down to 50 species
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

phylospace <- function(ape.phylo, species.niches, node.niches, subset.node, subset.to.root=FALSE, jitter.level=0)
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
	node.color[root.dist >= 0 & root.dist < one] <- "blue"
	node.color[root.dist >= one & root.dist < two] <- "cyan"
	node.color[root.dist >= two & root.dist < three] <- "green"
	node.color[root.dist >= three & root.dist <= four] <- "yellow"
	node.color[root.dist > four] <- "red"

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

	for(i in 1:dim(segments.to.plot)[1])
	{
		if(segments.to.plot$from[i] == "blue")
		{
			segments.to.plot$from.r[i] = 0
			segments.to.plot$from.g[i] = 0
			segments.to.plot$from.b[i] = 1
		}
		else if(segments.to.plot$from[i] == "cyan")
		{
			segments.to.plot$from.r[i] = 0
			segments.to.plot$from.g[i] = 1
			segments.to.plot$from.b[i] = 1
		}
		else if(segments.to.plot$from[i] == "green")
		{
			segments.to.plot$from.r[i] = 0
			segments.to.plot$from.g[i] = 1
			segments.to.plot$from.b[i] = 0
		}
		else if(segments.to.plot$from[i] == "yellow")
		{
			segments.to.plot$from.r[i] = 1
			segments.to.plot$from.g[i] = 1
			segments.to.plot$from.b[i] = 0
		}
		else if(segments.to.plot$from[i] == "red")
		{
			segments.to.plot$from.r[i] = 1
			segments.to.plot$from.g[i] = 0
			segments.to.plot$from.b[i] = 0
		}
	}

	for(i in 1:dim(segments.to.plot)[1])
	{
		if(segments.to.plot$to[i] == "blue")
		{
			segments.to.plot$to.r[i] = 0
			segments.to.plot$to.g[i] = 0
			segments.to.plot$to.b[i] = 1
		}
		else if(segments.to.plot$to[i] == "cyan")
		{
			segments.to.plot$to.r[i] = 0
			segments.to.plot$to.g[i] = 1
			segments.to.plot$to.b[i] = 1
		}
		else if(segments.to.plot$to[i] == "green")
		{
			segments.to.plot$to.r[i] = 0
			segments.to.plot$to.g[i] = 1
			segments.to.plot$to.b[i] = 0
		}
		else if(segments.to.plot$to[i] == "yellow")
		{
			segments.to.plot$to.r[i] = 1
			segments.to.plot$to.g[i] = 1
			segments.to.plot$to.b[i] = 0
		}
		else if(segments.to.plot$to[i] == "red")
		{
			segments.to.plot$to.r[i] = 1
			segments.to.plot$to.g[i] = 0
			segments.to.plot$to.b[i] = 0
		}
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

	output <- list(ape.phylo=ape.phylo, segments.to.plot=segments.to.plot, name.vector=name.vector, subset.node=subset.node, node.vector=node.vector)

	class(output) <- "phylospace"

	return(output)
}
