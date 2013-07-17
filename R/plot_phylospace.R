#' Plot a phylospace object.
#'
#' Plots a phylospace object and, if desired, available climate space, species labels,
#' and/or node labels. 
#'
#' @param phylospace.object Phylospace object. The results of a call to phylospace.
#' @param climate.points Two column data frame (or matrix?) of x, y points. Columns in
#' same order as species niche data frame. Here, each row refers to separate climate point.
#' @param species.labels Logical indicating whether species names should be plotted.
#' @param node.labels Logical indicating whether node names should be plotted.
#' @param label.adjust Numeric value by which species labels should be offset from their
#' points. Currently only works well if x & y axes are on similar scales.
#' @param x.label Optional label for x axis. Defaults to "trait1".
#' @param y.label Optional label for y axis. Defaults to "trait2".
#' @param x.limits Optional x limits for the plot.
#' @param y.limits Optional y limits for the plot.
#' @param ... Some arguments can be passed to graphical parameters, in particular "lwd". See: \code{\link{par}}
#' 
#' @details This function plots an object of class phylospace. Because currently the label
#' adjust argument calculates how much to offset the labels from the species' points as a
#' function of absolute slope (not slope as you see it plotted), this currently only works
#' well if the x and y axes are on similar scales. This will be updated in future
#' versions. 
#'
#' @return Returns a plotted phylospace object.
#'
#' @export
#'
#' @method plot phylospace
#'
#' @references Miller, E.T., A.E. Zanne, & R. E. Ricklefs. In press. Niche conservatism 
#' constrains Australian honeyeater assemblages in stressful environments. Ecology Letters. 
#'
#' @examples
#' #See examples under phylospace function

plot.phylospace <- function(phylospace.object, climate.points, species.labels=FALSE, node.labels=FALSE, label.adjust=0, x.label, y.label, x.limits, y.limits, ...)
{
	if(missing(x.label) & missing(y.label))
	{
		x.label <- "trait1"
		y.label <- "trait2"
	}
	else if(missing(x.label))
	{
		x.label <- "trait1"
	}
	else if(missing(y.label))
	{
		y.label <- "trait2"
	}

	##derive the x,y locations of all points and nodes and their respective labels
	slopes <- (phylospace.object$segments.to.plot[,6]-phylospace.object$segments.to.plot[,4])/(phylospace.object$segments.to.plot[,5]-phylospace.object$segments.to.plot[,3])
	x.name.location <- numeric(length=dim(phylospace.object$segments.to.plot)[1])
	y.name.location <- numeric(length=dim(phylospace.object$segments.to.plot)[1])
	phylospace.object$segments.to.plot$quadrant <- mapply(quadrant, x0=phylospace.object$segments.to.plot$x0, x1=phylospace.object$segments.to.plot$x1, y0=phylospace.object$segments.to.plot$y0, y1=phylospace.object$segments.to.plot$y1)

	##the general form of these logical statements follows x.location[quadrant in question] <- orig.x.location[quadrant in question] +/- (depending on quadrant) label.adjust/sqrt(1 + slope[quadrant in question]^2)
	##and y.location[quadrant in question] <- orig.x.location[quadrant] +/- slope[quadrant]*label.adjust/sqrt(1+slope[quadrant]^2)
		
	##quadrant 1
	x.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant1"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "quadrant1"] + 
		label.adjust/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant1"]^2)
	y.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant1"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "quadrant1"] + 
		(slopes[phylospace.object$segments.to.plot$quadrant == "quadrant1"]*label.adjust)/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant1"]^2)
		
	##quadrant 2. slope is negative, hence the addition sign for the y-values
	x.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant2"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "quadrant2"] + 
		label.adjust/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant2"]^2)
	y.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant2"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "quadrant2"] + 
		(slopes[phylospace.object$segments.to.plot$quadrant == "quadrant2"]*label.adjust)/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant2"]^2)
		
	##quadrant 3
	x.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant3"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "quadrant3"] - 
		label.adjust/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant3"]^2)
	y.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant3"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "quadrant3"] - 
		(slopes[phylospace.object$segments.to.plot$quadrant == "quadrant3"]*label.adjust)/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant3"]^2)
		
	##quadrant 4. slope is negative, hence the subtraction sign for the y-values
	x.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant4"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "quadrant4"] - 
		label.adjust/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant4"]^2)
	y.name.location[phylospace.object$segments.to.plot$quadrant == "quadrant4"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "quadrant4"] - 
		(slopes[phylospace.object$segments.to.plot$quadrant == "quadrant4"]*label.adjust)/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "quadrant4"]^2)

	#straight up line
	x.name.location[phylospace.object$segments.to.plot$quadrant == "straightUp"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "straightUp"]
	y.name.location[phylospace.object$segments.to.plot$quadrant == "straightUp"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "straightUp"] + 
		(slopes[phylospace.object$segments.to.plot$quadrant == "straightUp"]*label.adjust)/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "straightUp"]^2)
		
	#straight down line
	x.name.location[phylospace.object$segments.to.plot$quadrant == "straightDown"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "straightDown"]
	y.name.location[phylospace.object$segments.to.plot$quadrant == "straightDown"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "straightDown"] - 
		(slopes[phylospace.object$segments.to.plot$quadrant == "straightDown"]*label.adjust)/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "straightDown"]^2)

	#straight right line
	x.name.location[phylospace.object$segments.to.plot$quadrant == "straightRight"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "straightRight"] + 
		label.adjust/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "straightRight"]^2)
	y.name.location[phylospace.object$segments.to.plot$quadrant == "straightRight"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "straightRight"]
		
	#straight left line
	x.name.location[phylospace.object$segments.to.plot$quadrant == "straightRight"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "straightRight"] - 
		label.adjust/sqrt(1 + slopes[phylospace.object$segments.to.plot$quadrant == "straightRight"]^2)
	y.name.location[phylospace.object$segments.to.plot$quadrant == "straightRight"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "straightRight"]

	##for the root edge
	x.name.location[phylospace.object$segments.to.plot$quadrant == "stationary"] <- 
		phylospace.object$segments.to.plot[,5][phylospace.object$segments.to.plot$quadrant == "stationary"]
	y.name.location[phylospace.object$segments.to.plot$quadrant == "stationary"] <- 
		phylospace.object$segments.to.plot[,6][phylospace.object$segments.to.plot$quadrant == "stationary"]

	##set up an empty plot. if x and y limits aren't given, just let it figure out the appropriate limits
	if(missing(x.limits) & missing(y.limits))
	{
		plot(phylospace.object$segments.to.plot$y1[phylospace.object$segments.to.plot$descendant <= length(phylospace.object$ape.phylo$tip.label)] ~ 
		phylospace.object$segments.to.plot$x1[phylospace.object$segments.to.plot$descendant <= length(phylospace.object$ape.phylo$tip.label)], 
		xlab=x.label, ylab=y.label, type="n")
	}

	else if(missing(x.limits) | missing(y.limits))
	{
		warning("Need to specify either both x & y limits or no limits at all")
	}
	else
	{
		plot(phylospace.object$segments.to.plot$y1[phylospace.object$segments.to.plot$descendant <= length(phylospace.object$ape.phylo$tip.label)] ~ 
		phylospace.object$segments.to.plot$x1[phylospace.object$segments.to.plot$descendant <= length(phylospace.object$ape.phylo$tip.label)], 
		xlab=x.label, ylab=y.label, type="n", xlim=x.limits, ylim=y.limits)
	}

	if(missing(climate.points))
	{
		print("No climate points supplied")
	}
	else
	{
		points(climate.points, pch=20, cex=0.6, col="gray80")
	}

	##this nested for loop breaks considers an imaginary line segment between either x0 and x1 or y0 and y1, 
	##breaks it into 99 new points, adds the initial x or y point to the start of that vector of 
	##new points, and does this for each row of the segments to plot. this is just to get the appropriate
	##color scale between nodes. then it plots each line segment according to that color scale
	##because we added in the fake branch that just goes from blue to blue, this seems to cause this to crash
	##need to take that line out for this part

	temp.segments.to.plot <- phylospace.object$segments.to.plot[phylospace.object$segments.to.plot$descendant!=length(phylospace.object$ape.phylo$tip.label)+1, ]

	for(i in 1:dim(temp.segments.to.plot)[1])
	{
		breaks = 99
		x <- c()
		y <- c()
		temp.x <- c()
		temp.y <- c()
			
		for(j in 1:breaks)
		{
			temp.x[j] <- (temp.segments.to.plot[i,5]-temp.segments.to.plot[i,3])/breaks*j+temp.segments.to.plot[i,3]
			x <- c(temp.segments.to.plot[i,3],temp.x)
			temp.y[j] <- (temp.segments.to.plot[i,6]-temp.segments.to.plot[i,4])/breaks*j+temp.segments.to.plot[i,4]
			y <- c(temp.segments.to.plot[i,4],temp.y)
		}
		
		temp.quad <- quadrant(x0=temp.x[1], x1=temp.x[breaks], y0=temp.y[1], y1=temp.y[breaks])
		
		##calculate the euclidean distance between each successive x,y point
		
		xy <- data.frame(x,y)

		xydist <- as.matrix(dist(xy))[,1]

		##the way this works is if you want it to go from just blue (0,0,1) to cyan (0,1,1) you code it c(0,0),c(0,1),c(1,1) and if you wanted it to go from cyan to blue to cyan it would be c(0,0,0),c(1,0,1),c(1,1,1), etc. the first place in each argument refers to color 1, the second place color 2, etc.

		color.scale.lines(x,y,c(temp.segments.to.plot$from.r[i],temp.segments.to.plot$to.r[i]),c(temp.segments.to.plot$from.g[i],temp.segments.to.plot$to.g[i]),c(temp.segments.to.plot$from.b[i],temp.segments.to.plot$to.b[i]),colvar=xydist,...)
	}

	##add the species' points
	points(phylospace.object$segments.to.plot$y1[phylospace.object$segments.to.plot$descendant <= length(phylospace.object$ape.phylo$tip.label)] ~ 
		phylospace.object$segments.to.plot$x1[phylospace.object$segments.to.plot$descendant <= length(phylospace.object$ape.phylo$tip.label)], pch=20, col="red")

	##add text to the points, with offset equal to label.adjust argument. IMPORTANTLY, this will work best if x & y are on similar scales. need to revise this script so that it doesn'matter. to query graphical parameters after calling the blank plot: par("usr")
	if(species.labels == TRUE & node.labels == FALSE)
	{
		text(x=x.name.location[1:length(phylospace.object$name.vector)], y=y.name.location[1:length(phylospace.object$name.vector)], labels=phylospace.object$name.vector, cex=0.7)
	}

	else if(species.labels == TRUE & node.labels == TRUE)
	{
		text(x=x.name.location[1:length(phylospace.object$name.vector)], y=y.name.location[1:length(phylospace.object$name.vector)], labels=phylospace.object$name.vector, cex=0.7)
		text(x=x.name.location[(length(phylospace.object$name.vector)+1):(length(phylospace.object$name.vector)+length(phylospace.object$node.vector))], y=y.name.location[(length(phylospace.object$name.vector)+1):(length(phylospace.object$name.vector)+length(phylospace.object$node.vector))], labels=phylospace.object$node.vector, cex=0.7)
	}
	else if(species.labels == FALSE & node.labels == TRUE)
	{
		text(x=x.name.location[(length(phylospace.object$name.vector)+1):(length(phylospace.object$name.vector)+length(phylospace.object$node.vector))], y=y.name.location[(length(phylospace.object$name.vector)+1):(length(phylospace.object$name.vector)+length(phylospace.object$node.vector))], labels=phylospace.object$node.vector, cex=0.7)
	}	
}
