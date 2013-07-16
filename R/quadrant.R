quadrant <- function(x0, x1, y0, y1)
{
	if(x0 < x1 & y0 < y1)
	{
		quad <- "quadrant1"
	}
	else if(x0 < x1 & y0 > y1)
	{
		quad <- "quadrant2"
	}
	else if(x0 > x1 & y0 > y1)
	{
		quad <- "quadrant3"
	}
	else if(x0 > x1 & y0 < y1)
	{
		quad <- "quadrant4"
	}
	else if(x0 == x1 & y0 < y1)
	{
		quad <- "straightUp"
	}
	else if(x0 == x1 & y0 > y1)
	{
		quad <- "straightDown"
	}
	else if(x0 < x1 & y0 == y1)
	{
		quad <- "straightRight"
	}
	else if(x0 > x1 & y0 == y1)
	{
		quad <- "straightLeft"
	}
	else if(x0 == x1 & y0 == y1)
	{
		quad <- "stationary"
	}
	return(quad)
}
