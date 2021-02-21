
#' basic one-species continuous-time population model
#'
#' show plots of demographic parameters and/or time dynamics
#' 
#' @param timeMax maximum time for dynamics plot [time]
#' @param steps number of steps for dynamics plot
#' @param popMax maximum population for demographic parameter plot [number]
#' @param popSteps number of steps for demographic parameter plot
#' @param b0 \emph{per capita} birth rate at zero density [1/time]
#' @param bDD characteristic density for exponential decrease in per capita birth rate with increasing population density [number]
#' @param bAllee characteristic scale for Allee effect in birth rate [number]
#' @param d0 \emph{per capita} death rate at zero density [1/time]
#' @param dDD characteristic density for exponential increase in \emph{per capita} death rate with increasing population density [number]
#' @param dAllee characteristic scale for Allee effect in death rate [number]
#' @param N0 initial population size for dynamics plot [number]: if \code{N0}=0, simulations of time dynamics will not be run nor plotted
#' @param reportPcTotal whether to plot \emph{per capita} rates ("p"), total rates ("t"), both ("b"), or neither ("n")
#' @param reportDiff whether to plot the overall growth rate (birth-death) rather than birth and death separately
#' @param fontSize scaled font size
#' @param legendSize scaled legend size (base plots only)
#' @param title plot title
#' @param tlab label for time axis
#' @param plab label for population size axis
#' @param printPlots print plots (alternatively, return a list of plots)?
#' @param plotType "ggplot2" or "base"
#' @param logScale make y-axis logarithmic (for time dynamics only)?
#' @param \dots additional arguments passed down to plotting functions, including \code{bw} for black-and-white plotting
#' @details The basic model considered here is an exponential-density-dependence model, i.e.
#' \deqn{\frac{dN}{dt} = N (b_0 \exp(-N/b_{DD}) - d_0 \exp(N/d_{DD}))}
#' (more details to be added later)
#' @examples
#' bd()     ## basic plot
#' ## set initial population size to something other than 0
#' ## in order to get a dynamics plot
#' bd(N0=1) 
#' @export
bd <- function(N0=NULL, MaxTime=20, steps=100, popMax=100, b0=1, bDD=NULL, bAllee=NULL, d0=0.5, dDD=NULL, dAllee=NULL, reportPcTotal="b", popSteps=100, fontSize=1, legendSize=1, title="", tlab = "Time (years)", plab="Population size", showLog=TRUE, growMax=NULL){
	pop <- 1:popSteps*(popMax/popSteps)
	b <- rfun(b0, bDD, bAllee, pop, TRUE)
	d <- rfun(d0, dDD, dAllee, pop, FALSE)
	if (reportPcTotal == "p" | reportPcTotal == "b")
		bdplots(pop, b, d, reportTotal=FALSE, title, fontSize=fontSize, 
			legendSize=legendSize, plab=plab)
	if (reportPcTotal == "t" | reportPcTotal == "b")
		bdplots(pop, b, d, reportTotal=TRUE, title, fontSize=fontSize, 
			legendSize=legendSize, plab=plab)
	if(!is.null(N0)){
		sim <- bdsim(N0, MaxTime, steps, b0, bDD, bAllee, d0, dDD, dAllee)
		if(is.null(growMax)) growMax <- max(sim$N)
		# print(data.frame(sim$time, sim$N))
		plot(sim$time, sim$N,
			cex.lab=fontSize, cex.axis=fontSize,
			main=title, xlab = tlab, ylab = "Population",
			type = "l", ylim = c(0, growMax)
		)
		if (showLog){
			plot(sim$time, sim$N,
				cex.lab=fontSize, cex.axis=fontSize,
				main=title, xlab = tlab, ylab = "Population",
				type = "l", ylim = c(1, growMax), log="y"
			)
		}
	}
}

bdplots = function(pop, b, d, reportTotal=FALSE, title="Birth-death plot", fontSize, legendSize, plab){
	ylab <- "Per capita rate (1/t)"
	lpos <- "topright"
	if(reportTotal) {
		ylab <- "Total rate (pop/t)"
		lpos <- "bottomright"
		b <- b*pop
		d <- d*pop
	}
	par(cex=1.6)
	respPlot(pop, b, d, lpos, ylab, title, fontSize=fontSize,
	legendSize=legendSize, plab=plab)
}

respPlot <- function(pop, b, d, lpos, ylab, plab="Population size", title, logscale=FALSE, legendSize=1, fontSize=1){
	ymin = ifelse(logscale, min(c(b, d)), 0)
	ymax = max(c(b,min(d)))
	logPar <- ifelse(logscale, "y", "")

	plot(pop, b,
		cex.axis = fontSize, cex.lab=fontSize,
		ylim = c(ymin, ymax),
		xlab = plab,
		ylab = ylab,
		type = "l", lwd=2, col="blue", main=title, log=logPar
	)
	lines(pop, d, lty=2, lwd=2)
	legend(lpos, cex=legendSize,
		legend = c("Birth rate", "Death rate"),
		col = c("blue", "black"),
		lty = c(1, 2)
	)
}


## 2019 Mar 09 (Sat)
## Probably remove divOffset from the numerator? It's a divOffset!
rfun <- function(r0, DD, Allee, pop, birth=TRUE, divOffset=1/2, mmax=1000){
	mult <- 1 + 0*pop
	if (!is.null(DD)) mult <- mult*exp(pop/DD)
	if (!is.null(Allee))
		mult <- mult*exp((Allee+divOffset)/(pop+divOffset))
	mult <- mmax*mult/(mmax+mult)
	if (birth) {mult <- 1/mult}
	return(r0*mult)
}

ndot = function(time, vars, parms){
	ndot <- with(as.list(c(vars, parms)),
		rfun(b0, bDD, bAllee, exp(n), TRUE) 
			- rfun(d0, dDD, dAllee, exp(n), FALSE) 
	)
	list(c(ndot))
}

popSim = function (N0, MaxTime, steps, parms){
	sim <- as.data.frame(lsoda(
		y = c(n=log(N0)),
		times = (0:steps)*MaxTime/steps,
		func = ndot,
		parms
	))
	sim$N <- exp(sim$n)
	return(sim)
}


bdsim <- function(N0=1, MaxTime=20, steps=100, b0=1, bDD=NULL, bAllee=NULL, d0=0.5, dDD=NULL, dAllee=NULL){

	parms = list(b0=b0, bDD=bDD, bAllee=bAllee, d0=d0, dDD=dDD, dAllee=dAllee)
	return(popSim(N0, MaxTime, steps, parms))
}

args(bd)


