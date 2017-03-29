#' Plot Robust Autocovariance and Robust Autocorrelation Functions
#'
#'Plot method for objects of class "robacf".
#'
#' @param x an object of class "robacf".
#' @param type the type of plot to be drawn, default to histogram like vertical lines.
#' @param xlab the x label of the plot.
#' @param ylab the y label of the plot.
#' @param ylim numeric of length 2 giving the y limits for the plot.
#' @param main overall title for the plot.
#' @param max.mfrow positive integer; for multivariate x indicating how many rows and columns of plots should be put on one page, using par(mfrow = c(m,m))(see \code{\link[graphics]{par}}).
#' @param ask logical; if TRUE, the user is asked before a new page is started.
#' @param mar,oma,mgp,xpd,cex.main graphics parameters as in par(*), by default adjusted to use smaller than default margins for multivariate x only.
#' @param verbose logical. Should R report extra information on progress?
#' @param ... graphics parameters to be passed to the plotting routines.
#'
#' @method plot robacf
#'
#' @return None
#' @section Contributions: plot.acf (stats) - R Core
#' @importFrom grDevices dev.flush dev.hold dev.interactive
#' @importFrom graphics abline box frame mtext par title
#' @importFrom utils str
#' @examples
#' robacf(cbind(ldeaths,mdeaths))
plot.robacf <-  function (x, type = "h", xlab = "Lag", ylab = NULL,ylim = NULL, main = NULL, max.mfrow = 6,            ask = Npgs > 1 && dev.interactive(),
            mar = if(nser > 2) c(3,2,2,0.8) else par("mar"),
            oma = if(nser > 2) c(1,1.2,1,1) else par("oma"),
            mgp = if(nser > 2) c(1.5,0.6,0) else par("mgp"),
            xpd = par("xpd"),
            cex.main = if(nser > 2) 1 else par("cex.main"),
            verbose = getOption("verbose"),
            ...)
  {
    if((nser <- ncol(x$lag)) < 1L) stop("x$lag must have at least 1 column")
    if (is.null(ylab))
      ylab <- switch(x$type,
                     correlation = "Robust ACF",
                     covariance = "Robust ACF (cov)")
    if (is.null(snames <- x$snames))
      snames <- paste("Series ", if (nser == 1L) x$series else 1L:nser)
    clim0 <- c(0,0)
    Npgs <- 1L ## we will do [ Npgs x Npgs ] pages !
    nr <- nser
    if(nser > 1L) { ## at most m x m (m := max.mfrow)  panels per page
      sn.abbr <- if(nser > 2L) abbreviate(snames) else snames

      if(nser > max.mfrow) {
        ##  We need more than one page: The plots are laid out
        ##  such that we can manually paste the paper pages and get a
        ##  nice square layout with diagonal !
        ## NB: The same applies to pairs() where we'd want several pages
        Npgs <- ceiling(nser / max.mfrow)
        nr <- ceiling(nser / Npgs)  # <= max.mfrow
      }
      opar <- par(mfrow = rep(nr, 2L), mar = mar, oma = oma, mgp = mgp,
                  ask = ask, xpd = xpd, cex.main = cex.main)
      on.exit(par(opar))
      if(verbose) { # FIXME: message() can be suppressed but not str()
        message("par(*) : ", appendLF=FALSE, domain = NA)
        str(par("mfrow","cex", "cex.main","cex.axis","cex.lab","cex.sub"))
      }
    }

    if (is.null(ylim)) {
      ## Calculate a common scale
      ylim <- range(x$acf[, 1L:nser, 1L:nser], na.rm = TRUE)
      ylim <- range(c(-clim0, clim0, ylim))
    }

    for (I in 1L:Npgs) for (J in 1L:Npgs) {
      dev.hold()
      ## Page [ I , J ] : Now do   nr x nr  'panels' on this page
      iind <- (I-1)*nr + 1L:nr
      jind <- (J-1)*nr + 1L:nr
      if(verbose)
        message(gettextf("Page [%d,%d]: i =%s; j =%s", I, J, paste(iind,collapse=","), paste(jind,collapse=",")), domain = NA)
      for (i in iind) for (j in jind)
        if(max(i,j) > nser) {
          frame(); box(col = "light gray")
          ## the above is EXTREMELY UGLY; should have a version
          ## of frame() that really does advance a frame !!
        }
      else {
        clim <- clim0
        plot(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab,
             ylab = if(j==1) ylab else "", ylim = ylim, ...)
        abline(h = 0)
        title(if (!is.null(main)) main else
          if (i == j) snames[i]
          else paste(sn.abbr[i], "&", sn.abbr[j]),
          line = if(nser > 2) 1 else 2)
      }
      if(Npgs > 1) {                  # label the page
        mtext(paste("[",I,",",J,"]"), side=1, line = -0.2, adj=1,
              col = "dark gray", cex = 1, outer = TRUE)
      }
      dev.flush()
    }
    invisible()
  }
