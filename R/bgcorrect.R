#
# Spatial Transcriptomics using NanoString GeoMx DSP
#


#'
#' background correction using quantile-quantile approach
#'
#' @param x matrix of counts
#' @param bg vector of TRUE/FALSE corresponding to rows of x
#' @param bw.adjust bandwidth adjustment factor for density function
#' @param bg.quant subtract background quantile from all genes
#' @returns list where element 'y' contains corrected counts
#'
bgcorrect_qq <- function(x, bg, bw.adjust=1, bg.quant=0.5) {
  # subtract constant background noise
  bg.loq <- stats::quantile(x[bg], bg.quant)
  x <- pmax(0, x - bg.loq)

  # setup
  n <- 2^11
  xlog <- log2(x+1)

  # select bandwidth
  bwbg <- stats::bw.nrd0(xlog[bg]) * bw.adjust
  bw <- stats::bw.nrd0(xlog[!bg]) * bw.adjust
  bw <- max(bw, bwbg)

  xmin <- 0
  xmax <- max(xlog[!bg])
  xbgmax <- max(xlog[bg])

  # smooth ecdf functions for gene probes and neg probes
  d <- stats::density(xlog[!bg], bw=bw, n=n, from=0, to=xmax)
  dcdf <- cumsum(d$y) / sum(d$y)

  dbg <- stats::density(xlog[bg], bw=bw, n=n, from=0, to=xbgmax)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)

  # map gene probes to quantiles
  q <- stats::approx(d$x, dcdf, xout=xlog, ties="ordered")$y
  # map quantiles to background probes
  log2noise <- stats::approx(dbgcdf, dbg$x, yleft=0, yright=xbgmax, rule=2,
                             xout=q, ties="ordered")$y
  # convert from log2 to linear scale and round up to discrete counts
  noise <- ceiling(2^log2noise - 1)
  #noise <- ceiling(2^log2noise - 1)
  # ensure noise is not greater than x
  noise <- pmin(noise, x)
  y <- x - noise
  # TODO: do we need to ensure isotonic behavior?
  # ir <- stats::isoreg(x, y)
  # yf <- ir$yf[order(ir$ord)]
  # if (any(y != yf)) { print("PROBLEM") }

  res <- list(
    x = x,
    bg = bg,
    bw.adjust = bw.adjust,
    bg.quant = bg.quant,
    bg.loq = bg.loq,
    bw = bw,
    y = y,
    noise = noise,
    dx = d$x,
    dy = d$y,
    dcdf = dcdf,
    dbgx = dbg$x,
    dbgy = dbg$y,
    dbgcdf = dbgcdf
  )
}


#'
#' Background correction
#'
#' @param ds stgeomx dataset
#' @param method correction method to choose from
#' @param bw.adjust bandwidth adjustment factor
#' @param bg.quant background subtraction quantile
#' @returns stgeomx dataset
#' @export
#'
#' @examples
#' data(example_ds, package= "stgeomx")
#' bgcorrect(example_ds, "qq")
#' bgcorrect(example_ds, "bgsub")
#' bgcorrect(example_ds, "none")
#'
bgcorrect <- function(ds, method=c("qq", "bgsub", "none"),
                      bw.adjust=1, bg.quant=0.5) {
  # check method
  method <- match.arg(method)
  bg <- ds$meta$bg
  x <- ds$counts

  apply_bgcorrect_bgsub <- function(x) {
    xsub <- ceiling(stats::quantile(x[bg], bg.quant))
    y <- pmax(0, x - xsub)
    return(y)
  }
  apply_bgcorrect_qq <- function(x) {
    res <- bgcorrect_qq(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)
    return(res$y)
  }

  f <- switch(method,
              qq = apply_bgcorrect_qq,
              bgsub = apply_bgcorrect_bgsub,
              none = function(x) return(x))
  x <- apply(x, MARGIN=2, FUN=f)
  ds$x <- x
  return(ds)
}

