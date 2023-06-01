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
  # convert from log2 to linear scale and round to discrete counts
  noise <- round(2^log2noise - 1)
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
#    y = yf,
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
#' @param ds list dataset
#' @param method correction method to choose from
#' @param bw.adjust bandwidth adjustment factor
#' @param bg.quant background subtraction quantile
#' @returns matrix of background corrected counts-per-million
#' @export
bgcorrect <- function(ds, method=c("qq", "bgsub", "none"),
                      bw.adjust=1, bg.quant=0.5) {

  # check method
  method <- match.arg(method)
  bg <- ds$meta$bg
  x <- ds$counts

  apply_bgcorrect_bgsub <- function(x) {
    xsub <- stats::quantile(x[bg], bg.quant)
    y <- pmax(1, x - xsub)
    return(y)
  }
  apply_bgcorrect_qq <- function(x) {
    y <- bgcorrect_qq(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)$y
    return(y)
  }

  f <- switch(method,
              qq = apply_bgcorrect_qq,
              bgsub = apply_bgcorrect_bgsub,
              none = function(x) return(x))
  x <- apply(x, MARGIN=2, FUN=f)
  return(x)
}


#'
#' subtract constant background noise from each sample
#'
#' @param x count matrix
#' @param offset vector of bg noise levels per sample
#' @param min.count minimum count threshold
#' @returns matrix of background-subtracted counts
#'
obsolete_bgsub <- function(x, offset=0, min.count=1) {
  x <- sweep(x, 2, offset, FUN="-")
  x <- apply(x, 2, pmax, min.count)
  x
}

#'
#' background correction method assuming normal dist
#'
#' assumes normally distributed background (negative probes)
#' and gene probes
#'
#' @param x matrix of counts
#' @param bg vector of TRUE/FALSE corresponding to rows of x
#' @returns matrix of background-subtracted counts
#'
obsolete_bgcorrect_norm <- function(x, bg) {
  a <- stats::sd(x[bg])^2 / stats::sd(x[!bg])^2
  c <- a * mean(x[!bg]) - mean(x[bg])
  y <- (1-a)*x + c
  return(y)
}

