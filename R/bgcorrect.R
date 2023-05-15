#
# Spatial Transcriptomics using NanoString GeoMx DSP
#


#'
#' subtract background noise from each sample
#'
#' @param x count matrix
#' @param offset vector of bg noise levels per sample
#' @param min.count minimum count threshold
#' @returns matrix of background-subtracted counts
#'
background_subtract <- function(x, offset=0, min.count=1) {
  x <- sweep(x, 2, offset, FUN="-")
  x <- apply(x, 2, pmax, min.count)
  x
}

#'
#' background correction
#'
#' assumes normally distributed background (negative probes)
#' and gene probes
#'
#' @param x matrix of counts
#' @param bg vector of TRUE/FALSE corresponding to rows of x
#' @returns matrix of background-subtracted counts
#'
bgcorrect_norm <- function(x, bg) {
  a <- stats::sd(x[bg])^2 / stats::sd(x[!bg])^2
  c <- a * mean(x[!bg]) - mean(x[bg])
  y <- (1-a)*x + c
  return(y)
}


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
  # setup
  n <- 2^13
  xlog <- log2(x)

  # select bandwidth
  bwbg <- stats::bw.nrd0(xlog[bg]) * bw.adjust
  bw <- stats::bw.nrd0(xlog[!bg]) * bw.adjust
  bw <- max(bw, bwbg)

  # smooth ecdf functions for gene probes and neg probes
  d <- stats::density(xlog[!bg], bw=bw, n=n)
  dcdf <- cumsum(d$y) / sum(d$y)

  dbg <- stats::density(xlog[bg], bw=bw, n=n)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)

  q <- stats::approx(d$x, dcdf, xout=xlog, ties="ordered")$y
  q <- pmax(bg.quant, q)

  log2noise <- stats::approx(dbgcdf, dbg$x, xout=q, ties="ordered")$y
  noise <- 2^log2noise

  # ensure noise is not greater than (x - 1)
  noise <- pmin(noise, x - 1)
  y <- x - noise
  # ensure ymin == 1
  y <- y - (min(y) - 1)
  # ensure isotonic behavior
  ir <- stats::isoreg(x, y)
  yf <- ir$yf[order(ir$ord)]

  res <- list(
    y = y,
    yf = yf,
    noise = noise,
    dx = d$x,
    dy = d$y,
    dcdf = dcdf,
    dbgx = dbg$x,
    dbgy = dbg$y,
    dbgcdf = dbgcdf,
    bw.adjust = bw.adjust,
    bg.quant = bg.quant
  )
}


bgcorrect_kdeppv <- function(x, bg, bw.adjust=2, bg.quant=0.5) {
  # setup
  eps <- 1e-10
  n <- 2^13
  x <- log2(x)
  xmin <- min(x)
  xmax <- max(x)
  expressed <- !bg & (x > stats::quantile(x[bg], bg.quant))

  # select bandwidth
  bw.gene <- stats::bw.nrd0(x[expressed])
  bw.bg <- stats::bw.nrd0(x[bg])
  bw <- max(bw.gene, bw.bg) * bw.adjust

  d <- stats::density(x[expressed], bw=bw, from=xmin, to=xmax, n=n)
  dbg <- stats::density(x[bg], bw=bw, from=xmin, to=xmax, n=n)

  dcdf <- cumsum(d$y) / sum(d$y)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)



  # convert density to ecdf
  dcdf <- cumsum(d$y) / sum(d$y)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)

  # calc false omission rate (FOR) as FN / PN
  dppv <- (1 - dcdf + eps) / (1 - dcdf + 1 - dbgcdf + eps)
  dppv <- (dppv - min(dppv)) / (max(dppv) - min(dppv))

  # ensure isotonic behavior
  dppv.ir <- stats::isoreg(d$x, dppv)$yf
  # interpolate
  ppv <- stats::approx(d$x, dppv.ir, xout=x, ties="ordered")$y
  # transform
  y <- (2^x) * ppv
  y <- pmax(1, y)

  return(list(
    y = y,
    dx = d$x,
    dcdf = dcdf,
    dbgcdf = dbgcdf,
    dp = dppv.ir,
    bw = bw,
    bw.adjust = bw.adjust,
    bg.quant = bg.quant
  ))
}


bgcorrect_kdefor <- function(x, bg, bw.adjust=2, bg.quant=0.5) {
  # setup
  eps <- 1e-10
  n <- 2^13
  x <- log2(x)
  xmin <- min(x)
  xmax <- max(x)
  expressed <- !bg & (x > stats::quantile(x[bg], bg.quant))

  # select bandwidth
  bw.gene <- stats::bw.nrd0(x[expressed])
  bw.bg <- stats::bw.nrd0(x[bg])
  bw <- max(bw.gene, bw.bg) * bw.adjust

  # use kde to produce a smooth density estimate in log space
  d <- stats::density(x[expressed], bw=bw, from=xmin, to=xmax, n=n)
  dbg <- stats::density(x[bg], bw=bw, from=xmin, to=xmax, n=n)

  # convert density to ecdf
  dcdf <- cumsum(d$y) / sum(d$y)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)

  # calc false omission rate (FOR) as FN / PN
  dp <- dcdf / (dbgcdf + dcdf + eps)
  dp <- (dp - min(dp)) / (max(dp) - min(dp))

  # ensure isotonic behavior
  dp.ir <- stats::isoreg(d$x, dp)$yf

  # interpolate
  p <- stats::approx(d$x, dp.ir, xout=x, ties="ordered")$y
  #f <- splinefun(d$x, dp.ir, method="natural", ties="ordered")
  #p <- f(x)

  # transform
  y <- (2^x) * p
  y <- pmax(1, y)

  return(list(
    y = y,
    dx = d$x,
    dcdf = dcdf,
    dbgcdf = dbgcdf,
    dp = dp.ir,
    bw = bw,
    bw.adjust = bw.adjust,
    bg.quant = bg.quant
  ))
}


#'
#' background correction
#'
#' @param ds list dataset
#' @param method correction method to choose from
#' @param bw.adjust bandwidth adjustment factor
#' @param bg.quant background subtraction quantile
#' @returns matrix of background corrected counts-per-million
#' @export
st_geomx_bgcorrect <- function(ds, method=c("qq", "norm", "bgsub", "kdeppv", "kdefor", "none"),
                               bw.adjust=1, bg.quant=0.5) {
  bg <- ds$meta$bg
  x <- ds$counts

  apply_bgcorrect_qq <- function(x) {
    y <- bgcorrect_qq(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)$y
    return(y)
  }
  apply_bgcorrect_norm <- function(x) {
    y <- bgcorrect_norm(x, bg)
    y <- pmax(y, 1)
    return(y)
  }
  apply_bgcorrect_bgsub <- function(x) {
    xsub <- 2^stats::quantile(log2(x)[bg], bg.quant)
    y <- pmax(1, x - xsub)
    return(y)
  }
  apply_bgcorrect_kdeppv <- function(x) {
    y <- bgcorrect_kdeppv(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)$y
    return(y)
  }
  apply_bgcorrect_kdefor <- function(x) {
    y <- bgcorrect_kdefor(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)$y
    return(y)
  }

  if (method == "none") {
    return(x)
  } else if (method == "qq") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_qq)
  } else if (method == "norm") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_norm)
  } else if (method == "bgsub") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_bgsub)
  } else if (method == "kdeppv") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_kdeppv)
  } else if (method == "kdefor") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_kdefor)
  } else {
    stop("method not found")
  }
  return(x)
}
