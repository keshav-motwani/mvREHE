# https://rdrr.io/cran/bpnreg/src/R/lme4code.R
safeDeparse <- function(x, collapse = " "){

  paste(deparse(x, 500L), collapse = collapse)

}

mkReTrmsMV <-function(bars, fr, cov_re, drop.unused.levels=TRUE)
{
  if (!length(bars))
    stop("No random effects terms specified in formula",call.=FALSE)
  stopifnot(is.list(bars), vapply(bars, is.language, NA),
            inherits(fr, "data.frame"))
  names(bars) <- barnames(bars)
  term.names <- vapply(bars, safeDeparse, "")
  ## get component blocks
  blist <- lapply(bars, mkBlistMV, fr, cov_re, drop.unused.levels)
  nl <- vapply(blist, `[[`, 0L, "nl")   # no. of levels per term
  # (in lmer jss:  \ell_i)

  ## order terms stably by decreasing number of levels in the factor
  if (any(diff(nl) > 0)) {
    ord <- rev(order(nl))
    blist      <- blist     [ord]
    nl         <- nl        [ord]
    term.names <- term.names[ord]
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rbind, Ztlist)  ## eq. 7, JSS lmer paper
  names(Ztlist) <- term.names
  q <- nrow(Zt)

  ## Create and install Lambdat, Lind, etc.  This must be done after
  ## any potential reordering of the terms.
  cnms <- lapply(blist, `[[`, "cnms")   # list of column names of the
  # model matrix per term
  nc <- lengths(cnms)                   # no. of columns per term
  # (in lmer jss:  p_i)
  nth <- as.integer((nc * (nc+1))/2)    # no. of parameters per term
  # (in lmer jss:  ??)
  nb <- nc * nl                         # no. of random effects per term
  # (in lmer jss:  q_i)
  ## eq. 5, JSS lmer paper
  if (sum(nb) != q) {
    stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
                 sum(nb),q))
  }
  boff <- cumsum(c(0L, nb))             # offsets into b
  thoff <- cumsum(c(0L, nth))           # offsets into theta
  ### FIXME: should this be done with cBind and avoid the transpose
  ### operator?  In other words should Lambdat be generated directly
  ### instead of generating Lambda first then transposing?
  Lambdat <-
    t(do.call(sparseMatrix,
              do.call(rbind,
                      lapply(seq_along(blist), function(i)
                      {
                        mm <- matrix(seq_len(nb[i]), ncol = nc[i],
                                     byrow = TRUE)
                        dd <- diag(nc[i])
                        ltri <- lower.tri(dd, diag = TRUE)
                        ii <- row(dd)[ltri]
                        jj <- col(dd)[ltri]
                        ## unused: dd[cbind(ii, jj)] <- seq_along(ii)
                        data.frame(i = as.vector(mm[, ii]) + boff[i],
                                   j = as.vector(mm[, jj]) + boff[i],
                                   x = as.double(rep.int(seq_along(ii),
                                                         rep.int(nl[i], length(ii))) +
                                                   thoff[i]))
                      }))))
  thet <- numeric(sum(nth))
  ll <- list(Zt = drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x),
             Gp = unname(c(0L, cumsum(nb))))
  ## lower bounds on theta elements are 0 if on diagonal, else -Inf
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower) # initial values of theta are 0 off-diagonal, 1 on
  Lambdat@x[] <- ll$theta[ll$Lind]  # initialize elements of Lambdat
  ll$Lambdat <- Lambdat
  # massage the factor list
  fl <- lapply(blist, `[[`, "ff")
  # check for repeated factors
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else asgn <- seq_along(fl)
  names(fl) <- ufn
  ## DON'T need fl to be a data.frame ...
  ## fl <- do.call(data.frame, c(fl, check.names = FALSE))
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll
}
