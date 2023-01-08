mkBlistMV<-function(x,frloc, cov_re, drop.unused.levels=TRUE) {
  
  frloc <- factorize(x,frloc)
  ## try to evaluate grouping factor within model frame ...
  if (is.null(ff <- tryCatch(eval(substitute(makeFac(fac),
                                             list(fac = x[[3]])), frloc),
                             error=function(e) NULL)))
    stop("couldn't evaluate grouping factor ",
         deparse(x[[3]])," within model frame:",
         " try adding grouping factor to data ",
         "frame explicitly if possible",call.=FALSE)
  if (all(is.na(ff)))
    stop("Invalid grouping factor specification, ",
         deparse(x[[3]]),call.=FALSE)
  ## NB: *also* silently drops <NA> levels - and mkReTrms() and hence
  ##     predict.merMod() have relied on that property  :
  if (drop.unused.levels) ff <- factor(ff, exclude=NA)
  nl <- length(levels(ff))
  ## this section implements eq. 6 of the JSS lmer paper
  ## model matrix based on LHS of random effect term (X_i)
  ##    x[[2]] is the LHS (terms) of the a|b formula
  mm <- model.matrix(eval(substitute( ~ foo, list(foo = x[[2]]))), frloc)
  ## nc <- ncol(mm)
  ## nseq <- seq_len(nc)
  ## this is J^T (see p. 9 of JSS lmer paper)
  ## construct indicator matrix for groups by observations
  ## use fac2sparse() rather than as() to allow *not* dropping
  ## unused levels where desired
  sm <- fac2sparse(ff, to = "d",
                   drop.unused.levels = drop.unused.levels)
  ## looks like we don't have to filter NAs explicitly any more ...
  ## sm <- as(ff,"sparseMatrix")
  ## sm <- KhatriRao(sm[,!is.na(ff),drop=FALSE],t(mm[!is.na(ff),,drop=FALSE]))
  
  #browser()
  #### START source code modified by Eardi ----
  # Commented version gives annoying message
  #if (!is.null(cov_mat <- tryCatch(eval(get(as.character(x[[3]]),cov_re)),
  #                           error=function(e) NULL)))
  #{
  #  sm = chol(cov_mat)%*%sm;
  #  message("Covariance structure of ",as.character(x[[3]]), " succesfully applied") 
  #}
  #browser()
  cov_mat <- NULL
  try(cov_mat <- eval(get(as.character(x[[3]]),cov_re)), silent = TRUE)
  if(!is.null(cov_mat))
  {
    sm_labels = rownames(sm);
    sm = chol(cov_mat[sm_labels,sm_labels])%*%sm;
    message("Covariance structure of ",as.character(x[[3]]), " successfully applied") 
  }
  
  
  #### END source code modified by Eardi ----
  
  sm <- KhatriRao(sm, t(mm))
  dimnames(sm) <- list(
    rep(levels(ff),each=ncol(mm)),
    rownames(mm))

  # Alternative, apply premultiplication post KhatriRao
  # for (varbl in 1:ncol(mm))
  # {
  #   ind_rows = seq(from = varbl, to = length(levels(ff))*ncol(mm), by = ncol(mm));
  #   # This is t(t(chol(cov_re$id_genetic))) as A = L L^T, we need L^T (given by chol)
  #   sm[ind_rows,] =  chol(cov_re$id_genetic)%*%sm[ind_rows,];
  # }
  
  list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}