lFormulaMV <- function(formula, data=NULL, cov_re = NULL, REML = TRUE,
         subset, weights, na.action, offset, contrasts = NULL,
         control=lmerControl(), ...)
{
  control <- control$checkControl ## this is all we really need
  mf <- mc <- match.call()
  
  #browser()
  
  ignoreArgs <- c("start","verbose","devFunOnly","control")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(checkArgs, c(list("lmer"), l...))
  if (!is.null(list(...)[["family"]])) {
    ## lmer(...,family=...); warning issued within checkArgs
    mc[[1]] <- quote(lme4::glFormula)
    if (missing(control)) mc[["control"]] <- glmerControl()
    return(eval(mc, parent.frame()))
  }
  
  cstr <- "check.formula.LHS"
  checkCtrlLevels(cstr,control[[cstr]])
  denv <- checkFormulaData(formula, data,
                           checkLHS = control$check.formula.LHS == "stop")
  #mc$formula <- formula <- as.formula(formula,env=denv) ## substitute evaluated call
  formula <- as.formula(formula, env=denv)
  ## as.formula ONLY sets environment if not already explicitly set ...
  ## ?? environment(formula) <- denv
  # get rid of || terms so update() works as expected
  RHSForm(formula) <- expandDoubleVerts(RHSForm(formula))
  mc$formula <- formula
  
  ## (DRY! copied from glFormula)
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  fr.form <- subbars(formula) # substitute "|" by "+"
  environment(fr.form) <- environment(formula)
  ## model.frame.default looks for these objects in the environment
  ## of the *formula* (see 'extras', which is anything passed in ...),
  ## so they have to be put there ...
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x=.(i)))))
      assign(i,get(i,parent.frame()),environment(fr.form))
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  ## convert character vectors to factor (defensive)
  fr <- factorize(fr.form, fr, char.only=TRUE)
  ## store full, original formula & offset
  attr(fr,"formula") <- formula
  attr(fr,"offset") <- mf$offset
  n <- nrow(fr)
  ## random effects and terms modules
  reTrms <- mkReTrmsMV(findbars(RHSForm(formula)), fr, cov_re)
  wmsgNlev <- checkNlevels(reTrms$flist, n=n, control, allow.n=TRUE)
  wmsgZdims <- checkZdims(reTrms$Ztlist, n=n, control, allow.n=TRUE)
  if (anyNA(reTrms$Zt)) {
    stop("NA in Z (random-effects model matrix): ",
         "please use ",
         shQuote("na.action='na.omit'"),
         " or ",
         shQuote("na.action='na.exclude'"))
  }
  wmsgZrank <- checkZrank(reTrms$Zt, n=n, control, nonSmall = 1e6)
  
  ## fixed-effects model matrix X - remove random effect parts from formula:
  fixedform <- formula
  RHSForm(fixedform) <- nobars(RHSForm(fixedform))
  mf$formula <- fixedform
  ## re-evaluate model frame to extract predvars component
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr,"terms"), "predvars.fixed") <-
    attr(attr(fixedfr,"terms"), "predvars")
  ## so we don't have to fart around retrieving which vars we need
  ##  in model.frame(.,fixed.only=TRUE)
  attr(attr(fr,"terms"), "varnames.fixed") <- names(fixedfr)
  
  ## ran-effects model frame (for predvars)
  ## important to COPY formula (and its environment)?
  ranform <- formula
  RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr,"terms"), "predvars.random") <-
    attr(terms(ranfr), "predvars")
  
  ## FIXME: shouldn't we have this already in the full-frame predvars?
  X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
  ## backward compatibility (keep no longer than ~2015):
  if(is.null(rankX.chk <- control[["check.rankX"]]))
    rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
  X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
  if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
  X <- checkScaleX(X, kind=scaleX.chk)
  
  list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula,
       wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
}