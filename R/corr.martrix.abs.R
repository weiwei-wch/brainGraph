#corr.matrix.my
corr.matrix.weighted<-function (resids, densities, thresholds = NULL, what = c("resids", 
                                                                         "raw"), exclude.reg = NULL, type = c("pearson", 
                                                                                                              "spearman"), rand = FALSE) 
{
  Group <- Study.ID <- NULL
  what <- match.arg(what)
  type <- match.arg(type)
  stopifnot(inherits(resids, "brainGraph_resids"))
  if (isTRUE(rand)) {
    res.all <- as.matrix(resids$resids.all[, !c("Study.ID", 
                                                "Group")])
    corrs <- rcorr(res.all)
    r <- abs(corrs$r)
    N <- ncol(r)
    emax <- N * (N - 1)/2
    thresholds <- sort(r[lower.tri(r)])[emax - densities * 
                                          emax]
    r.thresh <- array(0, dim = c(dim(r), length(thresholds)), 
                      dimnames = list(rownames(r), colnames(r)))
    for (i in seq_along(thresholds)) r.thresh[, , i] <- ifelse(r > 
                                                                 thresholds[i], 1, 0)
    return(list(list(R = r, r.thresh = r.thresh)))
  }
  groups <- resids$groups
  out <- sapply(groups, function(x) NULL)
  if (what == "resids") {
    res.all <- resids$resids.all[, !"Study.ID"]
  }
  else if (what == "raw") {
    res.all <- dcast(resids$all.dat.tidy, "Study.ID + Group ~ region")
    setkey(res.all, Group, Study.ID)
    res.all <- res.all[, !"Study.ID"]
  }
  if (!is.null(exclude.reg)) 
    res.all <- res.all[, -exclude.reg, with = FALSE]
  for (g in groups) {
    corrs <- rcorr(as.matrix(res.all[g, !"Group"]), 
                   type = type)
    r <- abs(corrs$r)
    p <- corrs$P
    if (hasArg("densities")) {
      N <- ncol(r)
      emax <- N * (N - 1)/2
      thresholds <- sort(r[lower.tri(r)])[emax - densities * 
                                            emax]
    }
    r.thresh <- array(0, dim = c(dim(r), length(thresholds)), 
                      dimnames = list(rownames(r), colnames(r)))
    for (i in seq_along(thresholds)) r.thresh[, , i] <- ifelse(r > 
                                                                 thresholds[i], 1, 0)
    out[[g]] <- list(R = r, P = p, r.thresh = r.thresh, thresholds = thresholds, 
                     what = what, exclude.reg = exclude.reg, type = type)
    if (hasArg("densities")) 
      out[[g]] <- c(out[[g]], list(densities = densities))
  }
  return(out)
}