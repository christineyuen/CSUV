#' Graphical illustration of selection uncertainty
#' @export csuv.ci
#' @param csuv.fit fitted results from CSUV::csuv()
#' @param level significance level
#' @param type todo
#' @examples
#' \dontrun{
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' mod.0 = csuv(X, Y, intercept = FALSE, q = 0, method.names = NULL)
#' print(csuv.ci(mod.0, level = 0.1, type = "original"))
#' }
csuv.ci<-function(csuv.fit, level, type = "original"){
  mods = NULL
  B = csuv.fit$param$B
  mods = NULL
  if(type == "PR"){
    mods = csuv.fit$PR.mod.collection[,-1] # first is intercept
  } else if(type == "original"){
    mods = csuv.fit$mod.collection[,-1] # first is intercept
  } else if(type == "conditional"){
    mods = csuv.fit$mod.collection[,-1]
    mods[which(mods==0)] = NA
  } else if(type == "conditional.1"){
    mods = csuv.fit$mod.collection[,-1]
    pos.i = which(colMeans(mods>0)>colMeans(mods<0))
    neg.i = which(colMeans(mods>0)<colMeans(mods<0))
    if(length(pos.i)){ mods[,pos.i][which(mods[,pos.i]<0)] = NA }
    if(length(neg.i)){ mods[,neg.i][which(mods[,neg.i]>0)] = NA }
    mods[which(mods==0)] = NA
  }
  p = ncol(mods)

  li = ui = NULL
  if(type == "conditional"){
    m = nrow(mods)
    li = sapply(1:p, function(i){
      stats::quantile(mods[which(mods[,i]!=0),i], level/2, na.rm = TRUE)
      # num.non.zero = sum(!is.na(mods[,i]))
      # ci.margin = num.non.zero-(1-level)*m
      # if(ci.margin<=1){
      #   return (0)
      # } else{
      #   j = floor(ci.margin/2)
      #   return (sort(mods[,i],partial=j)[j])
      # }
    })
    ui = sapply(1:p, function(i){
      stats::quantile(mods[which(mods[,i]!=0),i], 1-level/2, na.rm = TRUE)
      # num.non.zero = sum(!is.na(mods[,i]))
      # ci.margin = num.non.zero-(1-level)*m
      # if(ci.margin<=1){
      #   return (0)
      # } else{
      #   j = floor(ci.margin/2)
      #   return (sort(mods[,i],decreasing = TRUE)[j])
      # }
    })
  } else{
    li = sapply(1:p, function(i)
      stats::quantile(mods[,i], level/2, na.rm = TRUE))
    ui = sapply(1:p, function(i)
      stats::quantile(mods[,i], 1-level/2, na.rm = TRUE))
  }

  names(li) = colnames(mods)
  li[which(is.na(li))] = 0
  ui[which(is.na(ui))] = 0
  return (rbind(li = li,
                ui = ui))
}

# ci.helper<-function(mods, alpha, type, b = NULL){
#   p = ncol(mods)
#
#   li = sapply(1:p, function(i) stats::quantile(mods[,i], alpha/2))
#   ui = sapply(1:p, function(i) stats::quantile(mods[,i], 1-alpha/2))
#   names(li) = colnames(mods)
#
#   if(type == "quantile"){
#     return (rbind(li = li,
#                   ui = ui))
#   } else{
#     if(is.null(b)){
#       b = sapply(1:p, function(i) median(mods[,i]))
#     }
#     basic.li = 2*b-ui
#     basic.ui = 2*b-li
#     if(type == "basic"){
#       return (rbind(li = basic.li,
#                     ui = basic.ui))
#     } else if(type == "basic2"){
#       i = which(b!=0)
#       li[i] = basic.li[i]
#       ui[i] = basic.ui[i]
#       return (rbind(li = li,
#                     ui = ui))
#     } else if(type == "basic3"){
#       i = which(b!=0 & (sign(ui)==sign(li) | sign(ui)==0 | sign(li)==0))
#       li[i] = basic.li[i]
#       ui[i] = basic.ui[i]
#       return (rbind(li = li,
#                     ui = ui))
#     }
#   }
# }

# # here the models does not include intercept
# get.ci.from.mod.collection<-function(mods, level = 0.1, csuv.sel = NULL){
#   ci.interval = HDCI::ci(Beta = NULL, Beta_bootstrap = mods, type = "quantile", a = level)
#   rownames(ci.interval) = c("li", "ui")
#   if(!is.null(csuv.sel)){
#     li = ci.interval["li",]
#     ui = ci.interval["ui",]
#     is.sel = sign(ui)==sign(li) & sign(li)!=0 & csuv.sel!=0
#     ci.interval["ui", which(!is.sel & ui<0)] = 0
#     ci.interval["li", which(!is.sel & li>0)] = 0
#   }
#   return (ci.interval)
# }
#
# ## ==== LPR helper ====
lm.PR<-function (x, y, intercept, est.b, ...)
{
  PR.mod = HDCI::PartRidge(x = x, y = y, intercept = intercept,
                           varset = est.b[-1] != 0, lambda2 = 1/length(y), ...)
  return (c(PR.mod$beta0, PR.mod$beta))
}

get.PR.mod.collection<-function(x, y, intercept, flds, csuv.mod.collection,
                                 variable.freq, num.core = 1){
  # # csuv.mod.collection includes intercept!!!
  # LPR.mod = do.call(rbind, lapply(1:length(flds),
  #                                 function(fld.i){
  #                                   fld = flds[[fld.i]]
  #                                   return (lm.PR(x = x[-fld,], y = y[-fld], intercept = intercept,
  #                                                  est.b = csuv.mod.collection[fld.i,]))
  #                                 }))
  # return (LPR.mod)
  result = NULL
  if(num.core>1){
    cl = makeCluster(min(detectCores()-1, num.core), outfile = "parallel_log.txt")
    clusterExport(cl=cl, list("lm.PR"), envir=environment())
    result = parLapply(cl, 1:length(flds), function(fld.i, param){
      fld = param$flds[[fld.i]]
      return (lm.PR(x = param$x[-fld,], y = param$y[-fld], intercept = param$intercept,
                     est.b = param$csuv.mod.collection[fld.i,]))
    }, param = list(x = x, y = y, intercept = intercept, flds = flds, csuv.mod.collection = csuv.mod.collection))
    stopCluster(cl)
    rm(cl)
  } else{
    result =  lapply(1:length(flds),
                     function(fld.i){
                       fld = flds[[fld.i]]
                       return (lm.PR(x = x[-fld,], y = y[-fld], intercept = intercept,
                                      est.b = csuv.mod.collection[fld.i,]))
                     })
  }
  return (do.call(rbind, result))
}
