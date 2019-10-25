###########################################################
## This file includes functions related to plotting for CSUV
## (include helper functions)
###########################################################

# @import graphics

## ------- plot helper functions ---------
#' Graphical illustration of selection uncertainty
#' @export plot.csuv
#' @param x fitted results from CSUV::csuv()
#' @param is.box.plot TRUE to get a box plot like uncertainty illustration. FALSE to get a "confidence interval like" uncertainty illustration. Default is TRUE
#' @param compare.method.fit (optional) fitted results from CSUV::lm.compare.methods()
#' @param cv.mod (optional) fitted results from cross validation
#' @param with.thr whether the selection by the CSUV should be show. Default is TRUE
#' @param to.shade whether to shade the graph by the relative frequency calculated by CSUV. Default is TRUE
#' @param level (only use when is.box.plot is FALSE) the significant level of the "confidence interval". Default is 0.1
#' @param ci.type todo
#' @param xlim range of x shown in the graph
#' @param ylim range of y shown in the graph
#' @param var.freq.thr minimim variable frequency to show
#' @param return.lim TRUE to return the xlim and ylim
#' @param ... additional argument for plot
#' @examples
#' \dontrun{
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' mod.0 = csuv(X, Y, intercept = FALSE, q = 0, method.names = NULL)
#' cv.mod = lm.cv(X, Y, intercept = FALSE, fit.percent = 0.5, num.repeat = 50)
#' compare.mod = lm.compare.method(X, Y, intercept = FALSE)
#' plot(mod.0, compare.method.fit = compare.mod, cv.mod = cv.mod$est.b)
#' }
plot.csuv<-function(x,
                   is.box.plot = TRUE,
                   compare.method.fit = NULL,
                   cv.mod = NULL,
                   with.thr = TRUE,
                   to.shade = TRUE,
                   ci.type = "original",
                   level = 0.1,
                   ylim = NULL, xlim = NULL,
                   var.freq.thr = 0.05,
                   return.lim = FALSE, ...){
  # plot.new()
  par(new = FALSE)
  return (csuv.plot.helper(new.fit = x, is.box.plot = is.box.plot,
                          compare.method.fit = compare.method.fit,
                          compare.method.names = rownames(compare.method.fit),
                          cv.fit = cv.mod,
                          print.compare.method.points = FALSE,
                          with.thr = with.thr, to.shade = to.shade,
                          ci.type = ci.type,
                          level = level,
                          var.freq.thr = var.freq.thr,
                          ylim = ylim, xlim = xlim,
                          return.lim, ...))
}

csuv.plot.helper<-function(new.fit,
                          is.box.plot = TRUE,
                          compare.method.fit = NULL,
                          compare.method.names = NULL,
                          cv.fit = NULL,
                          print.compare.method.points = FALSE,
                          with.thr = TRUE, to.shade = TRUE,
                          ci.type = "original",
                          level = 0.1,
                          xlim = NULL, ylim = NULL,
                          return.lim = FALSE,
                          var.freq.thr = 0.05,
                          xlab = NULL, ylab = NULL, ...){
  # old.par <- par(no.readonly = TRUE) # store the current setting
  # source("compare_method.r", local = T)
  shiny::req(new.fit)
  # get fitted coef, order and frequency
  mod = new.fit$mod.collection
  if(is.null(colnames(mod))){ # make sure mod has colnames
    colnames(mod) = paste0("X", 0:(ncol(mod)-1))
  }
  var.freq = new.fit$variable.freq

  ## ====== update var.order =======
  var.order = new.fit$variable.order
  if(!is.null(compare.method.fit) | !is.null(cv.fit)){ # if we need to compare csuv fit with other methods
    not.chosen.i = which(var.freq==0) #variables that never chosen by csuv
    if(length(not.chosen.i)){
      cv.dummy = rep(0, length(not.chosen.i))
      compare.dummy = rep(0, length(not.chosen.i))
      if(!is.null(cv.fit)){
        cv.dummy = cv.fit[-1][not.chosen.i]
      }
      if(!is.null(compare.method.fit)){
        compare.dummy = colMeans(compare.method.fit[,-1,drop=FALSE][,not.chosen.i, drop = FALSE]!=0)
      }
      new.order = order(-compare.dummy,-cv.dummy)
      # update the var.order for zero coefficients
      var.order = c(var.order[1:sum(var.freq!=0)], not.chosen.i[new.order])
    }
  }

  # initialize variables we needed for later use
  cv.est.b = compare.mod = thr = mod.ci = truncated.mod.ci = NULL
  legend.val = c()
  compare.method.sel = rep(0, length(var.order))

  ## ====== reorder things according to the var.order =======
  # reorder mod, var.freq, ci, m.est.b by the var.order
  if(!is.box.plot){
    mod.ci = csuv.ci(csuv.fit = new.fit, level = level, type = ci.type)
    mod.ci = mod.ci[,var.order]
  }

  var.freq = var.freq[var.order]
  mod = mod[,-1][,var.order] # remove the intercept
  m.est.b = new.fit$est.b["csuv.m",-1][var.order] # remove the intercept
  p = ncol(mod)

  # plot the legend
  # plot the space for legend
  # par(fig=c(0.75,1,0,1), new=TRUE, mar = c(3,0.5,3,0))
  # plot(1, type="n", xaxt='n', yaxt = "n", xlab = "", ylab = "", axes=FALSE)

  # reorder cv and comparing mod est.b
  if(length(compare.method.fit)){ ## get compare mod values (need to calculate here so that we can set the xlim)
    if(is.null(compare.method.names)){
      compare.method.names = rownames(compare.method.fit)
    }
    compare.mod = compare.method.fit[compare.method.names,-1,drop=FALSE][,var.order,drop=FALSE]
    compare.method.sel = colSums(compare.mod!=0)
  }
  if(!is.null(cv.fit)){ # get cv value (need to calculate here so that we can set the xlim)
    cv.est.b = cv.fit[-1][var.order]
  }

  ## ====== set the x and y lim, and keep the mod according to xlim =======
  if(is.null(xlim)){
    xlim = range(c(1, which(var.freq>var.freq.thr), which(cv.est.b!=0)), na.rm = TRUE) #which(compare.method.sel!=0),
    xlim[2] = xlim[2]+0.1
  }
  if(is.null(ylim)){
    ylim = range(c(mod, compare.mod), na.rm=TRUE)
    if(all(compare.method.sel==0)){
      ylim[1] = ylim[1]-0.75
    } else{
      ylim[1] = ylim[1]-1.25
    }
    ylim[2] = ylim[2]+0.5
  }
  truncated.mod = mod[,1:floor(xlim[2])]
  if(!is.box.plot){
    truncated.mod.ci = mod.ci[,1:floor(xlim[2])]
  }
  rm(mod)
  rm(mod.ci)
  truncated.p = ncol(truncated.mod)

  ## ====== plot empty graph =======
  if(is.null(ylab)){ ylab = "estimated coefficients" }
  if(is.null(xlab)){ xlab = "sorted coefficients" }
  graphics::plot(1, type="n", xaxt='n', yaxt = "n", xlab = xlab, ylab = ylab,
       ylim = ylim, xlim = xlim, ...)
  par(new=TRUE)

  ## ====== shade according the var.freq and create legend values =======
  plot.legend = plot.col = plot.lty = plot.lwd = plot.pt.bg = plot.pch = c()
  if(to.shade || with.thr){
    if(to.shade){
      thr = sapply(c(9:1)/10, function(i) max(c(which(var.freq>= i), 0)))+0.5
      thr = c(0.5, thr, max(var.order))
      for(i in 1:10){
        # get shade color and legend val
        if(thr[i]!=thr[i+1]){
          rect.col = rgb(i/10,i/10,i/10, alpha=0.5)
          rect(xleft = thr[i], ybottom = ylim[2]+1, xright = thr[i+1], ytop = ylim[1]-1, border = NA, col = rect.col)
          legend.val = rbind(legend.val,
                             c(idx = paste((10-i)/10, "-", (11-i)/10),
                               col = rect.col))
        }
      }
      # points legend param
      n.pt = nrow(legend.val)
      plot.legend = c(plot.legend, legend.val[,"idx"])
      plot.col = c(plot.col, legend.val[,"col"])
      plot.lty = c(plot.lty, rep(0, n.pt))
      plot.lwd = c(plot.lwd, rep(0, n.pt))
      plot.pt.bg = c(plot.pt.bg, legend.val[,"col"])
      plot.pch = c(plot.pch, rep(15, n.pt))
    }

    if(with.thr){
      # line legend param
      plot.legend = c(plot.legend, "csuv.m cutoff", "csuv.s cutoff") #, "csuv.ebic cutoff"
      plot.col = c(plot.col, "green", "blue") #, "red"
      plot.lty =c(plot.lty, 1, 2)
      plot.lwd = c(plot.lwd, 1, 1)
      plot.pt.bg=c(plot.pt.bg, NA, NA)
      plot.pch=c(plot.pch, NA, NA)
    }
    # print legend
    legend("topright",
           legend = plot.legend, col = plot.col,
           lty = plot.lty, lwd = plot.lwd,
           pt.bg = plot.pt.bg, pch = plot.pch,
           cex = 0.8, ncol = 1, xpd = NA) #round(nrow(legend.val)/3)
    # par(new = TRUE)
    # prepare another part for graph
    # par(fig=c(0,0.75,0,1), mar = c(3,3,3,0.5))
    par(new=TRUE)
  }

  ## ====== plot box =======
  if(is.box.plot){
    conditional.truncated.mod = truncated.mod
    if(ci.type == "conditional"){
      conditional.truncated.mod[which(conditional.truncated.mod==0)] = NA
    }
    boxplot(conditional.truncated.mod, col = "gold", outpch = "-", medcol="orangered",
            ylim = ylim, xlim = xlim, width = var.freq[1:truncated.p],
            cex.axis = 0.8, xaxt = 'n',
            xlab = NA,
            main = "")

    ## ====== plot ci =======
  } else{
    # ui = sapply(1:truncated.p, function(i) stats::quantile(truncated.mod[,i], probs = 1-level/2))
    # li = sapply(1:truncated.p, function(i) stats::quantile(truncated.mod[,i], probs = level/2))

    # truncated.mod.ci = ci.helper(truncated.mod, alpha = level,
    #                              b = m.est.b[1:truncated.p])

    # temperarily removed
    # plotrix::plotCI(x = 1:truncated.p, y = m.est.b[1:truncated.p],
    #                 ui= truncated.mod.ci["ui",], li= truncated.mod.ci["li",], ylim = ylim,
    #                 xlab = NA, ylab = NA,
    #                 xlim = xlim,
    #                 #xlab = "covariates", ylab = "estimated coefficients",
    #                 scol = "orange", lwd = 2.5, xaxt='n',
    #                 col = "orangered", pch = 19, cex = 0.6)
    plotrix::plotCI(x = 1:truncated.p, y = colMeans(truncated.mod.ci),
                    ui= truncated.mod.ci["ui",], li= truncated.mod.ci["li",], ylim = ylim,
                    xlab = NA, ylab = NA,
                    xlim = xlim,
                    #xlab = "covariates", ylab = "estimated coefficients",
                    scol = "orange", lwd = 2.5, xaxt='n', col = "orange", cex = 0.1, pch = 19)
  }
  points(y = m.est.b[1:truncated.p], x = 1:truncated.p,
         col = "orangered", pch = 19, cex = 1)
  axis(1, at=1:truncated.p, labels = colnames(truncated.mod), las = 2, cex = 0.5)
  abline(h = 0)
  y.i = ifelse(all(compare.method.sel==0), ylim[1]+0.5, ylim[1]+1)
  text(x = 1:truncated.p, y = y.i, label = round(var.freq[1:truncated.p]*100), cex = 0.8)

  ## ====== add comparing methods and counts =======
  methods = get.compare.methods()
  if(!is.null(compare.mod)){
    if(print.compare.method.points){
      col = rainbow(length(methods)+1)
      for(method.name in compare.method.names){
        j = which(names(methods) == method.name)
        est.b = compare.mod[method.name,]
        points(y = est.b[which(est.b!=0)], x = which(est.b!=0), pch = 21, bg = col[j+1], col = "black", cex = 1.2)
      }
    }
    text(x = which(compare.method.sel!=0), y = ylim[1]+0.5, label = compare.method.sel[which(compare.method.sel!=0)], col = "blue", cex = 0.8)
  }
  if(!is.null(cv.est.b)){
    points(y = cv.est.b[which(cv.est.b!=0)], x = which(cv.est.b!=0),
           col = "blue", bg = "white", pch = 21, cex = 1.2)
  }

  ## ====== draw threshold lines ======
  if(with.thr){
    abline(v = sum(var.freq>=0.5)+0.5, col = "green", lwd = 2)
    median.size = sum(new.fit$est.b["csuv.s",-1]!=0)
    abline(v = median.size+0.5, lty = 2, col = "blue")
  }
  if(return.lim){
    return (c(xlim = xlim, ylim = ylim))
  }
}

