###########################################################
## This file includes functions related to plotting for CSUV
## (include helper functions)
###########################################################

# @import graphics

## ------- plot helper functions ---------
#' Graphical illustration of selection uncertainty
#' @export
#' @param x fitted results from CSUV::csuv()
#' @param with.unconditional TRUE to get a unconditonal boxplot on the same graph. Default is FALSE
#' @param compare.method.fit (optional) fitted results from CSUV::lm.compare.methods()
#' @param cv.mod (optional) fitted results from cross validation
#' @param with.thr whether the selection by the CSUV should be show. Default is TRUE
#' @param to.shade whether to shade the graph by the relative frequency calculated by CSUV. Default is TRUE
#' @param level the significant level of the whiskers. Default is 0.1
#' @param var.freq.thr minimum variable frequency to show
#' @param ... additional argument for plot
#' @return a ggplot object
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
                    with.unconditional = FALSE,
                    compare.method.fit = NULL,
                    cv.mod = NULL,
                    with.thr = TRUE,
                    to.shade = TRUE,
                    level = 0.1,
                    var.freq.thr = 0.05, ...){
  # plot.new()
  par(new = FALSE)
  return (csuv.plot.helper(new.fit = x,
                           with.unconditional = with.unconditional,
                           compare.method.fit = compare.method.fit,
                           compare.method.names = rownames(compare.method.fit),
                           cv.fit = cv.mod,
                           print.compare.method.points = FALSE,
                           with.thr = with.thr, to.shade = to.shade,
                           level = level,
                           var.freq.thr = var.freq.thr, ...))
}

csuv.plot.helper<-function(new.fit,
                           with.unconditional = FALSE,
                           compare.method.fit = NULL,
                           compare.method.names = NULL,
                           cv.fit = NULL,
                           print.compare.method.points = FALSE,
                           with.thr = TRUE, to.shade = TRUE,
                           level = 0.1,
                           var.freq.thr = 0.05, ...){
  shiny::req(new.fit)

  # get fitted coef, order and frequency
  mod = new.fit$mod.collection
  if(is.null(colnames(mod))){ # make sure mod has colnames
    colnames(mod) = paste0("X", 0:(ncol(mod)-1))
  }

  var.freq = new.fit$variable.freq
  var.order = get.plot.var.order(new.fit$variable.order, var.freq, compare.method.fit, cv.fit, var.freq.thr = var.freq.thr)

  ## get compare mod values (need to calculate here so that we can set the xlim)
  if(length(compare.method.fit) && is.null(compare.method.names)){
    compare.method.names = rownames(compare.method.fit)
  }

  ## ======== convert to ggplot df  =======
  ggplot.df = get.df.for.gg.plot(mod = mod[,-1], csuv.mod = new.fit$est.b["csuv.m",-1],
                                 compare.mod = compare.method.fit[compare.method.names,-1,drop=FALSE],
                                 cv.mod = cv.fit[-1],
                                 tau = var.freq,
                                 var.order = var.order)

  ## ====== prepare info for plotting =======
  shade.info = get.shade.col.n.lab(var.freq)
  get.cond.boxplot.stat<-function(x){
    get.boxplot.stat(x, is.conditional = TRUE)
  }
  get.cond.errorbar.stat<-function(x){
    get.errorbar.stat(x, is.conditional = TRUE)
  }

  get.boxplot.stat<-function(x, is.conditional = FALSE){
    r = rep(0,6)
    cond.x = x
    if(is.conditional){ cond.x = x[which(x!=0)]}
    if(length(cond.x)){
      tau = ifelse(is.conditional, max(mean(x>0), mean(x<0))*0.8, 0.8)
      r = c(quantile(cond.x, probs = c(0.25, 0.25, 0.5, 0.75, 0.75)), tau)
    }
    names(r) = c("ymin", "lower", "middle", "upper", "ymax", "width")
    return (r)
  }
  get.errorbar.stat<-function(x, is.conditional = FALSE){
    r = rep(0,3)
    cond.x = x
    if(is.conditional){ cond.x = x[which(x!=0)]}
    if(length(cond.x)){
      tau = ifelse(is.conditional, max(mean(x>0), mean(x<0))*0.8, 0.8)
      r = c(quantile(cond.x, probs = c((level/2), (1-level/2))), tau)
    }
    names(r) = c("ymin", "ymax", "width")
    return (r)
  }

  y.min = min(ggplot.df[,"coefficients"])
  method.point.col = c(csuv = "red")
  method.point.sty = c(csuv = 16)
  if(!is.null(cv.fit)){
    method.point.col = c(method.point.col, cv = "blue")
    method.point.sty = c(method.point.sty, cv = 1)
  }
  if(with.thr){
    method.point.col = c(method.point.col, csuv_m = "green", csuv_s = "yellow")
  }
  ## ========= plot ============
  p = ggplot2::ggplot(ggplot.df, ggplot2::aes(x=variables, y =coefficients))

  if(to.shade){
    p = p + ggplot2::geom_tile(data = subset(ggplot.df, type == "tau"),
                               ggplot2::aes(fill = factor(pmin(coefficients%/%10*10, 90),
                                                          levels = shade.info$tau.factors, ordered = TRUE),
                                            y = 0, height = Inf),
                               alpha = 0.3)+
      ggplot2::scale_fill_manual(name = bquote(tau) , values = shade.info$col.rgb,
                        label = shade.info$col.labels)
  }
  p = p + ggplot2::geom_hline(yintercept = 0)


  # == add the box plot ==
  p = p + ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                                fun.data = get.cond.errorbar.stat, geom="errorbar", col = "black", lwd = 0.7)+
    ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                          fun.data = get.cond.boxplot.stat, geom="boxplot", varwidth = TRUE, fill="gold", color = "orangered")
  if(with.unconditional){
    p = p + ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                         fun.data = get.errorbar.stat, geom="errorbar", col = "brown", lwd = 0.7, alpha = 0.5)+
      ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                   fun.data = get.boxplot.stat, geom="boxplot", varwidth = TRUE, fill="palegreen", color = "royalblue", alpha = 0.5)
  }

  # == add coef estimation points ==
  p = p + ggplot2::geom_point(data = subset(ggplot.df, type%in%names(method.point.col)),
                              ggplot2::aes(color = type, shape = type)) +
    ggplot2::scale_shape_manual(name = "methods", values = method.point.sty) +
    ggplot2::scale_color_manual(name = "methods", values = method.point.col)

  # == add tau and selection frequency text ==
  p = p + ggplot2::geom_text(data = subset(ggplot.df, type == "tau"),
                             ggplot2::aes(y = y.min-1, label = coefficients)) +
    ggplot2::geom_text(label = "tau", x = 0.75, y =y.min-1, parse = TRUE)

  if(!is.null(compare.method.fit)){
    p = p + ggplot2::geom_text(data = subset(ggplot.df, type == "compare"),
                         ggplot2::aes(y = y.min-1.5, label = coefficients), col = "blue")
  }

  # geom_text(label = "compare", x = 0.75, y = -1.5) +


  # == add threshold ==
  cf = data.frame("type" = c("csuv_m", "csuv_s"), cutoff = c((sum(new.fit$est.b["csuv.m",]!=0)+0.5), (sum(new.fit$est.b["csuv.s",]!=0)+0.5)), dummy = 0)
  if(with.thr){
    p = p + ggplot2::geom_vline(data = cf, ggplot2::aes(xintercept = cutoff, linetype = type), color = c("chartreuse2", "blue"), size = c(1.5,1))+
      ggplot2::scale_linetype_manual(name = "cutoff",
                            values = c(csuv_m = "solid", csuv_s = "dashed"),
                            guide = ggplot2::guide_legend(override.aes = list(color = c("chartreuse2", "blue"),
                                                                              size = c(1.5,1))))
  }

  # update axis
  p = p + ggplot2::ylab("estimated coefficients") + ggplot2::xlab("sorted covariates")
  return (p)
}

