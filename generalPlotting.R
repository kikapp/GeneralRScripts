library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggthemes)
library(grid)
library(gridExtra)



# A theme for facet plots
theme_facet <- function() {
  return(theme_bw() + theme(strip.background = element_rect(color="black", fill = "white"),
                            panel.margin = unit(0.05,"null") ))
}

# A theme for cleaner looking heat maps
theme_heatmap <- function() {
  return(theme_bw() + theme(panel.background = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            axis.ticks = element_blank(),
                            panel.margin = rep(unit(0,"null"),4))
  )
}

# A blank theme, copied mostly or completely from a blog post 
# that I can't find anymore
blankground <- function() {
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        panel.margin = unit(0,"null"),
        plot.margin = rep(unit(0,"null"),4),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()
        
  )
}

## Taken from here: http://sape.inf.usi.ch/quick-reference/ggplot2/themes 
## and here :https://github.com/hadley/ggplot2/blob/ggplot2-0.9.2/R/theme-defaults.r
theme_complete_bw <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
  theme(
    axis.line =         element_blank(),
    axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
    axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
    axis.ticks =        element_line(colour = "black"),
    axis.title.x =      element_text(size = base_size, vjust = 0.5),
    axis.title.y =      element_text(size = base_size, angle = 0), #, vjust = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.margin = unit(0.1, "cm"),
    
    legend.background = element_rect(colour=NA), 
    legend.key =        element_rect(fill = NA, colour = "black", size = 0.25),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       element_text(size = base_size * 0.8),
    legend.title =      element_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
    
    panel.background =  element_rect(fill = NA, colour = "black", size = 0.25), 
    panel.border =      element_blank(),
    panel.grid.major =  element_blank(), #element_line(colour = "#aaaaaa", size = 0.05),
    panel.grid.minor =  element_blank(), #element_line(colour = "#aaaaaa", size = 0.05),
    panel.margin =      unit(0.25, "lines"),
    
    strip.background =  element_rect(fill = NA, colour = NA), 
    strip.text.x =      element_text(colour = "black", size = base_size * 0.8),
    strip.text.y =      element_text(colour = "black", size = base_size * 0.8, angle = -90),
    
    plot.background =   element_rect(colour = NA, fill = "white"),
    plot.title =        element_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
  )
}


# returns the legend from a ggplot object
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# For generating a residual QQ plot from a lm object

ggQQ_lims <- function(LM, title = "", lincol = "#999999", resids = F) # argument: a linear model
{
  if(!resids) { 
    resids <- LM$resid[!is.na(LM$resid)]
  }
  if(resids) { 
    resids <- LM[!is.na(LM)]
  }
  
  .df <- data.frame(y = sort(resids),
                    x = qnorm( seq(1, length(resids),1)/ (length(resids)+1) ))
  slope <- diff( quantile(.df$y, c(0.25, 0.75)))/diff(qnorm(c(0.25, 0.75)))
  int <- quantile(.df$y,0.5) - slope * quantile(.df$x,0.5)
  
  xlims <- c(min(.df$x), max(.df$x))
  ylims <- c(min(c(int + slope * min(.df$x), min(.df$y))), 
             max(c(int + slope * max(.df$x), max(.df$y))))
  p <- ggplot(data = .df, aes(y = y, x = x)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = slope, intercept = int, color=lincol, linetype = 2) +
    theme_bw() + 
    scale_y_continuous("Residuals", limits = ylims) +
    scale_x_continuous("Normal Theoretical Quantiles", limits = xlims) + 
    ggtitle(title)
  
  return(p)
}


# Plots bar charts en masse for categorical variables
# .groups = data as list object, list item names will be plotted on x axis
# .cat_vars = describes which variables to plot and how to order them
#  ex. .cat_vars <- list( list("VAR_1", "Pretty var 1 name to use in plot", c(ordered character array of VAR_1 levels for ordering in plot)),
#                        list("VAR_2", "Pretty var 2 name to use in plot", c(ordered character array of VAR_2 levels for ordering in plot)),
#                        etc.)
plot_categorical_summary <- function(.groups, .cat_vars, .output, w = 12, h = 8, pal = 'YlOrRd') {
  .var_cats <- names(.cat_vars)
  .df <- ldply(.groups, identity)
  .plots <- list()
  for (.var in .cat_vars) {
    print(.var[[1]])
    .to_plot <- ddply(.df, .(.id), function(..df, ..var) {
#       print(table(..df$.id, ..df[[..var]]))
      ..to_return <- data.frame(table(..df$.id, ..df[[..var]]) * 100/sum(table(..df$.id, ..df[[..var]]))) 
      #print(..to_return)
      return(..to_return)
    }, .var[[1]])
#     to_plot <- data.frame(table(.df[[.var[[1]]]], .df$.id))
    .plots[[.var[[1]]]] <- ggplot(data = .to_plot, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = 'identity', position = 'dodge', color = "black") +
      scale_x_discrete("") + 
      scale_y_continuous("%", breaks = seq(0,100,10)) +
      scale_fill_brewer(.var[[2]], breaks = rev(.var[[3]]), palette = pal) +
      coord_flip() +
      theme_bw()
  }
  
  pdf(.output, width = w, height = h)
  l_ply(.plots, print)
  dev.off()

}


# Generates boxplot plots en masse for continuous variables
# .groups = data as list object, list item names will be plotted on x axis
# .cat_vars = describes which variables to plot and how to order them
#  ex. .cat_vars <- list( list("VAR_1", "Pretty var 1 name to use in plot"),
#                        list("VAR_2", "Pretty var 2 name to use in plot"),
#                        etc.)
# .groups <- culled_groups
# .cont_vars <- cont.vars.to.summarize
# .var <- .cont_vars[[2]]

plot_continuous_summary <- function(.groups, .cont_vars, .output, w = 12, h = 8) {
  .var_cats <- names(.cont_vars)
  .df <- ldply(.groups, identity)
  .plots <- list()
  for (.var in .cont_vars) {
    print(.var[[1]])
#     .df[1,]
    .to_plot <- .df[ ,c(.var[[1]], ".id")]
    names(.to_plot) <- c('y', 'x')
    ktest <- kruskal.test(.to_plot, y ~ x)$p.value
    ktest <- ifelse(ktest < 0.001, "< 0.001", paste0("= ", round(ktest,3)))
    summary(.to_plot)
    .plots[[.var[[1]]]] <- ggplot(data = .to_plot, aes(x = x, y = y)) + geom_boxplot() +
      scale_x_discrete("") + 
      scale_y_continuous(.var[[1]]) +
      coord_flip() +
      ggtitle(paste0(.var[[2]], " - Kruskal Wallis p ", ktest)) +
      theme_bw()
    if (min(.to_plot$y, na.rm = T) > 0 ) {
      ktest <- kruskal.test(.to_plot, log(y) ~ x)$p.value
      ktest <- ifelse(ktest < 0.001, "< 0.001", paste0("= ", round(ktest,3)))
      .plots[[paste0(.var[[1]], " log")]] <- ggplot(data = .to_plot, aes(x = x, y = log(y))) + geom_boxplot() +
        scale_x_discrete("") + 
        scale_y_continuous(.var[[2]]) +
        coord_flip() +
        ggtitle(paste0(.var[[2]], " log_e - Kruskal Wallis p ", ktest)) +
        theme_bw()
    }
  }  
  pdf(.output, width = w, height = h)
  l_ply(.plots, print)
  dev.off()
  
}

# generates kaplan-meier plot using ggplot2
# sfit is object returned by survfit function
# requires survival and ggplot2

# fit <- survfit(Surv(time,status)~rx, data=colon)
# ggkm(fit, timeby=500, ystratalabs=c("Obs","Lev","Lev+5FU"))
# modified from:
# FROM http://pastebin.com/FjkWnCWm
ggkm <- function(sfit,
                 #table = FALSE,
                 returns = FALSE,
                 xlabs = "Time",
                 ylabs = "Survival Probability",
                 xlims = c(0,max(sfit$time)),
                 ylims = c(0,1),
                 ystratalabs = NULL,
                 ystrataname = NULL,
                 xbreaks = NULL,
                 timeby = 365.25,
                 main = "Kaplan-Meier Plot",
                 pval = TRUE,
                 subs = NULL,
                 pal = "Black",
                 strata_label = "Strata",
                 ...) {
  
  #############
  # libraries #
  #############
  
  require(ggplot2)
  require(survival)
  require(gridExtra)
  require(RColorBrewer)
  require(plyr)
  
  #################################
  # sorting the use of subsetting #
  #################################
  
  
  times <- seq(0, max(sfit$time), by = timeby)
  
  if(is.null(subs)){
    subs1 <- 1:length(levels(summary(sfit)$strata))
    subs2 <- 1:length(summary(sfit,censored=T)$strata)
    subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$strata)
  } else {
    for(i in 1:length(subs)){
      if(i==1){
        ssvar <- paste("(?=.*\\b=",subs[i],sep="")
      }
      if(i==length(subs)){
        ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],"\\b)",sep="")
      }
      if(!i %in% c(1, length(subs))){
        ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],sep="")
      }
      if(i==1 & i==length(subs)){
        ssvar <- paste("(?=.*\\b=",subs[i],"\\b)",sep="")
      }
    }
    subs1 <- which(regexpr(ssvar,levels(summary(sfit)$strata), perl=T)!=-1)
    subs2 <- which(regexpr(ssvar,summary(sfit,censored=T)$strata, perl=T)!=-1)
    subs3 <- which(regexpr(ssvar,summary(sfit,times = times,extend = TRUE)$strata, perl=T)!=-1)
  }
  
  if( !is.null(subs) ) pval <- FALSE
  
  ##################################
  # data manipulation pre-plotting #
  ##################################
  
  if(is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*","",names(sfit$strata))) #[subs1]
  if(is.null(ystrataname)) ystrataname <- "Strata"
  if(tolower(pal) == "black") {
    pall = rep("#000000", length(ystratalabs))
  } else {
    pall = brewer.pal(length(ystratalabs), pal)
  }
  m <- max(nchar(ystratalabs))
  times <- seq(0, max(sfit$time), by = timeby)
  if(is.null(xbreaks)) xbreaks <- times/timeby
  
  .df <- data.frame(                      # data to be used in the survival plot
    time = sfit$time[subs2]/timeby,
    n.risk = sfit$n.risk[subs2],
    n.event = sfit$n.event[subs2],
    surv = sfit$surv[subs2],
    strata = factor(summary(sfit, censored = T)$strata[subs2]),
    upper = sfit$upper[subs2],
    lower = sfit$lower[subs2]
  )
  
  
  
  strata_var <- strsplit(factorConvert(.df$strata, "character")[1], "=")[[1]][1]
  .df$strata_ch <- factorConvert(.df$strata, "character")
  .df$strata_ch <- gsub(pattern = paste0(strata_var, "="), replacement = "", x = .df$strata_ch)
  
  #SET LABELS IN DECREASING ORDER OF STRATA SIZE
  ystratalabs <- gsub(pattern = paste0(strata_var, "="), 
                      replacement = "", 
                      x = rownames(summary(sfit)$table)[order(summary(sfit)$table[,1], decreasing=TRUE)])  
  
  
  .df$strata <- factor(.df$strata_ch, levels =  ystratalabs)
  #levels(.df$strata) <- ystratalabs       # final changes to data for survival plot
  zeros <- data.frame(time = 0, surv = 1,
                      strata = factor(ystratalabs, levels=levels(.df$strata)),
                      upper = 1, lower = 1)
  .df <- rbind.fill(zeros, .df)
  d <- length(levels(.df$strata))
  
  ###################################
  # specifying plot parameteres etc #
  ###################################
  #ylims <- c(0.95, 1)
  p <- ggplot( .df, aes(x = time, y = surv, linetype = factor(strata), color = factor(strata))) +
    geom_step(size = 0.7) +
    scale_linetype_manual(strata_label, values = .df$strata) + 
    scale_color_manual(strata_label, values = pall) +
    theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    scale_x_continuous(xlabs, breaks = xbreaks, limits = xlims/timeby) +
    scale_y_continuous(ylabs, limits = ylims) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position = "bottom") +  
    #theme(legend.position = c(ifelse(m < 10, .28, .35),ifelse(d < 4, .25, .35))) +    # MOVE LEGEND HERE [first is x dim, second is y dim]
    theme(legend.key = element_rect(colour = NA)) +
    labs(linetype = ystrataname) +
    theme(plot.margin = unit(c(0, 1, .5,ifelse(m < 10, 1.5, 2.5)),"lines")) +
    ggtitle(main)
  
  ## Create a blank plot for place-holding
  ## .df <- data.frame()
  #   blank.pic <- ggplot(.df, aes(time, surv)) +
  #     geom_blank() + theme_bw() +
  #     theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
  #           axis.title.x = element_blank(),axis.title.y = element_blank(),
  #           axis.ticks = element_blank(),
  #           panel.grid.major = element_blank(),panel.border = element_blank())
  
  #####################
  # p-value placement #
  #####################a
  
  if(pval) {
    sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
    pvaltxt <- ifelse(pval < 0.001,"Log-rank p < 0.001",paste("Log-rank p =", round(pval, 3)))
    p <- p + annotate("text",x = 0.1 * max(sfit$time/timeby),y = mean(ylims[1]),label = pvaltxt)
  }
  
  if(returns) return(p)
}

theme_set(theme_complete_bw())

generalPlottingLoaded <- TRUE
