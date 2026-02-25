multitypeFig<- function(y, maxAll, Modeltype = ""){
  nstrain<- dim(y)[3]
  maxY<- maxAll
  plotlists<- list()

  for(i in 1:nstrain){
    Strain<- i
    sim.object<- y[,,i]

    spatdata<- sim.object

    ts_spatdata <- as.data.frame(t(spatdata))
    ts_spatdata$Time <- 1:ncol(spatdata)
    naming<- c(paste("u", 1:(ncol(ts_spatdata)-1), sep=""), "Time")
    colnames(ts_spatdata)<- naming

    Colors <- rep(c("blue", "red"), length.out = nrow(spatdata))
    Linetypes <- rep(c("dotted", "dashed", "dotdash", "longdash", "twodash"), length.out = nrow(spatdata))

    long_data <- reshape2::melt(ts_spatdata, id.vars = "Time")
    library(ggplot2)
    if(Strain==nstrain){
      rfigs<- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable, linetype = variable)) +
        geom_line() +
        scale_color_manual(values = Colors) +
        scale_linetype_manual(values = Linetypes) +
        ylim(0, maxY) +
        labs(x = "Time [month]", y = "Case counts", color = "Location") +
        guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
        theme(axis.title.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 16))
    }else{
      rfigs<- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable, linetype = variable)) +
        geom_line() +
        scale_color_manual(values = Colors) +
        scale_linetype_manual(values = Linetypes) +
        ylim(0, maxY) +
        labs(x = "Time [month]", y = "Case counts", color = "Location") +
        guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
        theme(axis.title.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16),
              legend.position = "none")
    }
    plotlists[[Strain]]<- rfigs
  }
  row_1<- cowplot::plot_grid(plotlist = plotlists[1:3], ncol = 3, labels = c("A", "B", "C"), label_size = 17)
  row_2<- cowplot::plot_grid(plotlist = plotlists[4:5], ncol = 3, labels = c("D", "E"), label_size = 17, rel_widths = c(1, 1.35, 0.65))
  finalplot<- cowplot::plot_grid(row_1, row_2, nrow = 2)
  print(finalplot)
  add_legend(0.58, -0.18, legend=substitute(paste(bold(Modeltype))),
             col="black",
             horiz=TRUE, bty='n', cex=5.0)
}


multitypeFig_allsim<- function(y_arraylist, Modeltype = "") {
  library(ggplot2)
  library(cowplot)
  library(reshape2)

  nmodel  <- length(y_arraylist)
  nstrain <- dim(y_arraylist[[1]][[1]])[3]

  maxVec<- numeric(nmodel)
  for(n in 1:nmodel){
    model_n<- y_arraylist[[n]][[1]]
    maxVec[n]<- max(model_n)
  }
  maxY<- max(maxVec)

  plotlists <- list()
  shared_legend <- NULL

  for(m in 1:nmodel){

    y<- y_arraylist[[m]][[1]]

    for(i in 1:nstrain){

      index<- (m - 1) * nstrain + i

      sim.object <- y[ , , i]
      ts_spatdata <- as.data.frame(t(sim.object))
      ts_spatdata$Time <- 1:ncol(sim.object)

      naming<- c(paste0("u", 1:(ncol(ts_spatdata)-1)), "Time")
      colnames(ts_spatdata) <- naming

      Colors<- rep(c("blue", "red"),
                       length.out = nrow(sim.object))
      Linetypes<- rep(c("dotted", "dashed",
                         "dotdash", "longdash",
                         "twodash"),
                       length.out = nrow(sim.object))

      long_data <- melt(ts_spatdata, id.vars = "Time")

      base_plot <- ggplot(long_data,
                          aes(Time, value,
                              color = variable,
                              linetype = variable)) +
        geom_line() +
        scale_color_manual(values = Colors, drop = FALSE) +
        scale_linetype_manual(values = Linetypes, drop = FALSE) +
        ylim(0, maxY) +
        labs(x = "Time [month]",
             y = "Case counts",
             color = "Location",
             linetype = "Location") +
        theme(
          axis.title = element_text(size = 18),
          axis.text  = element_text(size = 16)
        )
      # Titles for first row
      if (m == 1) {
        base_plot <- base_plot +
          ggtitle(paste0("Strain ", i)) +
          theme(plot.title = element_text(
            hjust = 0.5,
            size  = 18,
            face  = "bold"
          ))
      }
      # Extract legend
      if (m == nmodel && i == nstrain) {
        legend_plot  <- base_plot +
          theme(legend.position = "bottom",
                legend.title = element_text(size = 18),
                legend.text  = element_text(size = 16))
        shared_legend <- cowplot::get_legend(legend_plot)
      }

      plotlists[[index]] <- base_plot +
        theme(legend.position = "none")
    }
  }
  #left labels
  row_labels <- LETTERS[1:nmodel]
  row_list <- list()

  for (m in 1:nmodel) {

    row_plots <- plotlists[((m - 1) * nstrain + 1):(m * nstrain)]

    row_grid <- plot_grid(plotlist = row_plots,
                          ncol = nstrain)

    label_plot <- ggdraw() +
      draw_label(row_labels[m],
                 fontface = "bold",
                 size = 20,
                 x = 0.5, y = 0.5)

    row_list[[m]] <- plot_grid(label_plot,
                               row_grid,
                               ncol = 2,
                               rel_widths = c(0.06, 1))
  }
  main_grid <- plot_grid(plotlist = row_list,
                         ncol = 1)
  finalplot <- plot_grid(main_grid,
                         shared_legend,
                         ncol = 1,
                         rel_heights = c(1, 0.08))
  print(finalplot)
}


multitypeFig2 <- function(array.object, names = NULL){
  plotlists<- list()
  nstrain<- dim(array.object)[3]
  maxY<- max(array.object, na.rm = T)
  for(i in 1:nstrain){
    spatdata <- array.object[,,i]

  ts_spatdata <- as.data.frame(t(spatdata))
  ts_spatdata$Time <- seq.Date(from = as.Date("2011-01-01"), to = as.Date("2019-12-01"), by = "month")
  if(is.null(names)){
    colnames(ts_spatdata) <- c(paste("u", 1:(ncol(ts_spatdata) - 1), sep = ""), "Time")
  }else{
    colnames(ts_spatdata) <- c(names, "Time")
  }
  long_data <- reshape2::melt(ts_spatdata, id.vars = "Time")

  library(ggplot2)
  if(i==nstrain){
    a <- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
      geom_line() +
      ylim(0, maxY) +
      labs(x = "Time [month/year]", y = "Case counts", color = "Location") +
      guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
      scale_x_date(date_labels = "%b %Y", date_breaks = "1 years") +
      theme(axis.title.y = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 23))
  }else{
  a <- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
    geom_line() +
    ylim(0, maxY) +
    labs(x = "Time [month/year]", y = "Case counts", color = "Location") +
    guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
    scale_x_date(date_labels = "%b %Y", date_breaks = "1 years") +
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 23),
          legend.position = "none")
  }
  plotlists[[i]]<- a
  }
  row_1<- cowplot::plot_grid(plotlist = plotlists[1:3], ncol = 3, labels = c("A", "B", "C"), label_size = 20)
  row_2<- cowplot::plot_grid(plotlist = plotlists[4], ncol = 2, labels = c("D"), label_size = 20, rel_widths = c(0.82, 1))
  finalplot<- cowplot::plot_grid(row_1, row_2, nrow = 2)
  print(finalplot)
  #export ==> 33 x 13
}

mcmc.plot<- function(inf.object, Histograms=FALSE){
  if(Histograms){
     par(mfrow=c(3, 3))
      for (i in 1:ncol(inf.object)) {
        hist(inf.object[, i], main = colnames(inf.object)[i], xlab ="", col = "white", border = "black")
      }
  }
    par(mfrow=c(3, 3))
    for (i in 1:ncol(inf.object)) {
      plot(inf.object[, i], type = "l", main = colnames(inf.object)[i], xlab ="MCMC iterations", ylab = "", col = "purple")
      grid()
    }
}

inf.plot<- function(inf.object){

    rPosterior<- inf.object[, startsWith(colnames(inf.object), "r")]
    sPosterior<- inf.object[, startsWith(colnames(inf.object), "s")]
    inf.r<- colMeans(rPosterior)
    inf.s<- colMeans(sPosterior)
    uCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,2]
    lCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,1]
    uCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,2]
    lCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,1]

  par(mfrow=c(1,2))
  plot(0, type = "n", xlim = c(1,length(inf.r)), ylim = c(min(lCI.r, inf.r), max(uCI.r, inf.r)), ylab = "Trend component", xlab = "Time")
  polygon(c(1:length(inf.r), rev(1:length(inf.r))), c(lCI.r, rev(uCI.r)),
          col = "pink", border = NA)
  lines(1:length(inf.r), inf.r, col="red")
  grid()

  plot(0, type = "n", xlim = c(1,length(inf.s)), ylim = c(min(lCI.s, inf.s), max(uCI.s, inf.s)), ylab = "Seasonal component", xlab = "Season")
  polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
          col = "pink", border = NA)
  lines(1:length(inf.s), inf.s, col="red")
  grid()
  add_legend("topright", legend="Posterior means", lty=1, col="red",
             horiz=TRUE, bty='n', cex=1.1)
}

RecoverInf.plot<- function(inf.object, true_r, true_s, Modeltype=""){

    rPosterior<- inf.object[, startsWith(colnames(inf.object), "r")]
    sPosterior<- inf.object[, startsWith(colnames(inf.object), "s")]
    inf.r<- colMeans(rPosterior)
    inf.s<- colMeans(sPosterior)
    uCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,2]
    lCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,1]
    uCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,2]
    lCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,1]

  par(mfrow=c(1,2))
  plot(0, type = "n", xlim = c(1,length(inf.r)), ylim = c(min(lCI.r, inf.r), max(uCI.r, inf.r)), ylab = "Trend component", xlab = "Time")
  polygon(c(1:length(inf.r), rev(1:length(inf.r))), c(lCI.r, rev(uCI.r)),
          col = "pink", border = NA)
  lines(1:length(inf.r), inf.r, col="red")
  points(1:length(inf.r), true_r, pch = 19)
  grid()

  plot(0, type = "n", xlim = c(1,length(inf.s)), ylim = c(min(lCI.s, inf.s), max(uCI.s, inf.s)), ylab = "Seasonal component", xlab = "Season")
  polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
          col = "pink", border = NA)
  lines(1:length(inf.s), inf.s, col="red")
  points(1:length(inf.s), true_s, pch = 19)
  grid()
  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.4)
  add_legend("topleft", legend=substitute(paste(bold(Modeltype))),
             col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.5)

}

RecoverInfAllRS.plot<- function(allinf.object, true_r, true_s){

  allobject<- length(allinf.object)
  pdf("AllRS1.pdf", paper="a4", width=12,height=12, pointsize=12)

  par(mfrow=c(4,2))
  # Modeltype<- c("No outbreak", "Independent 1", "Independent 2", "Gaussian copula-dependent 1",
  #                "Gaussian copula-dependent 2", "Frank copula-dependent 1",
  #               "Frank copula-dependent 2", "General-dependent")
  alphabets<- c("A", "C", "E", "G")
  #alphabets<- c("B", "D", "F", "H")
  for(a in 1:allobject){
    inf.object<- allinf.object[[a]]

    if(is.data.frame(inf.object)){
      rPosterior<- inf.object[, startsWith(colnames(inf.object), "r")]
      sPosterior<- inf.object[, startsWith(colnames(inf.object), "s")]
      inf.r<- colMeans(rPosterior)
      inf.s<- colMeans(sPosterior)
      uCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,2]
      lCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,1]
      uCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,2]
      lCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,1]
    }else{
      inf.r<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
      uCI.r<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "r")[,1,]))[,2]
      lCI.r<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "r")[,1,]))[,1]
      inf.s<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
      uCI.s<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "s")[,1,]))[,2]
      lCI.s<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "s")[,1,]))[,1]
    }

    plot(0, type = "n", xlim = c(1,length(inf.r)), ylim = c(-0.5, 0.5), ylab = "Trend component", xlab = "Time [Month]", main = "", cex.main=2.0, cex.lab=1.5)
    polygon(c(1:length(inf.r), rev(1:length(inf.r))), c(lCI.r, rev(uCI.r)),
            col = "pink", border = NA)
    lines(1:length(inf.r), inf.r, col="red")
    points(1:length(inf.r), true_r, pch = 19)
    grid()
    legendary::labelFig(alphabets[a], adj = c(-0.15, 0.10), font=2, cex=2.0)

    plot(0, type = "n", xlim = c(1,length(inf.s)), ylim = c(-1.6, 1.6), ylab = "Seasonal component", xlab = "Season [Month]", main = "", cex.main=2.0, cex.lab=1.5)
    polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
            col = "pink", border = NA)
    lines(1:length(inf.s), inf.s, col="red")
    points(1:length(inf.s), true_s, pch = 19)
    grid()
  }
#  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
#             pch=c(19, NA), col=c("black", "red"),
#             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}

RecoverInfUA.plot <- function(inf.object, true_u, true_a_k,
                              Modeltype = "", burn.in = 100) {
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]

  #U violin plot
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Spatial comp.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #a_k violin plot
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercept", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  return(list(rfigs_u,rfigs_ak))

#  plotlists <- list(rfigs_u, rfigs_ak)
#  print(cowplot::plot_grid(plotlist = plotlists, ncol = 2,
#                           labels = c("A", ""),
#                           rel_widths = c(1.25, 1), label_size = 25))

  # Legends
  #  add_legend(0.85, 1.15, legend = "Truth",
  #             pch = 19, col = "black",
  #             horiz = TRUE, bty = 'n', cex = 1.8)
  #  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
  #             horiz = TRUE, bty = 'n', cex = 1.5)
}


AllViolinplots<- function(plotlists){
  library(ggplot2)
  library(cowplot)
  row_1<- cowplot::plot_grid(plotlist = plotlists[[1]], ncol = 5, labels = c("A", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1.60, 1))
  row_2<- cowplot::plot_grid(plotlist = plotlists[[2]], ncol = 5, labels = c("B", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1, 1))
  row_3<- cowplot::plot_grid(plotlist = plotlists[[3]], ncol = 5, labels = c("C", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1.60, 1))
  row_4<- cowplot::plot_grid(plotlist = plotlists[[4]], ncol = 5, labels = c("D", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1.60, 1))
  row_5<- cowplot::plot_grid(plotlist = plotlists[[5]], ncol = 5, labels = c("E", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1.60, 1))
  row_6<- cowplot::plot_grid(plotlist = plotlists[[6]], ncol = 5, labels = c("F", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1.60, 1))
  row_7<- cowplot::plot_grid(plotlist = plotlists[[7]], ncol = 5, labels = c("G", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1.60, 1))
  row_8<- cowplot::plot_grid(plotlist = plotlists[[8]], ncol = 5, labels = c("H", "", "", "", ""), label_size = 22, rel_widths = c(1.16, 1.16, 1.16, 1.60, 1))
  print(cowplot::plot_grid(row_1, row_2, row_3,row_4,row_5,row_6,row_7,row_8, nrow = 8))
  add_legend(0.75, 1.40, legend = "Truth",
             pch = 19, col = "black",
             horiz = TRUE, bty = 'n', cex = 5)
  #36 x 26 Export
}


RecoverInfUAB.plot <- function(inf.object, true_u, true_a_k, true_B,
                               Modeltype = "", burn.in = 100) {
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)
  n_B  <- length(true_B)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
    fullB.draws <- as.data.frame(inf.object$draws(variables = "B")[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
    fullB.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]
  B.draws  <- fullB.draws[thinning, , drop = FALSE]

  #U violin plot
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Spatial comp.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #a_k violin plot
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercept", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #B violin plot
  B_labels <- do.call(expression, lapply(1:n_B, function(i) {
    bquote(beta[.(i)])
  }))

  comp_B <- data.frame(
    value = as.vector(as.matrix(B.draws)),
    group = factor(rep(paste0("B", 1:n_B), each = nrow(B.draws)))
  )

  rfigs_B <- ggplot(comp_B, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("B", 1:n_B), levels = levels(comp_B$group)),
                                 y = true_B),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0.5, 2.5) +
    labs(x = "Regression coeff.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_B)) +
    scale_x_discrete(labels = B_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #combine plots
  return(list(rfigs_u, rfigs_ak, rfigs_B))
#  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B)
#  print(cowplot::plot_grid(plotlist = plotlists, ncol = 3,
                           #                           labels = c("A", "B", "C"),
#                           labels = c("H", "", ""),
#                           rel_widths = c(1.25, 1, 1), label_size = 25))

  # Legends
  #  add_legend(0.85, 1.15, legend = "Truth",
  #             pch = 19, col = "black",
  #             horiz = TRUE, bty = 'n', cex = 1.8)
  #  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
  #             horiz = TRUE, bty = 'n', cex = 1.5)
}

RecoverInfUABG.plot<- function(inf.object, true_u, true_a_k, true_B, true_G,
                               Modeltype = "", burn.in = 100){
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)
  n_B  <- length(true_B)
  n_G  <- length(true_G)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
    fullB.draws <- as.data.frame(inf.object$draws(variables = "B")[, 1, ])
    fullG.draws <- as.data.frame(inf.object$draws(variables = paste0("G", 1:n_G))[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
    fullB.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
    fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]
  B.draws  <- fullB.draws[thinning, , drop = FALSE]
  G.draws  <- fullG.draws[thinning, , drop = FALSE]

  #U violin plot
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Spatial comp.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #a_k violin plot
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercept", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # B violin plot
  B_labels <- do.call(expression, lapply(1:n_B, function(i) {
    bquote(beta[.(i)])
  }))

  comp_B <- data.frame(
    value = as.vector(as.matrix(B.draws)),
    group = factor(rep(paste0("B", 1:n_B), each = nrow(B.draws)))
  )

  rfigs_B <- ggplot(comp_B, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("B", 1:n_B), levels = levels(comp_B$group)),
                                 y = true_B),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0.5, 2.5) +
    labs(x = "Regression coeff.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_B)) +
    scale_x_discrete(labels = B_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #G violin plot
  G_index_pairs <- rep(list(c(0,1), c(1,0)), ncol(G.draws)/2)

  G_labels <- do.call(expression, lapply(G_index_pairs, function(idx) {
    bquote(gamma[.(idx[1])][.(idx[2])])
  }))

  comp_G <- data.frame(
    value = as.vector(as.matrix(G.draws)),
    group = factor(rep(paste0("G", 1:n_G), each = nrow(G.draws)))
  )

  rfigs_G <- ggplot(comp_G, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("G", 1:n_G),
                                            levels = levels(comp_G$group)),
                                 y = true_G),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0, 1) +
    labs(x = "Transition prob.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_G)) +
    scale_x_discrete(labels = G_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  return(list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G))

#  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G)
#  print(cowplot::plot_grid(plotlist = plotlists, ncol = 4,
                           #                           labels = c("A", "B", "C", "D","E"),
#                           labels = c("C", "", "", ""),
#                           rel_widths = c(1.25, 1, 1, 1.25), label_size = 25))

  # Legends
#  add_legend(0.85, 1.15, legend = "Truth",
#             pch = 19, col = "black",
#             horiz = TRUE, bty = 'n', cex = 2.0)
  #  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
  #             horiz = TRUE, bty = 'n', cex = 1.5)
}

RecoverInfUABGcop.plotFrank <- function(inf.object, true_u, true_a_k, true_B, true_G, true_cop,
                                        Modeltype = "", burn.in = 100) {
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)
  n_B  <- length(true_B)
  n_G  <- length(true_G)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
    fullB.draws <- as.data.frame(inf.object$draws(variables = "B")[, 1, ])
    fullG.draws <- as.data.frame(inf.object$draws(variables = paste0("G", 1:n_G))[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
    fullB.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
    fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
    fullcop.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "c")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]
  B.draws  <- fullB.draws[thinning, , drop = FALSE]
  G.draws  <- fullG.draws[thinning, , drop = FALSE]
  cop.draws<- fullcop.draws[thinning]

  #U violin plot
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Spatial comp.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #a_k violin plot
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercept", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #B violin plot
  B_labels <- do.call(expression, lapply(1:n_B, function(i) {
    bquote(beta[.(i)])
  }))

  comp_B <- data.frame(
    value = as.vector(as.matrix(B.draws)),
    group = factor(rep(paste0("B", 1:n_B), each = nrow(B.draws)))
  )

  rfigs_B <- ggplot(comp_B, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("B", 1:n_B), levels = levels(comp_B$group)),
                                 y = true_B),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0.5, 2.5) +
    labs(x = "Regression coeff.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_B)) +
    scale_x_discrete(labels = B_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #G violin plot
  G_index_pairs <- rep(list(c(0,1), c(1,0)), ncol(G.draws)/2)

  G_labels <- do.call(expression, lapply(G_index_pairs, function(idx) {
    bquote(gamma[.(idx[1])][.(idx[2])])
  }))

  comp_G <- data.frame(
    value = as.vector(as.matrix(G.draws)),
    group = factor(rep(paste0("G", 1:n_G), each = nrow(G.draws)))
  )

  rfigs_G <- ggplot(comp_G, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("G", 1:n_G),
                                            levels = levels(comp_G$group)),
                                 y = true_G),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0, 1) +
    labs(x = "Transition prob.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_G)) +
    scale_x_discrete(labels = G_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #Copula \psi violin plot
  comp_cop <- data.frame(
    value = as.vector(cop.draws),
    group = "copParam")

  rfigs_cop <- ggplot(comp_cop, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = "copParam",
                                 y = true_cop),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0, 20) +
    labs(x = "Copula param.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_G)) +
    scale_x_discrete(labels = expression(psi)) +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  return(list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G, rfigs_cop))

#  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G, rfigs_cop)
#  print(cowplot::plot_grid(plotlist = plotlists, ncol = 5,
#                           #                           labels = c("A", "B", "C", "D","E"),
#                           labels = c("G", "", "", "",""),
#                           rel_widths = c(1.25, 1, 1, 1.5, 0.7), label_size = 25))

  # Legends
  #  add_legend(0.85, 1.15, legend = "Truth",
  #             pch = 19, col = "black",
  #             horiz = TRUE, bty = 'n', cex = 2.0)
  #  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
  #             horiz = TRUE, bty = 'n', cex = 1.5)
}

RecoverInfUABGcop.plotGauss <- function(inf.object, true_u, true_a_k, true_B, true_G, true_cop,
                                        Modeltype = "", burn.in = 100) {
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)
  n_B  <- length(true_B)
  n_G  <- length(true_G)
  n_c <- length(true_cop)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
    fullB.draws <- as.data.frame(inf.object$draws(variables = "B")[, 1, ])
    fullG.draws <- as.data.frame(inf.object$draws(variables = paste0("G", 1:n_G))[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
    fullB.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
    fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
    fullcop.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "c")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]
  B.draws  <- fullB.draws[thinning, , drop = FALSE]
  G.draws  <- fullG.draws[thinning, , drop = FALSE]
  c.draws<- fullcop.draws[thinning, , drop = FALSE]

  #U violin plot
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Spatial comp.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #a_k violin plot
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercept", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #B violin plot
  B_labels <- do.call(expression, lapply(1:n_B, function(i) {
    bquote(beta[.(i)])
  }))

  comp_B <- data.frame(
    value = as.vector(as.matrix(B.draws)),
    group = factor(rep(paste0("B", 1:n_B), each = nrow(B.draws)))
  )

  rfigs_B <- ggplot(comp_B, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("B", 1:n_B), levels = levels(comp_B$group)),
                                 y = true_B),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0.5, 2.5) +
    labs(x = "Regression coeff.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_B)) +
    scale_x_discrete(labels = B_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #G violin plot
  G_index_pairs <- rep(list(c(0,1), c(1,0)), ncol(G.draws)/2)

  G_labels <- do.call(expression, lapply(G_index_pairs, function(idx) {
    bquote(gamma[.(idx[1])][.(idx[2])])
  }))

  comp_G <- data.frame(
    value = as.vector(as.matrix(G.draws)),
    group = factor(rep(paste0("G", 1:n_G), each = nrow(G.draws)))
  )

  rfigs_G <- ggplot(comp_G, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("G", 1:n_G),
                                            levels = levels(comp_G$group)),
                                 y = true_G),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0, 1) +
    labs(x = "Transition prob.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_G)) +
    scale_x_discrete(labels = G_labels) +   #labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  #Copula violin plot
  n <- (sqrt(8*n_c + 1) - 1)/2
  c_labels <- list()
  k <- 1
  for (i in 1:n) {
    for (j in i:n) {
      c_labels[[k]] <- bquote(rho[.(i) * .(j+1)])
      k <- k + 1
    }
  }
  copulaP <- data.frame(
    value = as.vector(as.matrix(c.draws)),
    group = factor(rep(paste0("c", 1:n_c), each = nrow(c.draws)))
  )

  rfigs_c <- ggplot(copulaP, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("c", 1:n_c), levels = levels(copulaP$group)),
                                 y = true_cop),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-1.1, 1.1) +
    labs(x = "Copula param.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_c)) +
    scale_x_discrete(labels = c_labels) +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  return(list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G, rfigs_c))
#  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G, rfigs_c)
#  print(cowplot::plot_grid(plotlist = plotlists, ncol = 5,
                           #                           labels = c("A", "B", "C", "D","E"),
#                           labels = c("E", "", "", "",""),
                           #                           rel_widths = c(1, 1, 1, 0.7, 1.4), label_size = 25)),
#                           rel_widths = c(1, 0.7, 0.7, 1.2, 1.3), label_size = 25))


  # Legends
  #  add_legend(0.85, 1.15, legend = "Truth",
  #             pch = 19, col = "black",
  #             horiz = TRUE, bty = 'n', cex = 1.8)
  #  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
  #             horiz = TRUE, bty = 'n', cex = 1.5)
}

Allmodels_RS_fig<- function(all.infobjects, time=108, burn.in=1000){
  rmat<- matrix(NA, nrow = 8, ncol = time)
  smat<- matrix(NA, nrow = 8, ncol = 12)
  for(i in 1:8){
    inf.object<- all.infobjects[[i]]
    Model<- i

    if(!is.data.frame(inf.object)){
      r.draws<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
      s.draws<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
    }else{
      r.draws<- colMeans(inf.object[-(1:burn.in), startsWith(colnames(inf.object), "r")])
      s.draws<- colMeans(inf.object[-(1:burn.in), startsWith(colnames(inf.object), "s")])
    }
    rmat[i, ]<- r.draws
    smat[i, ]<- s.draws
  }
  r_names<- c("A", "B", "C", "D", "E", "F", "G", "H")
  ts_rdata <- as.data.frame(t(rmat))
  ts_sdata <- as.data.frame(t(smat))
  ts_rdata$Time <- seq.Date(from = as.Date("2011-01-01"), to = as.Date("2019-12-01"), by = "month")
  ts_sdata$Time <- seq.Date(from = as.Date("2019-01-01"), to = as.Date("2019-12-01"), by = "month")

  colnames(ts_rdata) <- c(r_names, "Time")
  colnames(ts_sdata) <- c(r_names, "Time")

  long_data <- reshape2::melt(ts_rdata, id.vars = "Time")
  library(ggplot2)
  a<- (ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
         geom_line() +
         labs(x = "Time [month/year]", y = "Trend component", color = "Model") +
         scale_x_date(date_labels = "%b %Y", date_breaks = "1 years")+
         theme(axis.title.y = element_text(size=18),
               axis.title.x = element_text(size=18),
               axis.text.x = element_text(size=16),
               axis.text.y = element_text(size=16),
               legend.title = element_text(size = 18),
               legend.text = element_text(size = 16), legend.position = "none"))

  long_data2 <- reshape2::melt(ts_sdata, id.vars = "Time")
  library(ggplot2)
  b<- (ggplot2::ggplot(data = long_data2, mapping = aes(x = Time, y = value, color = variable)) +
         geom_line() +
         labs(x = "Time [month]", y = "Seasonal component", color = "Model") +
         scale_x_date(date_labels = "%b", date_breaks = "1 months") +
         theme(axis.title.y = element_text(size=18),
               axis.title.x = element_text(size=18),
               axis.text.x = element_text(size=18),
               axis.text.y = element_text(size=16),
               legend.title = element_text(size = 22),
               legend.text = element_text(size = 20)))
  plotlists<- list(a, b)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 2, labels = c("I", "II"), label_size = 20, rel_widths = c(1.08, 1.15)))
  #export==> 23 x 9
}


relativemedian_maps<- function(all.infobjects, burn.in=1000){
  plotlists<- list()
  maxVec<- numeric(length(all.infobjects))
  minVec<- numeric(length(all.infobjects))
  for(m in 1:length(all.infobjects)){
    model_m<- all.infobjects[[m]]
    maxVec[m]<- max(colMeans(model_m[-(1:burn.in),startsWith(colnames(model_m), "u")]))
    minVec[m]<- min(colMeans(model_m[-(1:burn.in),startsWith(colnames(model_m), "u")]))
  }
  minU<- min(minVec)
  maxU<- max(maxVec)

  for(i in 1:8){
    inf.object<- all.infobjects[[i]]
    Model<- i

    fullu.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]

    #get required regions
    poly <- cshapes::cshp(date=as.Date("2019-12-31"), useGW=TRUE)
    #newpoly<- poly[poly$country_name %in% c("Austria","Belgium","Cyprus","Czech Republic",
    #                                      "Denmark","Estonia","Finland","France","German Federal Republic",
    #                                      "Greece","Hungary","Iceland","Ireland","Italy/Sardinia","Latvia",
    #                                      "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
    #                                      "Portugal","Rumania","Slovakia","Slovenia","Spain","Sweden","United Kingdom"), ]


    newpoly<- poly[poly$country_name %in% c("Albania", "Austria", "Belarus (Byelorussia)", "Belgium", "Bosnia-Herzegovina",
                                            "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland",
                                            "France", "German Federal Republic", "Greece", "Hungary", "Iceland",
                                            "Ireland", "Italy/Sardinia", "Kosovo", "Latvia", "Lithuania", "Luxembourg", "Malta",
                                            "Macedonia (FYROM/North Macedonia)", "Moldova", "Montenegro", "Netherlands", "Norway",
                                            "Poland", "Portugal", "Rumania", "Russia (Soviet Union)", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
                                            "Switzerland", "Ukraine", "United Kingdom"), ]

    Europe<- c("Austria","Belgium","Cyprus","Czech Republic",
               "Denmark","Estonia","Finland","France","German Federal Republic",
               "Greece","Hungary","Ireland","Italy/Sardinia","Latvia",
               "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
               "Portugal","Slovakia","Slovenia","Spain","Sweden","United Kingdom",
               "Albania", "Belarus (Byelorussia)", "Bosnia-Herzegovina", "Bulgaria",
               "Croatia", "Iceland", "Kosovo", "Rumania", "Russia (Soviet Union)",
               "Moldova", "Montenegro", "Macedonia (FYROM/North Macedonia)",
               "Serbia", "Switzerland", "Ukraine")

    sortedpoly <- newpoly[order(newpoly$country_name), ]
    #sortedpoly$country_name<- c("Austria","Belgium","Cyprus","Czechia",
    #                          "Denmark","Estonia","Finland","France","Germany",
    #                          "Greece","Hungary","Iceland","Ireland","Italy","Latvia",
    #                          "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
    #                          "Portugal","Romania","Slovakia","Slovenia","Spain","Sweden","United Kingdom")

    u.draws<- colMeans(fullu.draws)
    #expUi<- c(exp(u.draws), rep(NA, 15))
    expUi<- c(u.draws, rep(NA, 15))
    uidf<- data.frame(expUi = expUi, countryname = Europe)
    sorteduidf<- uidf[order(uidf$countryname), ]
    rrDF<- data.frame(rr=sorteduidf$expUi, gwcode=sortedpoly$gwcode, index=c(1:nrow(uidf)), name=sorteduidf$countryname)

    library(sf)
    sortedpoly<- st_make_valid(newpoly)
    bbox <- st_bbox(c(xmin = -29, ymin = 30, xmax = 40.17875, ymax = 71.15471), crs = st_crs(sortedpoly))
    sortedpoly<- st_crop(sortedpoly, bbox)
    shapefile <- ggplot2::fortify(sortedpoly, region = 'gwcode')

    shp_rrDF<- sp::merge(shapefile, rrDF,
                         by.x="gwcode",
                         by.y="gwcode",
                         all.x=F)

    library(ggplot2)
    library(viridis)
    if(Model == 8){
      rfig<- (ggplot(data = shp_rrDF) +
                geom_sf(aes(fill = rr)) +
                #geom_sf_text(aes(label = country_name), size = 3) +
                #scale_fill_viridis(option = "turbo", direction = 1, alpha=1, begin=0.6, end=1, na.value = "lightgrey") +  # Reverse the color scale
                scale_fill_gradient2(low = "lightblue", mid = "blue", high = "red", midpoint = 0, limits = c(minU, maxU), oob = scales::squish, na.value = "lightgrey", labels = ~ format(round(exp(.),1), nsmall=1)) +
                coord_sf() + theme_void() +
                ggtitle("") + theme(plot.title = element_text(hjust = 0.55)) +
                labs(fill = expression(Exp(u[i]))) +
                theme(legend.title = element_text(size = 18),
                      legend.text = element_text(size = 16)))
    }else{
      rfig<- (ggplot(data = shp_rrDF) +
                geom_sf(aes(fill = rr)) +
                #geom_sf_text(aes(label = country_name), size = 3) +
                #scale_fill_viridis(option = "turbo", direction = 1, alpha=1, begin=0.6, end=1, na.value = "lightgrey") +  # Reverse the color scale
                scale_fill_gradient2(low = "lightblue", mid = "blue", high = "red", midpoint = 0, limits = c(minU, maxU), oob = scales::squish, na.value = "lightgrey") +
                coord_sf() + theme_void() +
                ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) +
                labs(fill = "Median relative risk") +
                theme(legend.position = "none")) #+
      #theme_minimal()
    }
    plotlists[[Model]]<- rfig
  }
  row_1<- cowplot::plot_grid(plotlist = plotlists[1:4], ncol = 4, labels = c("A", "B", "C", "D"), label_size = 17)
  row_2<- cowplot::plot_grid(plotlist = plotlists[5:8], ncol = 4, labels = c("E", "F", "G", "H"), label_size = 17, rel_widths = c(1.16, 1.16, 1.16, 1.60))
  print(cowplot::plot_grid(row_1, row_2, nrow = 2))
  #export ==> 23 x 9
}


heat_maps<- function(outbreakProb_arraylist, location_names, pdfname=NULL){
  countries<- location_names
  time<- ncol(outbreakProb_arraylist[[1]][,,1])
  ndept<- nrow(outbreakProb_arraylist[[1]][,,1])
  nstrain<- dim(outbreakProb_arraylist[[1]])[3]
  n_models<- length(outbreakProb_arraylist)
  strain_names<- c("NEIMENI_B", "NEIMENI_W", "NEIMENI_Y", "NEIMENI_C")

  pdf(paste0(pdfname,".pdf"), paper="special", width=24,height=46, pointsize=14)
  par(mfrow=c(n_models+1,4))

  for(i in 1:n_models){
    outProb_model_i<- outbreakProb_arraylist[[i]]
    for(k in 1:nstrain){
      mean_Xit<- outProb_model_i[,,k]
      mean_Xit<- mean_Xit[ndept:1, ]
      par(mar = c(4, 7.5, 4, 1))
      if(i==1){
        image(x=1:time, y=1:ndept, t(mean_Xit), main = strain_names[k], axes=F, ylab = "", xlab = "Time [month/year]", cex.lab=1.8, zlim=c(0,1), cex.main=2.5)
      }else{
       image(x=1:time, y=1:ndept, t(mean_Xit), main = "", axes=F, ylab = "", xlab = "Time [month/year]", cex.lab=1.8, zlim=c(0,1))
      }
      #custom Y-axis
      axis(2, at=seq(1, length(countries), length.out=length(countries)), labels=rev(countries), lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.40)
      #custom X-axis
      years<- 2011:2019
      axis(1, at = seq(1, time, by = 12), labels = years, cex.axis = 1.8)
      if(k==1) {legendary::labelFig(LETTERS[i], adj = c(-0.15, 0.05), font=2, cex=2.5)}
    }
  }

  par(mar = c(5, 19, 7, 17))
  #c(bottom, left, top, right) ==> specification order

  zseq <- seq(0, 1, length.out = 101)
  xseq <- c(0, 1)
  zmat <- matrix(zseq[-1], nrow = 1)

  image(x = xseq, y = zseq, z = zmat,
        axes = FALSE, xlab = "", ylab = "", main = "")

  axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, cex.axis=2.0)
  mtext("Posterior probability of outbreak", side = 4, line = 5.0, cex = 2.0)
  box()
  dev.off()
}


# Super-impose true values
simulation_heat_maps<- function(Truth_arraylist, matrix_arraylist, Outbreaktype=""){
  time<- ncol(matrix_arraylist[[1]][,,1])
  ndept<- nrow(matrix_arraylist[[1]][,,1])
  nstrain<- dim(matrix_arraylist[[1]])[3]
  n_models<- length(matrix_arraylist)
  pdf(paste0("Alloutbreaksperstrain",Outbreaktype,".pdf"), paper="special", width=36,height=46, pointsize=14)
  par(mfrow=c(n_models+1,5))

  for(m in 1:n_models){
    matrix_array<- matrix_arraylist[[m]]
    Truth_array<- Truth_arraylist[[m]]
  for(i in 1:nstrain){
    Truth<- Truth_array[,,i]
    smallTruth<- Truth[c(1,3,5,7,9), ]
    bigTruth<- Truth[c(2,4,6,8), ]
    bigsmallTruth<- Truth[c(2,4,6,8,1,3,5,7,9), ]

    X_it<- matrix_array[,,i]
    smallxit<- X_it[c(1,3,5,7,9), ]
    bigxit<- X_it[c(2,4,6,8), ]
    bigsmallxit<- X_it[c(2,4,6,8,1,3,5,7,9), ]
    par(mar = c(4, 7.5, 4, 1))
    if(m==1){
      image(x=1:time, y=1:ndept, t(bigsmallxit), zlim = c(0,1),  main =paste0("Strain ", i), axes=F, ylab="spatial location", xlab="Time [month]", cex.lab=1.80, cex.main=3.0)
    }else{
      image(x=1:time, y=1:ndept, t(bigsmallxit), zlim = c(0,1),  main ="", axes=F, ylab="spatial location", xlab="Time [month]", cex.lab=1.80, cex.main=2.5)
    }
    abline(h=4.5, col="black", lty=2)
    #custom Y-axis
    axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    #custom X-axis
    axis(1, cex.axis = 1.8)
    # Build grid of row/col indices
    grid <- expand.grid(x = seq(nrow(bigsmallTruth)),
                        y = seq(ncol(bigsmallTruth)))
    out <- transform(grid, z = bigsmallTruth[as.matrix(grid)])

    # Select true outbreak cells
    outbreakcell <- out$z == 1

    # Overlay truth as colored dots at cell centers
    points(out$y[outbreakcell],   # time axis
           out$x[outbreakcell],   # location axis
           pch = 16, cex = 2.0, col = "magenta")
    if(i==1) {legendary::labelFig(LETTERS[m], adj = c(-0.15, 0.05), font=2, cex=3.0)}
    }
  }
    par(mar = c(3, 21, 3, 21))
    #c(bottom, left, top, right) ==> specification order

    zseq <- seq(0, 1, length.out = 101)
    xseq <- c(0, 1)
    zmat <- matrix(zseq[-1], nrow = 1)

    image(x = xseq, y = zseq, z = zmat,
          axes = FALSE, xlab = "", ylab = "", main = "")

    axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, cex.axis=2.0)
    mtext("Posterior probability of outbreak", side = 4, line = 5.0, cex = 2.0)

    add_legend(-0.60, -0.90, legend=c("Truth", "Small cities", "Large cities"), lty=c(NA, 1, 1),
               pch=c(16, NA, NA), col=c("magenta", "blue", "red"),
               horiz=TRUE, bty='n', cex=5.5)

  #  add_legend(0.50, -0.4, legend=substitute(paste(bold(Outbreaktype))),
  #             col="black",
  #             horiz=TRUE, bty='n', cex=3.0)
#  add_legend(0.58, -0.30, legend=substitute(paste(bold(Outbreaktype))),
#             col="black",
#             horiz=TRUE, bty='n', cex=5.0)

  dev.off()
}
