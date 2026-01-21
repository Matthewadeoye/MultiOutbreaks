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


multitypeFig2 <- function(array.object, names = NULL){
  plotlists<- list()
  nstrain<- dim(array.object)[3]
  maxY<- max(array.object, na.rm = T)
  for(i in 1:nstrain){
    spatdata <- array.object[,,i]

  ts_spatdata <- as.data.frame(t(spatdata))
  ts_spatdata$Time <- seq.Date(from = as.Date("2013-01-01"), to = as.Date("2019-12-01"), by = "month")
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
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 16),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16))
  }else{
  a <- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
    geom_line() +
    ylim(0, maxY) +
    labs(x = "Time [month/year]", y = "Case counts", color = "Location") +
    guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
    scale_x_date(date_labels = "%b %Y", date_breaks = "1 years") +
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.position = "none")
  }
  plotlists[[i]]<- a
  }
  row_1<- cowplot::plot_grid(plotlist = plotlists[1:3], ncol = 3, labels = c("A", "B", "C"), label_size = 17)
  row_2<- cowplot::plot_grid(plotlist = plotlists[4], ncol = 2, labels = c("D"), label_size = 17, rel_widths = c(0.9, 1))
  finalplot<- cowplot::plot_grid(row_1, row_2, nrow = 2)
  print(finalplot)
  #export ==> 23 x 9
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
  pdf("AllRS2.pdf", paper="a4", width=12,height=12, pointsize=12)

  par(mfrow=c(4,2))
  # Modeltype<- c("No outbreak", "Independent 1", "Independent 2", "Gaussian copula-dependent 1",
  #                "Gaussian copula-dependent 2", "Frank copula-dependent 1",
  #               "Frank copula-dependent 2", "General-dependent")
  #alphabets<- c("A", "C", "E", "G")
  alphabets<- c("B", "D", "F", "H")
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
  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.8)
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

  plotlists <- list(rfigs_u, rfigs_ak)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 2,
                           labels = c("A", ""),
                           rel_widths = c(1.25, 1), label_size = 25))

  # Legends
  #  add_legend(0.85, 1.15, legend = "Truth",
  #             pch = 19, col = "black",
  #             horiz = TRUE, bty = 'n', cex = 1.8)
  #  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
  #             horiz = TRUE, bty = 'n', cex = 1.5)
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
  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 3,
                           #                           labels = c("A", "B", "C"),
                           labels = c("H", "", ""),
                           rel_widths = c(1.25, 1, 1), label_size = 25))

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

  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 4,
                           #                           labels = c("A", "B", "C", "D","E"),
                           labels = c("C", "", "", ""),
                           rel_widths = c(1.25, 1, 1, 1.25), label_size = 25))

  # Legends
  add_legend(0.85, 1.15, legend = "Truth",
             pch = 19, col = "black",
             horiz = TRUE, bty = 'n', cex = 2.0)
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

  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G, rfigs_cop)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 5,
                           #                           labels = c("A", "B", "C", "D","E"),
                           labels = c("G", "", "", "",""),
                           rel_widths = c(1.25, 1, 1, 1.5, 0.7), label_size = 25))

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

  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G, rfigs_c)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 5,
                           #                           labels = c("A", "B", "C", "D","E"),
                           labels = c("E", "", "", "",""),
                           #                           rel_widths = c(1, 1, 1, 0.7, 1.4), label_size = 25)),
                           rel_widths = c(1, 0.7, 0.7, 1.2, 1.3), label_size = 25))


  # Legends
  #  add_legend(0.85, 1.15, legend = "Truth",
  #             pch = 19, col = "black",
  #             horiz = TRUE, bty = 'n', cex = 1.8)
  #  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
  #             horiz = TRUE, bty = 'n', cex = 1.5)
}
