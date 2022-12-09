plot.Dirichlet <- function(a1, a2=NULL, a3=NULL, main='Dirichlet Distribution', ...) {

  if (length(a1) == 3) {
    a2 <- a1[2]
    a3 <- a1[3]
    a1 <- a1[1]
  }

  x1 <- x2 <- seq(0.01, .99, by=.01)

  ddirichlet <- function(x1, x2){
    term1 <- gamma(a1+a2+a3)/(gamma(a1)*gamma(a2)*gamma(a3))
    term2 <- x1^(a1-1)*x2^(a2-1)*(1-x1-x2)^(a3-1)
    term3 <- (x1 + x2 < 1)
    term1*term2*term3
  }

  z <- outer(x1, x2, ddirichlet)
  z[z<=0] <- NA

  p <- persp(x1, x2, z,
             col = c("lightblue"),
             theta = 80,
             phi = 30,
             r = 50,
             d = 0.1,
             expand = 0.5,
             ltheta = 10,
             lphi = 180,
             shade = 0.75,
             ticktype = "detailed",
             nticks = 5,
             border=T,
             zlim = if(length(na.omit(unique(as.vector(z)))) < 2) {
               c(0,2.1)
             } else {
               range(z, na.rm = TRUE)
             }, main=main, ...
  )

  text(trans3d(1.05,-0.05,0,p), expression(alpha[1]))
  text(trans3d(0.01,1.1,0,p), expression(alpha[2]))
  text(trans3d(-0.1,-0.01,0,p), expression(alpha[3]))
}

plot.signal <- function(signals, ylabel='Mean signal', type='enhancer.signal', fun=colMeans, smoothness=25) {

  half.window <- nrow(signals[[1]][[type]])/2
  enhancer.means <- as.data.frame(do.call(rbind, lapply(signals, function(mod) fun(t(as.matrix(mod[[type]]))) )))

  colnames(enhancer.means) <- c(-(half.window-1):half.window)
  enhancer.means$Data <- rownames(enhancer.means)
  enhancer.means <- melt(enhancer.means, id='Data')

  ggplot(enhancer.means, aes(x=variable, y=value, color=Data, group=Data)) +
    geom_line(size=.4, alpha=0.6) + geom_point(size=.8) +
    geom_smooth(size=1, method='gam', formula=y ~ s(x, k=smoothness, bs='cs'))  +
    xlab('Position') + ylab(ylabel) +
    scale_colour_brewer(palette="Paired") +
    scale_x_discrete(breaks = seq.int(-half.window, half.window, 100)) + theme_minimal()

}

plot.region.overlap.heatmap <- function(loci.gr,
                                        labels,
                                        annotation.grlist,
                                        resize.loci=50) {

  library(RColorBrewer)
  library(gplots)

  if(is.null(names(annotation.grlist)))
    stop('Annotation list must have names.')

  if(length(loci.gr)!=length(labels))
    stop('Loci must have the same length as labels.')

  loci.gr.split <- split(loci.gr, labels)

  overlapMatrix <- matrix(0, length(annotation.grlist), length(loci.gr.split))
  for (annot in seq_along(annotation.grlist)) {
    for (clust in seq_along(loci.gr.split)) {
      if(resize.loci)
        lgs <- resize(loci.gr.split[[clust]], resize.loci,'center')
      else
        lgs <- loci.gr.split[[clust]]

      #TODO: countOverlap(GRangesList, GRangesList) or something would better without loops
      overlapMatrix[annot, clust] <- sum(countOverlaps(annotation.grlist[[annot]], lgs))
    }
  }

  row.names(overlapMatrix) <- names(annotation.grlist)
  colnames(overlapMatrix) <- paste0('cluster', seq_along(loci.gr.split))

  heatmap.2(overlapMatrix,
            cexRow=0.55,
            cexCol=1,
            margins=c(3,6),
            col=brewer.pal(9, 'YlGn'))
}

plot.circular.ideogram <- function(gr, labels, ...) {

  gr <- keepStandardChromosomes(gr, species = 'Homo sapiens', style = 'UCSC')
  base <- GRanges(seqnames=seqlevels(gr), ranges=IRanges(start = 1, end=seqlengths(gr)), seqlengths = seqlengths(gr))
  values(gr)$Label <- factor(labels)
  values(gr)$y <- jitter(rep(1, length(gr)))

  gg <- ggbio() +
    circle(base, geom='ideo', fill='gray70') +
    circle(base, geom='scale', size=1) +
    circle(base, geom='text', aes(label=seqnames), vjust=-1, size=1) +
    circle(gr,
           geom='point',
           aes(colour=Label, y=y),
           size=1,
           alpha=0.6) +
    facet_wrap(~Label) +
    ggtitle('Circular ideogram of clusters')

  if (length(unique(labels)) <= 12)
    gg <- gg + scale_color_brewer(palette = 'Paired')

  if (length(list(...)) > 0)
    for (l in list(...))
      gg <- gg + l

  gg
}


#' DO NOT USE THIS This functions plots the means of the clusters before or after clustering'
#' @param data
#' @param labels
#' @param ...
#' @param fun
#' @param rotatex
#' @param title
#' @param smooth
#' @param xbreaks
#' @param facet
#' @param smoothness
#' @param scales
#' @param ncol
#'
#' @return
#' @export
#'
#' @examples
plot.clusters <- function(data, labels, ..., fun = mean, rotatex = T,
                          title='Cluster mean signals',smooth='gam', xbreaks=10,
                          facet='grid', smoothness=25, scales=NULL, ncol=NULL) {

  linesize<-ifelse(smooth, .2, .6)

  if (missing(data)) {
    S <- nrow(fitted(labels)[[1]])
    M <- length(fitted(labels))
    data <- melt(fitted(labels),
                 varnames = c('Bin', 'Component'),
                 value.name='Value')
    names(data)[which(names(data) == 'L1')] <- 'Datatype'
    title <- 'Dirichlet parameters'

    if (is.null(scales)) {
      if (!any(diff(sapply(fitted(labels), nrow)))) scales <- 'fixed' else scales <- 'free_x'
    }

    } else {

    if(!is.vector(labels)) labels <- mixture(labels, assign=T)
    M <- length(data)
    S <- ncol(data[[1]])

    if (is.null(scales)) {
      if (!any(diff(sapply(data, ncol)))) scales <- 'fixed' else scales <- 'free_x'
    }

    data <- ldply(data, function(d) {
      d <- as.data.frame(d)
      d$id <- rownames(d)
      d$Component <- factor(labels)
      d <- melt(d,measure=1:(ncol(d)-2), variable.name = 'Bin', value.name = 'Value')
      d$Bin <- as.numeric(d$Bin)
      d <- ddply(d, c('Component', 'Bin'), here(summarize), Value=fun(Value))
      d}, .id = 'Datatype')
  }

  gg <- ggplot(data, aes(Bin, Value, group=Component)) +
    geom_point(size=1, shape=1, aes(color=Component)) +
    labs(x='Position', y='Clusters', title=title) +
    guides(color=F) +
    scale_x_discrete(breaks=seq(0, S, xbreaks)) +
    theme_minimal()

  if (is.character(smooth) && smooth == 'gam') {
    #gg <- gg + geom_smooth(aes(group=labels), method='loess', span=.1)
    gg <- gg + geom_smooth(size=.6, method='gam', formula=y~s(x, k=smoothness, bs='cs'))
  }
  else if (is.character(smooth) && smooth == 'auto')
    gg <- gg + geom_smooth(size=.6)
  else
    gg <- gg + geom_line(size=linesize, aes(color=Component))

  if(facet == 'grid')
    gg <- gg + facet_grid(Component~Datatype, scales=scales)
  else if(facet == 'wrap') {
    if (scales == 'fixed') scales <- 'free_y' else scales <- 'free'
    if(is.null(ncol)) ncol<-M
    gg <- gg + facet_wrap(~Component+Datatype, scales=scales, ncol=ncol)
  } else {
    if (scales == 'fixed') scales <- 'free_y' else scales <- 'free'
    if(is.null(ncol)) ncol<-M
    gg <- gg + facet_wrap(~Datatype, scales=scales, ncol=ncol)
  }

  if (length(list(...)) > 0)
    for (l in list(...))
      gg <- gg + l

  if (rotatex)
    gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  gg
}

#' plot.EM plots EM diagnostics
#'
#' @param fit
#' @param smoothness.scale
#' @param EM_lambda_optim_message
#' @param skip
#' @return
#' @export
#'
#' @examples
# plot.EM <- function(fit, smoothness.scale='free', k=10) {
#
#   EM.diagnostics <- fit@EM.diagnostics
#   K <- ncol(mixture(fit))
#   #plot hkm values
#   comp.labels=paste0("Cluster ", seq(1,k,1))
#   names(comp.labels)=seq(1,k,1)
#
#   hkm.plot <- ggplot(EM.diagnostics,
#                      aes(seq_along(hkm), hkm, color=EM.iter)) +
#     geom_line(aes(group=Datatype)) +
#     geom_point(shape=1) +
#     facet_wrap(~Datatype+Component, scale=smoothness.scale, ncol=K, labeller=labeller(Component=comp.labels)) +
#     labs(title='Numerical optimization iterations vs. regularization term', x='Numerical optimization iterations', y=expression(h[k]^(m))) +
# #     scale_color_gradientn('EM iterations', colours=RColorBrewer::brewer.pal(n=9, name='YlOrRd'))
# #     scale_color_brewer(palette='Blues') +
#     scale_color_continuous('EM iterations') +
#     theme_minimal()
#
#
#   #plot number of num.opt. iterations
#   nop <- ggplot(EM.diagnostics, aes(EM.iter, NO.iter.count)) +
#     geom_line(aes(color=factor(Datatype))) +
#     geom_point(aes(color=factor(Datatype))) +
#     #geom_smooth(method='gam', formula=y~s(x, k=k, bs='cs')) +
#     facet_grid(Datatype~Component, labeller=labeller(Component=comp.labels)) +
#     guides(color=F) +
#     labs(title='EM iterations vs. Numerical optimization iterations', x='EM iterations', y='Numerical optimization iterations') +
#     theme_minimal()
#
#   #create a grid layout and print ggplots on it
#   grid.newpage()
#   pushViewport(viewport(layout = grid.layout(2, 1)))
#   vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#   print(hkm.plot, vp = vplayout(1, 1))
#   print(nop, vp = vplayout(2, 1))
#
# }

plot.EM <- function(fit,  smoothness.scale='free', skip=1, plot.det=TRUE, plot.gradient=TRUE) {


  K=ncol(fit@fit$Estimate[[1]])
  M=length(fit@Data)



  comp.labels=paste0("Cluster ", seq(1,K,1))
  names(comp.labels)=seq(1,K,1)


  nll.plot <- ggplot(fit@EM.diagnostics,
                     aes(EM.iter, nll, )) +
    geom_line(aes(group=Datatype)) +
    geom_point(shape=1) +
    facet_wrap(~Datatype+Component, scale=smoothness.scale, ncol=K, labeller=labeller(Component=comp.labels)) +
    labs(title='EM iterations vs. BFGS neg.log.lik', x='EM iterations',
         y="nll" ) +
    #     scale_color_gradientn('EM iterations', colours=RColorBrewer::brewer.pal(n=9, name='YlOrRd'))
    #     scale_color_brewer(palette='Blues') +
    scale_color_continuous('EM iterations') +
    theme_minimal()

  if(plot.det==TRUE){
  detHes.plot <- ggplot(fit@EM.diagnostics,
                        aes(EM.iter, detH, )) +
    geom_line(aes(group=Datatype)) +
    geom_point(shape=1) +
    facet_wrap(~Datatype+Component, scale=smoothness.scale, ncol=K, labeller=labeller(Component=comp.labels)) +
    labs(title='EM iterations vs. BFGS Hessian determinant', x='EM iterations',
         y="Hessian determinant" ) +
    #     scale_color_gradientn('EM iterations', colours=RColorBrewer::brewer.pal(n=9, name='YlOrRd'))
    #     scale_color_brewer(palette='Blues') +
    scale_color_continuous('EM iterations') +
    theme_minimal()
  }


  hkm.plot <- ggplot(fit@EM.diagnostics,
                     aes(seq_along(hkm), hkm, color=EM.iter)) +
    geom_line(aes(group=Datatype)) +
    geom_point(shape=1) +
    facet_wrap(~Datatype+Component, scale=smoothness.scale, ncol=K, labeller=labeller(Component=comp.labels)) +
    labs(title='Numerical optimization iterations vs. regularization term', x='Numerical optimization iterations', y=expression(h[k]^(m))) +
    #     scale_color_gradientn('EM iterations', colours=RColorBrewer::brewer.pal(n=9, name='YlOrRd'))
    #     scale_color_brewer(palette='Blues') +
    scale_color_continuous('EM iterations') +
    theme_minimal()

  if(plot.gradient==TRUE){
  gradient.plot <- ggplot(fit@EM.diagnostics,
                          aes(seq_along(gradient), gradient, color=EM.iter)) +
    geom_line(aes(group=Datatype)) +
    geom_point(shape=1) +
    facet_wrap(~Datatype+Component, scale=smoothness.scale, ncol=K, labeller=labeller(Component=comp.labels)) +
    labs(title='Numerical optimization iterations vs. gradient norm', x='Numerical optimization iterations', y="gradient norm") +
    #     scale_color_gradientn('EM iterations', colours=RColorBrewer::brewer.pal(n=9, name='YlOrRd'))
    #     scale_color_brewer(palette='Blues') +
    scale_color_continuous('EM iterations') +
    theme_minimal()

  }

  #plot number of num.opt. iterations
  nop <- ggplot(fit@EM.diagnostics, aes(EM.iter, NO.iter.count)) +
    geom_line(aes(color=factor(Datatype))) +
    geom_point(aes(color=factor(Datatype))) +
    #geom_smooth(method='gam', formula=y~s(x, k=k, bs='cs')) +
    facet_grid(Datatype~Component, labeller=labeller(Component=comp.labels)) +
    guides(color=F) +
    labs(title='EM iterations vs. Numerical optimization iterations', x='EM iterations', y='Numerical optimization iterations') +
    theme_minimal()

  #create a grid layout and print ggplots on it
  grid.newpage()
  if(plot.det==TRUE && plot.gradient==TRUE){
    pushViewport(viewport(layout = grid.layout(5, 1)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    print(hkm.plot, vp = vplayout(1, 1))
    print(gradient.plot, vp = vplayout(2, 1))
    print(detHes.plot, vp = vplayout(3, 1))
    print(nll.plot, vp = vplayout(4, 1))
    print(nop, vp = vplayout(5, 1))
  }
  if(plot.det==FALSE && plot.gradient==FALSE){
    pushViewport(viewport(layout = grid.layout(3, 1)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    print(hkm.plot, vp = vplayout(1, 1))
    print(nll.plot, vp = vplayout(2, 1))
    print(nop, vp = vplayout(3, 1))
    }

  if(plot.det==FALSE && plot.gradient==TRUE){
    pushViewport(viewport(layout = grid.layout(4, 1)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    print(hkm.plot, vp = vplayout(1, 1))
    print(gradient.plot, vp = vplayout(2, 1))
    print(nll.plot, vp = vplayout(3, 1))
    print(nop, vp = vplayout(4, 1))
  }



}

plot.label.distribution <- function(fits,
                                    groupColumn,
                                    alphaLines = 1/10,
                                    useSplines = T,
                                    showPoints = T,
                                    scale='uniminmax', ...) {

  labels <- as.data.frame(sapply(fits, mixture, assign=T)) #get hard labels
  labels <- as.data.frame(lapply(labels, factor))
  N <- nrow(labels)
  L <- ncol(labels)

  if (missing(groupColumn)) groupColumn <- L
  labels[,(L+1):(2*L)] <- lapply(labels, function(l)jitter(as.numeric(l))) #add jitter
  colnames(labels)[(L+1):(2*L)] <- as.character(seq_len(L))

  if ('useSplines' %in% names(formals(ggparcoord))) {

    ggp <- ggparcoord(labels, L+seq(L), scale=scale, groupColumn=groupColumn,
                      alphaLines = alphaLines,
                      showPoints = showPoints,
                      useSplines = useSplines,
                      ...)
    if(useSplines)
      ggp <- ggp + scale_x_continuous(breaks=seq(L))
  } else {
    ggp <- ggparcoord(labels, L+seq(L), scale=scale, groupColumn=groupColumn,
                      alphaLines = alphaLines,
                      showPoints = showPoints,
                      ...)
  }

  ggp <- ggp +
    #     scale_color_brewer(palette = 'Paired') +
    xlab('Number of clusters') +
    ylab('Memberships of genomic loci') +
    guides(colour = F) +
    theme_minimal() +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank())
  ggp
}

plot.clustering.comparison <- function(fits1, fits2, num, method='exact') {

  K <- length(fits1)

  labels1 <- data.frame(lapply(fits1, mixture, assign=T))
  labels1 <- data.frame(lapply(labels1, factor))
  names(labels1) <- seq(K)

  labels2 <- data.frame(lapply(fits2, mixture, assign=T))
  labels2 <- data.frame(lapply(labels2, factor))
  names(labels2) <- seq(K)

  for (k in seq(K)) {
    #let's first match clusters
    m <- matchClasses(table(labels1[, k], labels2[, k]), method=method)
    labels1[, k] <- m[labels1[, k]]
  }

  l1 <- labels1[, num]
  l2 <- labels2[, num]
  count <- as.data.frame(table(l1, l2))
  count$Freq <- sqrt(count$Freq / max(count$Freq))
  ggplot(count, aes(l1, l2, height = Freq, width = Freq)) +
    geom_tile(colour = "white", fill='darkgray') +
    theme_minimal() +
    theme(aspect.ratio = 1)

}





#' DO NOT USE THIS This functions draws the heatmaps after clustering
#'
#' @param cl can be DMN object or just the data
#' @param labels
#' @param show.param
#'
#' @return
#' @export
#'
#' @examples
plot.heatmap <- function(cl, labels, show.param=T) {

  #cl=fits[[2]]
  #labels=mixture(fits[[2]],assign=TRUE)


  if (class(cl) == "DMN")
    data <- cl@Data
  else {
    data <- cl
    #show.param <- F
  }


  if (missing(labels))
    labels <- mixture(cl, assign=T)

  if (!any(diff(sapply(data, ncol)))) scales <- 'fixed' else scales <- 'free_x'

  S <- sapply(data, ncol)
  K <- length(unique(labels))
  M <- length(data)
  N <- length(labels)

  datatypes <- names(data)

  data <- ldply(data, function(d) {
    d <- as.data.frame(d)
    d$id <- factor(1:nrow(d))
    d$Component <- factor(labels)
    d <- melt(d,measure=1:(ncol(d)-2), variable.name = 'Bin', value.name = 'Reads')
    d$Bin <- as.numeric(d$Bin)
    #      d <- ddply(d, c('Component', 'Bin'), here(summarize), Reads=sqrt(Reads))
    d}, .id = 'Datatype')

  data$Reads <- sqrt(data$Reads)

  title.margin <- 0.3
  param.height <- 0.7
  data.height <- 1.3

  if (show.param) {
    params <- fitted(cl)

    params <- melt(params,
                   varnames = c('Bin', 'Component'),
                   value.name='Reads')
    names(params)[which(names(params) == 'L1')] <- 'Datatype'

    data.height <- as.vector((table(labels)/N)*(data.height*K))
    heights <- as.vector(t(data.frame(param.height, data.height)))

    heights[1] <- heights[1] + title.margin
    heights[2:K*2] <- heights[2:K*2] - title.margin/((K*2)-1)

    heights <- unit(heights, "null")
  } else {
    heights <- as.vector((table(labels)/N)*(K))
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(ifelse(show.param, K*2, K),
                                             M,
                                             heights=heights)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

  for (k in seq_len(K)) {
    comp.label <- sort(unique(labels))[k]

    for (m in seq_len(M)) {

      if (show.param) {
        margin <- c(-2, 2, -10, -2)

        if (k %% K == 1) margin[1] <- margin[1] + 5
        if (m %% M == 1) margin[4] <- margin[4] + 2

        margin <- unit(margin,"mm")


        ggparams <- ggplot(subset(params, Component==comp.label & Datatype==datatypes[m]),
                           aes(Bin, Reads, group=Component)) +
          #geom_point(size=1, shape=1, aes(color=Component)) +
          geom_area(fill='brown') +
          labs(x='', y='') +
          guides(color=F) +
          scale_x_discrete(breaks=seq(0, S[m], 10), expand=c(0,0)) +
          theme_minimal() +
          # theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
          theme(plot.margin=margin,
                axis.ticks.margin=unit(c(0,0),"mm"),
                panel.margin=unit(c(0,0,0,0),"mm"))

        if (k %% K == 1) ggparams <- ggparams + ggtitle(datatypes[m])

        print(ggparams, vp = vplayout((k-1)*2+1, m))

      }
      margin <- c(0, 2, -8, -2)

       if (m %% M == 1) {
         margin[4] <- margin[4] + 2
         yaxis <- paste('Cluster', comp.label)
       } else {
         yaxis  <- ''
       }

      margin <- unit(margin,"mm")

      gg <- ggplot(subset(data, Component==comp.label & Datatype==datatypes[m]),
                   aes(Bin, id, group=Component)) +
        geom_raster(aes(fill=Reads)) +
        labs(x='', y=yaxis) +
        guides(fill=F) +
        scale_x_continuous(expand= c(0,0))+
        scale_fill_gradientn(colours = RColorBrewer::brewer.pal(name = 'YlOrRd', 9)) +
        theme_minimal() +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
        theme(plot.margin=margin,
              axis.ticks.margin=unit(c(0,0,0,0),"mm"),
              panel.margin=unit(c(0,0,0,0),"mm"),
              plot.background=element_blank(),
              #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )
#         scale_x_discrete(expand = c(0, 0), breaks=seq(0, S[m], 10)) + scale_y_discrete(expand = c(0, 0))


      if (show.param)
        print(gg, vp = vplayout((k-1)*2+2, m))
      else
        print(gg, vp = vplayout(k, m))
    }
  }
}


heatmapdmn <-
  function(count, fit1, fitN, ntaxa=30, ..., transform=sqrt,
           lblwidth=NULL, col=.gradient, order.rows=TRUE,
           plot.weights=FALSE, simple.labels=FALSE)
  {
    if(is.matrix(count)) count <- list(Datatype=count)

    count <- do.call(cbind, count)
    p1 <- do.call(rbind, fitted(fit1, scale=TRUE))
    pN <- do.call(rbind, fitted(fitN, scale=TRUE))

    datatypes <- names(fitted(fitN))

    diff <- rowSums(abs(pN - as.vector(p1)))
    if (order.rows)
      taxa <- rev(head(order(diff, decreasing=TRUE), ntaxa))
    else
      taxa <- rev(head(1:length(diff), ntaxa))

    pN <- pN[taxa, , drop=F]

    if (is.null(lblwidth)) lblwidth <- .2 * nrow(count)

    cl <- mixture(fitN, assign=TRUE)
    #     thetas <- round(mixturewt(fitN)$theta)
    pis <- round(mixturewt(fitN)$pi, 2)

    ncl <- length(unique(cl))
    nms <- names(cl)
    grp <- factor(cl, levels=as.character(seq(1, ncl)))
    idx <- split(nms, grp)
    ## 2 * ncl + 1 (for labels) panels
    mwd <- .15 * length(cl) / ncl    # 'm's take up 15% of total width
    wd <- c(unlist(Map(c, lapply(idx, length), mwd), use.names=FALSE), lblwidth)
    layout(matrix(seq(1, 2 * ncl + 1), nrow=1), widths=wd)
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)
    topmargin <- ifelse(plot.weights, 4, 1)
    par(mar=c(1, 0, topmargin, 0))
    for (i in seq_along(idx)) {
      image(transform(count[idx[[i]], taxa, drop=FALSE]),
            col=col, xaxt="n", yaxt="n")
      if (plot.weights)
        #             axis(3, c(0.5), labels=paste0('Group', i, ' \n\u03b8:',
        #                                           thetas[i], ' \u03c0:', pis[i]))
        axis(3, c(0.5), labels=paste0('Group', i, ' \u03c0:', pis[i]))

      image(t(transform(pN[, i, drop=FALSE])),
            col=col, xaxt="n", yaxt="n")
    }
    if (simple.labels) {
      xat <- cumsum(sapply(fitted(fitN), nrow)) - sapply(fitted(fitN), nrow) + 1
      xat <- xat[xat < ntaxa]
      xat <- xat / ntaxa

      axis(4, 1-xat, labels=datatypes[seq_along(xat)], las=1)
    } else {
      xat <- (seq_len(nrow(pN)) - 1) / (nrow(pN) - 1)
      axis(4, xat, labels=rownames(pN), las=1)
    }
  }



# show.param=TRUE
# sqroot=FALSE
#
# figure_path="figures/"
# figure_name=paste0("shifted-flipped-data2")
# evaluation=BIC
# title.margin=0.1
# param.height=0.4
# data.height=1.3
# param="average"
# figure_width=1100
# figure_height=1400
# figure_res=160
# legend_down=0.25
# legend_downshift=0.05
# DNAlabels=c(paste0("-",window/2000, "kb"),"0",paste0(window/2000, "kb"))
#
# data=data$data
plot.heatmap.before.clustering<- function(
  show.param=TRUE,
  sqroot=FALSE,
  bin_size=40,
  figure_path,
  figure_name,
  evaluation=BIC,
  title.margin=0.3,
  param.height=0.3,
  data.height=1.3,
  param="DirichletParameters",
  figure_width=1200,
  pdf_width=7,
  pdf_height=8,
  figure_height=800,
  figure_res=80,
  legend_down=0.15,
  legend_downshift=0.075,
  DNAlabels=c("-1kb","0","1kb"), data=NULL) {

  library("plyr")
  library("reshape2")
  library("egg")
  library("gridExtra")
  library("grid")

  i=1
  for(d in data){ #sig is window x N
    rownames(data[[i]]) <- paste0( seq(nrow(d)))
    colnames(data[[i]]) <- paste0('bin', seq(ncol(d)))
    i=i+1
  }

  S <- sapply(data, ncol) #40 window/bin_size

  M <- length(data) #9

  datatypes <- names(data)

  dataframe <- ldply(data, function(d) {
    d <- as.data.frame(d)
    d$id <- factor(1:nrow(d))
    d <- melt(d, measure = 1:(ncol(d)-1), variable.name = "Bin",
              value.name = "Reads")
    d$Bin <- as.numeric(d$Bin)
    d
  }, .id = "Datatype")

  #50*9*1000
  #[1] 450000
  #data is data.frame':	S[1]*K*M*N=2245500 obs. of  5 variables:
  #Datatype: modifications
  #id str of "1"-"998"
  #Component, cluster index as str
  #Bin, vector of 1-50
  #Reads: bin count 0-74

  if(sqroot==TRUE){
    dataframe$Reads <- sqrt(dataframe$Reads)
  }

  #grid.newpage()
  #grid.payout descrives a subdivision of a rectangular region
  #nrow=K x 2=10 or K=5
  #ncol=M=9
  #heights: heights of the rows in the layout

  #viewport create viewports which describe rectangular regions on a
  #graphics device and define a number of coordinate systems within those regions
  #layout: A Grid layout object which splits the viewport into subregions.
  #pushViewport(viewport(layout = grid.layout(ifelse(show.param,
  #                                                  K * 2, K), M,
  #                                           heights = heights)))
  #print(paste0(figure_path,"/",figure_name))
  png(paste0(figure_path,"/",figure_name, ".png"),
      width=figure_width,
      height=figure_height,
      res = figure_res)
  #postscript(paste0(figure_path,figure_name), paper="special", width=10, height=5)

  #postscript(paste0(figure_path,figure_name), paper="special", width=1000, height=500,
  #           units="px", pointsize=12,res = 80)

  gg_list<-list() #of length 2*K*M

  #layout.pos.row: A numeric vector giving the rows occupied by this viewport in its parent's layout.
  #layout.pos.col: A numeric vector giving the columns occupied by this viewport in its parent's layout.
  #vplayout <- function(x, y) viewport(layout.pos.row = x,
  #                                    layout.pos.col = y)

  tmp=subset(dataframe,  Datatype == datatypes[1])
  profile_length=length(unique(dataframe$Bin))
  bin=paste0("bin",1:profile_length) #20, 2000/100



  #over clusters

  for (m in seq_len(M)) { #over datatypes

    margin <- c(0, 0, 0, 0)
    margin <- unit(margin, "mm")
    #parameters for certain datatype and cluster
    #These are dirichlet parameters
    if(sqroot==TRUE){
      params=data.frame(Reads=sqrt(colMeans(data[[datatypes[m]]])), Bin=1:profile_length)
    }else{
      params=data.frame(Reads=colMeans(data[[datatypes[m]]]), Bin=1:profile_length)
    }



    #ggparams <- ggplot(subset(params, Bin %in% seq(995,1005,1)), aes(Bin,Reads))

    ggparams <- ggplot(params, aes(Bin, Reads))
    ggparams <- ggparams + geom_area(fill = "brown")
    #ggparams <- ggparams + geom_point(fill = "black")

    ggparams <- ggparams + labs(x = "", y = "") + guides(color = F)

    #ggparams <- ggparams +xlim(995,1005)
    #int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
    #ggparams <- ggparams + scale_x_continuous(breaks = int_breaks)
    ggparams <- ggparams + scale_x_continuous(breaks = NULL, expand = c(0, 0,0.01,0))
    #plot(ggparams)
    #ggparams <- ggparams + scale_x_discrete(name="Bin", breaks = seq(995,1005,1), labels= seq(995,1005,1))
    #,
    #                                          expand = c(0, 0,0.01,0))


    ggparams <- ggparams + theme_minimal()
    #plot.margin: margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    ggparams <- ggparams + theme(plot.margin = margin,
                                 axis.ticks.margin = unit(c(0, 0), "mm"),
                                 panel.margin = unit(c(0, 0, 0, 0), "mm"),
                                 plot.background = element_blank())

    if(bin_size==1){
      ggparams <- ggparams + geom_vline(xintercept=profile_length/2+1, linetype="dashed",
                                        color="black", size=0.5)
    }else{
      ggparams <- ggparams + geom_vline(xintercept=profile_length/2+0.5, linetype="dashed",
                                        color="black", size=0.5)
    }


    ggparams <- ggparams + ggtitle(datatypes[m])
    ggparams <- ggparams + theme(plot.title = element_text(hjust = 0.5))
    #print(ggparams, vp = vplayout((k - 1) * 2 + 1, m))


    #add cluster name
    margin <- c(0, 0, 0, 0)
    margin <- unit(margin, "mm")
    #Plot the data

    data_subset=subset(dataframe,  Datatype == datatypes[m])
    if(sort(data_subset$Reads)[0.95*length(data_subset$Reads)] ==1){
      data_subset$Reads[which(data_subset$Reads>sort(data_subset$Reads)[0.95*length(data_subset$Reads)])]= 5
    }


    #data_subset$Reads[which(data_subset$Reads>sort(data_subset$Reads)[0.95*length(data_subset$Reads)])]= sort(data_subset$Reads)[0.95*length(data_subset$Reads)]
    gg <- ggplot(data_subset,
                 aes(Bin, id))

    gg <- gg+ geom_raster(aes(fill = Reads))
    gg <- gg + labs(x = "", y="") #, y = yaxis)
    #gg <- gg + guides(fill = F) #uncomment if you want the figure without legend
    gg <- gg + scale_y_discrete(breaks = NULL)
    #gg <- gg + scale_x_continuous(breaks = c(1,profile_length/2-1,profile_length),
    #                              labels= DNAlabels,expand = c(0, 0,0.01,0)


    if(bin_size==1){
      gg <- gg + scale_x_continuous(breaks = c(1,profile_length/2+1,profile_length),
                                    labels= DNAlabels,expand = c(0, 0,0.01,0) )
    }else{
      gg <- gg + scale_x_continuous(breaks = c(1,profile_length/2+0.5,profile_length),
                                    labels= DNAlabels,expand = c(0, 0,0.01,0))
    }


    #col_fun = try( colorRamp2( c( min(dataframe$Reads) , quantile(max(dataframe$Reads), 0.95) ),
    #                           c("yellow", "red") ) )

    myPalette=RColorBrewer::brewer.pal(name = "YlOrRd", 9)

    #gg <- gg + scale_fill_gradientn(colours = myPalette)

    #print( seq(sqrt( min(data_subset$Reads)), sqrt(max(data_subset$Reads)) ,
    #     length = 5) )

    #ceiling(max(data_subset$Reads)/10)*10
    #sort(data_subset$Reads)[0.95*length(data_subset$Reads)]

    #seq(min(data_subset$Reads),
    #    ceiling(sort(data_subset$Reads)[0.95*length(data_subset$Reads)]/10)*10,
    #    10)

    #rint(seq(min(data_subset$Reads),
    #          ceiling(sort(data_subset$Reads)[0.95*length(data_subset$Reads)]/10)*10,
    #          10))


    #print(sqrt(seq(min(data_subset$Reads),
    #               ceiling(sort(data_subset$Reads)[0.95*length(data_subset$Reads)]/10)*10,
    #               10))^2)
    my_breaks=sqrt(seq(min(data_subset$Reads),
                       ceiling(max(data_subset$Reads)/10)*10,
                       10))^2
    if(length(my_breaks)>5){
      my_breaks=sqrt(seq(min(data_subset$Reads),
                         ceiling(max(data_subset$Reads)/50)*50, 50))^2
      if(length(my_breaks)==2){
        my_breaks=c(0,10,30,50)
      }
      if(length(my_breaks)>5){
        my_breaks=sqrt(seq(min(data_subset$Reads),
                           ceiling(max(data_subset$Reads)/100)*100, 100))^2

        if(length(my_breaks)>5){
          if( (length(my_breaks) %% 2)==0 ){ #even
            my_breaks=my_breaks[c(1,seq(2,length(my_breaks),2))]
          }else{
            my_breaks=my_breaks[seq(1,length(my_breaks),2)]
          }
        }


      }
    }else if(length(my_breaks)<3){
      my_breaks=sqrt(seq(min(data_subset$Reads),
                         ceiling(max(data_subset$Reads)/5)*5,
                         5))^2
    }

    if(length(my_breaks)==2){
      my_breaks=c(0,1,2,3,4,5)
    }
    if(300 %in% my_breaks ){
      my_breaks=c(0,100,300)
    }
    if(400 %in% my_breaks ){
      my_breaks=c(0,200,400)
    }
    if(500 %in% my_breaks ){
      my_breaks=c(0,200,500)
    }
    if(700 %in% my_breaks ){
      my_breaks=c(0,100,500,700)
    }


    if(1000 %in% my_breaks){
      my_breaks=c(0,500,1000)
    }

    print(my_breaks)

    #print(unique(round_any( ( seq(sqrt( min(data_subset$Reads)), sqrt(max(data_subset$Reads)) ,
    #                              length = 5) ), 5)^2))
    #my_breaks <- unique(round_any( ( seq(sqrt( min(data_subset$Reads)), sqrt(max(data_subset$Reads)) ,
    #                                length = 5)), 5)^2)
    gg <- gg + scale_fill_gradientn(colours = myPalette,
                                    trans="sqrt",
                                    breaks=my_breaks)


    gg <- gg + theme_minimal()
    #gg <- gg + theme(axis.ticks.y = element_blank())
    gg <- gg + theme(plot.margin = margin)
    #gg <- gg + geom_vline(xintercept=profile_length/2-1,
    #                      linetype="dashed",
    #                      color="black", size=0.5)
    if(bin_size==1){
      gg <- gg + geom_vline(xintercept=profile_length/2+1,
                            linetype="dashed",
                            color="black", size=0.5)
    }else{
      gg <- gg + geom_vline(xintercept=profile_length/2+0.5,
                            linetype="dashed",
                            color="black", size=0.5)
    }




    gg <- gg + theme(axis.ticks.margin = unit(c(0, 0, 0, 0), "mm"))
    gg <- gg + theme(panel.margin = unit(c(0, 0, 0, 0), "mm"))
    gg <- gg + theme(plot.background = element_blank())
    gg <- gg + theme(panel.grid.major = element_blank())
    gg <- gg + theme(panel.grid.minor = element_blank())
    gg <- gg + labs(fill="")

    gg <- gg + theme(plot.margin=margin(t=0,r=0,b=9, l=0,unit="mm"),
                     legend.position = c(0.5,-0.15),
                     legend.direction = "horizontal",
                     legend.key.size = unit(5, "mm"),
                     legend.key.width = unit(4.5,"mm"))
    #plot(gg)






    #reduce plot margins

    #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    ggparams <- ggparams + theme(plot.margin=margin(t=0,r=0,b=-3, l=0,
                                                    unit="mm"))

    gg <- gg + theme(plot.margin=margin(t=-3,r=0,b=9, l=0,
                                        unit="mm"))



    gg_list[[m]] <-ggarrange(ggparams, gg, ncol=1,
                             heights=c(0.25,1))

    #plot(gg_list[[m]])



  }

  #ggarrange(plots=gg_list, ncol=9)
  grid.arrange(grobs=gg_list, ncol=length(data), top="", bottom="",
               left="", right="")

  #,
  #heights=heights, top="", bottom="",
  #           left="", right="")

  dev.off()
}




# show.param=TRUE
# sqroot=FALSE
# concat=FALSE
# labels=data.tmp$labels
# figure_path=paste0("../",results_path,name,"/",analysisName,"/figures/")
# figure_name=paste0(as.character(coverage_values[1,]),sep="-", collapse="")
# param="average"
# figure_width=500
# figure_height=500
# figure_res=80
# legend_down=0.25
# legend_downshift=0.05
# DNAlabels=c("-1.8kb","0","1.8kb")
# cluster_nro=K.simul
# data=data$data
#
# fits=NULL
# M.orig=NULL
# evaluation=BIC
# title.margin=0.3
# param.height=0.3
# data.height=1.3
# N_start=NULL
# shift=FALSE
# flip=FALSE


#' Title
#'
#' @param fits To visualize DMN results with unconcatenated data, object given by dmn(), do not give labels
#' @param show.param
#' @param concat Whether to visualize results fromclustering of concatenated data, then give both data and fits, no labels
#' @param labels
#' @param sqroot
#' @param figure_path
#' @param figure_name
#' @param evaluation
#' @param title.margin
#' @param param.height
#' @param data.height
#' @param param
#' @param figure_width
#' @param figure_height
#' @param figure_res
#' @param legend_down
#' @param legend_downshift
#' @param DNAlabels
#' @param cluster_nro The number of clusters visualized, if not specified, selected by BIC
#' @param data To visualize data before clustering(give labels) or visualize clustering results of concatenated data (do not give labels)
#' @param N_start Visualize just small set of regions
#'
#' @return
#' @export
#'
#' @examples
plot.heatmap.final<- function( #data.bool=FALSE,
  fits=NULL,
  show.param=TRUE,
  concat=FALSE,
  M.orig=NULL,
  labels=NULL,
  #path_to_orig_file=NULL,
  sqroot=FALSE,
  #path_to_file="",
  figure_path,
  figure_name,
  evaluation=BIC,
  title.margin=0.3,
  param.height=0.3,
  data.height=1.3,
  param="DirichletParameters",
  figure_width=1200,
  pdf_width=7,
  pdf_height=8,
  figure_height=800,
  figure_res=80,
  legend_down=0.15,
  legend_downshift=0.075,
  DNAlabels=c("-1kb","0","1kb"),
  cluster_nro=NULL,
  data=NULL,
  N_start=NULL,
  shift=FALSE, flip=FALSE
) {

  library("plyr")
  library("reshape2")
  library("egg")
  library("gridExtra")
  library("grid")
  library("scales")

  if(is.null(cluster_nro) ){
    cluster_nro=as.integer(which.min(sapply(fits, evaluation))) #4
  }

  if(is.null(data)==TRUE){ #data is NULL, extract data from fits
    cl=fits[[cluster_nro]]
    if (class(cl) == "DMN"){
      data <- cl@Data
    }else{
      stop("fits element needs to be of class DMN")
    }
  }else{#data to be visualized is given as data

    #data=data
  }

  data_orig=data #
  #figure_name_final=paste0(figure_name,"-",param,"-",cluster_nro,"-clusters.pdf")
  figure_name_final=paste0(figure_name,"-",param,"-",cluster_nro,"-clusters.png")

  if(is.null(labels) && !is.null(fits)){
    labels <- mixture(fits[[cluster_nro]], assign = T)
  }

  #if (!any(diff(sapply(data, ncol)))) #True if same number of bins
  #  scales <- "fixed"
  #else scales <- "free_x"

  f_cols <-function(long.matrix, num){
    long.matrix=t(long.matrix)
    t(lapply(split(long.matrix,
                   rep(seq(num), each=(ncol(long.matrix)/num)*nrow(long.matrix))),
             function(x) matrix(x, nrow=nrow(long.matrix))
    ))
  }
  f_rows<-function(long.matrix, num){
    lapply(split(long.matrix,
                 rep(seq(num), each=(ncol(long.matrix)/num)*nrow(long.matrix))),
           function(x) matrix(x, nrow=nrow(long.matrix))
    )
  }


  Lx <- sapply(data, ncol) #note this is different from Lx_orig if shifting applied
  K <- length(unique(labels))
  M <- length(data)
  N <- length(labels)
  datatypes <- names(data)

  if(is.null(N_start)){
    data <- plyr::ldply(data, function(d) {
      d <- as.data.frame(d)
      d$id <- factor(1:nrow(d))
      d$Component <- factor(labels)
      d <- reshape2::melt(d, measure = 1:(ncol(d) - 2), variable.name = "Bin",
                          value.name = "Reads")
      d$Bin <- as.numeric(d$Bin)
      d
    }, .id = "Datatype")
  }else{
    data <- plyr::ldply(data, function(d) {

      d <- d[N_start,]
      d <- as.data.frame(d)
      d$id <- factor(1:nrow(d))
      d$Component <- factor(labels[N_start])
      d <- reshape2::melt(d, measure = 1:(ncol(d) - 2), variable.name = "Bin",
                          value.name = "Reads")
      d$Bin <- as.numeric(d$Bin)
      d
    }, .id = "Datatype")

  }
  if(sqroot==TRUE){
    data$Reads <- sqrt(data$Reads)
  }


  if (show.param) {
    if(param=="DirichletParameters"){

      params=fitted(fits[[cluster_nro]])

      if(concat==TRUE && flip==TRUE && shift==FALSE){
        params_new<-list()
        for(m in seq_len(M)){ #if concat==TRUE, M is the orig number of chromatin features
          params_new[[names(data_orig)[m]]]=matrix(0.0, nrow=Lx[[m]], ncol=K)
          for(k in seq_len(K)){
            params_new[[names(data_orig)[m]]][,k]=params[[1]][ ( (m-1)*Lx[[m]]+1 ):( m*Lx[[m]] ),k] #Tähän jäin
          }
        }
        params=params_new

      }

      if(concat==TRUE && shift==FALSE && flip==FALSE){

        #reformulate params/unconcatenate params
        params_tmp=params
        params_tmp$Data=params_tmp$Data[1:(nrow(params[[1]])/M.orig),]
        params2=rep(params_tmp, M.orig)

        tmp=f_cols(params$Data, M.orig)

        names(params2)=datatypes

        for(ii in 1:M.orig){
          va=t(tmp[[ii]])
          params2[[ii]][1:(nrow(params[[1]])/M.orig),]=va

        }
        params=params2

      }
      if(concat==FALSE && shift==TRUE){

        params_new<-list()
        S=length(params[[1]][,1])-unique(Lx)+1 #number of shift
        for(m in seq_len(M)){
          #print(m)
          params_new[[names(data_orig)[m]]]=matrix(0.0, nrow=Lx[[m]], ncol=K)
          for(k in seq_len(K)){
            params_new[[names(data_orig)[m]]][,k]=params[[names(data_orig)[m]]][floor(S/2+1):( length(params[[1]][,1])- floor(S/2)),k]
          }
        }

        params=params_new



      }
      if(concat==TRUE && shift==TRUE){
        params_new<-list()
        S=length(params[[1]][,1])-2*unique(Lx)+1 #number of shift
        for(m in seq_len(M)){ #if concat==TRUE, M is the orig number of chromatin features
          params_new[[names(data_orig)[m]]]=matrix(0.0, nrow=Lx[[m]], ncol=K)
          for(k in seq_len(K)){
            params_new[[names(data_orig)[m]]][,k]=params[[1]][ (floor(S/2+1) + (m-1)*Lx[[m]] ):( floor(S/2)+ncol(data_orig[[m]])+(m-1)*Lx[[m]]),k] #Tähän jäin
          }
        }
        params=params_new

      }



      # if(concat==TRUE){
      #  # if(shift==FALSE){ #model without shift, shift==FALSE (is.null(unshifted_data)==TRUE && is.null(fits)==TRUE
      #     params_tmp=params
      #     params_tmp$Data=params_tmp$Data[1:50,]
      #     params2=c(params_tmp, params_tmp)
      #     tmp=f_cols(params$Data, M.orig)
      #
      #     names(params2)=datatypes
      #
      #     for(ii in 1:M.orig){
      #       va=t(tmp[[ii]])
      #       params2[[ii]][1:50,]=va
      #
      #     }
      #     params=params2
      #   # }else{ #shift==TRUE
      #   #   params_new<-list()
      #   #   S=length(params[[1]][,1])-2*unique(Lx)+1 #number of shift
      #   #   for(m in seq_len(M)){ #if concat==TRUE, M is the orig number of chromatin features
      #   #     params_new[[names(unshifted_data)[m]]]=matrix(0.0, nrow=Lx[[m]], ncol=K)
      #   #     for(k in seq_len(K)){
      #   #       params_new[[names(unshifted_data)[m]]][,k]=params[[1]][ (floor(S/2+1) + (m-1)*Lx[[m]] ):( floor(S/2)+ncol(unshifted_data[[m]])+(m-1)*Lx[[m]]),k] #Tähän jäin
      #   #     }
      #   #   }
      #   #   params=params_new
      #   # }
      #
      # }else{ #concat is false
      #   #the length of params need to be same the length of data (may differ if shifting applied)
      #   if( max(data$Bin)!=  length(params[[1]][,1]) ){ #50 vs 70 this is considered when shifting
      #     params_new<-list()
      #     S=length(params[[1]][,1])-unique(Lx)+1 #number of shift
      #     for(m in seq_len(M)){
      #       params_new[[names(data_orig)[m]]]=matrix(0.0, nrow=Lx[[m]], ncol=K)
      #       for(k in seq_len(K)){
      #         params_new[[names(data_orig)[m]]][,k]=params[[names(data_orig)[m]]][floor(S/2+1):( length(params[[1]][,1])- floor(S/2)),k]
      #       }
      #     }
      #
      #    params=params_new
      #
      #   }
      # }
      #

      params <- melt(params, varnames = c("Bin", "Component"),
                     value.name = "Reads")
      names(params)[which(names(params) == "L1")] <- "Datatype"
      data.height <- as.vector((table(labels)/N) * (data.height *
                                                      K))
      #1.7368 1.3208 1.1544 0.9880

      heights <- as.vector(t(data.frame(param.height, data.height)))
      #0.7000 1.7368 0.7000 1.3208 0.7000 1.1544 0.7000 0.9880
      heights[1] <- heights[1] + title.margin
      heights[2:K * 2] <- heights[2:K * 2] - title.margin/((K * 2) - 1)
      #1.0000000 3.1000074 0.7000000 1.3205466 0.7000000 0.8674381 0.7000000 0.5278117
      heights <- unit(heights, "null")
    }
    if(param=="average"){
      #params should be a list over modifications, elements are #bins x cluster_nro matrices
      if(sqroot==TRUE){
        params <-list()

        for(m in seq_len(M)){
          cm_matrix<-matrix(NA, nrow=ncol(data_orig[[m]]), ncol=K )
          for(k in seq_len(K)){
            cm_matrix[,k]=sqrt( colMeans( data_orig[[m]][which((labels==k) ==TRUE),]) )
          }
          params[[ names(data_orig)[m] ]] =cm_matrix

        }
        params <- melt(params, varnames = c("Bin", "Component"),
                       value.name = "Reads")
        names(params)[which(names(params) == "L1")] <- "Datatype"
        data.height <- as.vector((table(labels)/N) * (data.height *
                                                        K))
        heights <- as.vector(t(data.frame(param.height, data.height)))
        heights[1] <- heights[1] + title.margin
        heights[2:K * 2] <- heights[2:K * 2] - title.margin/((K * 2) - 1)
        heights <- unit(heights, "null")
      }else{ #sqroot==FALSE

        params <-list()

        #if(concat==FALSE){
        for(m in seq_len(M)){
          cm_matrix<-matrix(NA, nrow=ncol(data_orig[[m]]), ncol=K )
          for(k in seq_len(K)){
            cm_matrix[,k]=colMeans( data_orig[[m]][which((labels==k) ==TRUE),])
          }
          params[[ names(data_orig)[m] ]] =cm_matrix

        }
        #}
        # else{ #concat==TRUE
        #   #data_orig is N x M.orig'unique(Lx) matrix
        #
        #   for(m in seq_len(M.orig)){
        #     cm_matrix<-matrix(NA, nrow=ncol(data_orig[[1]])/M.orig, ncol=K )
        #     for(k in seq_len(K)){
        #       cm_matrix[,k]=colMeans( data_orig[[1]][which((labels==k) ==TRUE), ( 1+(m-1)*unique(Lx) ):( m*unique(Lx)  )  ])
        #     }
        #     params[[ m ]] =cm_matrix
        #
        #   }
        #
        #   names(params)=datatypes
        #
        # }

        params <- melt(params, varnames = c("Bin", "Component"),
                       value.name = "Reads")
        names(params)[which(names(params) == "L1")] <- "Datatype"
        data.height <- as.vector((table(labels)/N) * (data.height *
                                                        K))
        heights <- as.vector(t(data.frame(param.height, data.height)))
        heights[1] <- heights[1] + title.margin
        #heights[2:K * 2] <- heights[2:K * 2] - title.margin/((K * 2) - 1) is this to add the legend??
        heights <- unit(heights, "null")

      } #sqroot==FALSE
    }


  }else{ #do not show param
    heights <- as.vector((table(labels)/N) * (K))
  }

  tmp=subset(data,  Datatype == datatypes[1])
  profile_length=length(unique(data$Bin)) #50
  bin=paste0("bin",1:profile_length) #20, 2000/100


  gAB_grobList<- list()
  i=1
  fontsize=14
  for (k in seq_len(K)) {
    #print(k)

    comp.label <- sort(unique(labels))[k]
    for (m in seq_len(M)) {
      #print(m)
      if (show.param) {
        #margin <- c(-2, 2, -10, -2) #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        if(m < M){
          margin <- c(2, 1, -6, 0)
        }else{
          margin <- c(2, 4, -6, 0)
        }

        #if (k%%K == 1)
        #  margin[1] <- margin[1] + 5
        #if (m%%M != 1)
        #  margin[4] <- margin[4] + 2
        margin <- unit(margin, "mm") #3mm   2mm   -10mm 0mm
        ggparams <- ggplot(subset(params, Component ==
                                    comp.label & Datatype == datatypes[m]),
                           aes(Bin, Reads, group = Component))

        ggparams <- ggparams +  geom_area(fill = "brown")
        ggparams <- ggparams + labs(x = "", y = "") + guides(color = F)
        #ggparams <- ggparams +  scale_x_discrete(breaks = seq(0, Lx[m], 10),
        #                   expand = c(0, 0))
        ggparams <- ggparams +  scale_x_discrete(breaks = NULL,
                                                 expand = c(0, 0,0.01,0))


        ggparams <- ggparams + theme_minimal()
        ggparams <- ggparams + theme(plot.margin = margin,
                                     axis.ticks.margin = unit(c(0, 0), "mm"),
                                     panel.margin = unit(c(0, 0, 0, 0), "mm"))

        #if(bin_size==1){
        #  ggparams <- ggparams + geom_vline(xintercept=profile_length/2+1, linetype="dashed",
        #                                    color="black", size=0.5)
        #}else{
        ggparams <- ggparams + geom_vline(xintercept=profile_length/2+0.5, linetype="dashed",
                                          color="black", size=1)
        #}

        #####Try this
        #ggparams <- ggparams + theme(plot.margin=margin(t=0,r=0,b=-3, l=0,
        #                                                unit="mm"))



        if (k%%K == 1) {
          ggparams <- ggparams + ggtitle(datatypes[m])
          ggparams <- ggparams + theme(plot.title = element_text(hjust = 0.5))
        }
        #print(ggparams, vp = vplayout((k - 1) * 2 + 1,
        #                              m))

        ggparams <- ggparams + theme(axis.text.y= element_text(size=fontsize),
                                     text=element_text(size=fontsize),
                                     panel.grid.minor=element_blank())
        ggparams <- ggparams + scale_y_continuous(breaks = scales::breaks_extended(Q = c(1, 2,4, 3), w=c(0.6, 0.05, 0.5, 0.05)))



      }

      if(m < M){
        margin <- c(0, 1, 2, 0)
      }else{
        margin <- c(0, 4, 2, 0)
      }





      if (m%%M == 1 | M==1) {
        #3mm   2mm   -10mm 0mm
        #if (k%%K == 1)
        ##  margin[1] <- margin[1] + 5
        ##if (m%%M == 1)
        #  margin[4] <- margin[4] + 2


        #margin[4] <- margin[4] + 5
        yaxis <- paste("Cluster", comp.label)
      }else {
        #margin[4] <- margin[4] + 5
        yaxis <- ""
      }
      margin <- unit(margin, "mm")

      data_subset=subset(data, Component == comp.label &
                           Datatype == datatypes[m])


      #What this does?
      #if(sort(data_subset$Reads)[0.95*length(data_subset$Reads)] ==1){
      #  data_subset$Reads[which(data_subset$Reads>sort(data_subset$Reads)[0.95*length(data_subset$Reads)])]= 5
      #}

      gg <- ggplot(data_subset, aes(Bin, id, group = Component))
      gg <- gg + geom_raster(aes(fill = Reads))
      gg <- gg + labs(x = "", y = yaxis) #+ guides(fill = F) this removes legend
      #gg <- gg + scale_x_continuous(expand = c(0, 0,0.01,0))
      gg <- gg + scale_x_continuous(breaks = c(1,profile_length/2+0.5,profile_length),
                                    labels= DNAlabels,expand = c(0, 0,0.01,0))
      gg <- gg + geom_vline(xintercept=profile_length/2+0.5,
                            linetype="dashed",
                            color="black", size=1)


      my_breaks=sqrt(seq(min(data_subset$Reads),
                         ceiling(max(data_subset$Reads)/10)*10,
                         10))^2
      if(length(my_breaks)>5){
        my_breaks=sqrt(seq(min(data_subset$Reads),
                           ceiling(max(data_subset$Reads)/50)*50, 50))^2
        if(length(my_breaks)==2){
          my_breaks=c(0,10,30,50)
        }
        if(length(my_breaks)>5){
          my_breaks=sqrt(seq(min(data_subset$Reads),
                             ceiling(max(data_subset$Reads)/100)*100, 100))^2

          if(length(my_breaks)>5){
            if( (length(my_breaks) %% 2)==0 ){ #even
              my_breaks=my_breaks[c(1,seq(2,length(my_breaks),2))]
            }else{
              my_breaks=my_breaks[seq(1,length(my_breaks),2)]
            }
          }


        }
      }else if(length(my_breaks)<=3){

        my_breaks=sqrt(seq(min(data_subset$Reads),
                           ceiling(max(data_subset$Reads)/5)*5,
                           5))^2
      }



      if(length(my_breaks)==2){
        my_breaks=c(0,1,2,3,4,5)
      }
      if(300 %in% my_breaks ){
        my_breaks=c(0,100,300)
      }
      if(400 %in% my_breaks ){
        my_breaks=c(0,200,400)
      }
      if(500 %in% my_breaks ){
        my_breaks=c(0,200,500)
      }
      if(700 %in% my_breaks ){
        my_breaks=c(0,100,500,700)
      }


      if(1000 %in% my_breaks){
        my_breaks=c(0,500,1000)
      }

      print(my_breaks)



      gg <- gg + scale_fill_gradientn(colours = RColorBrewer::brewer.pal(name = "YlOrRd", 9),
                                      trans="sqrt", breaks=my_breaks)
      gg <- gg + theme_minimal()
      gg <- gg + theme(axis.ticks.y = element_blank(),
                       axis.text.y = element_blank())
      gg <- gg + labs(fill="")

      gg <- gg + theme(plot.margin = margin,
                       axis.ticks.margin = unit(c(0, 0, 0, 0), "mm"),
                       panel.margin = unit(c(0, 0, 0, 0), "mm"),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       #legend.position = c(0.5,-legend_down-(k-1)*legend_downshift) ,
                       legend.position="bottom",
                       legend.justification = "center",
                       legend.direction = "horizontal",
                       legend.box.just="top",
                       #legend.justification=c("top","center"),
                       legend.margin=margin(0,0,0,0),
                       legend.box.margin=margin(-10,-5,-5,-5),
                       legend.key.size = unit(6, "mm"),
                       legend.key.width = unit(6,"mm"),
                       legend.text=element_text(size=fontsize),
                       axis.text.x= element_text(size=fontsize),
                       text=element_text(size=fontsize))

      #Try this
      #gg <- gg + theme(plot.margin=margin(t=-3,r=0,b=9, l=0,
      #                                    unit="mm"))




      if (show.param){
        #print(gg, vp = vplayout((k - 1) * 2 + 2, m))
        #print(ggarrange(ggparams, gg, ncol=1,
        #                heights=c(0.25,1) ), vp = vplayout((k - 1)  + 2, m))

        #gA=ggplot_gtable(ggplot_build(ggparams)) #"gtable" "gTree"  "grob"   "gDesc"
        #gB=ggplot_gtable(ggplot_build(gg))
        gA=ggplotGrob(ggparams) #"gtable" "gTree"  "grob"   "gDesc"
        gB=ggplotGrob(gg)

        #gAB_grobList[[i]] <- ggarrange(ggparams, gg, ncol=1, nrow=2, align="v")

        #maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
        #gA$widths[2:3] <- as.list(maxWidth)
        #gB$widths[2:3] <- as.list(maxWidth)

        gAB_grobList[[i]]=arrangeGrob(cowplot::plot_grid(ggparams, gg, align="v", axis="l", ncol=1, rel_heights=as.numeric( heights[c( (k-1)*2+1 ,(k-1)*2+2)])))

        #gAB_grobList[[i]]=ggplotGrob ( as_ggplot( arrangeGrob(gA,gB,nrow=2,heights=heights[c( (k-1)*2+1 ,(k-1)*2+2)])))
        #grid.arrange(gAB)
        i=i+1
      }else{
        #print(gg, vp = vplayout(k, m))
        gAB_grobList[[i]]=gB
        i=i+1
      }
    }
  }

  if(K==1){
    heights=heights[-which(is.na(heights)==TRUE)]
  }

  print(figure_height)
   png(paste0(figure_path,figure_name_final),
       width=figure_width,
       height=figure_height,
       res = figure_res)
  #pdf(paste0(figure_path,figure_name_final), width=pdf_width, height=pdf_height)
  #     width=figure_width,
  #     height=figure_height,
  #     res = figure_res)
  # grid.newpage()
  #egg::ggarrange(plots=gAB_grobList, ncol = M, nrow = K)
  #ggpubr::ggarrange(plots=gAB_grobList, ncol = M, nrow = K)
  if(K==1){
    grid.arrange(arrangeGrob(grobs=gAB_grobList, nrow=cluster_nro, heights=heights[1]+heights[2]))
  }else{
    grid.arrange(arrangeGrob(grobs=gAB_grobList, nrow=cluster_nro, heights=heights[seq(1,length(heights),2)]+heights[seq(2,length(heights),2)]))
  }

  dev.off()
  dev.off()


}


#param="average" #"DirichletParameters"



