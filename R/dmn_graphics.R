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
#'
#' @return
#' @export
#'
#' @examples
plot.EM <- function(fit, smoothness.scale='free', k=10) {

  EM.diagnostics <- fit@EM.diagnostics
  K <- ncol(mixture(fit))
  #plot hkm values
  comp.labels=paste0("Cluster ", seq(1,k,1))
  names(comp.labels)=seq(1,k,1)
  
  hkm.plot <- ggplot(EM.diagnostics,
                     aes(seq_along(hkm), hkm, color=EM.iter)) +
    geom_line(aes(group=Datatype)) +
    geom_point(shape=1) +
    facet_wrap(~Datatype+Component, scale=smoothness.scale, ncol=K, labeller=labeller(Component=comp.labels)) +
    labs(title='Numerical optimization iterations vs. regularization term', x='Numerical optimization iterations', y=expression(h[k]^(m))) +
#     scale_color_gradientn('EM iterations', colours=RColorBrewer::brewer.pal(n=9, name='YlOrRd'))
#     scale_color_brewer(palette='Blues') +
    scale_color_continuous('EM iterations') +
    theme_minimal()
  
  
  #plot number of num.opt. iterations
  nop <- ggplot(EM.diagnostics, aes(EM.iter, NO.iter.count)) +
    geom_line(aes(color=factor(Datatype))) +
    geom_point(aes(color=factor(Datatype))) +
    #geom_smooth(method='gam', formula=y~s(x, k=k, bs='cs')) +
    facet_grid(Datatype~Component, labeller=labeller(Component=comp.labels)) +
    guides(color=F) +
    labs(title='EM iterations vs. Numerical optimization iterations', x='EM iterations', y='Numerical optimization iterations') +
    theme_minimal()

  #create a grid layout and print ggplots on it
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  print(hkm.plot, vp = vplayout(1, 1))
  print(nop, vp = vplayout(2, 1))

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