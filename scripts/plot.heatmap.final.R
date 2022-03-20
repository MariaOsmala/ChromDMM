
library("plyr")
library("reshape2")
library("egg")
library("gridExtra")
library("grid")
library("scales")

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

