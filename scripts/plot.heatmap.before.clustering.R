library("plyr")
library("reshape2")
library("egg")
library("gridExtra")
library("grid")
#library("ChromDMM")
#library("EnrichedHeatmap")
#library("circlize")

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


