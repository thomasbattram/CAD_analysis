######Forest plot function



forestPlot<-function(result,sigP=NULL,coll=NULL, columns)
{
  ### result must be a list format 
  ### sigP: threshold for significant P  
  
  library('RColorBrewer')  
  
  #metabolite names
  rname=result[[1]][["Metabolite"]]
  cname=names(result)
  
  #x-axis range  
  ##This finds the range between each data.frame under the "2.5 %" and "97.5 %" columns 
  xrange1=range(sapply(result,function(x){x[,'2.5 %']}), na.rm=T)
  xrange2=range(sapply(result,function(x){x[,'97.5 %']}), na.rm=T)
  xrange=range(c(xrange1,xrange2))
  
  #show the plots in four columns - needs to be four because the metabolite list is so long
  #picks out the metabolite that is a quarter of the way down the list
  if (columns == 4) {
    dividename=rname[round(median(1:length(rname))/2)]
    #match() is similar to %in% and so this function converts the name of the metab found above into a number (i.e. 1/4 of the total number of metabs)
    ind=match(dividename,rname)
    #these variables allow us to split the forest plot into four columns each containing a quarter of the values
    ind1=1:ind;
    ind2=(ind+1):(length(rname)*0.5) ####CHANGED
    ind3=(ind2[length(ind2)] + 1):(length(rname)*0.75)
    ind4=(ind3[length(ind3)] + 1):length(rname)
  } else if (columns == 2) {
    dividename=rname[round(median(1:length(rname)))]
    #match() is similar to %in% and so this function converts the name of the metab found above into a number (i.e. 1/4 of the total number of metabs)
    ind=match(dividename,rname)
    #these variables allow us to split the forest plot into four columns each containing a quarter of the values
    ind1=1:ind;
    ind2=(ind+1):(length(rname))
  } else if (columns == 1) {
    dividename=rname[length(rname)]
    #match() is similar to %in% and so this function converts the name of the metab found above into a number (i.e. 1/4 of the total number of metabs)
    ind=match(dividename,rname)
    #these variables allow us to split the forest plot into four columns each containing a quarter of the values
    ind1=1:ind;
  }
  
  #This selects the colour pallette based on whether we entered a value for "coll" (it is set to NULL)
  #and will give you the same number of colours as there are data.frames in the list
  if (is.null(coll)) coll=c(brewer.pal(8,'Set1'))[1:length(result)]
  #NOT SURE ABOUT THIS
  pchh=rep(21,length(result))
  #background colour
  col.bgg=coll
  #setting mfrow as c(1,4) gives you 4 columns in the plot 
  par(mfrow=c(1,columns),mar=c(5,1,2,1))
  #This for loop produces each column
  for (s in 1:columns)
  {
    if (s==1) scutind=ind1
    if (s==2) scutind=ind2
    if (s==3) scutind=ind3
    if (s==4) scutind=ind4
    #These produce matrices containing each of the variables from the regression analysis
    va=sapply(result,function(x){x[scutind,'Estimate']})
    vlb=sapply(result,function(x){x[scutind,'2.5 %']})
    vub=sapply(result,function(x){x[scutind,'97.5 %']})
    vp=sapply(result,function(x){x[scutind,'Pr(>|t|)']})
    #This inputs the row names of the matrices 
    rownames(va)=rownames(vlb)=rownames(vub)=rownames(vp)=mnames=rname[scutind]
    
    #the background colour for the points - UNSURE ABOUT THIS BIT
    if (!is.null(sigP))
    {  
      bgg=vp
      for (i in 1:dim(bgg)[1])
      {
        for (j in 1:dim(bgg)[2])
        {
          bgg[i,j]=if (is.na(vp[i,j])) c(NA) else if (vp[i,j]<sigP) coll[j] else 'white' 
        }
      }
    } else
    {
      bgg=matrix(NA,nrow = dim(vp)[1],ncol=dim(vp)[2])
      bgg=apply(bgg,1,function(x){x[1:length(x)]=col.bgg})
      bgg=t(bgg)
    }
    
    
    #plot - UNSURE ON MOST OF THIS
    cexall=1
    cexadj=0.8 
    #Object that can be used to place the metab names in the correct place
    yv0=seq(1,by=2.5,length.out=dim(va)[1]);
    yv0=rev(yv0); 
    yrange=c(min(yv0)-1.25,max(yv0)+1.25)
    xrange=range(pretty(xrange))
    offset=seq(2.5,0,length.out=length(result)+2)
    offset=offset[-c(1,length(offset))]
    for (j in 1:length(result))
    {    
      if (j==1)     
      {
        
        #plot	
        plot(va[,j],yv0,col="transparent",xlab=NA,ylab=NA,xaxt="n",yaxt="n",xlim=c(xrange[1]-abs(xrange[1])*0.3,xrange[2]),ylim=yrange,bty="n",yaxs='i',xpd=T); 
        
        #the rectangular bars
        for (k in 1:length(yv0))
          rect(par('usr')[1],yv0[k]-1.25,par('usr')[2],yv0[k]+1.25,col=rep(c('white','grey92'),length(yv0))[k], border='transparent')
        
        #
        cut=pretty(xrange)
        ind=which.min(abs(cut))
        mid=cut[ind]
        cutd=cut[-which(cut==mid)];
        #abline() sets a straight line in the plot - here they set the lines going down through the
        #the plots and the middle line at y=0
        abline(v=cutd,lty=2,col="grey",ylim=c(min(yv0)-1.25,max(yv0)+1.25));
        abline(v=mid,col="black", ylim=c(min(yv0)-1.25,max(yv0)+1.25));	
        #sets the axes
        axis(side=1,line=-1.2,at=cut,labels=cut,tick=F,cex.axis=cexadj);
        axis(side=2,line=-1,at=yv0,labels=mnames,cex.axis=0.8,las=1,tick=F,hadj=0);
        #sets the line at the bottom of the graph
        lines(c(min(cut),par('usr')[2]),c(min(yv0)-1.25,min(yv0)-1.25)) 
        yv=yv0-1.25+offset[j]; 
        #segments() draws line segments between pairs of points - here it is set so that it shows the 95% confident intervals
        segments(vlb[,j],yv,vub[,j],yv,col=coll[j],lwd=1.5,xpd=NA)
        cnum=length(result)
        if (cnum<=2) cex.p=0.9 else if (cnum>2&cnum<=4) cex.p=0.7 else cex.p=0.6
        #points() adds points to the graph
        points(va[,j],yv,pch=pchh[j],bg=bgg[,j],col=coll[j],xpd=NA,lwd=1.2,cex=cex.p); 
      } else
      {
        yv=yv0-1.25+offset[j];
        segments(vlb[,j],yv,vub[,j],yv,col=coll[j],lwd=1.5,xpd=NA) 
        points(va[,j],yv,pch=pchh[j],bg=bgg[,j],col=coll[j],cex=cex.p,xpd=NA,lwd=1.2);          
      }
    }
  }
  #mtext() allows to write text in the margins of the plot
  mtext(text = 'Beta coefficients (95%CI)',line = -4,side=1,cex=1,outer = T)
  
  if(!is.null(sigP)) mtext(text = paste('Closed symblols: P < ',format(sigP, digits = 3, scientific = TRUE),';   Open symbols: P >= ', format(sigP, digits = 3, scientific = TRUE),sep=''),line = -2.8,side=1,cex=1,outer = T)
  #gives the plot a legend
  if(!is.null(names(result))) {
    x <- length(cname)
    if (columns == 4) {
    legend('bottomleft',legend=cname,col=coll,pch=pchh,pt.bg = coll,lty=1,lwd=1,border='transparent',bg='transparent', box.col = 'transparent',horiz = T,cex=1.1,xpd=NA,inset=c(-2.0, -0.075),text.width=c(rep(max(strwidth(cname))*1.8,length(result))))
    } else if (columns == 2) {
      legend('bottomleft',legend=cname,col=coll,pch=pchh,pt.bg = coll,lty=1,lwd=1,border='transparent',bg='transparent', box.col = 'transparent',horiz = T,cex=1.1,xpd=NA,inset=c(-1,-0.12),text.width=c(rep(max(strwidth(cname))*1.8,length(result))))
    } else {
      legend('bottom', legend=cname, col=coll,pch=pchh,pt.bg = coll,lty=1,lwd=1,border='transparent',bg='transparent', box.col = 'transparent',horiz = T,cex=1.1,xpd=NA, inset=c(0, -0.5), text.width=c(rep(max(strwidth(cname))*1.8,length(result))))
    }
  }
}














