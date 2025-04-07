pdf(file='BafA1_CRISPRa.pdf',width=4.5,height=4.5);
gstable=read.table('BafA1_CRISPRa.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("FPR2","FCGR2A","PLAUR","SELPLG","TREM1","ITGAX","CD19","NLRP3","HLA-DRB5","TLR4")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='0,1_vs_2,3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(856.0415878971099,1336.8753140843855,841.0936070979322,855.7168611173627),c(1037.9767165288054,1040.7162525946273,1178.7156888203417,928.947781086202),c(723.5339219572625,635.5773117364471,704.860135525732,718.7479182126818),c(681.4679962620728,683.8279453499471,693.0137466933667,777.0614285582391),c(637.2987742821236,692.9789275869902,611.7813661285765,653.6537671292691))
targetgene="FPR2"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(615.214163292149,547.3951192703955,508.54854916082223,579.0667190128587),c(939.1217911451096,892.6367218497484,810.6314643861359,747.2266093116749),c(781.3745697881483,841.0584583318691,719.2450362507469,659.0782797195535),c(674.1064592654145,678.836500493378,834.3242420508664,747.2266093116749),c(853.9382916123504,781.1611200530417,760.7073971640252,786.5543255912368))
targetgene="FCGR2A"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(452.20870122328904,536.5803220811628,514.4717435770049,511.2603116343039),c(635.1954779973642,665.5259808758609,601.6273185579777,671.2834330476934),c(499.5328676303774,628.9220519276885,556.7802751211664,557.3686686517212),c(723.5339219572625,790.3121022900848,632.089461269774,691.6253552612599),c(564.7350524579214,697.9703724435591,512.7794023152384,435.3171353703224))
targetgene="PLAUR"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(605.7493300107313,595.6457528838955,561.8572989064658,523.4654649624438),c(740.3602922353383,684.6598528260419,740.3993020228277,661.7905360146957),c(749.825125516756,707.1213546806022,638.0126556859567,600.7647693739963),c(746.6701810896168,576.5118809337145,598.2426360344448,568.21769383229),c(666.7449222687563,673.8450556368091,649.8590445183219,591.2718723409986))
targetgene="SELPLG"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(637.2987742821236,687.1555752543263,550.8570807049838,682.1324582282622),c(813.9756622019203,748.7167284853434,810.6314643861359,637.3802293584159),c(655.1767927025792,565.6970837444817,529.702814932903,644.1608700962714),c(658.3317371297185,700.4660948718436,599.9349772962112,560.0809249468634),c(660.435033414478,494.15304080032666,521.2411086240706,522.1093368148727))
targetgene="TREM1"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(658.3317371297185,680.5003154455677,595.7041241417951,751.2949937543882),c(621.5240521464275,638.9049416408265,583.0115646785466,535.6706182905837),c(695.1394221130095,668.8536107802402,571.1651758461813,645.5169982438425),c(641.5053668516425,650.551646306154,780.1693216743396,629.2434604729893),c(523.7207749051115,623.0986995950248,477.24023581814265,450.23454499360446))
targetgene="ITGAX"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(717.224033102984,693.8108350630849,625.3200962227082,697.0498678515443),c(699.3460146825283,731.2466714873522,693.0137466933667,717.3917900651107),c(781.3745697881483,853.5370704732916,660.859262719804,604.8331538167096),c(626.7822928583262,470.0277239935767,582.1653940476633,560.0809249468634),c(598.3877930140732,589.8224005512317,603.3196598197442,602.1208975215674))
targetgene="CD19"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(826.5954399104771,743.7252836287745,693.0137466933667,836.7310670513673),c(902.3141061618186,841.0584583318691,743.7839845463606,772.9930441155258),c(793.9943474967051,692.1470201108953,707.3986474183816,673.9956893428356),c(691.9844776858702,741.22956120049,664.2439452433368,709.2550211796841),c(782.426217930528,750.3805434375331,632.089461269774,650.9415108341269))
targetgene="NLRP3"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(638.3504224245033,879.3262022322311,779.3231510434563,652.297638981698),c(697.2427183977688,614.7796248340765,759.0150559022587,697.0498678515443),c(669.8998666958956,600.6371977404644,579.6268821550136,566.8615656847189),c(808.7174214900215,941.7192629393431,852.9399959302975,816.3891448378008),c(677.2614036925538,708.7851696327918,617.7045605447591,682.1324582282622))
targetgene="HLA-DRB5"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(743.5152366624775,790.3121022900848,775.0922978890402,736.3775841311061),c(732.9987552386801,821.9245863816882,869.8634085479621,800.1156070669477),c(687.7778851163512,634.7454042603523,560.1649576446994,568.21769383229),c(800.3042363509836,760.3634331506711,639.7049969477231,645.5169982438425),c(769.8064402219711,736.2381163439211,727.7067425595792,724.1724308029662))
targetgene="TLR4"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("SPI1","TGFBR3","LPAR2","IL15RA","CSF1R","AIF1","SORL1","MERTK","TGFBR2","CD52")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='0,1_vs_2,3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(531.0823119017697,438.41523990197334,863.0940435008963,645.5169982438425),c(563.6834043155417,471.69153894576635,661.7054333506871,726.8846870981084),c(600.4910892988327,660.5345360192919,836.0165833126329,767.5685315252414),c(609.9559225802503,579.8395108380938,779.3231510434563,915.386499610491),c(454.31199750804853,397.65177357332686,716.7065243580972,549.2318997662946))
targetgene="SPI1"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(770.8580883643508,722.0956892503091,772.5537859963905,716.0356619175396),c(558.425163603643,517.4464501309818,769.1691034728575,612.9699227021362),c(614.1625151497693,569.024713648861,573.703687738831,490.9183894207374),c(623.627348431187,612.283902405792,660.859262719804,591.2718723409986),c(443.79551608425106,504.9678379895594,720.9373775125133,574.9983345701455))
targetgene="TGFBR3"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(715.1207368182245,763.6910630550503,676.0903340757021,669.9273049001223),c(726.6888663844016,765.35487800724,796.246563661121,810.9646322475164),c(626.7822928583262,805.2864368597916,737.860790130178,732.3091996883928),c(709.8624961063258,596.4776603599903,1052.63626481874,787.9104537388079),c(553.1669228917442,712.944707013266,682.0135284918847,650.9415108341269))
targetgene="LPAR2"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(534.2372563289089,526.5974323680249,504.31769600640604,576.3544627177166),c(573.1482375969593,585.6628631707575,635.474143793307,599.4086412264252),c(883.3844395989832,767.0186929594296,773.3999566272737,893.6884492493534),c(677.2614036925538,681.3322229216626,934.1723764950877,969.6316255133349),c(782.426217930528,813.6055116207399,808.9391231243694,762.144018934957))
targetgene="IL15RA"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1042.1833090983243,910.9386863238345,1396.1815409573321,1338.4984816526735),c(912.830587585616,985.8103591723689,1226.1012441498026,1174.4069757965706),c(871.8163100328061,846.0499031884381,924.018328924489,922.1671403483465),c(922.2954208670336,757.0358032462917,689.6290641698338,762.144018934957),c(933.8635504332108,991.6337115050327,966.3268604686506,1026.589007711321))
targetgene="CSF1R"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(938.0701430027298,920.0896685608776,892.7100155818093,847.5800922319361),c(736.1536996658193,686.3236677782315,792.0157105067049,623.818947882705),c(652.02184827544,626.426329499404,643.0896794712561,691.6253552612599),c(573.1482375969593,559.873731411818,833.4780714199832,682.1324582282622),c(663.5899778416172,651.3835537822488,611.7813661285765,791.9788381815212))
targetgene="AIF1"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(677.2614036925538,680.5003154455677,639.7049969477231,687.5569708185466),c(697.2427183977688,521.605987511456,588.088588463846,650.9415108341269),c(614.1625151497693,559.873731411818,543.2415450270347,505.8357990440195),c(459.5702382199472,506.63165294174905,721.7835481433966,503.1235427488773),c(557.3735154612632,664.694073399766,547.4723981814509,633.3118449157026))
targetgene="SORL1"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(525.824071189871,643.8963864973954,673.5518221830524,569.5738219798611),c(691.9844776858702,653.8792762105332,646.474361994789,625.175076030276),c(1028.5118832473877,1007.4399535508344,1135.560986645297,1166.270206911144),c(130.404369655088,183.01964474086174,198.00392762667622,185.7895562172404),c(675.1581074077943,759.5315256745762,531.3951561946694,698.4059959991154))
targetgene="MERTK"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(700.3976628249081,868.5114050429984,833.4780714199832,943.8651907094841),c(805.5624770628823,803.622621907602,742.9378139154774,659.0782797195535),c(589.9746078750352,552.3865641269645,739.5531313919445,695.6937397039732),c(600.4910892988327,577.3437884098093,750.5533495934264,684.8447145234044),c(574.199885739339,597.3095678360851,550.8570807049838,531.6022338478704))
targetgene="TGFBR2"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(756.1350143710345,668.0217033041454,672.7056515521691,726.8846870981084),c(777.1679772186293,757.0358032462917,639.7049969477231,692.981483408831),c(755.0833662286547,655.5430911627229,881.7097973803274,892.3323211017823),c(829.7503843376164,769.5144153877142,786.0925160905222,899.1129618396378),c(408.0394792433399,404.3070333820855,368.930395065089,408.1945724189004))
targetgene="CD52"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd1_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-baf-10-rd2_S20_L001_R1_001-clipped-trimmed-aligned-counts")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("BafA1_CRISPRa_summary.Rnw");
library(tools);

texi2dvi("BafA1_CRISPRa_summary.tex",pdf=TRUE);

