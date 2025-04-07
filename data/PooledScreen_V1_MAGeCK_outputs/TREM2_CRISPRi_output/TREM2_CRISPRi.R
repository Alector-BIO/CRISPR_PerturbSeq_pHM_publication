pdf(file='TREM2_CRISPRi.pdf',width=4.5,height=4.5);
gstable=read.table('TREM2_CRISPRi.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("TREM2","TYROBP","SPI1","FCGR2A","SIGLEC9","GPR141","PTPRC","LRP1","CD226","CD48")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2_vs_0,3 neg.'


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
targetmat=list(c(1877.538369023071,1811.761038767669,385.66669354020854,464.2082168420685),c(2871.146545501122,2767.9017719085987,170.67401352268396,267.15656561114963),c(3760.6535846868046,3252.553833137731,206.50612685893807,302.209022801265),c(2314.3913558749737,1844.0711761829446,413.9552040688302,430.1031233597941),c(3176.0141618573434,2418.473619121175,286.65690669003277,347.68248077763093))
targetgene="TREM2"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(2079.2343225270347,2436.423695462995,215.93563036847863,309.78793246399266),c(1381.1990180041007,1511.3964279812194,284.7710059881246,233.05147212887522),c(1381.1990180041007,1519.7731302740685,298.91526125243547,283.2617486444459),c(1344.949514839794,1784.2375883768789,178.21761633031642,237.78829066807998),c(1738.1172030065063,2307.1831458018933,498.82073565469517,503.0501288635477))
targetgene="TYROBP"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(711.0479466844799,609.1059238657488,245.16709124805433,269.9986567346725),c(517.7172631415102,605.5159085973849,256.482495459503,274.7354752738773),c(750.085873169118,546.8789925474404,259.31134651236516,323.998388081607),c(747.2974498487868,624.6626566953258,300.80116195434357,362.8403001030862),c(738.9321798877928,585.1724887433224,216.87858071943268,274.7354752738773))
targetgene="SPI1"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(802.1364418153022,587.565832255565,612.9177281201358,513.4711296497982),c(1050.3061173247872,758.689893380913,710.9845646193576,727.5753276218543),c(1116.298802572628,879.5537407491656,808.1084507676252,701.9965075101485),c(950.8523522329712,819.7201529431,680.8101533888278,627.1547745907129),c(973.1597387956215,816.130137674736,820.3668053300279,653.6809584102598))
targetgene="FCGR2A"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(787.2648507735353,728.7730994778801,672.3236002302413,529.5763126830944),c(765.8869386509954,935.7973132868674,671.3806498792873,533.3657675144583),c(843.9627916202716,798.1800613329164,614.803628822044,498.3133103243429),c(952.711301113192,922.6339239695329,841.1117130510172,898.1007950332264),c(673.8689690800627,843.6535880655263,955.2087055164578,653.6809584102598))
targetgene="SIGLEC9"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(716.6247933251425,767.0665956737622,693.0685079512305,412.1032129108159),c(896.9428347065661,934.6006415307461,805.2795997147631,667.891414027874),c(722.201639965805,1106.9213744122153,762.8468339218306,824.2064258216318),c(673.8689690800627,613.892610890234,694.0114583021846,564.6287698732099),c(647.8436847569706,665.3494964034505,512.964990919006,470.8397627969552))
targetgene="GPR141"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(759.3806175702223,588.7625040116864,631.776735139217,539.0499497615041),c(693.3879323223817,716.806381916667,646.8639407544817,588.3128625692337),c(773.3227341718788,787.4100155278245,647.8068911054359,608.2075004338939),c(918.320746829106,806.5567636257655,868.4572732286847,519.155311896844),c(716.6247933251425,737.1498017707294,624.2331323315844,509.6816748184344))
targetgene="PTPRC"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(839.3154194197194,832.8835422604344,741.1589758498873,878.2061571685664),c(823.5143539378421,855.6203056267394,798.6789472580847,1138.7311768248294),c(946.204980032419,1063.8411911918479,729.8435716384387,783.4697863844707),c(824.4438283779525,737.1498017707294,662.8940967207008,607.2601367260529),c(1289.1810484331681,1419.2527027598783,764.7327346237387,803.3644242491308))
targetgene="LRP1"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(271.4065365122459,309.93798483542025,315.8883675696085,142.10455617614343),c(307.6560396765527,236.94100771202014,273.45560177667596,84.3153699978451),c(751.9448220493389,807.7534353818868,722.2999688308062,381.78757425990534),c(707.3300489240381,716.806381916667,644.9780400525736,631.8915931299177),c(772.3932597317684,800.573404845159,880.7156277910875,823.2590621137908))
targetgene="CD226"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(725.9195377262467,910.6672064083198,916.5477411273415,981.4688013232305),c(752.8742964894493,800.573404845159,615.7465791729979,564.6287698732099),c(837.4564705394986,920.2405804572903,942.9503509540551,893.3639764940216),c(1133.0293424946158,1096.1513286071236,839.225812349109,756.9436025649239),c(883.9301925450201,1036.3177408010579,966.5241097279065,939.7847981782285))
targetgene="CD48"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("SELPLG","IRF8","TYRO3","DPYD","CD37","ADGRG1","PLCB3","NFATC2","GRN","CD300A")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='1,2_vs_0,3 pos.'


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
targetmat=list(c(808.6427628960752,987.2541988000839,1029.7017832418283,1208.83609120506),c(746.3679754086763,877.160397236923,1102.3089602652904,998.5213480643678),c(748.2269242888972,840.0635727971624,1058.933244121404,1103.6787196347138),c(702.682676723486,825.7035117237066,857.1418690172361,774.9435130139021),c(748.2269242888972,807.7534353818868,958.03755656932,977.6793464918667))
targetgene="SELPLG"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(546.5309707849335,744.3298323074572,893.9169327044442,972.942527952662),c(738.0027054476824,601.9258933290208,800.5648479599928,856.4167918882243),c(791.9122229740874,1047.0877866061496,1022.1581804341957,1057.257897950507),c(633.9015681553141,630.6460154759324,815.6520535752577,1007.9949851427773),c(667.3626479992896,720.3963971850309,950.4939537616875,925.5743425606141))
targetgene="IRF8"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(685.0226623613878,775.4432979666113,839.225812349109,978.6267101997076),c(720.3426910855842,744.3298323074572,776.0481388351874,776.8382404295841),c(827.2322516982838,946.5673590919592,1003.2991734151146,1088.5209003092586),c(858.8343826620385,796.983389576795,1006.1280244679768,1227.783365361879),c(654.3500058377435,670.1361834279357,669.4947491773792,999.4687117722087))
targetgene="TYRO3"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(747.2974498487868,767.0665956737622,841.1117130510172,965.3636182899343),c(798.4185440548604,902.2905041154706,1152.2853288658553,1145.362722779716),c(1024.2808330016953,903.487175871592,971.2388614826767,917.9954328978865),c(958.2881477538546,877.160397236923,1312.5868885280447,1230.625456485402),c(828.1617261383942,844.8502598216476,906.1752872668469,813.7854250353813))
targetgene="DPYD"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(742.6500776482345,767.0665956737622,892.0310320025361,960.6267997507296),c(586.4983717096821,840.0635727971624,1003.2991734151146,816.6275161589042),c(883.9301925450201,848.4402750900115,1045.7319392080472,1115.9944478366463),c(755.6627198097806,664.1528246473292,794.9071458542685,865.8904289666339),c(798.4185440548604,786.2133437717032,937.2926488483308,902.8376135724312))
targetgene="CD37"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(848.6101638208237,811.3434506502508,956.1516558674119,1065.7841713210757),c(745.4385009685659,668.9395116718144,670.4376995283332,816.6275161589042),c(994.5376509181614,901.0938323593493,1065.5338965780822,1247.6780032265392),c(614.382604912995,706.0363361115751,987.2690174488957,761.6804211041288),c(909.0260024280018,974.0908094827495,943.8933013050091,1025.9948955917555))
targetgene="ADGRG1"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(678.5163412806148,919.043908701169,778.8769898880495,857.3641555960653),c(401.5329581277063,433.1951757159156,654.4075435621143,835.5747903157234),c(714.7658444449216,650.9894353299948,790.1923940994982,661.2598680729874),c(780.7585296927622,708.4296796238178,884.4874291949037,735.154237284582),c(720.3426910855842,769.4599391860048,848.6553158586496,988.1003472781173))
targetgene="PLCB3"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(877.4238714642471,822.1134964553426,698.7262100569549,1045.8895334564156),c(825.3733028180629,954.9440613848084,1062.7050455252202,1321.572372438134),c(656.2089547179645,794.5900460645523,719.4711177779441,1134.9417219934655),c(759.3806175702223,653.3827788422374,756.2461814651522,899.9955224489083),c(697.1058300828234,793.393374308431,859.0277697191442,749.3646929021963))
targetgene="NFATC2"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(898.801783586787,805.3600918696442,1103.2519106162445,1477.8873842318915),c(801.2069673751918,686.8895880136341,757.1891318161063,1029.7843504231193),c(843.9627916202716,737.1498017707294,644.0350897016197,904.7323409881132),c(701.7532022833756,745.5265040635785,839.225812349109,787.2592412158345),c(946.204980032419,920.2405804572903,797.7359969071306,839.3642451470871))
targetgene="GRN"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(732.4258588070198,801.7700766012803,799.6218976090387,746.5226017786734),c(853.2575360213759,710.8230231360604,882.6015284929956,810.9433339118584),c(827.2322516982838,685.6929162575128,902.4034858630307,991.889802109481),c(1009.4092419599283,1148.8048858764612,1005.1850741170227,1114.0997204209643),c(510.28146762062676,510.978839863801,743.9878269027495,826.1011532373137))
targetgene="CD300A"
collabel=c("TRE-KRAB-7d-dox-rd2-trem2-low_S14_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-low_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd1-trem2-high_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-KRAB-7d-dox-rd2-trem2-high_S15_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("TREM2_CRISPRi_summary.Rnw");
library(tools);

texi2dvi("TREM2_CRISPRi_summary.tex",pdf=TRUE);

