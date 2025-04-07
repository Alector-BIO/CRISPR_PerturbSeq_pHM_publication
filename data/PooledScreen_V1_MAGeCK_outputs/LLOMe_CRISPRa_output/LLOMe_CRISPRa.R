pdf(file='LLOMe_CRISPRa.pdf',width=4.5,height=4.5);
gstable=read.table('LLOMe_CRISPRa.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("LRRK2","NLRP3","FCGR2A","MMP9","CD300A","LTBR","IL1R2","ICAM3","CD86","EPHA1")
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
targetmat=list(c(805.005234870796,663.8215998314738,520.4467774040182,664.6541326372894),c(762.459756532701,729.3186643481793,694.2421814107752,770.3945628295854),c(725.5123674496186,730.2037598146212,818.247550756137,789.8162744975582),c(835.2349168478635,796.5859197977686,575.8734197629299,664.6541326372894),c(773.655935042726,900.1420893714785,742.153346839665,643.0744530062085))
targetgene="LRRK2"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(880.0196308879633,791.2753469991168,729.9406968283794,677.6019404159379),c(960.6321161601431,894.8315165728268,859.5826738712575,738.0250433829642),c(845.3114775068859,736.399428079715,700.8182237245444,681.917876342154),c(736.7085459596436,788.6200605997909,686.7267044807533,725.0772356043157),c(832.9956811458584,798.3561107306525,776.9124276410164,625.8107093013439))
targetgene="NLRP3"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(654.9764428364613,582.3928169188131,581.5100274604464,605.3100136518171),c(999.8187409452306,949.7074354922286,801.3377276635877,718.6033317149914),c(831.876063294856,894.8315165728268,713.03087373583,708.8924758810051),c(717.6750424926012,722.2379006166435,863.3404123362685,704.5765399547889),c(909.1296950140282,831.1046429890052,746.8505199209287,718.6033317149914))
targetgene="FCGR2A"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(841.9526239538784,819.5984019252597,817.3081161398842,981.8754232141774),c(802.765999168791,810.7474472608401,730.8801314446321,767.1576108849233),c(837.4741525498684,1055.0337959988224,930.9797047064658,863.187185243233),c(531.8184792261865,793.9306333984428,431.2004888600079,469.3580319760079),c(695.2826854725512,636.3836403717729,777.8518622572691,650.6273408770868))
targetgene="MMP9"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(651.6175892834538,694.7999411569426,661.3619698419293,555.6767505003312),c(613.5505823493688,576.1971486537193,475.35391582388667,543.8079267032367),c(780.373642148741,677.0980318281033,623.7845851918197,638.7585170799924),c(547.4931291402214,579.7375305194871,444.35257348754624,493.0956795701968),c(671.7707106014988,682.4086046267552,654.7859275281601,762.8416749587071))
targetgene="CD300A"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1228.22078254974,1168.326015703394,1033.3780778780144,1006.6920547899205),c(947.1967019481132,831.1046429890052,759.0631699322144,880.4509289480976),c(932.6416698850808,921.3843805660857,825.7630276861589,892.3197527451921),c(697.5219211745562,677.0980318281033,818.247550756137,857.7922653354627),c(1212.5461326357051,923.1545714989696,928.1614008577076,920.3733362655971))
targetgene="LTBR"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(597.8759324353339,706.3061822206881,578.6917236116882,771.4735468111395),c(691.9238319195437,637.2687358382149,804.1560315123459,864.266169224787),c(951.6751733521232,792.1604424655588,656.6647967606656,910.6624804316108),c(638.1821750714238,813.402733660166,501.65808507896344,548.1238626294529),c(752.3831958736786,777.9989150024874,802.2771622798404,690.5497481945863))
targetgene="IL1R2"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(548.612746991224,649.6600723684024,546.750946659095,523.30723105371),c(890.0961915469858,833.7599293883311,778.7912968735219,661.4171806926272),c(794.9286742117735,776.2287240696035,802.2771622798404,865.345153206341),c(789.330584956761,680.6384136938711,759.0631699322144,749.8938671800586),c(876.6607773349558,800.1263016635364,685.7872698645006,649.5483568955328))
targetgene="ICAM3"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(697.5219211745562,756.7566238078801,701.7576583407972,653.8642928217489),c(746.7851066186661,722.2379006166435,757.1843006997088,557.8347184634393),c(839.7133882518734,760.2970056736481,726.1829583633684,754.2098031062748),c(719.9142781946061,690.3744638247329,588.0860697742155,554.5977665187771),c(626.9859965613988,612.4860627778398,590.9043736229738,612.8629015226954))
targetgene="CD86"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(667.2922391974887,700.1105139555945,715.8491775845882,633.3635971722222),c(850.9095667618984,902.7973757708045,779.7307314897746,951.6638717306644),c(539.6558041832039,708.0763731535721,664.1802736906875,523.30723105371),c(913.6081664180383,1182.4875431664655,685.7872698645006,893.398736726746),c(727.7516031516236,751.4460510092284,799.4588584310821,784.4213545897879))
targetgene="EPHA1"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("CD40","PTPRE","FCER1G","TGFBR2","IRF8","CTSB","LDLR","ITGB3","CD68","TMEM106B")
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
targetmat=list(c(754.6224315756835,689.4893683582908,740.2744776071595,679.7599083790459),c(626.9859965613988,583.277912385255,658.5436659931711,688.3917802314783),c(644.8998821774388,616.9115401100497,593.722677471732,567.5455742974257),c(573.2443397132789,605.4052990463041,420.8667080812277,468.2790479944539),c(1293.1586179078852,845.2661704520767,1844.1101517041295,1689.6889151136286))
targetgene="CD40"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(660.5745320914738,669.1321726301256,890.584016207598,914.978416357827),c(819.5602669338284,919.6141896332018,801.3377276635877,820.0278259810713),c(1054.680015644353,874.4743208446615,915.948750846422,966.769647472421),c(732.2300745556336,839.0705021869829,740.2744776071595,702.4185719916808),c(694.1630676215486,724.8931870159695,692.3633121782698,713.2084118072212))
targetgene="PTPRE"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(849.7899489108959,721.3528051502016,927.2219662414549,1043.3775101627577),c(799.4071456157835,681.5235091603131,906.5544046838946,1094.0897572957977),c(661.6941499424762,760.2970056736481,675.4534890857204,758.525739032491),c(868.8234523779383,904.5675667036884,934.7374431714768,894.4777207083001),c(1030.048422922298,769.1479603380677,969.4965239728282,952.7428557122183))
targetgene="FCER1G"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(745.6654887676636,924.0396669654116,1069.0765932956187,1182.566443783229),c(857.6272738679133,855.0022205829383,835.1573738486863,787.65830653445),c(628.1056144124013,587.7033897174648,690.4844429457643,730.4721555120859),c(639.3017929224263,614.2562537107239,741.2139122234122,713.2084118072212),c(611.3113466473638,635.498544905331,596.5409813204902,826.5017298703956))
targetgene="TGFBR2"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(819.5602669338284,881.5550845761973,832.339069999928,937.6370799704617),c(953.9144090541282,890.4060392406169,931.9191393227186,984.0333911772856),c(609.0721109453589,671.7874590294515,639.7549736681162,599.9150937440469),c(910.2493128650308,1060.3443687974743,958.2233085777953,1078.9839815540413),c(797.1679099137785,827.5642611232374,1004.2556047741796,1183.645427764783))
targetgene="IRF8"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(484.7945294840816,551.4144755933443,650.0887544468965,708.8924758810051),c(435.53134403997166,482.377029210871,411.47236191870036,484.4638077177645),c(622.5075251573888,662.0514088985899,636.936669819358,762.8416749587071),c(606.8328752433539,737.284523546157,651.0281890631492,740.1830113460722),c(824.0387383378385,808.9772563279562,1025.8626009479926,1022.876814513231))
targetgene="CTSB"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(822.9191204868359,803.6666835293044,812.6109430586205,920.3733362655971),c(798.287527764781,829.3344520561213,657.6042313769184,814.6329060733011),c(547.4931291402214,581.5077214523711,605.9353274830177,554.5977665187771),c(726.6319853006212,664.7066952979158,844.5517200112137,988.3493271035018),c(747.9047244696686,839.0705021869829,1043.7118586567947,917.136384320935))
targetgene="LDLR"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(705.3592461315736,716.0422323515498,782.5490353385328,715.3663797703293),c(641.5410286244313,716.0422323515498,791.0039468848075,771.4735468111395),c(749.024342320671,696.5701320898265,812.6109430586205,694.8656841208025),c(784.852113552751,624.8773993080274,671.6957506207094,608.5469655964793),c(684.0865069625262,682.4086046267552,806.9743353611041,971.085583398637))
targetgene="ITGB3"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(695.2826854725512,771.8032467373936,866.1587161850267,813.553922091747),c(753.5028137246811,757.6417192743221,791.9433815010602,771.4735468111395),c(625.8663787103964,738.1696190125989,928.1614008577076,809.2379861655309),c(657.2156785384662,659.396122499264,771.2758199435,698.1026360654646),c(996.4598873922231,985.1112541499072,1074.7132009931352,1007.7710387714744))
targetgene="CD68"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(901.2923700570108,875.3594163111035,836.096808464939,836.2125857043819),c(747.9047244696686,797.4710152642106,984.527477832872,978.6384712695153),c(635.9429393694188,713.3869459522239,676.3929237019731,686.2338122683702),c(787.091349254756,793.9306333984428,917.8276200789275,717.5243477334374),c(793.809056360771,859.4276979151481,789.125077652302,728.3141875489778))
targetgene="TMEM106B"
collabel=c("TRE-VPH-7d-dox-untreated-rd2_S18_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-untreated-rd1_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd1_S13_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-llome-100-rd2_S19_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("LLOMe_CRISPRa_summary.Rnw");
library(tools);

texi2dvi("LLOMe_CRISPRa_summary.tex",pdf=TRUE);

