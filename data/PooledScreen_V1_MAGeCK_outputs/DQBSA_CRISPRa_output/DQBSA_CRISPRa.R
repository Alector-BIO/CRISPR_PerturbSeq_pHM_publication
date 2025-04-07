pdf(file='DQBSA_CRISPRa.pdf',width=4.5,height=4.5);
gstable=read.table('DQBSA_CRISPRa.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ICAM3","CD44","LILRA5","ITGAX","IL15RA","IFNGR1","CD200R1","TGFBI","LSR","TGFB1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='4,5_vs_1,3 neg.'


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
targetmat=list(c(1008.4480291208142,891.1425973118014,546.2458310458205,501.8713334462032),c(2195.067088121631,2061.941185349086,1118.7114619818403,1045.6276067925367),c(1482.3829686016209,1351.3227910056005,876.178312997496,881.0796286134536),c(709.1206989224099,677.2266342568608,379.0946067457994,371.72884161365573),c(2788.970521054973,2782.994504719642,1739.2467260498925,1673.90170529449))
targetgene="ICAM3"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1331.531496636552,1122.7979329127616,761.4666884778737,860.1371586633885),c(110.46603852560156,78.26193770302706,85.214349643148,85.26577051097938),c(425.23485798027264,445.5712986559007,176.98364925884584,277.48772683836273),c(943.1186514981466,990.2743850689691,1131.82136192694,1003.7426668924064),c(2141.61577915763,1554.803829033471,1105.6015620367407,1112.1947434195292))
targetgene="CD44"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(391.97626573600553,352.70046591497527,334.3024486000421,392.67131156372085),c(792.2671795330779,924.5343573984263,908.9530628602453,1203.4440767733845),c(840.9672610336119,725.227289381384,929.7104044399864,953.6303280833221),c(425.23485798027264,494.615446283131,99.41674125033933,109.20002188248236),c(353.96644602827166,526.9637138670488,390.0195233667158,440.5398143067268))
targetgene="LILRA5"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1622.5441787738896,1612.1959166823574,1796.0562924786577,1603.5948418906999),c(880.1648876072123,891.1425973118014,734.1543969255827,828.7234537382908),c(1010.8236428525476,890.0991048090943,258.9205239157189,270.7562186401275),c(405.04214126053904,329.7436308554207,237.07069067388608,216.15620769888633),c(298.1395233325375,259.82963317404983,429.3492232020149,356.76993450646637))
targetgene="ITGAX"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(848.094102228812,1048.7099652205625,1800.4262591270242,1898.2853119023305),c(441.86415410240625,434.0928811261234,332.1174652758589,320.1206120938524),c(1822.0957322394922,1385.7580435949324,2087.751566257126,1669.4140331623332),c(1881.4860755328264,2031.6799027705824,694.8246970902836,836.2029072918855),c(1705.6906593845574,1996.2011576785435,3461.0135855063186,2935.685519785913))
targetgene="IL15RA"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1254.3240503552177,1212.538288145566,1681.3446679590354,1400.9016505882842),c(1294.709483794685,1206.2773331293238,703.5646303870168,711.2960329468543),c(2034.7131612296287,2093.245960430297,1384.186935870109,1277.4906669539719),c(1721.1321486408242,1627.8483042229627,2000.3522332897946,2278.9894977802996),c(232.81014570986994,226.43787308742495,191.18604086603716,192.22195632738334))
targetgene="IFNGR1"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(2216.4476117072313,1994.1141726731294,1088.1216954432743,1115.9344701963266),c(73.64402568373437,75.13146019490597,105.97169122288918,83.02193444490098),c(1608.2904963834892,1683.1534068664353,1560.0780934668633,1765.8989840037045),c(1351.7242133562856,1188.5379605833043,854.3284797556632,791.3261859703174),c(400.2909137970723,394.4401660232564,237.07069067388608,338.8192459778391))
targetgene="CD200R1"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1154.5482736224162,1282.4522858269368,1342.6722527106267,1144.3563936999865),c(691.3035959344097,583.3123090132283,260.01301557781056,281.2274536151601),c(2120.23525557203,1932.5481150134149,1603.777759950529,1733.7373337232475),c(630.7254457752089,578.0948464996932,664.2349305517178,672.4028744681619),c(1774.583457604825,1329.409448448753,1567.7255351015046,1296.1893008379584))
targetgene="TGFBI"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(2214.0719979754977,2166.290435619789,1343.7647443727183,1269.2632680450176),c(24.94394418320035,45.913670119109206,15.294883269282973,29.169868859019264),c(646.1669350314758,732.5317369003333,716.6745303321164,635.0056067001885),c(1522.7684020410882,1540.1949339955725,1303.3425528753276,1366.4961642417486),c(1207.999582586417,1221.929720669929,1433.3490606642329,1279.7345030200502))
targetgene="LSR"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1835.1616077640258,1719.6756444611813,1846.3109089348732,1930.4469621827875),c(2233.0769078293647,1987.8532176568872,1841.9409422865067,1783.1017271769724),c(1250.7606297576176,1166.6246180264566,741.8018385602242,848.9179783329964),c(454.93002962693976,498.7894162939591,233.79321568761117,355.27404379574745),c(501.2544973957404,534.268161385998,400.9444399876322,431.5644700424132))
targetgene="TGFB1"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("SCARB1","FCGR1A","FCGR2C","AXL","TREM2","SQSTM1","SIGLEC9","CSF1R","IL15RA","CCR2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='4,5_vs_1,3 pos.'


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
targetmat=list(c(1077.3408273410819,1254.277988253847,2069.179208001568,2087.515486808276),c(1548.900153090155,1318.9745234216828,2748.7090218225685,2887.817017042907),c(1054.7724968896148,1193.7554230968394,1959.930041792404,2220.649760062261),c(867.0990120826789,744.0101544301106,3534.2105268664586,3262.5376400780005),c(1787.6493331293584,1461.9329962925453,3526.563085231817,2768.145760185392))
targetgene="SCARB1"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1241.258174830684,945.4042074525669,4286.937282047599,3372.4856073158426),c(270.81996541760384,336.0045858716628,467.58643137522233,531.0412023052224),c(1413.4901703813532,1076.8842627936524,1910.76791699828,1895.2935304808925),c(1028.6407458405479,1004.8832801068675,1658.402343055111,1919.975727207755),c(377.7225833456053,396.5271510286704,611.795330771319,612.5672460394045))
targetgene="FCGR1A"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(592.715626067475,470.61511872086936,1266.1978363642118,1395.666033100768),c(59.39034329333417,86.60987772468327,146.3938827202799,126.40276505575014),c(574.8985230794748,545.7465789157753,1852.865858907423,1993.274372032983),c(649.7303556290758,586.4427865213494,690.454730441917,876.5919564812968),c(462.05687082213984,480.00655124523263,989.7974458550267,976.8166340994655))
targetgene="FCGR2C"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(465.6202914197399,480.00655124523263,869.6233630249462,1004.4906122477659),c(1101.0969646584156,1432.7152062167486,1885.6406087701723,1922.2195632738335),c(129.4709483794685,104.34925027070274,133.2839827751802,92.74522406457406),c(753.0695529594773,820.1851071277235,1579.7429433845127,1713.5428091285419),c(635.4766732386756,771.1409595004933,1477.0487271478985,1565.4496287673671))
targetgene="AXL"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(60.57815015920085,80.34892270844111,88.49182462942292,96.4849508413714),c(1573.8440972733554,1528.7165164657952,3167.1333284036673,3517.587006255579),c(586.7765917381416,661.5742467162554,1208.295778273355,1225.134492078809),c(921.7381279125464,936.0127749282036,1189.723420017797,1706.063355574947),c(549.9545788962744,523.8332363589278,1030.2196373524175,1555.726339147694))
targetgene="TREM2"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1312.5265867826852,1270.9738682971595,1229.0531198530962,1310.4002625897883),c(1873.1714274717597,1638.283229250033,7025.813878911343,5744.968274516075),c(800.5818275941446,858.7943297278836,1267.2903280263035,1104.7152898659347),c(1731.8224104336243,1641.4137067581541,2037.4969498009104,1769.638710780502),c(927.6771622418797,745.0536469328176,1293.5101279165028,1143.608448344627))
targetgene="SQSTM1"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(981.1284712058805,871.3162397603679,1679.159684634852,1492.8989292974986),c(1118.9140676464158,1027.840115166422,1985.0573500205116,1400.9016505882842),c(823.1501580456115,774.2714370086144,905.6755878739704,875.8440111259373),c(96.21235613520136,79.30543020573408,73.19694136013995,71.80275411450896),c(1481.1951617357543,1500.5422188927055,2270.1976738264298,2148.099060592393))
targetgene="SIGLEC9"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1965.820363009361,1893.9388924132547,2833.9233714657166,3490.6609734626386),c(261.3175104906703,277.5690057200693,430.44171486410653,441.2877596620863),c(1420.6170115765533,1965.9398751000397,1847.4034005969647,1971.5839567275584),c(2475.389508466168,1370.105656054327,2017.8320998832608,1433.0633008687412),c(1583.346552200289,1693.5883318935055,2292.047507068263,2520.575847561408))
targetgene="CSF1R"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(848.094102228812,1048.7099652205625,1800.4262591270242,1898.2853119023305),c(441.86415410240625,434.0928811261234,332.1174652758589,320.1206120938524),c(1822.0957322394922,1385.7580435949324,2087.751566257126,1669.4140331623332),c(1881.4860755328264,2031.6799027705824,694.8246970902836,836.2029072918855),c(1705.6906593845574,1996.2011576785435,3461.0135855063186,2935.685519785913))
targetgene="IL15RA"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(642.6035144338757,657.4002767054272,771.2991134366985,857.1453772419507),c(1513.2659471141546,1614.2829016877715,1532.7658019145722,1357.520819977435),c(2679.692289395238,2437.5984863236163,2907.1203128258567,2418.1073338771607),c(1152.1726598906828,1026.796622663715,3038.2193122768535,3001.504711057546),c(1662.9296122133567,1438.9761612329908,1517.4709186452892,1587.887989428151))
targetgene="CCR2"
collabel=c("TRE-VPH-undiff-DQBSA-Low-rep2_S11_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-Low-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep2_S12_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-undiff-DQBSA-High-rep1_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("DQBSA_CRISPRa_summary.Rnw");
library(tools);

texi2dvi("DQBSA_CRISPRa_summary.tex",pdf=TRUE);

