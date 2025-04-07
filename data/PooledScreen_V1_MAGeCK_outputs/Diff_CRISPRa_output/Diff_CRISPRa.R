pdf(file='Diff_CRISPRa.pdf',width=4.5,height=4.5);
gstable=read.table('Diff_CRISPRa.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("NFATC1","SLC24A4","FLT1","NLRP3","PILRA","VLDLR","CLDN7","CD300A","Non","TYRO3")
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
targetmat=list(c(1416.460260099681,934.0600411851273,893.4832505482183,800.0149950351924),c(1444.1217735965772,988.9470950218189,997.3520032809886,1067.6891851033583),c(1221.8758203284117,839.967948893656,825.6505957023275,820.0654961638941),c(1384.0295201378028,822.3256815890051,819.2912843105252,911.2952762994862),c(1280.0603832011932,911.517144073629,836.2494480219979,781.969544019361))
targetgene="NFATC1"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1096.922086945881,999.7284805968834,866.9861197490421,972.4493047420259),c(901.3838018816152,1045.7944007812496,887.123939156416,847.1336726876411),c(475.9688012052127,853.6897123528288,813.99185815069,724.8256158025616),c(1713.106146221567,856.6300902369373,853.2076117334706,874.2018492113882),c(1724.552289737524,968.3644498330596,928.4594632031308,813.0478207688485))
targetgene="SLC24A4"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1175.137400971587,906.6165142667815,843.6686446457672,915.3053765252265),c(1104.5528492898522,1275.143875741711,1046.1067239514725,964.4291042905453),c(1324.8911119720249,983.066339253602,969.7949872498455,988.4897056449872),c(1446.0294641825699,1029.132259437968,1076.843395678517,1098.7674618528458),c(1551.9062917051724,1029.132259437968,921.0402665793614,966.4341544034155))
targetgene="FLT1"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1048.2759770030636,881.1332392711746,830.9500218621627,732.8458162540422),c(1084.5220981369275,1191.833169025304,689.9852860105458,858.161448308427),c(640.984036893593,825.2660594731135,640.1706801080948,615.5503846511381),c(767.8454608621166,836.0474450481779,561.7391729425334,671.6917878115024),c(648.6147992375644,745.8758566021845,568.0984843343357,531.3382799105915))
targetgene="NLRP3"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1907.6905859928363,1206.5350584458465,1015.3700522244284,1154.9088650132103),c(720.1531962122957,713.5316998769913,644.4102210359629,755.903892552049),c(1184.6758539015514,1042.854022897141,910.441414259691,1020.5705074509096),c(1330.6141837300033,887.0139950393916,777.9557602638105,741.8685417619579),c(804.0915819959805,1085.9795651973989,933.758889362966,990.4947557578573))
targetgene="PILRA"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1123.6297551497805,894.8550027303476,889.24370962035,969.4417295727206),c(1345.875708417946,793.9020287092898,715.4225315777549,796.0048948094521),c(652.4301804095501,664.5254018085166,670.907351835139,648.6337115134957),c(1083.568252843931,880.1531133098051,768.416793176107,730.840766141172),c(953.8452929964182,716.4720777610997,674.0870075310402,544.3711056442475))
targetgene="VLDLR"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1172.275865092598,989.9272209831884,875.4652016047785,919.3154767509668),c(1159.8758762836444,1000.7086065582529,954.9565940023068,1203.030067722094),c(1487.044811781416,1115.3833440384835,913.6210699555921,982.4745553063767),c(1114.0913022198165,964.4439459875816,1054.585805807209,1025.583132733085),c(1550.952446412176,797.8225325547677,891.3634800842842,802.0200451480626))
targetgene="CLDN7"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(730.6454944352563,825.2660594731135,777.9557602638105,675.7018880372427),c(1421.2294865646631,646.8831345038657,565.9787138704016,584.4721079016506),c(1343.0141725389567,848.7890825459814,760.9975965523378,955.4063787826296),c(611.4148328107041,606.6979700877165,501.3257147204119,570.4367571115596),c(1067.352882862992,775.2796354432694,729.2010395933264,751.8937923263087))
targetgene="CD300A"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(933.8145418434934,1064.4167940472698,1091.6817889260556,1238.1184446973216),c(529.3841376130121,708.6310700701438,674.0870075310402,641.6160361184501),c(1011.0760105762033,932.0997892623883,907.2617585637898,873.1993241549532),c(882.3068960216868,837.0275710095474,911.5012994916581,941.3710279925385),c(802.1838914099877,742.9354787180761,762.0574817843047,768.9367182857051),c(785.9685214290486,960.5234421421036,843.6686446457672,814.0503458252836),c(912.8299453975721,803.7032883229847,790.674383047415,812.0452957124133),c(1079.7528716719453,840.9480748550254,840.4889889498661,814.0503458252836),c(907.1068736395937,711.5714479542522,797.0336944392172,787.9846943579715),c(604.7379157597292,678.2471652676895,781.1354159597116,900.2675006787003),c(519.8456846830479,912.4972700349985,922.1001518113285,766.9316681728349),c(859.4146089897728,874.2723575415882,798.0935796711843,755.903892552049),c(943.3529947734576,1002.6688584809918,960.256020162142,840.1159972925956),c(1164.6451027486266,897.795380614456,874.4053163728115,839.1134722361605),c(856.5530731107835,814.4846738980491,600.954926525314,705.777639730295),c(986.2760329582964,704.7105662246659,746.1592033047991,725.8281408589967),c(718.2455056263029,656.6843941175606,819.2912843105252,865.1791237034726),c(705.8455168173494,763.5181239068354,862.746578821174,787.9846943579715),c(373.90735485459595,977.185583485385,809.7523172228218,993.5023309271626),c(1202.7989144684834,1006.5893623264698,929.5193484350978,834.1008469539852),c(706.7993621103459,919.3581517645849,839.429103717899,952.3988036133244),c(738.2762567792277,896.8152546530865,993.1124623531205,893.2498252836548),c(731.5993397282527,915.4376479191069,896.6629062441194,860.1664984212972),c(687.7224562504175,666.4856537312555,632.7514834843254,548.3812058699879),c(679.1378486134497,902.6960104213035,926.3396927391966,935.355877653928),c(965.2914365123752,896.8152546530865,996.2921180490216,960.419004064805),c(800.2762008239948,765.4783758295744,706.9434497220185,657.6564370214113),c(847.9684654738157,724.3130854520557,691.0451712425129,695.7523891659444),c(768.7993061551131,641.0023787356488,718.602187273656,705.777639730295),c(964.3375912193787,898.7755065758255,883.9442834605148,742.871066818393),c(670.553240976482,803.7032883229847,854.2674969654377,738.8609665926527),c(1153.1989592326695,932.0997892623883,919.9803813473943,898.2624505658301),c(420.6457742114204,749.7963604476626,686.8056303146446,791.9947945837118),c(1306.7680514050928,1027.172007515229,982.5136100334499,1070.6967602726636),c(722.0608867982886,722.3528335293167,767.35690794414,805.0276203173679),c(848.9223107668122,1007.5694882878392,974.0345281777136,842.1210474054658),c(784.0608308430558,934.0600411851273,996.2921180490216,983.4770803628118),c(644.7994180655787,687.068298920015,715.4225315777549,877.2094243806935),c(1330.6141837300033,1086.9596911587682,1054.585805807209,800.0149950351924),c(1119.814373977795,819.3853037048965,839.429103717899,802.0200451480626),c(1069.2605734489848,795.8622806320287,716.4824168097218,723.8230907461265),c(1212.3373673984474,860.5505940824153,925.2798075072296,998.5149562093379),c(1122.6759098567843,871.3319796574797,870.1657754449433,819.0629711074589),c(1332.5218743159962,968.3644498330596,987.8130361932851,936.3584027103631),c(1042.5529052450852,977.185583485385,887.123939156416,832.095796841115),c(1140.7989704237161,1003.6489844423613,958.136249698208,904.2776009044406),c(868.9530619197369,920.3382777259544,803.3930058310195,867.1841738163428),c(1062.58365639801,877.2127354256967,950.7170530744387,878.2119494371286),c(1001.5375576462391,1000.7086065582529,903.0222176359216,932.3483024846228),c(884.2145866076796,734.1143450657506,603.0746969892482,726.8306659154317),c(893.7530395376439,756.6572421772489,947.5373973785375,612.5428094818328),c(429.23038184838816,844.8685787005033,868.0460049810092,710.7902650124705),c(497.9072429441303,794.8821546706592,668.7875813712049,838.1109471797255),c(1113.13745692682,868.3916017733712,830.9500218621627,888.2372000014793),c(605.6917610527255,763.5181239068354,691.0451712425129,687.7321887144637),c(285.19974260592903,762.5379979454659,695.284712170381,698.7599643352496),c(946.2145306524468,631.2011191219538,759.9377113203707,709.7877399560355),c(1194.2143068315156,1035.9931411675545,929.5193484350978,964.4291042905453),c(589.4763910717865,897.795380614456,861.6866935892069,868.1866988727778),c(842.2453937158373,797.8225325547677,845.7884151097013,871.1942740420831),c(565.6302587468759,993.8477248286664,850.0279560375694,883.224574719304),c(953.8452929964182,747.8361085249235,819.2912843105252,800.0149950351924),c(856.5530731107835,763.5181239068354,729.2010395933264,703.772589617425),c(567.5379493328688,955.6228123352562,785.3749568875797,1077.714435667709),c(1180.8604727295658,757.6373681386185,793.8540387433161,787.9846943579715),c(1012.9837011621961,969.344575794429,883.9442834605148,793.999844696582),c(1054.9528940540386,803.7032883229847,848.9680708056025,789.9897444708416),c(490.27648060015895,652.7638902720827,815.051743382657,719.8129905203862),c(862.276144868762,689.0285508427539,769.4766784080741,671.6917878115024),c(1068.3067281559884,923.2786556100629,742.979547608898,788.9872194144066),c(1086.4297887229204,1002.6688584809918,847.9081855736354,859.163973364862),c(1062.58365639801,906.6165142667815,839.429103717899,827.0831715589396),c(772.6146873270987,641.9825046970183,721.7818429695571,770.9417683985752),c(891.845348951651,754.69699025451,806.5726615269207,905.2801259608757),c(471.19957474023056,670.4061575767336,715.4225315777549,763.9240930035296),c(711.5685885753279,672.3664094994725,739.7998919129968,717.8079404075161),c(991.0452594232785,1001.6887325196224,987.8130361932851,980.4695051935065),c(693.445528008396,953.6625604125171,766.297022712173,882.2220496628689),c(970.0606629773573,791.9417767865507,800.2133501351184,704.77511467386),c(518.8918393900515,784.1007690955948,755.6981703925026,818.0604460510239),c(1003.445248232232,1084.0193132746597,1030.208445471967,957.4114288954997),c(1263.845013220254,839.967948893656,794.9139239752832,829.0882216718097),c(1016.7990823341818,724.3130854520557,736.6202362170957,748.8862171570034),c(615.2302139826897,717.4522037224692,701.6440235621833,772.9468185114454),c(788.8300573080378,769.3988796750524,718.602187273656,633.5958356669695),c(177.41522449733378,764.4982498682049,659.2486142835015,687.7321887144637),c(796.4608196520091,709.6111960315133,909.3815290277239,787.9846943579715),c(1025.3836899711496,735.0944710271201,788.5546125834809,629.5857354412292),c(915.6914812765615,1000.7086065582529,850.0279560375694,931.3457774281877),c(948.1222212384397,923.2786556100629,1056.7055762711432,1120.8230130944175),c(459.7534312242736,733.1342191043811,736.6202362170957,637.6059358927098),c(1034.9221429011138,857.6102161983068,850.0279560375694,837.1084221232903),c(612.3686781037005,983.066339253602,915.7408404195262,834.1008469539852),c(569.4456399188616,554.7512941351333,680.4463189228424,711.7927900689056),c(280.43051614094696,809.5840440912016,780.0755307277445,997.5124311529029),c(1131.2605174937519,1032.0726373220766,857.4471526613388,949.3912284440191),c(321.4458637397929,555.7314200965028,666.6678109072708,665.676637472892),c(545.5995075939512,911.517144073629,929.5193484350978,717.8079404075161),c(1216.1527485704332,877.2127354256967,878.6448573006796,931.3457774281877),c(443.5380612433345,950.7221825284087,906.2018733318228,826.0806465025045))
targetgene="Non"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1106.4605398758451,993.8477248286664,837.309333253965,890.2422501143495),c(1062.58365639801,1002.6688584809918,922.1001518113285,868.1866988727778),c(1353.5064707619174,776.2597614046389,895.6030210121523,934.353352597493),c(1213.291212691444,770.3790056364219,716.4824168097218,865.1791237034726),c(562.7687228678867,1083.0391873132903,829.8901366301956,735.8533914233475))
targetgene="TYRO3"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("FCGR2C","TGFBR3","P2RY13","CD86","CLEC2B","C1QC","IL10","PTPRC","CLEC2D","GRN")
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
targetmat=list(c(478.8303370842019,741.9553527567066,709.0632201859526,792.9973196401469),c(409.1996306954634,661.5850239244081,734.5004657531616,744.8761169312631),c(361.50736604564247,747.8361085249235,784.3150716556127,665.676637472892),c(451.1688235873058,774.2995094818998,789.614497815448,767.9341932292699),c(557.0456511099082,816.4449258207881,869.1058902129762,807.032670430238))
targetgene="FCGR2C"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1281.0142284941896,912.4972700349985,949.6571678424716,886.2321498886092),c(401.56886835149203,602.7774662422385,870.1657754449433,879.2144744935637),c(817.4454160979304,731.1739671816422,730.2609248252935,702.7700645609899),c(774.5223779130915,817.4250517821575,735.5603509851287,855.1538731391217),c(251.8151573510544,592.9762066285435,854.2674969654377,767.9341932292699))
targetgene="TGFBR3"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(477.8764917912055,771.3591315977914,903.0222176359216,985.482130475682),c(529.3841376130121,847.8089565846119,1023.8491340801647,1125.835638376593),c(581.8456287278151,782.1405171728558,808.6924319908547,795.0023697530171),c(560.8610322818939,947.7818046443002,865.9262345170752,1069.6942352162284),c(568.4917946258653,1066.377045970009,1065.1846581268794,1098.7674618528458))
targetgene="P2RY13"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(658.1532521675285,671.386283538103,658.1887290515344,788.9872194144066),c(433.9996083133703,634.1414970060623,702.7039087941503,828.0856966153747),c(567.5379493328688,724.3130854520557,782.1953011916786,826.0806465025045),c(475.9688012052127,717.4522037224692,798.0935796711843,792.9973196401469),c(342.43046018571414,716.4720777610997,824.5907104703604,822.0705462767642))
targetgene="CD86"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(397.7534871795064,1046.774526742619,923.1600370432955,1034.6058582410008),c(309.99972022383594,822.3256815890051,837.309333253965,781.969544019361),c(710.6147432823316,637.0818748901708,913.6210699555921,886.2321498886092),c(1030.1529164361316,1078.1385575064428,895.6030210121523,1200.0224925527887),c(562.7687228678867,809.5840440912016,767.35690794414,743.8735918748281))
targetgene="CLEC2B"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(307.13818434484665,837.0275710095474,783.2551864236457,774.9518686243155),c(557.9994964029046,851.7294604300898,854.2674969654377,889.2397250579145),c(651.4763351165536,599.83708835813,700.5841383302162,820.0654961638941),c(624.7686669126539,676.2869133449504,696.3445974023481,715.8028902946459),c(495.045707065141,947.7818046443002,1125.598116349001,903.2750758480055))
targetgene="C1QC"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(270.89206321098277,917.3978998418459,760.9975965523378,688.7347137708988),c(985.3221876653,974.2452056012766,986.7531509613182,984.4796054192469),c(477.8764917912055,1010.5098661719478,1052.4660353432748,1047.6386839746567),c(883.2607413146833,965.4240719489511,930.5792336670648,1094.7573616271054),c(454.03035946629507,750.776486409032,818.2313990785581,887.2346749450443))
targetgene="IL10"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(464.52265768925565,725.2932114134252,797.0336944392172,762.9215679470946),c(523.6610658550336,872.3121056188492,817.1715138465911,932.3483024846228),c(421.59961950441686,830.166689279961,993.1124623531205,841.1185223490306),c(933.8145418434934,841.9282008163949,817.1715138465911,742.871066818393),c(1162.7374121626337,986.9868430990799,890.3035948523171,1027.5881828459553))
targetgene="PTPRC"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(391.07657012853144,906.6165142667815,880.7646277646137,873.1993241549532),c(335.7535431347392,545.9301604828079,647.589876731864,712.7953151253406),c(717.2916603333065,771.3591315977914,894.5431357801854,857.158923251992),c(509.35338646008734,657.6645200789301,721.7818429695571,654.6488618521062),c(686.7686109574211,894.8550027303476,775.8359897998764,867.1841738163428))
targetgene="CLEC2D"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(485.5072541351769,722.3528335293167,711.1829906498866,701.7675395045547),c(381.53811719856725,609.6383479718249,588.2363037417095,687.7321887144637),c(433.04576302037384,860.5505940824153,857.4471526613388,886.2321498886092),c(436.8611441923595,849.7692085073509,777.9557602638105,752.8963173827437),c(486.46109942817327,756.6572421772489,557.4996320146653,674.6993629808077))
targetgene="GRN"
collabel=c("TRE-VPH-7d-dox-undiff-rd1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-undiff-rd2_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-diff-rd1_S4_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("Diff_CRISPRa_summary.Rnw");
library(tools);

texi2dvi("Diff_CRISPRa_summary.tex",pdf=TRUE);

