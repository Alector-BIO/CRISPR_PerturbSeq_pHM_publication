pdf(file='Bodipy_CRISPRi.pdf',width=4.5,height=4.5);
gstable=read.table('Bodipy_CRISPRi.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("CD300A","CD69","IL1R1","ILDR1","SPI1","GBA","SELPLG","VLDLR","UNC13D","PLEK")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='0,2_vs_1,3 neg.'


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
targetmat=list(c(702.5899150001461,796.790543277285,974.5365918917008,883.0892629125481),c(854.1501431974041,944.6068321190686,1027.9106450831255,638.3896682730934),c(1030.6095517413544,798.8154239463506,776.6844981648683,551.5607798526418),c(1633.602745354731,1436.6528347019917,1256.1307345912858,1203.764134919898),c(810.8472208553303,734.0192425362537,756.4391676439831,678.8440367417129))
targetgene="CD300A"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(640.883250662691,705.6709131693364,673.6173609676346,590.0417644935237),c(584.5894516179952,704.6584728348035,688.3412377100965,784.4200715256711),c(416.7906275424595,450.53594886707987,431.5936370134161,518.0132547811036),c(838.9941203776783,667.1981804570913,805.2120093533883,1016.2926712848318),c(2667.460016271741,876.7733297053734,758.2796522367908,1268.8858012352366))
targetgene="CD69"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(754.5534218106345,598.3522377088634,747.2367446799443,535.7737092307415),c(572.6811479739249,654.0364561081653,411.3483064925309,455.8516642073712),c(848.7372779046449,835.26327598953,667.1756648928074,651.2166631533873),c(1026.2792595071471,936.5073094428064,866.8682432124477,771.5930766453772),c(1291.5096588523486,1183.5427510688007,1050.9167024932221,1246.1918872162548))
targetgene="IL1R1"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(953.7468645841736,807.9273869571455,610.1206425157674,779.4866119563274),c(978.646044930866,914.2336220830856,1030.671371972337,988.6652976965063),c(1042.5178553854248,741.1063248779831,779.4452250540799,705.4847184161697),c(772.9571638060158,764.3924525722366,464.7223596839555,774.5531523869835),c(729.6542414639422,885.8852927161682,813.4941900210232,953.1443887972306))
targetgene="ILDR1"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(768.6268715718085,538.6182579714304,476.68550953720586,601.882067459949),c(558.607698212751,548.7426613167579,544.7834394710924,430.1976744467832),c(600.8280474962728,573.0412293455443,285.27511188520043,359.1558566482318),c(751.305702634979,545.7053403131597,738.0343217159057,762.7128494205582),c(546.6993945686806,454.5857102052109,593.5562811804976,667.9904256891565))
targetgene="SPI1"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(618.1492164331023,769.4546542449004,679.1388147460577,620.6292138234555),c(1011.1232366874214,852.474761676587,696.6234183777314,875.1957276015979),c(643.0483967797948,674.2852627988207,533.7405319142459,413.4239119110141),c(473.08442658715535,560.8919453311511,441.71630227385873,435.131134016127),c(897.4530655394778,798.8154239463506,772.0832866828489,821.9143642526843))
targetgene="GBA"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(744.810264283668,856.524523014718,703.0651144525584,520.97333052271),c(1231.9681406319974,1081.2862772809917,1273.6153382229593,1619.1614306586496),c(1346.7208848384926,1063.0623512594018,937.7269000355459,960.051232194312),c(844.4069856704375,804.8900659535472,783.1261942396953,719.2984052103325),c(869.3061660171298,937.5197497773393,968.0948958168736,815.007520855603))
targetgene="SELPLG"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1019.783821155836,1079.261396611926,1132.8182668731667,1076.4808780308267),c(1067.417035732117,1005.3532521910344,967.1746535204699,887.0360305680231),c(614.9014972574467,655.0488964426982,842.941943505947,586.0949968380486),c(867.1410199000262,824.1264323096697,623.9242769618254,426.2509067913081),c(729.6542414639422,921.320704424815,899.0767235865833,849.5417378410099))
targetgene="VLDLR"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(659.2869926580723,568.9914680074132,451.83896753430133,620.6292138234555),c(991.6369216334881,1049.900626910476,946.9293229995847,1001.4922925768003),c(1207.068960285305,876.7733297053734,746.3165023835405,950.1843130556243),c(497.98360693384774,669.2230611261568,522.6976243573995,551.5607798526418),c(899.6182116565815,828.1761936478007,914.7208426254491,758.7660817650832))
targetgene="UNC13D"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1051.1784398538396,942.581951450003,643.2493651863068,1148.5093877432469),c(936.4256956473441,802.8651852844816,746.3165023835405,778.4999200424586),c(1798.1538502546111,934.4824287737409,815.3346746138309,1333.0207756367065),c(1400.8495377660847,1187.5925124069317,1314.1059992647297,1655.669031471794),c(940.7559878815515,977.004922824117,901.8374504757949,868.2888842045165))
targetgene="PLEK"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("CLEC7A","CD226","ITGAX","MRC1","CLDN7","ADAM10","C9orf72","SIRPA","SORL1","CD200R")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='0,2_vs_1,3 pos.'


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
targetmat=list(c(350.7536709707971,598.3522377088634,774.8440135720605,601.882067459949),c(403.79975083983743,603.4144393815271,825.4573398742735,583.1349210964423),c(929.9302572960331,788.691020601023,568.7097391775931,849.5417378410099),c(564.0205635055102,665.1732997880258,514.4154436897646,612.7356785125054),c(772.9571638060158,704.6584728348035,797.8500709821573,671.9371933446315))
targetgene="CLEC7A"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(291.21215275044574,574.0536696800771,1297.54163792946,527.8801739197913),c(71.44982186442164,328.0306683886154,138.95658675698473,138.13686794162763),c(562.9379904469583,718.8326375182622,806.1322516497921,709.4314860716447),c(854.1501431974041,905.1216590722908,972.696107298893,990.6386815242438),c(1388.9412341220145,759.3302508995728,874.2301815836787,1641.8553446776311))
targetgene="CD226"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1129.1237000695721,1226.0652451191768,1037.1130680471642,1032.0797419067321),c(1142.1145767721944,854.4996423456525,681.8995416352694,1220.537897455667),c(3074.507486287234,1755.5715400798124,2560.114068595573,8259.598010995463),c(894.2053463638223,764.3924525722366,954.2912613708156,874.2090356877292),c(828.1683897921598,818.0517903024731,774.8440135720605,751.8592383680018))
targetgene="ITGAX"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(690.6816113560758,771.479534913966,911.9601157362375,847.5683540132724),c(1137.7842845379869,771.479534913966,1032.5118565651449,1823.4066568294847),c(762.1314332204975,772.4919752484988,766.5618329044256,746.9257787986579),c(908.2787961249962,716.8077568491967,787.7274057217147,737.0588596599703),c(594.3326091449618,655.0488964426982,756.4391676439831,500.2528003314658))
targetgene="MRC1"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(680.9384538291092,969.9178404823876,982.8187725593357,764.6862332482958),c(908.2787961249962,890.9474943888321,741.7152909015211,919.5968637256925),c(632.2226661942763,742.1187652125158,857.665820248409,1112.988478843971),c(882.297042719752,939.5446304464048,912.8803580326413,1059.7071154950577),c(793.5260519185009,798.8154239463506,964.4139266312582,796.2603744920964))
targetgene="CLDN7"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(921.2696728276184,701.6211518312052,767.4820752008295,997.5455249213252),c(533.7085178660585,752.2431685578434,849.3836395807741,964.9846917636559),c(1013.288382804525,1046.8633059068777,787.7274057217147,1291.5797152542184),c(767.5442985132566,915.2460624176184,1213.7995889567076,597.9352998044739),c(736.1496798152532,748.1934072197124,819.0156437994464,536.7604011446102))
targetgene="ADAM10"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1274.188489915519,903.0967784032252,1147.542143615629,1934.9028430966555),c(608.4060589061357,655.0488964426982,729.7521410482708,777.5132281285897),c(488.24044940688117,608.476641054191,530.0595627286305,528.8668658336601),c(584.5894516179952,638.849851090174,426.99242553139675,530.8402496613976),c(806.516928621123,877.7857700399062,1234.0449194775929,1179.0968370731787))
targetgene="C9orf72"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1161.6008918261275,906.1340994068236,1205.5174082890728,1133.7090090352153),c(1296.922524145108,1149.1197796946867,803.3715247605805,1254.085422527205),c(705.8376341758016,800.8403046154161,834.6597628383122,875.1957276015979),c(834.6638281434709,965.8680791442566,1015.947495229875,849.5417378410099),c(687.4338921804202,742.1187652125158,1038.033310343568,985.7052219549))
targetgene="SIRPA"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(835.7464012020228,929.4202271010771,760.1201368295986,1112.988478843971),c(1326.1519967260076,830.2010743168662,761.9606214224063,1026.1595904235196),c(1064.1693165564616,1114.696808320573,963.4936843348544,1001.4922925768003),c(1518.8500011482356,1129.8834133385644,1107.0514825738585,1947.7298379769495),c(581.3417324423397,749.2058475542451,1072.0822753105112,791.3269149227525))
targetgene="SORL1"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(748.0579834593235,753.2556088923762,641.408880593499,887.0360305680231),c(1053.3435859709432,694.5340694894759,761.0403791260024,715.3516375548573),c(492.57074164108855,622.6508057376496,977.2973187809124,562.4143909051982),c(822.7555244994006,753.2556088923762,650.6113035575378,624.5759814789307),c(510.97448363646987,505.2077269318491,618.4028231834022,452.89158846576487))
targetgene="CD200R"
collabel=c("TRE-PMA-no-serum-Unsorted-KRAB-rep2_S2_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-KRAB-rep1_S1_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep1_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-KRAB-rep2_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("Bodipy_CRISPRi_summary.Rnw");
library(tools);

texi2dvi("Bodipy_CRISPRi_summary.tex",pdf=TRUE);

