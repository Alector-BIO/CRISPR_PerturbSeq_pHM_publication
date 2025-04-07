pdf(file='Bodipy_CRISPRa.pdf',width=4.5,height=4.5);
gstable=read.table('Bodipy_CRISPRa.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("S100A4","LAPTM5","ADAM10","LPAR2","CR1","MERTK","NFATC1","ABCA1","TLR4","HEXB")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='2,3_vs_0,1 neg.'


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
targetmat=list(c(619.5864615383676,512.8654438880907,447.96597165437396,710.3039639397579),c(284.4160418499153,201.54712180865317,160.22140592017755,280.94112006572516),c(794.8327095469014,770.1979297687817,865.413580276333,651.9954295864942),c(574.5778623230613,634.33357533527,608.1873775745515,642.4540330559602),c(7927.259240517168,862.8736152432963,921.0005986568028,875.688170469015))
targetgene="S100A4"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(708.6460301984422,761.200290402324,786.9377896215523,753.7703259121909),c(658.8492821304435,822.3842380942366,451.2357962649898,641.3938778859008),c(800.5784881701319,725.2097329364931,849.0644572232537,711.3641191098172),c(714.3918088216727,476.8748864222597,815.2762695802231,664.7172916272062),c(1220.977957436505,1034.728527142639,742.2501866098021,894.770963530083))
targetgene="LAPTM5"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(617.6712019972908,648.7297983216024,478.4843346867887,486.6112230572371),c(587.0270493400609,585.7463227563983,759.6892511997534,698.6422570691051),c(543.9337096658313,764.799346148907,693.202817450564,521.5963436691953),c(2256.1757393885537,529.0611947477146,683.3933436187165,850.2444463875908),c(578.408381405215,540.7581259241095,555.8701838046976,653.0555847565536))
targetgene="ADAM10"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(888.680427059668,690.118939407308,682.3034020818445,960.5005840737622),c(1016.0451865412799,911.4608678221681,926.4503063411627,866.1467739384809),c(849.4176064675921,1052.7238058755545,875.2230541081807,918.0943772713886),c(1514.9702969918046,748.6035952892831,783.6679650109364,690.1610157086304),c(784.2987820709785,900.6637005824188,560.2299499521854,685.920395028393))
targetgene="LPAR2"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(767.0614462012867,939.3535498581871,617.9968514063991,797.2366878846237),c(1033.2825224109718,1074.318140355053,1014.7355708277911,1034.711445977916),c(622.459350849983,625.3359359688122,440.33638089627027,296.8434476166152),c(900.1719843061293,926.7568547451463,833.8052757070463,720.9055156403513),c(597.5609768159836,534.4597783675891,546.06070997285,706.0633432595205))
targetgene="CR1"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(582.2389004873687,895.2651169625442,709.5519405036434,639.273567545782),c(610.0101638329833,797.190847868155,535.1612946041304,373.17461986088773),c(1253.537369634812,894.3653530258983,729.1708881673386,888.410032509727),c(440.50969444768026,385.99872882103665,464.3150947074533,359.3926026501163),c(609.0525340624449,667.6248409911636,543.8808268991061,547.0400677506195))
targetgene="MERTK"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(904.9601331588213,736.0069001762423,571.129365320905,705.0031880894612),c(946.1382132919741,767.4986379588444,626.7163837013748,899.0115842103204),c(905.9177629293598,1050.924278002263,622.3566175538868,572.4837918320437),c(1007.426518606434,793.5917921215719,924.2704232674188,670.018067477503),c(834.0955301389772,557.8536407203793,781.4880819371925,649.8751192463754))
targetgene="NFATC1"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(724.9257362975956,663.1260213079347,789.1176726952962,541.7392919003228),c(719.1799576743649,653.2286180048312,637.6157990700943,733.6273776810634),c(469.23858756383333,628.9349917153953,473.0346270024289,358.33244748005694),c(719.1799576743649,646.0305065116651,534.0713530672584,500.3932402680085),c(787.1716713825939,646.9302704483108,694.2927589874361,596.8673607434085))
targetgene="ABCA1"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(655.0187630482898,411.1921190471183,793.477438842784,645.6344985661382),c(490.3064425156789,511.9656799514449,274.6652672917329,338.1894992489295),c(593.7304577338299,490.3713454719463,542.7908853622341,710.3039639397579),c(681.832396623366,599.2427818060848,937.3497217098823,830.1014981564633),c(552.5523776006772,633.4338113986242,700.8324082086677,516.2955678188987))
targetgene="TLR4"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(905.9177629293598,822.3842380942366,910.1011832880833,566.1228608116876),c(654.0611332777513,457.97984375269846,640.8856236807102,631.8524813553668),c(770.8919652834404,537.1590701775265,668.1341621025091,595.8072055733492),c(779.5106332182863,811.5870708544873,689.9329928399482,602.1681365937052),c(1010.2994079180493,1030.22970745941,674.6738113237408,578.8447228523997))
targetgene="HEXB"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("CLEC2D","C3AR1","LPL","C1QC","GRN","MS4A4A","CD52","LRP1B","Non","FPR2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='2,3_vs_0,1 pos.'


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
targetmat=list(c(302.6110074901456,671.2238967377467,437.06655628565437,504.6338609482459),c(203.97514112468673,188.05066275896655,303.00374725040376,341.3699647591075),c(432.84865628337275,547.05647348063,584.2086637633685,482.3706023769998),c(566.9168241587537,513.7652078247364,576.5790730052647,579.904878022459),c(483.60303412190984,755.8017067824494,698.6525251349238,718.7852053002326))
targetgene="CLEC2D"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(798.6632286290551,470.5765388657393,869.773346423821,1170.4113077455115),c(524.7811142550626,745.0045395427001,619.086792943271,716.6648949601139),c(685.6629157055197,636.1331032085616,838.1650418545343,621.2509296547734),c(412.73843110206565,739.6059559228254,509.00269771920347,693.3414812188084),c(347.61960670545204,638.8323950184988,526.4417623091548,614.8899986344172))
targetgene="C3AR1"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(311.2296754249915,292.4232794098762,328.0724025984588,488.7315333973558),c(518.0777058612936,682.9208279141417,400.0085440320079,705.0031880894612),c(452.0012516941415,647.8300343849565,857.7839895182294,905.3725152306764),c(299.73811817853027,508.36662420486175,546.06070997285,502.5135506081272),c(814.9429347282085,1011.334664789849,814.1863280433512,1094.080135501239))
targetgene="LPL"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(441.4673242182187,672.1236606743925,565.6796576365452,562.9423953015096),c(503.713259303217,461.5788994992816,403.2783686426237,597.9275159134678),c(588.9423088811377,561.4526964669624,597.287962205832,783.4546706738523),c(378.263759362682,505.6673323949245,769.498725031601,640.3337227158414),c(482.64540435137144,615.4385326657087,478.4843346867887,747.4093948918348))
targetgene="C1QC"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(502.75562953267854,609.1401851091883,561.3198914890573,491.9119989075338),c(366.77220211622074,501.1685127116956,492.65357466612414,496.1526195877712),c(307.3991563428377,489.4715815353005,368.4002394627211,426.18237836385475),c(429.9757669717575,505.6673323949245,372.76000561020896,530.0775850296701),c(295.90759909637654,308.6190302695002,466.4949777811972,484.4909127171184))
targetgene="GRN"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(693.3239538698272,637.0328671452073,826.1756849489427,1254.1635661801993),c(799.6208583995935,634.33357533527,734.6205958516983,633.9727916954854),c(607.137274521368,516.4644996346738,729.1708881673386,857.6655325780062),c(234.61929378191667,355.40675497508033,337.88187643030636,307.4449993172086),c(586.0694195695224,762.9998182756156,706.2821158930275,665.7774467972656))
targetgene="MS4A4A"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(378.263759362682,545.2569456073385,433.79673167503853,430.4229990440921),c(773.7648545950558,821.4844741575907,907.9213002143393,900.0717393803798),c(609.0525340624449,967.2462318942061,856.6940479813575,506.7541712883646),c(734.5020340029799,900.6637005824188,653.9649221231737,621.2509296547734),c(430.9333967422959,400.394951807369,617.9968514063991,673.198532987681))
targetgene="CD52"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(666.510320294751,534.4597783675891,1056.1533492289254,788.7554465241491),c(601.3914958981374,519.1637914446111,467.5849193180691,461.167498975813),c(759.4004080369792,516.4644996346738,821.8159188014548,499.3330850979492),c(880.0617591248221,578.548211263232,1166.237444452993,800.4171533948017),c(604.2643852097527,586.646086693044,704.1022328192836,552.3408436009162))
targetgene="LRP1B"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(949.9687323741279,708.1142181402233,783.6679650109364,783.4546706738523),c(401.2468738556044,696.4172869638284,446.87603011750195,442.08470591474486),c(775.6801141361326,663.1260213079347,812.0064449696073,572.4837918320437),c(791.0021904647476,709.9137460135149,779.3081988634485,819.4999464558699),c(842.7141980738231,666.7250770545178,858.8739310551014,522.6564988392547),c(584.1541600284455,411.1921190471183,392.3789532739042,503.57370577818654),c(552.5523776006772,502.06827664834134,561.3198914890573,651.9954295864942),c(477.85725549867925,220.44216447821438,354.23099948338574,558.7017746212722),c(778.553003447748,760.3005264656782,828.3555680226866,592.6267400631712),c(626.2898699321368,795.3913199948634,728.0809466304667,914.9139117612106),c(586.0694195695224,845.7781004470266,791.2975557690401,851.3046015576501),c(699.0697324930578,613.6390047924172,501.3731069610998,638.2134123757228),c(765.1461866602099,847.5776283203182,653.9649221231737,753.7703259121909),c(703.85788134575,696.4172869638284,706.2821158930275,817.3796361157513),c(561.1710455355231,573.1496276433575,539.5210607516183,826.9210326462853),c(700.0273622635963,897.9644087724815,547.1506515097219,1178.8925491059863),c(687.5781752465966,807.0882511712584,794.567380379656,755.8906362523095),c(548.7218585185235,837.6802250172146,635.4359159963503,583.085343532637),c(613.8406829151371,902.4632284557104,686.6631682293323,756.9507914223689),c(587.9846791105992,771.0976937054274,537.3411776778744,860.8459980881842),c(516.1624463202166,753.1024149725121,603.8276114270636,827.9811878163447),c(587.0270493400609,824.1837659675281,697.5625835980519,738.9281535313601),c(775.6801141361326,800.789903614738,965.688201668553,778.1538948235557),c(450.08599215306464,471.4763028023851,458.8653870230935,417.70113700338),c(513.2895570086014,446.28291257630343,707.3720574298994,656.2360502667316),c(599.4762363570605,528.1614308110687,524.2618792354109,701.8227225792832),c(642.5695760312901,579.4479751998779,617.9968514063991,445.26517142492287),c(684.7052859349812,546.1567095439842,658.3246882706615,831.1616533265227),c(776.6377439066711,565.0517522135455,674.6738113237408,543.8596022404415),c(678.0018775412123,799.8901396780923,495.92339927674,774.9734293133777),c(347.61960670545204,713.512801760098,451.2357962649898,541.7392919003228),c(1111.8081635951235,1087.8145994047395,1046.343875397078,1037.891911488094),c(494.1369615978326,591.1449063762728,476.3044516130448,642.4540330559602),c(1001.6807399832034,983.44198275383,785.8478480846803,1081.358273460527),c(666.510320294751,732.4078444296592,652.8749805863017,548.1002229206788),c(764.1885568896714,750.4031231625747,828.3555680226866,641.3938778859008),c(567.8744539292921,638.8323950184988,591.8382545214721,648.8149640763162),c(703.85788134575,688.3194115340164,673.5838697868688,755.8906362523095),c(1354.0884955413476,831.3818774606942,1085.5817707244682,837.5225843468787),c(569.7897134703691,746.8040674159917,609.2773191114235,874.6280152989556),c(682.7900263939044,751.3028870992205,367.31029792584917,539.6189815602041),c(756.527518725364,552.4550571005046,769.498725031601,600.0478262535865),c(1005.5112590653571,1191.2874521190035,1057.2432907657974,1066.5161010796962),c(906.8753926998983,1126.504448680508,1253.4327674027495,732.567222511004),c(926.027988110667,978.0433991339553,877.4029371819246,823.7405671361073),c(1176.926987991737,1093.2131830246142,1070.3225892082608,1070.7567217599335),c(596.6033470454453,728.8087886830762,577.6690145421367,557.6416194512129),c(929.8585071928208,686.5198836607249,1012.5556877540472,712.4242742798766),c(876.2312400426683,1015.8334844730778,606.0074945008075,892.6506531899644),c(2016.768296753945,871.871254609754,862.1437556657172,657.2962054367908),c(500.8403699916017,471.4763028023851,374.93988868395286,626.55170550507),c(572.6626027819843,533.5600144309434,730.2608297042106,558.7017746212722),c(543.9337096658313,796.2910839315091,670.314045176253,711.3641191098172),c(907.8330224704367,913.2603956954596,912.2810663618272,786.6351361840303),c(615.7559424562139,597.4432539327933,626.7163837013748,609.5892227841206),c(319.84834335983743,485.87252578871744,536.2512361410023,702.8828777493425),c(720.1375874449034,377.9008533912247,515.5423469404352,544.9197574105008),c(801.5361179406704,831.3818774606942,880.6727617925405,794.0562223744457),c(358.1535341813748,389.5977845676197,728.0809466304667,594.7470504032898),c(743.1207019378259,773.7969855153648,709.5519405036434,750.5898604020128),c(535.3150417309854,341.91029592539377,651.7850390494298,585.2056538727558),c(820.6887133514391,704.5151623936403,608.1873775745515,664.7172916272062),c(752.6969996432102,788.1932085016972,628.8962667751186,601.1079814236458),c(458.70466008791055,812.486834791133,707.3720574298994,617.010308974536),c(836.9684194505925,825.0835299041738,753.1496019785217,677.4391536679184),c(706.7307706573653,658.6272016247059,802.1969711377596,643.5141882260194),c(654.0611332777513,654.128381941477,731.3507712410825,661.5368261170282),c(494.1369615978326,741.405483796117,504.6429315717157,573.543947002103),c(635.866167637521,791.7922642482803,789.1176726952962,453.7464127853976),c(1281.3086329804264,787.2934445650515,1162.9676198423772,705.0031880894612),c(732.586774461903,678.4220082309129,744.430069683546,840.7030498570567),c(807.281896563901,682.9208279141417,736.8004789254422,592.6267400631712),c(558.2981562239078,605.5411293626053,403.2783686426237,521.5963436691953),c(566.9168241587537,638.8323950184988,576.5790730052647,577.7845676823404),c(472.1114768754486,583.047030946461,480.6642177605326,720.9055156403513),c(646.4000951134439,858.3747955600675,518.812171551051,576.724412512281),c(722.0528469859803,462.47866343592733,709.5519405036434,698.6422570691051),c(718.2223279038265,884.4679497227949,644.155448291326,554.4611539410349),c(716.3070683627496,645.1307425750192,731.3507712410825,956.2599633935248),c(554.467637141754,607.3406572358967,684.4832851555884,519.4760333290767),c(735.4596637735183,1060.8216813053664,667.0442205656371,735.7476880211821),c(733.5444042324415,819.6849462842993,722.6312389461069,709.2438087696985),c(753.6546294137486,608.2404211725426,473.0346270024289,521.5963436691953),c(570.7473432409074,566.851280086837,504.6429315717157,494.0323092476525),c(574.5778623230613,509.26638814150755,505.7328731085876,869.3272394486589),c(259.51766781591596,593.8441981862102,493.7435162029961,423.0019128536767),c(611.9254233740602,667.6248409911636,671.4039867131249,578.8447228523997),c(906.8753926998983,837.6802250172146,792.387497305912,778.1538948235557),c(927.9432476517438,685.620119724079,1082.3119461138524,957.3201185635842),c(763.230927119133,740.5057198594712,859.9638725919733,926.5756186318632),c(400.28924408506595,610.0399490458341,542.7908853622341,816.3194809456919),c(832.1802705979003,798.0906118048007,564.5897160996732,662.5969812870876),c(292.07708001422276,285.22516791671006,326.98246106158683,559.7619297913316),c(339.95856854114453,723.4102050632015,617.9968514063991,435.7237748943888),c(348.57723647599045,697.3170509004741,491.5636331292522,626.55170550507),c(639.6966867196749,569.5505718967744,753.1496019785217,692.2813260487491),c(245.15322125783945,515.5647356980279,498.10328235048394,441.0245507446855),c(649.2729844250591,986.1412745637673,1135.7190814205783,689.100860538571),c(770.8919652834404,835.8806971439232,665.9542790287652,710.3039639397579),c(601.3914958981374,1082.416015784865,788.0277311584242,947.77872203305))
targetgene="Non"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(949.9687323741279,835.8806971439232,858.8739310551014,741.0484638714787),c(869.5278316488993,723.4102050632015,728.0809466304667,805.7179292450985),c(792.9174500058244,546.1567095439842,1125.9096075887305,906.4326704007358),c(640.6543164902132,729.7085526197219,830.5354510964305,761.1914121026063),c(478.81488526921765,523.6626111278399,577.6690145421367,556.5814642811536))
targetgene="FPR2"
collabel=c("TRE-PMA-no-serum-Unsorted-VPH-rep2_S4_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-Unsorted-VPH-rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts","TRE-PMA-no-serum-BODIPY-high-VPH-rep1_S7_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("Bodipy_CRISPRa_summary.Rnw");
library(tools);

texi2dvi("Bodipy_CRISPRa_summary.tex",pdf=TRUE);

