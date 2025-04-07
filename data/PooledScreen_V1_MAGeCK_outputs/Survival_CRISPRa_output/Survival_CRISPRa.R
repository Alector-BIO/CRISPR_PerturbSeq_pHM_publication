pdf(file='Survival_CRISPRa.pdf',width=4.5,height=4.5);
gstable=read.table('Survival_CRISPRa.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("CD86","TGFBR3","LPL","GRN","P2RY13","SIRPA","GBA","IFNGR2","CLEC2D","LRP1B")
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
targetmat=list(c(935.8321028120745,884.7089149597,730.6958568485135,746.9857759976118),c(801.8900315998002,828.0875444022793,780.1646877776088,642.4293253889995),c(839.781538587483,923.1305592665212,784.2870903550335,695.2465014696387),c(766.6421181228859,720.9113787043042,806.9603045308689,514.1590406217327),c(776.3352943290373,808.8767222488686,788.4094929324581,617.637589677679))
targetgene="CD86"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(908.514969867466,962.5632994761536,912.0815702551967,1031.5517859014637),c(779.8600856767288,715.8558991902487,578.1669614838027,487.2115018050801),c(964.030433593606,745.1776803717702,758.5220742461296,780.4007241302611),c(905.8713763566974,847.2983665556899,873.9493464140189,771.7775117089323),c(763.1173267751944,660.245624535639,584.3505653499395,404.2130822497899))
targetgene="TGFBR3"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(726.1070176244344,744.1665844689591,780.1646877776088,482.8998955944157),c(927.0201244428459,724.9557623155484,586.4117666386519,627.338703651674),c(792.1968553936488,783.5993246785914,826.5417167736358,462.4197660937597),c(787.7908662090344,648.112473701906,585.3811659942957,526.01595770106),c(1030.1202713628202,969.6409707958312,1174.884734566016,1089.7584697454336))
targetgene="LPL"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(775.4540964921144,733.0445295380372,741.001863292075,432.23852261910866),c(777.2164921659602,748.2109680802034,654.431409166158,514.1590406217327),c(682.9283236152145,676.4231589806163,651.3396072330896,639.1956207310011),c(726.9882154613573,823.0320648882238,807.9909051752251,561.5867089390415),c(769.2857116336545,693.6117893284048,658.5538117435826,540.0286778857193))
targetgene="GRN"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(839.781538587483,831.1208321107125,783.2564897106773,701.7139107856353),c(854.7619018151717,866.5091887091005,816.2357103300743,679.0779781796471),c(856.5242974890174,776.5216533589139,773.981083911472,657.519947126325),c(901.4653871720831,934.2526141974432,916.2039728326213,680.1558797323131),c(820.3951861751802,948.4079568367985,1003.8050276028944,817.0493769209087))
targetgene="P2RY13"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(842.4251320982517,776.5216533589139,741.001863292075,729.7393511549541),c(843.3063299351745,944.3635732255541,904.8673657447035,907.5931073448617),c(920.851739584386,989.862888852053,966.7034044060729,1060.6551278234488),c(977.2484011474488,928.1860387805767,978.0400114939905,769.6217086036),c(839.781538587483,986.8296011436197,705.9614413839657,541.1065794383854))
targetgene="SIRPA"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(720.8198306028974,709.7893237733822,846.1231290164028,643.5072269416655),c(852.1183083044031,952.4523404480427,738.9406620033627,690.9348952589743),c(747.2557657105831,541.947403906742,577.1363608394464,430.08271951377645),c(815.107999153643,867.5202846119116,852.3067328825397,704.9476154436336),c(929.6637179536145,904.9308330159217,897.6531612342105,658.597848678991))
targetgene="GBA"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(793.0780532305716,752.2553516914478,704.9308407396096,545.4181856490499),c(808.0584164582602,837.187407527579,804.8991032421566,619.7933927830112),c(841.5439342613288,808.8767222488686,925.4793786318266,613.3259834670146),c(845.9499234459431,871.5646682231559,878.0717489914435,606.8585741510179),c(917.3269482366945,884.7089149597,782.2258890663212,837.5295064215647))
targetgene="IFNGR2"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(819.5139883382574,930.2082305861989,898.6837618785667,528.1717608063922),c(735.8001938305858,749.2220639830146,662.6762143210073,487.2115018050801),c(809.8208121321059,917.0639838496547,878.0717489914435,836.4516048688986),c(778.9788878398059,739.1111049549037,677.1046233419935,634.8840145203367),c(896.1782001505459,881.6756272512667,846.1231290164028,759.920594629605))
targetgene="CLEC2D"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(815.9891969905659,683.500830300294,753.3690710243488,833.2179002109003),c(912.9209590520802,764.3885025251808,660.615013032295,697.4023045749709),c(828.3259667074859,650.1346655075281,570.9527569733095,718.960335628293),c(808.9396142951831,886.7311067653222,762.6444768235542,886.0350762915394),c(980.7731924951403,906.9530248215439,784.2870903550335,655.3641440209927))
targetgene="LRP1B"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("AXL","CLDN7","NFATC1","LTBR","FLT1","CSF1R","PILRA","TYRO3","FPR3","C1QA")
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
targetmat=list(c(756.9489419167345,902.9086412102996,935.7853850753881,887.1129778442056),c(762.2361289382716,791.6880919010802,1048.120855310209,750.21948065561),c(884.7226282705487,899.8753535018664,1007.927430180319,1131.7966302994118),c(871.5046607167059,977.72973801832,1315.046422198453,1431.4532619405893),c(823.0387796859487,838.1985034303901,1083.1612772183182,1080.0573557714386))
targetgene="AXL"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(934.0697071382288,877.6312436400224,921.356976054402,1223.4182622760309),c(1037.169854058203,928.1860387805767,1071.8246701304006,1116.7060085620863),c(890.0098152920859,949.4190527396095,964.6422031173605,1299.9492725153243),c(878.5542434120888,1108.16110948095,1007.927430180319,1163.0557753267287),c(949.0500703659173,763.3774066223698,877.0411483470874,1060.6551278234488))
targetgene="CLDN7"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(837.1379450767145,748.2109680802034,940.938388297169,1158.7441691160643),c(901.4653871720831,917.0639838496547,843.0313270833343,1080.0573557714386),c(809.8208121321059,815.9543935685462,890.4389567237174,1150.1209566947355),c(778.097690002883,912.0085043355994,746.1548665138558,1057.4214231654505),c(925.2577287690002,772.4772697476695,744.0936652251435,1036.9412936647943))
targetgene="NFATC1"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(869.7422650428603,846.2872706528788,912.0815702551967,1233.1193762500257),c(904.9901785197745,941.3302855171208,1072.8552707747567,1066.0446355867793),c(852.1183083044031,711.8115155790043,880.1329502801558,827.8283924475697),c(882.0790347597803,784.6104205814025,929.6017812092513,874.1781592122123),c(964.030433593606,1053.5619307291513,1316.0770228428091,1048.7982107441217))
targetgene="LTBR"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(913.8021568890031,946.3857650311762,973.9176089165659,1026.1622781381334),c(1006.3279297659031,1061.65069795164,1237.7513738717414,1132.8745318520778),c(952.5748617136088,976.7186421155088,1071.8246701304006,1264.3785212773428),c(912.0397612151573,1176.9156308721037,1008.9580308246751,1254.677407303348),c(916.4457503997717,1035.3622044785518,1011.0192321133875,1456.24499765191))
targetgene="FLT1"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(871.5046607167059,1182.9822062889702,1300.6180131774668,1333.3642206479738),c(838.0191429136373,998.9627519773527,1126.4465042812767,1139.3419411680745),c(857.4054953259403,946.3857650311762,1103.7732901054414,728.6614496022879),c(821.2763840121031,831.1208321107125,1030.6006443561544,854.7759312642224),c(931.4261136274603,1004.0182314914082,983.1930147157713,942.085957030177))
targetgene="CSF1R"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(897.9405958243916,1420.5897434495755,1198.5885493862077,1595.2942979458376),c(823.0387796859487,909.9863125299772,835.8171225728412,752.3752837609422),c(963.1492357566831,934.2526141974432,1035.7536475779352,1005.6821486374773),c(845.9499234459431,805.8434345404353,988.3460179375521,1291.3260600939955),c(886.4850239443945,1128.3830275371718,1029.5700437117982,931.306941503516))
targetgene="PILRA"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(954.3372573874545,946.3857650311762,949.1831934520181,983.046216031489),c(874.1482542274745,811.9100099573019,1030.6006443561544,1088.6805681927674),c(903.2277828459288,803.8212427348132,940.938388297169,1128.5629256414134),c(838.9003407505602,775.5105574561028,751.3078697356366,935.6185477141804),c(897.0593979874687,853.3649419725564,969.7952063391413,695.2465014696387))
targetgene="TYRO3"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(848.5935169567116,1132.427411148416,1343.9032402404252,1321.5073035686464),c(820.3951861751802,990.8739847548641,942.9995895858813,746.9857759976118),c(808.9396142951831,844.2650788472566,955.3667973181551,734.0509573656185),c(929.6637179536145,1025.251245450441,1154.2727216788928,1040.1749983227928),c(886.4850239443945,878.6423395428335,1022.3558392013051,1020.7727703748028))
targetgene="FPR3"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(809.8208121321059,925.1527510721435,1005.8662288916067,1029.3959827961316),c(941.1192898336117,1060.639602048829,839.9395251502658,975.5009051628264),c(887.3662217813173,1124.3386439259273,844.0619277276904,1165.211578432061),c(907.6337720305431,1146.5827537877713,1257.3327861145083,1074.6678480081082),c(671.4727517352173,927.1749428777656,844.0619277276904,960.4102834255008))
targetgene="C1QA"
collabel=c("mCMV-VPH-lib-plasmid-1_S10_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-T0_S5_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd2-unsorted_S7_L001_R1_001-clipped-trimmed-aligned-counts","TRE-VPH-7d-dox-rd1-unsorted_S6_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("Survival_CRISPRa_summary.Rnw");
library(tools);

texi2dvi("Survival_CRISPRa_summary.tex",pdf=TRUE);

