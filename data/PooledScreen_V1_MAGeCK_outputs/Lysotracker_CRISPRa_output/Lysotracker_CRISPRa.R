pdf(file='Lysotracker_CRISPRa.pdf',width=4.5,height=4.5);
gstable=read.table('Lysotracker_CRISPRa.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ICAM1","MS4A4A","Non","PLAUR","HLA-DRB5","SORL1","C9orf72","PILRB","IL15","ITGAX")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='0,3_vs_1,2 neg.'


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
targetmat=list(c(889.9628361959963,1209.3265622547867,909.4507322536075,1020.5713364103167),c(1242.7408973908055,1661.8360112030648,1164.8581515539047,886.2525996379871),c(1535.3863345183179,1334.8040732076936,1228.9604058096656,1285.4428640641627),c(1110.449124442752,1111.5133844253555,903.4411459171299,813.4443124156029),c(856.2885667183099,1218.2186693301896,976.5577796776072,1163.6772802612097))
targetgene="ICAM1"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1817.6087834741652,1873.2705572181992,1479.3598364962315,974.1246704236232),c(1116.8632710099305,957.3835284517062,1127.7990358122931,1088.3583624449502),c(1143.321625599541,1102.6212773499526,890.4203755214285,943.9971032971195),c(541.1936166056734,579.9629836957189,446.71258434483366,465.72197516387143),c(830.6319804495964,919.8390763555609,899.4347550261449,818.4655736033536))
targetgene="MS4A4A"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1309.287668025281,1425.701167756256,1123.792644921308,1134.8050284316437),c(699.1419758224403,762.7451846901106,559.8931270151614,713.0190886605902),c(1034.2811339575092,953.4314808626383,798.2733850287723,893.7844914196131),c(1158.5552236965896,1075.944956123744,853.3612597798168,1030.613858785818),c(1075.1713183232712,1132.261634267962,706.126394536116,856.1250325114834),c(767.2922830987103,875.3785409785467,483.77170008644543,689.1680980187747),c(793.7506376883209,802.2656605807899,752.1998897824441,981.6565622052492),c(799.363015934602,937.6232905063665,808.2893622562349,804.6571053370393),c(850.6761884720288,925.7671477391627,849.3548688888317,1257.8259275315343),c(1031.07406067392,881.3066123621486,913.4571231445926,1166.1879108550852),c(529.9688601131113,496.9699843252923,553.8835406786839,716.7850345514032),c(934.8618621662447,961.3355760407742,1335.1297644207696,1137.315659025519),c(988.5803396663633,889.2107075402844,893.4251686896673,940.2311574063065),c(893.9716778004827,1000.8560519314534,961.5338138364133,1039.4010658643815),c(873.9274697780503,851.6662554441391,803.2813736425036,813.4443124156029),c(1694.136462055982,1278.4873950634756,1243.9843716508597,2252.03564270616),c(1019.8493041813579,1122.3815152952923,1220.9476240276956,859.8909784022964),c(734.4197819419212,659.0039354770774,739.1791193867426,900.0610679043014),c(1156.951687054795,1024.568337465861,1043.664827101607,1308.0385394090406),c(989.3821079872606,994.9279805478516,895.4283641351599,908.8482749828651),c(173.18195731381547,171.91407012445498,186.29717643080505,136.82936736620482),c(338.346231418658,269.72724795388626,350.55920296119234,313.8288242344147),c(1294.0540699282324,1179.6862053367772,1003.6009181917564,1334.4001606447314),c(815.3983823525479,618.4954476891312,869.386823343757,861.146293699234),c(970.9414366066229,970.227683116177,885.4123869076972,981.6565622052492),c(867.5133232108719,801.277648683523,771.2302465146231,852.3590866206704),c(546.8059948519544,657.0279116825435,776.2382351283544,836.0399877604808),c(853.883261755618,831.9060174987994,911.4539276991001,975.379985720561),c(882.7469213079206,741.008922950237,784.2510169103246,677.8702603463358),c(749.6533800389698,679.752185319684,816.302144038205,883.7419690441118),c(964.5272900394446,795.3495772999211,864.3788347300257,972.8693551266856),c(895.5752144422772,1201.422467076651,1041.6616316561144,1026.847912895005),c(789.7417960838345,692.5963399841548,770.2286487918768,799.6358441492887),c(1433.5617577643616,1276.5113712689417,1372.1888801623813,1392.1446643038637),c(756.8692949270454,780.5293988409163,679.0832560219668,524.7217941199414),c(1524.1615780257557,1349.6242516666985,1371.187282439635,1378.3361960375494),c(861.900944964591,954.4194927599053,1420.265570854202,1024.3372823011296),c(933.2583255244501,1019.6282779795262,920.4683072038164,782.0614299921615),c(1113.6561977263411,922.8031120473619,1090.7399200706814,1177.485748527524),c(789.7417960838345,933.6712429172986,903.4411459171299,888.7632302318625),c(920.4300323900934,888.2226956430175,1018.6248840329503,1016.8053905195037),c(1176.9958950772275,1103.6092892472197,1182.8869105633376,1117.2306142745165),c(1341.3584008611726,1029.5083969521959,1099.7542995753977,1110.9540377898281),c(1599.5278001901013,1838.6901408138547,1375.1936733306202,1915.6111431268673),c(1175.392358435433,1146.0938008296998,980.5641705685922,1120.9965601653294),c(1064.7483301516063,988.9999091642496,1635.6090812446487,1634.4205166128318),c(911.6105808602232,926.7551596364298,864.3788347300257,901.3163832012391),c(1180.2029683608166,1016.6642422877252,1129.8022312577855,1185.01764030915),c(1126.4844908606979,924.7791358418958,1027.6392635376667,1010.5288140348154),c(2095.020622504629,2225.990804542512,826.3181212656676,1206.35800035709),c(901.9893610094556,807.2057200671248,1006.6057113599951,794.6145829615381),c(873.125701457153,967.263647424376,947.5114457179656,1093.379623632701),c(838.6496636585694,775.5893393545814,826.3181212656676,798.3805288523511),c(1144.9251622413356,1063.1008014592733,813.2973508699662,956.550256266496),c(884.3504579497152,1161.9019911859716,1173.8725310586212,888.7632302318625),c(679.8995361209052,689.6323042923539,733.169533050265,740.6360251932188),c(732.8162453001266,707.4165184431596,1055.6839997745622,928.9333197338675),c(996.5980228753363,1124.3575390898263,1029.6424589831593,883.7419690441118),c(801.7683208972938,848.7022197523381,649.0353243395789,843.5718795421068),c(1020.6510725022551,1032.4724326439969,968.544997895637,955.2949409695584),c(845.0638102257477,1142.141753240632,967.5434001728909,981.6565622052492),c(1158.5552236965896,877.3545647730806,1125.7958403668006,1104.6774613051398),c(1160.1587603383844,1087.801098890948,973.5529865093685,1142.3369202132697),c(723.9967937702564,695.5603756759558,670.0688765172505,861.146293699234),c(1108.8455878009574,1147.081812726967,1034.6504475968904,1053.2095341306958),c(1034.2811339575092,930.7072072254977,1004.6025159145025,1110.9540377898281),c(946.888386979704,944.5393737872355,785.2526146330708,893.7844914196131),c(736.0233185837158,826.9659580124645,736.1743262185039,716.7850345514032),c(744.842770113586,613.5553882027963,818.3053394836975,692.9340439095877),c(1233.119677540038,1364.444430125703,1311.0914190748592,987.9331386899376),c(965.3290583603418,804.2416843753239,1112.7750699710991,1197.5707932785265),c(986.9768030245688,829.9299937042655,981.5657682913385,908.8482749828651),c(617.3616070909163,666.9080306552133,591.9442541430419,779.5507993982861),c(936.4653988080393,833.8820412933334,990.5801477960549,1105.9327766020774),c(835.4425903749802,798.313612991722,850.356466611578,947.7630491879324),c(927.645947278169,979.1197901915798,910.4523299763538,805.912420633977),c(1352.5831573537348,1355.5523230503002,1150.8357834354572,1340.6767371294197),c(1072.766013360579,1189.566324309447,1089.738322347935,1363.2724124742974),c(1094.4137580248062,1041.3645397193998,1226.957210364173,766.9976464289095),c(819.4072239570344,1147.081812726967,825.3165235429213,756.9551240534083),c(406.496538694928,437.6892704892733,364.58157107964,455.67945278837016),c(1026.2634507485361,665.9200187579463,986.5737569050698,956.550256266496),c(994.9944862335417,1148.069824624234,1014.6184931419652,967.848093938935),c(724.7985620911537,732.1168158748342,856.3660529480555,736.8700793024058),c(721.5914888075645,982.0838258833808,929.4826867085328,690.4234133157124),c(618.1633754118136,741.996934847504,732.1679353275189,547.3174694648193),c(1089.6031480994225,865.4984220058768,1033.6488498741442,994.2097151746258),c(1210.6701645549138,1441.5093581125277,1146.829392544472,1207.613315654028),c(1360.6008405627076,1005.7961114177883,969.5465956183833,1171.2091720428357),c(1236.3267508236272,1074.9569442264772,970.5481933411296,1179.9963791213993),c(817.8036873152398,987.0238853697157,1005.6041136372488,992.9543998776882),c(952.5007652259851,970.227683116177,879.4028005712196,992.9543998776882),c(608.542155561046,586.8790669765878,557.889931569669,466.9772904608091),c(971.7432049275202,1170.7940982613745,795.2685918605334,750.67854756872),c(1221.894921047476,847.7142078550711,786.2542123558171,548.5727847617569),c(1023.8581457858443,999.8680400341865,1001.5977227462638,824.7421500880419),c(600.5244723520731,409.03692546853085,656.0465083988028,814.6996277125406),c(1161.762296980179,805.2296962725909,824.3149258201751,1044.4223270521322),c(887.5575312333043,1169.8060863641076,1114.7782654165917,1016.8053905195037),c(929.2494839199636,1066.0648371510742,1229.962003532412,962.8268327511844))
targetgene="Non"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1067.9554034351954,1362.468406331169,835.3325007703841,832.2740418696678),c(1239.5338241072163,923.7911239446288,1006.6057113599951,1237.7408827805316),c(853.883261755618,863.5223982113429,924.4746980948015,1083.3371012571995),c(1043.9023538082765,939.5993143009006,944.5066525497267,831.0187265727302),c(1015.0386942559741,1079.897003712812,991.5817455188012,1156.1453884795837))
targetgene="PLAUR"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(893.1699094795854,821.0378866288626,856.3660529480555,896.2951220134885),c(943.6813136961149,1005.7961114177883,716.1423717635786,834.7846724635432),c(1494.4961501525559,1377.288584790174,1084.7303337342037,957.8055715634337),c(824.2178338824181,957.3835284517062,835.3325007703841,969.1034092358726),c(1142.5198572786437,934.6592548145657,1221.949221750442,1105.9327766020774))
targetgene="HLA-DRB5"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(990.1838763081579,1043.3405635139338,895.4283641351599,775.7848535074731),c(1263.5868737341352,1351.6002754612323,895.4283641351599,913.8695361706157),c(1054.3253419799414,1048.2806230002686,992.5833432415474,954.0396256726208),c(732.0144769792294,952.4434689653713,777.2398328511007,839.8059336512938),c(943.6813136961149,870.4384814922117,848.3532711660854,871.1888160747352))
targetgene="SORL1"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(994.1927179126444,688.644292395087,715.1407740408323,807.1677359309147),c(1100.0261362710871,1026.544361260395,1151.8373811582035,923.9120585461169),c(1039.0917438828928,1023.580325568594,861.3740415617868,841.0612489482314),c(850.6761884720288,756.8171133065088,952.5194343316969,1118.485929571454),c(1355.790230637324,1175.7341577477093,857.3676506708018,911.3589055767403))
targetgene="C9orf72"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1073.5677816814766,1306.151728186951,978.5609751230998,1265.3578193131602),c(577.2731910460516,738.0448872584361,769.2270510691305,710.5084580667149),c(776.1117346285805,714.3326017240285,770.2286487918768,745.6572863809694),c(1472.0466371674315,1549.202654914629,1148.8325879899646,1002.9969222531894),c(1266.7939470177243,1320.9719066459559,1089.738322347935,1375.8255654436741))
targetgene="PILRB"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1062.3430251889145,1109.5373606308215,1121.7894494758154,1217.6558380295291),c(941.276008733423,1138.189705651564,1007.6073090827414,981.6565622052492),c(1460.8218806748694,1666.7760706893998,1076.7175519522336,1179.9963791213993),c(245.34110619457192,226.25472447413904,267.4265919732524,214.65891577633965),c(929.2494839199636,1007.7721352123224,916.4619163128314,1041.9116964582568))
targetgene="IL15"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1227.507299293757,1594.65120218891,1056.6855974973082,990.4437692838128),c(1018.2457675395632,1342.7081683858294,931.4858821540254,1018.0607058164413),c(1100.8279045919844,1024.568337465861,1080.7239428432185,1231.4643062958432),c(908.403507576634,772.6253036627804,942.5034571042343,728.0828722238422),c(836.2443586958775,779.5413869436493,827.3197189884139,921.4014279522416))
targetgene="ITGAX"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetgenelist=c("FCGR2C","Non","TGFBR2","IL1R1","LTBR","SPI1","PTPRE","TLR7","CD48","CD80")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='0,3_vs_1,2 pos.'


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
targetmat=list(c(1144.1233939204383,781.5174107381833,947.5114457179656,982.9118775021869),c(957.3113751513689,1023.580325568594,1020.6280794784428,872.4441313716729),c(923.6371056736825,850.6782435468721,1049.6744134380845,1384.6127725222377),c(586.8944108968191,628.375566661801,679.0832560219668,1004.2522375501271),c(757.6710632479427,772.6253036627804,818.3053394836975,701.7212509881513))
targetgene="FCGR2C"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1309.287668025281,1425.701167756256,1123.792644921308,1134.8050284316437),c(699.1419758224403,762.7451846901106,559.8931270151614,713.0190886605902),c(1034.2811339575092,953.4314808626383,798.2733850287723,893.7844914196131),c(1158.5552236965896,1075.944956123744,853.3612597798168,1030.613858785818),c(1075.1713183232712,1132.261634267962,706.126394536116,856.1250325114834),c(767.2922830987103,875.3785409785467,483.77170008644543,689.1680980187747),c(793.7506376883209,802.2656605807899,752.1998897824441,981.6565622052492),c(799.363015934602,937.6232905063665,808.2893622562349,804.6571053370393),c(850.6761884720288,925.7671477391627,849.3548688888317,1257.8259275315343),c(1031.07406067392,881.3066123621486,913.4571231445926,1166.1879108550852),c(529.9688601131113,496.9699843252923,553.8835406786839,716.7850345514032),c(934.8618621662447,961.3355760407742,1335.1297644207696,1137.315659025519),c(988.5803396663633,889.2107075402844,893.4251686896673,940.2311574063065),c(893.9716778004827,1000.8560519314534,961.5338138364133,1039.4010658643815),c(873.9274697780503,851.6662554441391,803.2813736425036,813.4443124156029),c(1694.136462055982,1278.4873950634756,1243.9843716508597,2252.03564270616),c(1019.8493041813579,1122.3815152952923,1220.9476240276956,859.8909784022964),c(734.4197819419212,659.0039354770774,739.1791193867426,900.0610679043014),c(1156.951687054795,1024.568337465861,1043.664827101607,1308.0385394090406),c(989.3821079872606,994.9279805478516,895.4283641351599,908.8482749828651),c(173.18195731381547,171.91407012445498,186.29717643080505,136.82936736620482),c(338.346231418658,269.72724795388626,350.55920296119234,313.8288242344147),c(1294.0540699282324,1179.6862053367772,1003.6009181917564,1334.4001606447314),c(815.3983823525479,618.4954476891312,869.386823343757,861.146293699234),c(970.9414366066229,970.227683116177,885.4123869076972,981.6565622052492),c(867.5133232108719,801.277648683523,771.2302465146231,852.3590866206704),c(546.8059948519544,657.0279116825435,776.2382351283544,836.0399877604808),c(853.883261755618,831.9060174987994,911.4539276991001,975.379985720561),c(882.7469213079206,741.008922950237,784.2510169103246,677.8702603463358),c(749.6533800389698,679.752185319684,816.302144038205,883.7419690441118),c(964.5272900394446,795.3495772999211,864.3788347300257,972.8693551266856),c(895.5752144422772,1201.422467076651,1041.6616316561144,1026.847912895005),c(789.7417960838345,692.5963399841548,770.2286487918768,799.6358441492887),c(1433.5617577643616,1276.5113712689417,1372.1888801623813,1392.1446643038637),c(756.8692949270454,780.5293988409163,679.0832560219668,524.7217941199414),c(1524.1615780257557,1349.6242516666985,1371.187282439635,1378.3361960375494),c(861.900944964591,954.4194927599053,1420.265570854202,1024.3372823011296),c(933.2583255244501,1019.6282779795262,920.4683072038164,782.0614299921615),c(1113.6561977263411,922.8031120473619,1090.7399200706814,1177.485748527524),c(789.7417960838345,933.6712429172986,903.4411459171299,888.7632302318625),c(920.4300323900934,888.2226956430175,1018.6248840329503,1016.8053905195037),c(1176.9958950772275,1103.6092892472197,1182.8869105633376,1117.2306142745165),c(1341.3584008611726,1029.5083969521959,1099.7542995753977,1110.9540377898281),c(1599.5278001901013,1838.6901408138547,1375.1936733306202,1915.6111431268673),c(1175.392358435433,1146.0938008296998,980.5641705685922,1120.9965601653294),c(1064.7483301516063,988.9999091642496,1635.6090812446487,1634.4205166128318),c(911.6105808602232,926.7551596364298,864.3788347300257,901.3163832012391),c(1180.2029683608166,1016.6642422877252,1129.8022312577855,1185.01764030915),c(1126.4844908606979,924.7791358418958,1027.6392635376667,1010.5288140348154),c(2095.020622504629,2225.990804542512,826.3181212656676,1206.35800035709),c(901.9893610094556,807.2057200671248,1006.6057113599951,794.6145829615381),c(873.125701457153,967.263647424376,947.5114457179656,1093.379623632701),c(838.6496636585694,775.5893393545814,826.3181212656676,798.3805288523511),c(1144.9251622413356,1063.1008014592733,813.2973508699662,956.550256266496),c(884.3504579497152,1161.9019911859716,1173.8725310586212,888.7632302318625),c(679.8995361209052,689.6323042923539,733.169533050265,740.6360251932188),c(732.8162453001266,707.4165184431596,1055.6839997745622,928.9333197338675),c(996.5980228753363,1124.3575390898263,1029.6424589831593,883.7419690441118),c(801.7683208972938,848.7022197523381,649.0353243395789,843.5718795421068),c(1020.6510725022551,1032.4724326439969,968.544997895637,955.2949409695584),c(845.0638102257477,1142.141753240632,967.5434001728909,981.6565622052492),c(1158.5552236965896,877.3545647730806,1125.7958403668006,1104.6774613051398),c(1160.1587603383844,1087.801098890948,973.5529865093685,1142.3369202132697),c(723.9967937702564,695.5603756759558,670.0688765172505,861.146293699234),c(1108.8455878009574,1147.081812726967,1034.6504475968904,1053.2095341306958),c(1034.2811339575092,930.7072072254977,1004.6025159145025,1110.9540377898281),c(946.888386979704,944.5393737872355,785.2526146330708,893.7844914196131),c(736.0233185837158,826.9659580124645,736.1743262185039,716.7850345514032),c(744.842770113586,613.5553882027963,818.3053394836975,692.9340439095877),c(1233.119677540038,1364.444430125703,1311.0914190748592,987.9331386899376),c(965.3290583603418,804.2416843753239,1112.7750699710991,1197.5707932785265),c(986.9768030245688,829.9299937042655,981.5657682913385,908.8482749828651),c(617.3616070909163,666.9080306552133,591.9442541430419,779.5507993982861),c(936.4653988080393,833.8820412933334,990.5801477960549,1105.9327766020774),c(835.4425903749802,798.313612991722,850.356466611578,947.7630491879324),c(927.645947278169,979.1197901915798,910.4523299763538,805.912420633977),c(1352.5831573537348,1355.5523230503002,1150.8357834354572,1340.6767371294197),c(1072.766013360579,1189.566324309447,1089.738322347935,1363.2724124742974),c(1094.4137580248062,1041.3645397193998,1226.957210364173,766.9976464289095),c(819.4072239570344,1147.081812726967,825.3165235429213,756.9551240534083),c(406.496538694928,437.6892704892733,364.58157107964,455.67945278837016),c(1026.2634507485361,665.9200187579463,986.5737569050698,956.550256266496),c(994.9944862335417,1148.069824624234,1014.6184931419652,967.848093938935),c(724.7985620911537,732.1168158748342,856.3660529480555,736.8700793024058),c(721.5914888075645,982.0838258833808,929.4826867085328,690.4234133157124),c(618.1633754118136,741.996934847504,732.1679353275189,547.3174694648193),c(1089.6031480994225,865.4984220058768,1033.6488498741442,994.2097151746258),c(1210.6701645549138,1441.5093581125277,1146.829392544472,1207.613315654028),c(1360.6008405627076,1005.7961114177883,969.5465956183833,1171.2091720428357),c(1236.3267508236272,1074.9569442264772,970.5481933411296,1179.9963791213993),c(817.8036873152398,987.0238853697157,1005.6041136372488,992.9543998776882),c(952.5007652259851,970.227683116177,879.4028005712196,992.9543998776882),c(608.542155561046,586.8790669765878,557.889931569669,466.9772904608091),c(971.7432049275202,1170.7940982613745,795.2685918605334,750.67854756872),c(1221.894921047476,847.7142078550711,786.2542123558171,548.5727847617569),c(1023.8581457858443,999.8680400341865,1001.5977227462638,824.7421500880419),c(600.5244723520731,409.03692546853085,656.0465083988028,814.6996277125406),c(1161.762296980179,805.2296962725909,824.3149258201751,1044.4223270521322),c(887.5575312333043,1169.8060863641076,1114.7782654165917,1016.8053905195037),c(929.2494839199636,1066.0648371510742,1229.962003532412,962.8268327511844))
targetgene="Non"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1442.3812092942317,1235.0148715837283,1932.082007177543,2259.567534487786),c(1031.07406067392,955.4075046571722,1048.6728157153382,991.6990845807505),c(743.2392334717914,626.399542867267,810.2925577017274,607.5726037178268),c(746.4463067553806,736.0688634639021,788.2574078013096,785.8273758829745),c(1118.466807651725,1061.1247776647394,1212.9348422457253,906.3376443889897))
targetgene="TGFBR2"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(829.8302121286991,818.0738509370617,794.2669941377872,728.0828722238422),c(425.73897839646304,358.6483187079147,456.7285615722963,495.84954229037527),c(815.3983823525479,781.5174107381833,978.5609751230998,1040.6563811613191),c(1096.819062987498,1170.7940982613745,1553.478067979455,1304.2725935182275),c(988.5803396663633,1180.6742172340444,1089.738322347935,1094.6349389296386))
targetgene="IL1R1"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(903.5928976512502,925.7671477391627,1008.6089068054877,1077.0605247725114),c(982.9679614200822,962.3235879380411,1185.8917037315764,1276.655656985599),c(738.4286235464077,881.3066123621486,994.5865386870399,1039.4010658643815),c(993.3909495917471,1073.9689323292102,1083.7287360114574,1144.847550807145),c(2285.0397145572874,2451.257517119384,1489.3758137236944,2936.182479537184))
targetgene="LTBR"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1378.2397436224483,1466.2096555442024,1173.8725310586212,1107.188091899015),c(825.8213705242127,692.5963399841548,1041.6616316561144,1077.0605247725114),c(1055.928878621736,1026.544361260395,1085.73193145695,1089.6136777418878),c(1066.3518667934009,932.6832310200316,1079.7223451204725,1110.9540377898281),c(709.5649639941051,811.1577676561927,772.2318442373694,1046.9329576460075))
targetgene="SPI1"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1083.189001532244,1100.6452535554188,1234.9699921461433,1131.0390825408306),c(927.645947278169,1040.3765278221329,1111.7734722483528,1326.8682688631054),c(981.3644247782877,1095.7051940690837,1311.0914190748592,1382.1021419283625),c(1188.2206515697894,1114.4774201171565,1080.7239428432185,1309.2938547059782),c(706.3578907105159,823.0139104233966,816.302144038205,710.5084580667149))
targetgene="PTPRE"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(532.3741650758031,552.2986505722433,656.0465083988028,715.5297192544656),c(911.6105808602232,884.2706480539495,951.5178366089506,779.5507993982861),c(881.1433846661259,941.5753380954345,1003.6009181917564,1123.5071907592046),c(736.825086904613,745.9489824365719,914.4587208673388,942.7417880001818),c(871.5221648153585,906.9949216910901,862.3756392845331,758.210439350346))
targetgene="TLR7"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1031.07406067392,1242.918966761864,946.5098479952193,1001.7416069562518),c(920.4300323900934,920.8270882528278,881.4059960167122,792.1039523676627),c(725.600330412051,623.4355071754661,932.4874798767715,916.3801667644909),c(1031.07406067392,1156.9619316996368,1053.6808043290696,1335.655475941669),c(1023.8581457858443,1230.0748120973935,1060.6919883882933,862.4016089961716))
targetgene="CD48"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
targetmat=list(c(1181.8065050026112,1148.069824624234,840.3404893841154,1064.5073718031347),c(1286.8381550401566,1190.5543362067142,1049.6744134380845,1088.3583624449502),c(1020.6510725022551,1047.2926111030017,985.5721591823236,1104.6774613051398),c(919.6282640691961,739.032899155703,1172.870933335875,1073.2945788816983),c(950.8972285841905,1206.3625265629857,951.5178366089506,1233.9749368897187))
targetgene="CD80"
collabel=c("PMA-THP1-LysoTracker-VPH-Low-Rep1_S3_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-Low-Rep2_S7_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep1_S4_L001_R1_001-clipped-trimmed-aligned-counts","PMA-THP1-LysoTracker-VPH-High-Rep2_S8_L001_R1_001-clipped-trimmed-aligned-counts")

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
Sweave("Lysotracker_CRISPRa_summary.Rnw");
library(tools);

texi2dvi("Lysotracker_CRISPRa_summary.tex",pdf=TRUE);

