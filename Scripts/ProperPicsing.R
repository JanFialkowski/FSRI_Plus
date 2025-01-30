library(data.table)
library(Matrix)
library(graphics)
library(glue)
library("RColorBrewer")
library(viridis)

risk_measures<-function(lost_bt, num_banks, num_trials){
  
  el <- Matrix(0, 1, num_banks)
  var <- Matrix(0, 1, num_banks)
  es <- Matrix(0, 1, num_banks)
  
  for (j in 1:num_banks){
    
    el[,j] <- mean(lost_bt[j,])
    var[,j] <- quantile(lost_bt[j,],0.95)
    es[,j] <- mean(lost_bt[j,which(lost_bt[j,] >= var[j])])
  }
  
  return(list(el_b=el, var_b=var, es_b=es ))
  
}
CalculateRiskMeasures <- function(Data,banks=4,trials=1000){
  Risks <- risk_measures(lost_bt = t(Data), num_banks = banks, num_trials = trials )
  Risks <-  rbind(as.vector(Risks$el_b),
                  as.vector(Risks$var_b),
                  as.vector(Risks$es_b))
  return(Risks)
}
ScenariowiseLosses <- function(Folder="./data/",Trials=c("CovidShock","CovidShock_Direct")){
  Equities <- readRDS(paste0(Folder,"Equities_",Trials[1],".rds"))
  
  Direct_risks <- readRDS(paste0(Folder,"TotalLosses_No_Cascade_",Trials[2],".rds"))
  MIB_WO_risks <- readRDS(paste0(Folder,"TotalLosses_With_Cascade_",Trials[2],".rds"))
  MSC_W_risks <- readRDS(paste0(Folder,"TotalLosses_No_Cascade_",Trials[1],".rds"))-Direct_risks
  MIB_W_risks <- readRDS(paste0(Folder,"TotalLosses_With_Cascade_",Trials[1],".rds"))
  
  Direct_risks <- t(t(Direct_risks)/Equities)
  MIB_WO_risks <- t(t(MIB_WO_risks)/Equities)
  MIB_W_risks <- t(t(MIB_W_risks)/Equities)
  MSC_W_risks <- t(t(MSC_W_risks)/Equities)
  return(list(Direct_risks=Direct_risks,MIB_WO_risks=MIB_WO_risks,MIB_W_risks=MIB_W_risks,MSC_W_risks=MSC_W_risks))
}

GetPlotParams <- function(){
  output <- list()
  output$cyan<-brewer.pal(12,"Set3")[1]
  output$lightyellow<-brewer.pal(12,"Set3")[2]
  output$lightviolet<-brewer.pal(12,"Set3")[3]
  output$red<-brewer.pal(12,"Set3")[4]
  output$blue<-brewer.pal(12,"Set3")[5]
  output$orange<-brewer.pal(12,"Set3")[6]
  output$green<-brewer.pal(12,"Set3")[7]
  output$pink<-brewer.pal(12,"Set3")[8]
  output$gray<-brewer.pal(12,"Set3")[9]
  output$violet<-brewer.pal(12,"Set3")[10]
  output$lightgreen<-brewer.pal(12,"Set3")[11]
  output$yellow<-brewer.pal(12,"Set3")[12]
  output$lightred <- brewer.pal(6,"Reds")[2]
  output$darkred <- brewer.pal(6,"Reds")[5]
  output$darkviolet <- brewer.pal(6,"Purples")[5]
  output$darkgreen <- brewer.pal(6,"Greens")[5]
  output$cex_points <- 2.2
  output$lwd <- 4
  output$cex_xylab <- 2
  output$cex_ax <- 2
  output$cex.names <- 1.5
  return(output)
}
PlotInterbankScatter <- function(MIB_WO,MIB_W,logplot=F,title="",ylabel=expression("IB"["W"]), xlabel=expression("IB"["WO"]),legendoffset=c(0,0),textoffset=c(0,0)){
  library(viridis)
  nb <- dim(MIB_WO)[2]
  limits <- c(max(10^-5,min(MIB_WO[MIB_WO>0])),max(MIB_W))
  if(!logplot){
    limits <- c(min(MIB_WO),max(MIB_W))
  }
  EvilBanks <- which(colSums(MIB_WO)==0)
  NB <- nb-length(EvilBanks)
  c <- viridis(nb)
  
  plotlog = ""
  if(logplot){
    plotlog="xy"
  }
  
  plot(NULL,xlim=limits,ylim=limits,log=plotlog,
       ylab=ylabel, xlab=xlabel,main=title)
  for(i in (1:nb)[! (1:nb) %in% EvilBanks]){
    points(MIB_WO[,i],MIB_W[,i],col=c[i],pch=20)
  }
  legend("topright",title="Bank ID", legend=seq(nb)[! (1:nb) %in% EvilBanks],
         col=c[! (1:nb) %in% EvilBanks],pch=rep(20,nb)[! (1:nb) %in% EvilBanks],
         y.intersp=rep(0.8,nb),bty="n",cex=par("cex.lab"),inset = legendoffset)
  
  length <- NB*dim(MIB_WO)[1]
  if(!logplot){
    BigFit <- lm((MIB_W[,! (1:nb) %in% EvilBanks][1:length]~MIB_WO[,! (1:nb) %in% EvilBanks][1:length]))
  }else{
    xvals = log(MIB_WO[,! (1:nb) %in% EvilBanks][1:length])
    yvals = log(MIB_W[,! (1:nb) %in% EvilBanks][1:length])
    BigFit <- lm(yvals[xvals>-Inf]~xvals[xvals>-Inf])
  }
  if(logplot){
    text(limits[1]+textoffset[1],limits[2]+textoffset[2],adj=c(0,1),cex = par("cex.lab"),labels=glue("y = {format(exp(BigFit$coefficients[1]),digits=3)}*x^{format(BigFit$coefficients[2],digits=3)} 
                                                             \n R\U00B2={format(summary(BigFit)$r.squared,digits=4)}")
    )
    curve(exp(BigFit$coefficients[1])*x^BigFit$coefficients[2],add=T,from=limits[1],to=limits[2])
  }else{
    text(0.05+textoffset[1],0+textoffset[2],adj=c(0,-1.5),offset=0,cex = par("cex.lab"),labels=glue("y = x*{format(BigFit$coefficients[2],digits=3)}+{format(BigFit$coefficients[1],digits=3)}"))
    text(0.05+textoffset[1],0+textoffset[2],adj=c(0,0.5),cex = par("cex.lab"),labels=glue("R\U00B2 = {format(summary(BigFit)$r.squared,digits=4)}"))
    abline(a=min(BigFit$coefficients[1],-10,na.rm=T),b=min(BigFit$coefficients[2],-10,na.rm=T))
  }
}
IBAmplifierBoxPlot <- function(){
  Risks <- ScenariowiseLosses()
  MIB_W_Banks <- list()
  MIB_W_Scenarios <- list()
  MIB_WO_Banks <- list()
  MIB_WO_Scenarios <- list()
  Amplifiers_Scenarios <- list()
  Amplifiers_Banks <- list()
  EvilBanks <- which(colSums(Risks$MIB_WO_risks)==0)
  nb <- dim(Risks$MIB_WO_risks)[2]
  c <- viridis(nb)
  for(i in (1:nb)[-EvilBanks]){
    MIB_W_Banks[[i]] <- Risks$MIB_W_risks[1:1000,i]
    MIB_WO_Banks[[i]] <- Risks$MIB_WO_risks[1:1000,i]
  }
  for(i in 1:1000){
    MIB_W_Scenarios[[i]] <- Risks$MIB_W_risks[i,(1:nb)[-EvilBanks]]
    MIB_WO_Scenarios[[i]] <- Risks$MIB_WO_risks[i,(1:nb)[-EvilBanks]]
  }
  for(i in 1:nb){
    Amplifiers_Banks[[i]] <- Risks$MIB_W_risks[1:1000,i]/Risks$MIB_WO_risks[1:1000,i]
  }
  for(i in 1:1000){
    Amplifiers_Scenarios[[i]] <- Risks$MIB_W_risks[i,(1:nb)[-EvilBanks]]/Risks$MIB_WO_risks[i,(1:nb)[-EvilBanks]]
  }
  boxplot(Amplifiers_Banks,notch=T,whisklty="solid",col=c,xlab="",ylab="",outcol=c,pch=20)
  title(xlab="Bank ID",line=2.5)
  title(ylab=expression("IB"^"W"*"/IB"^"WO"),line=2)
  text(1,max(unlist(Amplifiers_Banks)),adj=c(0.55,0.1),labels="b)",cex = par("cex.lab"))
  ScenarioPoints <- boxplot(Amplifiers_Scenarios,plot=F,range=0)
  plot(ScenarioPoints$stats[3,],pch=20,cex=0.1,col=rgb(160/255, 32/255, 240/255,1),ylim=c(min(ScenarioPoints$stats[4,]),max(ScenarioPoints$stats[4,])),xlab="",ylab="")
  polygon(x=c(seq(1000),rev(seq(1000))),y=c(ScenarioPoints$stats[2,],rev(ScenarioPoints$stats[4,])),col=rgb(160/255, 32/255, 240/255,0.75),border=NA)
  title(xlab="Synthetic COVID shock scenario",line=2.5)
  title(ylab=expression("IB"^"W"*"/IB"^"WO"),line=2)
  text(0,max(ScenarioPoints$stats[4,]),adj=c(1,0.6),labels="c)",cex = par("cex.lab"))
}

Figure2FSRI <- function(){
  Folder <- "./data/"
  Trial <- "SingleFirm_FSRI"
  OG_Data <- readRDS(paste0(Folder, "TotalLosses_No_Cascade_",Trial,".rds"))
  Cascade_Data <- readRDS(paste0(Folder, "TotalLosses_With_Cascade_",Trial,".rds"))
  Cascade_Data <- Cascade_Data + OG_Data
  Equities <- readRDS(paste0(Folder, "Equities_",Trial,".rds"))
  OG_Data <- colSums(t(OG_Data)/sum(Equities))
  Cascade_Data <- colSums(t(Cascade_Data)/sum(Equities))
  OID <- order(OG_Data,decreasing=T)
  TargetFolder <- "./pics/"
  SubPlotSize <- 0.55
  PointSampling <- 1
  PointCoordinates <- seq(1,length(OID),length.out=length(OID)%/%PointSampling)
  LOG1 <- F
  LOG2 <- F
  xvsy <- F
  logax1 <- "x"
  logax2 <- ""
  SubPlotPosition <- c(0.4,0.2)
  if(!xvsy){logax2 <- "x"}
  if(LOG1){logax1 <- "xy"
  SubPlotPosition <- c(0.5,0.1)}
  if(LOG2){logax2 <- "xy"}
  
  TextCEX <- 1.5
  SubCEX <- 0.55 + TextCEX
  jpeg(paste0(TargetFolder,"Figure2.jpeg"),
      width=8,height=4.5,units="in",res=600)
  par(cex.axis=TextCEX,cex.lab=TextCEX,cex.main=TextCEX,cex.sub=TextCEX,cex=1)
  par(mar = c(3.1,4.1,1,1)+0.1) #c(bottom, left, top, right)
  ylim <- range(Cascade_Data)
  ylim <- c(0,ylim[2]*1.1)
  if(LOG1){ylim[1] <- 10^-5}
  #par(mgp=c(2,1,0))
  plot(OG_Data[OID][PointCoordinates],log=logax1,col="#22A884FF",ylim=ylim,pch=17,xaxt="n",ylab="",xlab="")
  title(xlab="Rank",line=2)
  title(ylab="Estimated financial system losses",line=3)
  axis(1,at=c(1,10,100,1000,10000),labels=c(expression(10^0),"",expression(10^2),"",expression(10^4)))
  points(Cascade_Data[OID][PointCoordinates],col="#482576ff",pch=20)
  abline(h=axTicks(2),col="lightgrey",lwd=0.5,lty=2)
  abline(v=axTicks(1),col="lightgrey",lwd=0.5,lty=2)
  legend("topright",legend=c("Direct and SC losses, FSRI",expression("Direct, SC and interbank losses, FSRI+"))# "IB"^"W"*", FSRI+"))
         ,col=c("#22A884FF","#482576ff"),pch=c(17,20),bty="n",cex=TextCEX,pt.cex=1.5)
  par(new=T)
  
  ylim <- range(Cascade_Data)
  ylim <- c(0,ylim[2]*1.1)
  if(LOG2){ylim[1] <- 10^-5}
  layout(matrix(c(0,0,0,0,1,0,0,0,0),nrow=3),widths = c(SubPlotPosition[1],SubPlotSize,1-SubPlotPosition[1]-SubPlotSize),heights = c(SubPlotPosition[2],SubPlotSize,1-SubPlotPosition[2]-SubPlotSize))
  if(!xvsy){
    par(mar=c(5.1,4.8,4.1,2.1))
    plot(Cascade_Data[OID][PointCoordinates]/OG_Data[OID][PointCoordinates],xlab="Rank",ylab=expression("FSRI"^"+"*"/FSRI"),log=logax2,pch=20,xaxt="n",xlim=c(1,10^4),bg="white",cex.axis=SubCEX,cex.lab=SubCEX,cex.main=SubCEX,cex.sub=SubCEX,cex=1)
    #title(xlab="Rank",line=3,cex.sub=SubCEX,cex.main=SubCEX,cex.axis=SubCEX)
    axis(1,at=c(1,10,100,1000,10000),labels=c(expression(10^0),"",expression(10^2),"",expression(10^4)),cex.axis=SubCEX,padj=0.3)}
  if(xvsy){plot(OG_Data[OID][PointCoordinates],Cascade_Data[OID][PointCoordinates],xlab="FSRI",ylab="FSRI+",xlim=ylim,ylim=ylim,log=logax2,pch=20)}
  if(xvsy){abline(0,1)}
  dev.off()
}

Figure2_FSRI_CCDF <- function(){
  Folder <- "./data/"
  Trial <- "SingleFirm_FSRI"
  OG_Data <- readRDS(paste0(Folder, "TotalLosses_No_Cascade_",Trial,".rds"))
  Cascade_Data <- readRDS(paste0(Folder, "TotalLosses_With_Cascade_",Trial,".rds"))
  Cascade_Data <- Cascade_Data + OG_Data
  Equities <- readRDS(paste0(Folder, "Equities_",Trial,".rds"))
  OG_Data <- colSums(t(OG_Data)/sum(Equities))
  Cascade_Data <- colSums(t(Cascade_Data)/sum(Equities))
  OID <- order(OG_Data,decreasing=T)
  
  TextCEX <- 1.5
  PPar <- GetPlotParams()
  OGCDF <- ecdf(OG_Data)
  CascadeCDF <- ecdf(Cascade_Data)
  from <- 10^-5
  ylim <- c(5*10^-6,1)
  
  TargetFolder <- "./pics/"
  jpeg(paste0(TargetFolder,"SingleFirmFSRI_Survival.jpeg"), width=8,height=8,units="in",res=300)
  par(cex.axis=TextCEX,cex.lab=TextCEX,cex.main=TextCEX,cex.sub=TextCEX,cex=1)
  par(mar = c(4.1,5.4,1,1)+0.1) #c(bottom, left, top, right)
  plot(function(v){1-OGCDF(v)},pch=17,type="p",col="#22A884FF",log="xy",ylim=ylim,from=from,ylab="",xlab="")
  par(new=T)
  plot(function(v){1-CascadeCDF(v)},pch=20,type="p",col="#482576ff",log="xy",ylim=ylim,from=from,ylab="1-CDF",xlab="Relative equity lost in the system")
  legend("topright",legend=c("Direct and SC losses, FSRI",expression("Direct, SC and interbank losses, FSRI+"))# "IB"^"W"*", FSRI+"))
         ,col=c("#22A884FF","#482576ff"),pch=c(17,20),bty="n",cex=TextCEX,pt.cex=1.5)
  abline(h=axTicks(2),col="lightgrey",lwd=0.5,lty=2)
  abline(v=axTicks(1),col="lightgrey",lwd=0.5,lty=2)
  fit <- fit_power_law(Cascade_Data[10^-4<Cascade_Data&Cascade_Data<10^-2])
  dev.off()
}

Figure2_SI <- function(){
  Folder <- "./data/"
  Trial <- "SingleFirm_FSRI"
  OG_Data <- readRDS(paste0(Folder, "TotalLosses_No_Cascade_",Trial,".rds"))
  Cascade_Data <- readRDS(paste0(Folder, "TotalLosses_With_Cascade_",Trial,".rds"))
  Cascade_Data <- Cascade_Data + OG_Data
  Equities <- readRDS(paste0(Folder, "Equities_",Trial,".rds"))
  OG_FSRI <- colSums(t(OG_Data)/sum(Equities))
  Cascade_FSRI <- colSums(t(Cascade_Data)/sum(Equities))
  OID <- order(OG_FSRI,decreasing=T)
  Amplifications <- Cascade_FSRI[OID]/OG_FSRI[OID]
  
  banks <- dim(OG_Data)[2]
  
  count_non_zeros <- function(v){
    return(sum(v!=0))
  }
  
  select_color <- function(n){
    viridis(banks)[n]
  }
  
  jpeg("./pics/Figure2_SI.jpeg",width=8,height = 9,units="in",res=600)
  layout(matrix(c(1,2)),heights=c(1,1.2))
  par(mar=c(1,4,1,1)+0.1)
  plot(Amplifications,xlim=c(1,2e4),log="x",col=select_color(apply(FUN=which.max,X=t(OG_Data[OID,])/Equities,MARGIN =2)),ylab="FSRI+/FSRI",xlab="Rank",pch=20,xaxt="n")
  legend("topleft",col=viridis(19),legend=seq(19),cex=1.2,pch=20,ncol=3,bty="n")
  
  par(mar = c(5,4,1,1)+0.1)
  plot(Amplifications,xlim=c(1,2e4),log="x",col=ifelse(1==(apply(FUN=count_non_zeros,X=t(OG_Data[OID,])/Equities,MARGIN =2)),"red","black"),ylab="FSRI+/FSRI",xlab="Rank",pch=20)
  legend("topleft",legend=c("Only 1 bank shocked","Multiple banks shocked"),col = c("red","black"),pch=20,bty="n")
  dev.off()
}

Figure2_DebtRanks <- function(){
  Equities <- readRDS(paste0("./data/Equities_SingleFirm_FSRI.rds"))
  nb <- length(Equities)
  DRs <- readRDS("./data/DebtRank.rds")
  DebtRanks <- list()
  for(i in 1:nb){
    DebtRanks[[i]] <- sum(Equities*DRs[i,])/sum(Equities)
  }
  DebtRanksAdded <- list()
  for(i in 1:nb){
    DebtRanksAdded[[i]] <- sum(Equities*DRs[i,]-Equities[i])/sum(Equities)
  }
  Order <- order(unlist(DebtRanks),decreasing=T)
  
  TextCEX <- 1.5
  SubCEX <- 0.55 + TextCEX
  pdf("./pics/Fig2_SI_DebtRanks.pdf",width=8,height=4.5)#,units="in")
  par(cex.axis=TextCEX,cex.lab=TextCEX,cex.main=TextCEX,cex.sub=TextCEX)
  par(mar = c(3.1,4.1,1,1)+0.1) #c(bottom, left, top, right)
  plot(unlist(DebtRanks)[Order],ylab="",xlab="",
       pch=20,col="black",log="")
  points(unlist(DebtRanksAdded)[Order],pch=20,col="blue")
  legend("topright",col=c("black","blue"),pch=20,bty="n",
         legend=c("DebtRank","Additional interbank losses"),cex=TextCEX)
  title(xlab="Rank",line=2)
  title(ylab="Estimated financial system losses",line=3)
  abline(h=axTicks(2),col="lightgrey",lwd=0.5,lty=2)
  dev.off()
}

Figure3 <- function(banks=4){
  # Data ingestion for Original Plot, all cascades ####
  Folder <- "./data/"
  Target <- "Fig3_VaR"
  
  Risks <- ScenariowiseLosses()
  Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks+Risks$MIB_W_risks)
  OG_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks)
  
  # Data ingestion for direct interbank contagion ####
  
  Direct_Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MIB_WO_risks)
  Direct_OG_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  Direct_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  # Data Setup for the VaR ####
  VaR <- array(dim=c(2,banks,3))
  VaR[1,,1] <- Direct_Cascade_risks[2,]
  VaR[2,,1] <- Cascade_risks[2,]
  VaR[1,,2] <- Direct_OG_risks[2,]
  VaR[2,,2] <- OG_risks[2,]
  VaR[1,,3] <- Direct_risks[2,]
  VaR[2,,3] <- Direct_risks[2,]
  # Plotting setup ####
  
  PPar <- GetPlotParams()
  TargetFolder = "./pics/"
  ylim = c(0,ceiling(max(Cascade_risks)))
  PPar <- GetPlotParams()
  pdf(file = paste0(TargetFolder, Target,".pdf"),  width = 8, height = 4.5,# horizontal = FALSE,
      paper = "special", bg="white")
  par(mar = c(3.6,3.9,1,1)+0.1) #c(bottom, left, top, right)
  par(mgp=c(2.5,1,0))
  
  PPar$cex_ax <- 1.2
  PPar$cex.names <- 1.2
  PPar$cex_points <- 1.2
  PPar$cex_xylab <- 1.2
  d <- -30
  colors <- viridis(3,begin=(0.1),end=0.9)
  ylim <- c(0,max(ceiling(10*VaR)/10))
  # Plotting the Value at Risk ####
  mp <- barplot( VaR[,,1],  
                 col = c(colors[1],colors[1]), xlab = "Bank ID", ylab=" Bank equity losses, VaR",
                 cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
                 border = NA,cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),xaxt="n")
  abline(h=seq(0,1,0.1),col="lightgrey",lwd=0.5,lty=2)
  axis(1,at=colMeans(mp),label=seq(banks),cex.axis=PPar$cex_ax,gap.axis=0.,tick=F)
  barplot( VaR[,,1], col = c(colors[1],colors[1]),add=T, beside=T,border=NA,density=c(-d,d),yaxt="n",xaxt="n")
  barplot( VaR[,,2], col = c(colors[2],colors[2]),add=T, beside=T,border=NA,density=c(-d,d),yaxt="n",xaxt="n")
  barplot( VaR[,,3], col = c(colors[3],colors[3]),add=T, beside=T,border=NA,density=c(-d,d),yaxt="n",xaxt="n")
  legend("topright",title="",legend=c("DI","", expression("IB"^"WO"),"DI","SC",expression("IB"^"W"))
         ,fill=c(colors[3],0,colors[1],colors[c(3,2,1)]),density=c(-d,-d,-d,d,d,d),border = c("black","#00000000","black","black","black","black")
         ,cex=PPar$cex_ax,bty="n",ncol=2,inset = c(0.,-0.1))
  dev.off()
}

Figure3_SI_OtherRisks <- function(banks=4){
  # Data ingestion for Original Plot, all cascades ####
  Folder <- "./data/"
  Target <- "Fig3_OtherRisks"
  
  Risks <- ScenariowiseLosses()
  Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks+Risks$MIB_W_risks)
  OG_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks)
  
  # Data ingestion for direct interbank contagion ####
  
  Direct_Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MIB_WO_risks)
  Direct_OG_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  Direct_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  # Data Setup for the 3 different cases ####
  EL <- array(dim=c(2,banks,3))
  EL[1,,1] <- Direct_Cascade_risks[1,]
  EL[2,,1] <- Cascade_risks[1,]
  EL[1,,2] <- Direct_OG_risks[1,]
  EL[2,,2] <- OG_risks[1,]
  EL[1,,3] <- Direct_risks[1,]
  EL[2,,3] <- Direct_risks[1,]
  
  ES <- array(dim=c(2,banks,3))
  ES[1,,1] <- Direct_Cascade_risks[3,]
  ES[2,,1] <- Cascade_risks[3,]
  ES[1,,2] <- Direct_OG_risks[3,]
  ES[2,,2] <- OG_risks[3,]
  ES[1,,3] <- Direct_risks[3,]
  ES[2,,3] <- Direct_risks[3,]
  # Plotting setup ####
  
  TargetFolder = "./pics/"
  PPar <- GetPlotParams()
  par(mar = c(4.1,5.4,1,1)+0.1) #c(bottom, left, top, right)
  
  PPar$cex_ax <- 1.2
  PPar$cex.names <- 1.2
  PPar$cex_points <- 1.2
  PPar$cex_xylab <- 1.2
  d <- -30
  colors <- viridis(3,begin=(0.1),end=0.9)
  
  pdf(file = paste0(TargetFolder, Target,"_SI.pdf"),  width = 8, height = 9,# horizontal = FALSE,
      paper = "special", bg="white")
  layout((matrix(c(1,2))),heights=c(1,1.1))
  ylim <- c(0,1)
  # Plotting for EL ####
  par(mar = c(0.5,5.4,1,1)+0.1) #c(bottom, left, top, right)
  mp <- barplot( EL[,,1],
           col = c(colors[1],colors[1]), xlab = "", ylab="Bank equity losses, EL",
           cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
           border = "white",cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),yaxt="n",xaxt="n")
  par(new=TRUE)
  barplot( EL[,,2],
           col = c(colors[2],colors[2]), xlab = "", ylab="",
           cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
           border = "white",cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),yaxt="n",xaxt="n")
  par(new=TRUE)
  barplot( EL[,,3],  
           col = c(colors[3],colors[3]), xlab = "", ylab="",
           cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
           border = "white",cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),yaxt="n",xaxt="n")
  text(0,ylim[2],labels="a)",adj=c(0,1),cex=PPar$cex_ax)
  axis(2,at=seq(0,1,0.2),cex.axis=PPar$cex_ax,labels=T,tick=T)
  axis(1,at=colMeans(mp),label=F,tick=F)
  abline(h=seq(0,1,0.2),col="lightgrey",lwd=0.5,lty=2)
  legend("topright",title="",legend=c("DI","", expression("IB"^"WO"),"DI","SC",expression("IB"^"W"))
         ,fill=c(colors[3],0,colors[1],colors[c(3,2,1)]),density=c(-d,-d,-d,d,d,d),border = c("black","#00000000","black","black","black","black")
         ,cex=PPar$cex_ax,bty="n",ncol=2,inset = c(0.2,0.))
  
  
  # Plotting for ES ####
  par(mar = c(4.1,5.4,1,1)+0.1) #c(bottom, left, top, right)
  barplot( ES[,,1],
           col = c(colors[1],colors[1]), xlab = "Bank ID", ylab="Bank equity losses, ES",
           cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
           border = "white",cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),yaxt="n",xaxt="n")
  par(new=TRUE)
  barplot( ES[,,2],
           col = c(colors[2],colors[2]), xlab = "Bank ID", ylab="",
           cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
           border = "white",cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),yaxt="n",xaxt="n")
  par(new=TRUE)
  barplot( ES[,,3],
           col = c(colors[3],colors[3]), xlab = "Bank ID", ylab="",
           cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
           border = "white",cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),yaxt="n",xaxt="n")
  text(0,ylim[2],labels="b)",adj=c(0,1),cex=PPar$cex_ax)
  abline(h=seq(0,1,0.2),col="lightgrey",lwd=0.5,lty=2)
  legend("topright",title="",legend=c("DI","", expression("IB"^"WO"),"DI","SC",expression("IB"^"W"))
         ,fill=c(colors[3],0,colors[1],colors[c(3,2,1)]),density=c(-d,-d,-d,d,d,d),border = c("black","#00000000","black","black","black","black")
         ,cex=PPar$cex_ax,bty="n",ncol=2,inset = c(0.2,0.))
  axis(2,at=seq(0,1,0.2),cex.axis=PPar$cex_ax,labels=T,tick=T)
  axis(1,at=colMeans(mp),cex.axis=PPar$cex_ax,tick=F,gap.axis=0,labels=seq(banks))
  dev.off()
}

Figure3_SystemSI <- function(banks=4){
  # Data ingestion for Original Plot, all cascades ####
  Folder <- "./data/"
  Target <- "Fig3_System_SI"
  Equities <- readRDS(paste0(Folder, "Equities_CovidShock.rds"))
  
  Risks <- ScenariowiseLosses()
  Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks+Risks$MIB_W_risks)
  OG_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks)
  
  # Data ingestion for direct interbank contagion ####
  
  Direct_Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MIB_WO_risks)
  Direct_OG_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  Direct_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  Rescale <- Equities/sum(Equities)
  VaR <- array(dim=c(2,banks,3))
  VaR[1,,1] <- Direct_Cascade_risks[2,]*Rescale
  VaR[2,,1] <- Cascade_risks[2,]*Rescale
  VaR[1,,2] <- Direct_OG_risks[2,]*Rescale
  VaR[2,,2] <- OG_risks[2,]*Rescale
  VaR[1,,3] <- Direct_risks[2,]*Rescale
  VaR[2,,3] <- Direct_risks[2,]*Rescale
  # Plotting setup ####
  
  TargetFolder = "./pics/"
  PPar <- GetPlotParams()
  pdf(file = paste0(TargetFolder, Target,".pdf"),  width = 8, height = 4.5)
  par(mar = c(3.6,3.9,1,1)+0.1) #c(bottom, left, top, right)
  par(mgp=c(2.5,1,0))
  
  PPar$cex_ax <- 1.2
  PPar$cex.names <- 1.2
  PPar$cex_points <- 1.2
  PPar$cex_xylab <- 1.2
  d <- -30
  colors <- viridis(3,begin=(0.1),end=0.9)
  ylim <- c(0,max(ceiling(100*VaR)/100))
  yat = seq(0,ylim[2],0.005)
  ylabel = character()
  for(tick in yat){
    if(tick%%0.01 ==0){
      ylabel <- append(ylabel,as.character(tick))
    } else {
      ylabel <- append(ylabel,"")
    }
  }
  
  # Plotting the Value at Risk ####
  mp <- barplot( VaR[,,1],  
                 col = c(colors[1],colors[1]), xlab = "Bank ID", ylab=" Financial system losses, VaR",
                 cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, beside=TRUE, names.arg = 1:banks,
                 border = NA,cex.names = PPar$cex.names,ylim=ylim,density=c(-d,d),xaxt="n",yaxt="n")
  axis(1,at=colMeans(mp),label=seq(banks),cex.axis=PPar$cex_ax,gap.axis=0.,tick=F)
  axis(2,at=yat,label=ylabel,cex.axis=PPar$cex_ax,gap.axis=0.,tick=T)
  barplot( VaR[,,1], col = c(colors[1],colors[1]),add=T, beside=T,border=NA,density=c(-d,d),yaxt="n",xaxt="n")
  barplot( VaR[,,2], col = c(colors[2],colors[2]),add=T, beside=T,border=NA,density=c(-d,d),yaxt="n",xaxt="n")
  barplot( VaR[,,3], col = c(colors[3],colors[3]),add=T, beside=T,border=NA,density=c(-d,d),yaxt="n",xaxt="n")
  abline(h=seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/6),col="lightgrey",lwd=0.5,lty=2)
  legend("topright",title="",legend=c("DI","", expression("IB"^"WO"),"DI","SC",expression("IB"^"W"))
         ,fill=c(colors[3],0,colors[1],colors[c(3,2,1)]),density=c(-d,-d,-d,d,d,d),border = c("black","#00000000","black","black","black","black")
         ,cex=PPar$cex_ax,bty="n",ncol=2,inset = c(0.,-0.1))
  dev.off()
}

Figure4 <- function(){
  Risks <- ScenariowiseLosses()
  PlotWidth <- 8
  WRatio <- 1.
  TextCEX <- 1.5
  pdf("./pics/Figure4.pdf",width = PlotWidth,height = 2*PlotWidth/(1+WRatio),pointsize=10)
  layout(matrix(c(1,3,2,3),nrow=2),widths=c(1,WRatio),heights=c(1,1))
  par(mar=c(3.5,4,0.5,0.5),cex.axis=TextCEX,cex.lab=TextCEX,cex.main=TextCEX,cex.sub=TextCEX)
  PlotInterbankScatter(Risks$MIB_WO_risks,Risks$MIB_W_risks,xlabel="",ylabel="")
  title(xlab=expression("IB"^"WO"),line=2.5)
  title(ylab=expression("IB"^"W"),line=2)
  text(0,0.3,adj=c(0,3),labels="a)",cex=TextCEX)
  TestCEX <- 1.25
  par(mar=c(3.5,4,0.5,0.5),cex.axis=TextCEX,cex.lab=TextCEX,cex.main=TextCEX,cex.sub=TextCEX)
  IBAmplifierBoxPlot()
  dev.off()
}

Figure4aAggregated <- function(){
  Folder=""
  Trial1="CovidShock"
  Trial2="CovidShock_Direct"
  Target="Figure4aAggregated.pdf"
  ### Pure IB Amplifier vs. IB Amplifier
  PPar <- GetPlotParams()
  Folder <- paste0("./data/",Folder)
  
  Risks <- ScenariowiseLosses()
  Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks+Risks$MIB_W_risks)
  OG_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MSC_W_risks)
  
  # Data ingestion for direct interbank contagion ####
  
  Direct_Cascade_risks <- CalculateRiskMeasures(Risks$Direct_risks+Risks$MIB_WO_risks)
  Direct_OG_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  Direct_risks <- CalculateRiskMeasures(Risks$Direct_risks)
  
  SC_amp <- (Cascade_risks-OG_risks) #IB amplifier with previous SC
  IB_amp <- (Direct_Cascade_risks-Direct_OG_risks) #IB Amp without previous SC
  
  SC_amp <- SC_amp
  IB_amp <- IB_amp
  TargetFolder <- "./pics/"
  
  EvilBanks <- which(colSums(Risks$MIB_W_risks)==0)
  
  ELFit <- lm(SC_amp[1,]~IB_amp[1,])
  VaRFit <- lm(SC_amp[2,]~IB_amp[2,])
  ESFit <- lm(SC_amp[3,]~IB_amp[3,])
  pdf(file = paste0(TargetFolder, Target),  width = 12, height = 5.5,# horizontal = FALSE,
      paper = "special", bg="white")
  
  par(mar = c(4.1,5.4,1,1)+0.1) #c(bottom, left, top, right)
  par(mfrow = c(1,3))
  ylim = c(0,ceiling(max(SC_amp)*10)/10)
  
  plot( t(rbind(IB_amp[1,],SC_amp[1,])),  
        xlab = "", ylab=expression("MIB"^"W"),
        cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, ylim=ylim,xlim=ylim,cex=4)
  text( t(rbind(IB_amp[1,],SC_amp[1,])),labels=(c(1:4))[],cex=PPar$cex.names)
  lines(ylim,ylim,type="l",lty="dotted") #The diagonal
  lines(c(0,ylim[2]),c(ELFit$coefficients[1],(ELFit$coefficients[1]+ylim[2]*ELFit$coefficients[2])),type="l")
  legend("topleft",title="Expected Loss",legend="",bty="n",cex=PPar$cex_ax)
  text(ylim[2]/4,0,adj=c(0,0),labels=glue("y = {format(ELFit$coefficients[1],digits=3)} + {format(ELFit$coefficients[2],digits=3)} * x \nR\U00B2 = {format(summary(ELFit)$r.squared,digits=4)}"),cex=PPar$cex_ax)
  
  par(mar = c(4.1,0.65,1,1)+0.1)
  plot( t(rbind(IB_amp[2,],SC_amp[2,])),  
        ylab = "", xlab=expression("MIB"^"WO"),
        cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, ylim=ylim,xlim=ylim,yaxt="n",cex=4)
  text( t(rbind(IB_amp[2,],SC_amp[2,])),labels=(c(1:4)),cex=PPar$cex.names)
  lines(ylim,ylim,type="l",lty="dotted")
  lines(c(0,ylim[2]),c(VaRFit$coefficients[1],(VaRFit$coefficients[1]+ylim[2]*VaRFit$coefficients[2])),type="l")
  legend("topleft",title="Value at Risk",legend="",bty="n",cex=PPar$cex_ax)
  text(ylim[2]/4,0,adj=c(0,0),labels=glue("y = {format(VaRFit$coefficients[1],digits=3)} + {format(VaRFit$coefficients[2],digits=3)} * x \nR\U00B2 = {format(summary(VaRFit)$r.squared,digits=4)}"),cex=PPar$cex_ax)
  par(mar = c(4.1,0.65,1,1)+0.1)
  
  plot( t(rbind(IB_amp[3,],SC_amp[3,])),  
        xlab = "", ylab="",
        cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax, ylim=ylim,xlim=ylim,yaxt="n",cex=4)
  text( t(rbind(IB_amp[3,],SC_amp[3,])),labels=(c(1:4))[-EvilBanks],cex=PPar$cex.names)
  lines(ylim,ylim,type="l",lty="dotted")
  lines(c(0,ylim[2]),c(ESFit$coefficients[1],(ESFit$coefficients[1]+ylim[2]*ESFit$coefficients[2])),type="l")
  legend("topleft",title="Expected Shortfall",legend="",bty="n",cex=PPar$cex_ax)
  text(ylim[2]/4,0,adj=c(0,0),labels=glue("y = {format(ESFit$coefficients[1],digits=3)} + {format(ESFit$coefficients[2],digits=3)} * x \nR\U00B2 = {format(summary(ESFit)$r.squared,digits=4)}"),cex=PPar$cex_ax)
  dev.off()
}

Figure4aLogged <- function(){
  Risks <- ScenariowiseLosses()
  PlotWidth <- 8
  WRatio <- 1-9/16
  TextCEX <- 1.2
  jpeg("./pics/Figure4_SI_Log.jpeg",width = PlotWidth,height = PlotWidth/(1+WRatio),units="in",res=300)
  par(mar=c(5,4,1,1)+0.1)
  par(cex.axis=TextCEX,cex.lab=TextCEX)#,cex.main=TextCEX,cex.sub=TextCEX)
  PlotInterbankScatter(Risks$MIB_WO_risks,Risks$MIB_W_risks,xlabel="",ylabel="",logplot=T,legendoffset=c(0,0.125))
  title(xlab=expression("IB"^"WO"),line=3)
  title(ylab=expression("IB"^"W"),line=2.5)
  dev.off()
}

Figure4bc_SI <- function(Folder="", Trial1="CovidShock",Trial2="CovidShock_Direct",Target="ScenariowiseComparison"){
  Folder <- paste0("./data/",Folder)
  Equities <- readRDS(paste0(Folder,"Equities_",Trial1,".rds"))
  
  #### B(L(S(phi)))-L(S(Phi)) vs. B(L(phi))-L(Phi)
  
  IB_risks_with <- readRDS(paste0(Folder,"TotalLosses_With_Cascade_",Trial1,".rds"))
  IB_risks_without <- readRDS(paste0(Folder,"TotalLosses_With_Cascade_",Trial2,".rds"))
  
  IB_risks_with <- t(t(IB_risks_with)/Equities)
  IB_risks_without <- t(t(IB_risks_without)/Equities)
  EvilBanks <- which(colSums(IB_risks_with)==0)
  
  TextCEX <- 1.2
  
  bla <- matrix(0,ncol=1000,nrow=length(Equities))
  for(i in 1:1000){
    bla[,i] <- IB_risks_with[i,]/IB_risks_without[i,]
  }
  pdf("./pics/Figure4b_SI.pdf",width=8,height = 4.5)
  par(mar=c(5,4,1,1)+0.1)
  par(cex.axis=TextCEX,cex.lab=TextCEX)#,cex.main=TextCEX,cex.sub=TextCEX)
  boxplot(bla[,1:100],ylim=c(1,5),xlab="",ylab="")
  title(ylab=expression("MIB"^"W"*" / MIB"^"WO"),line=2.5)
  title(xlab = "Scenario",line=3)
  dev.off()
  
  pdf("./pics/Figure4c_SI.pdf",width=8,height = 4.5)
  par(mar=c(5,4,1,1)+0.1)
  par(cex.axis=TextCEX,cex.lab=TextCEX)
  plot(c(1:1000),bla[1,1:1000],pch=20,xlab="",ylab="",col="black")
  points(c(1:1000),bla[4,1:1000],pch=20,col="purple")
  abline(h=mean(bla[4,1:1000]),col="purple")
  abline(h=mean(bla[1,1:1000]),col="black")
  text(500,axTicks(2)[5],adj=c(0,1),label=paste0("Means: ",format(mean(bla[1,1:1000]),digits=3),", ",format(mean(bla[4,1:1000]),digits=3)))
  title(ylab=expression("MIB"^"W"*" / MIB"^"WO"),line=2.5)
  title(xlab = "Scenario",line=3)
  legend("topleft",legend=c("Largest bank","Second largest bank"),col=c("black","purple"),pch=c(20,20),bty="n")
  dev.off()
}

SI_Figure_BankingShock <- function(){
  Risks <- ScenariowiseLosses()
  PlotWidth <- 8
  WRatio <- 1-9/16
  TextCEX <- 1.2
  jpeg("./pics/SI_Bankingshock.jpeg",width = PlotWidth,height = PlotWidth/(1+WRatio),units="in",res=100)
  par(mar=c(5,4,1,1)+0.1)
  par(cex.axis=TextCEX,cex.lab=TextCEX)#,cex.main=TextCEX,cex.sub=TextCEX)
  PlotInterbankScatter(Risks$Direct_risks,Risks$Direct_risks+Risks$MSC_W_risks,xlabel="",ylabel="",logplot=F,legendoffset=c(0,0.125),textoffset=c(0.1,0))
  title(xlab=expression("Direct"),line=3)
  title(ylab=expression("Direct + MSC"^"W"),line=2.5)
  dev.off()
}

SI_Figure_BankingShock_Logged <- function(){
  Risks <- ScenariowiseLosses()
  PlotWidth <- 8
  WRatio <- 1-9/16
  TextCEX <- 1.2
  jpeg("./pics/SI_Bankingshock_logged.jpeg",width = PlotWidth,height = PlotWidth/(1+WRatio),units="in",res=100)
  par(mar=c(5,4,1,1)+0.1)
  par(cex.axis=TextCEX,cex.lab=TextCEX)#,cex.main=TextCEX,cex.sub=TextCEX)
  PlotInterbankScatter(Risks$Direct_risks,Risks$Direct_risks+Risks$MSC_W_risks,xlabel="",ylabel="",logplot=T,legendoffset=c(0,0.125))
  title(xlab=expression("Direct"),line=3)
  title(ylab=expression("Direct + MSC"^"W"),line=2.5)
  dev.off()
}

Figure4_AmpCCDF <- function(){
  Risks <- ScenariowiseLosses()
  MIB_W_Banks <- list()
  MIB_W_Scenarios <- list()
  MIB_WO_Banks <- list()
  MIB_WO_Scenarios <- list()
  Amplifiers_Scenarios <- list()
  Amplifiers_Banks <- list()
  EvilBanks <- which(colSums(Risks$MIB_WO_risks)==0)
  nb <- dim(Risks$Direct_risks)[2]
  c <- viridis(nb)
  for(i in (1:nb)[-EvilBanks]){
    MIB_W_Banks[[i]] <- Risks$MIB_W_risks[1:1000,i]
    MIB_WO_Banks[[i]] <- Risks$MIB_WO_risks[1:1000,i]
  }
  for(i in 1:1000){
    MIB_W_Scenarios[[i]] <- Risks$MIB_W_risks[i,(1:nb)[-EvilBanks]]
    MIB_WO_Scenarios[[i]] <- Risks$MIB_WO_risks[i,(1:nb)[-EvilBanks]]
  }
  for(i in 1:nb){
    Amplifiers_Banks[[i]] <- Risks$MIB_W_risks[1:1000,i]/Risks$MIB_WO_risks[1:1000,i]
  }
  for(i in 1:1000){
    Amplifiers_Scenarios[[i]] <- Risks$MIB_W_risks[i,(1:nb)[-EvilBanks]]/Risks$MIB_WO_risks[i,(1:nb)[-EvilBanks]]
    Amplifiers_Scenarios[[i]][is.nan(Amplifiers_Scenarios[[i]])] <- 1
  }
  GCDF <- ecdf(unlist(Amplifiers_Scenarios))
  from <- 5
  ylim <- c(10^-05,1)
  TextCEX <- 1.5
  pdf("./pics/Figure4_Amps.pdf",width=9,height=5,pointsize=10)
  par(cex.axis=TextCEX,cex.lab=TextCEX,cex.main=TextCEX,cex.sub=TextCEX)
  plot(function(v){1-GCDF(v)},pch=20,type="p",col="#482576ff",log="xy",ylim=ylim,from=from,ylab="1-CDF",xlab=expression("Amplification Factor IB"^"W"*" / IB"^"WO"))
  dev.off()
}

Figure4_FactorBreaking <- function(){
  #### Plot the Factor breaking figures, 4 in total. Random Uniform shocks, random shocks, random IB and all random.
  IB_Scatter_clean <- function(MIB_WO,MIB_W,logplot=F,title="",ylabel=expression("IB"["W"]), xlabel=expression("IB"["WO"]),legendoffset=c(0,0),textoffset=c(0,0),limits = c(0,max(MIB_W))){
    library(viridis)
    EvilBanks <- which(colSums(MIB_WO)==0)
    NB <- dim(MIB_WO)[2]-length(EvilBanks)
    nb <- dim(MIB_W)[2]
    c <- viridis(nb)
    
    plotlog = ""
    if(logplot){
      plotlog="xy"
    }
    
    plot(NULL,xlim=limits,ylim=limits,log=plotlog,
         ylab=ylabel, xlab=xlabel,main=title,xaxt="n",yaxt="n")
    for(i in (1:nb)[! (1:nb) %in% EvilBanks]){
      points(MIB_WO[,i],MIB_W[,i],col=c[i],pch=20)
    }
    legend("topright",title="Bank ID", legend=seq(nb)[! (1:nb) %in% EvilBanks],
           col=c[! (1:nb) %in% EvilBanks],pch=rep(20,nb)[! (1:nb) %in% EvilBanks],
           y.intersp=rep(0.8,nb),bty="n",cex=par("cex.lab"),inset = legendoffset,ncol=2)
    
    length <- NB*dim(MIB_WO)[1]
    if(!logplot){
      BigFit <- lm((MIB_W[,! (1:nb) %in% EvilBanks][1:length]~MIB_WO[,! (1:nb) %in% EvilBanks][1:length]))
    }else{
      xvals = log(MIB_WO[,! (1:nb) %in% EvilBanks][1:length])
      yvals = log(MIB_W[,! (1:nb) %in% EvilBanks][1:length])
      BigFit <- lm(yvals[xvals>-Inf]~xvals[xvals>-Inf])
    }
    if(logplot){
      text(limits[1]+textoffset[1],limits[2]+textoffset[2],adj=c(0,1),cex = par("cex.lab"),labels=glue("y = {format(exp(BigFit$coefficients[1]),digits=3)}*x^{format(BigFit$coefficients[2],digits=3)} 
                                                             \n R\U00B2={format(summary(BigFit)$r.squared,digits=4)}")
      )
      curve(exp(BigFit$coefficients[1])*x^BigFit$coefficients[2],add=T,from=limits[1],to=limits[2])
    }else{
      text(0.05+textoffset[1],0+textoffset[2],adj=c(0,-1.5),offset=0,cex = par("cex.lab"),labels=glue("y = x*{format(BigFit$coefficients[2],digits=2)}+{format(BigFit$coefficients[1],digits=2)}"))
      text(0.05+textoffset[1],0+textoffset[2],adj=c(0,0.5),cex = par("cex.lab"),labels=glue("R\U00B2 = {format(summary(BigFit)$r.squared,digits=2)}"))
      abline(a=BigFit$coefficients[1],b=BigFit$coefficients[2])
    }
  }
  Dummy <- ScenariowiseLosses()
  OG <- list()
  OG$BigSim$IB <- Dummy$MIB_WO_risks
  OG$BigSim$SC <- Dummy$MIB_W_risks
  RS <- readRDS("./data/HeteroShocks.rds")
  RIB <- readRDS("./data/RandomSmallIB.rds")
  AllR <- readRDS("./data/RandomIBHetShocks.rds")
  limits = c(0,max(max(OG$BigSim$SC),max(RS$BigSim$SC),max(RIB$BigSim$SC),max(AllR$BigSim$SC)))
  
  jpeg("./pics/Figure_4_SI_Factor.jpeg",width=8,height = 8,units="in",res=300)
  layout(matrix(c(1,3,2,4),nrow = 2),widths = c(1.15,1),heights = c(1,1.15))
  PPar <- GetPlotParams()
  TextCEX = 1.5
  res = 0.5*10^-1
  tickloc = seq(limits[1],floor(limits[2]/res)*res,length.out=6)
  toffset = c(0.0,0)
  loffset = c(0,0.05*3)
  
  par(mar = c(0.5,4,1,1)+0.1,cex.axis=TextCEX,cex.lab=TextCEX,cex.main=TextCEX,cex.sub=TextCEX)
  IB_Scatter_clean(OG$BigSim$IB,OG$BigSim$SC,limits = limits,xlab="",ylab="",textoffset = toffset,legendoffset = loffset)
  axis(2,at=tickloc,labels=T,tick=T)
  axis(1,at=tickloc,label=F,tick=T,tcl=-0.25)
  title(ylab=expression("IB"^"W"),line=2)
  text(0,0.3,adj=c(0,3),labels="a)",cex=TextCEX)
  
  par(mar = c(0.5,0.5,1,1)+0.1)
  IB_Scatter_clean(RS$BigSim$IB,RS$BigSim$SC,limits = limits,xlab="",ylab="",textoffset = toffset,legendoffset = loffset)
  axis(2,at=tickloc,labels=F,tick=T,tcl = -0.25)
  axis(1,at=tickloc,label=F,tick=T,tcl=-0.25)
  text(0,0.3,adj=c(0,3),labels="b)",cex=TextCEX)
  
  par(mar = c(3.5,4,1,1)+0.1)
  IB_Scatter_clean(RIB$BigSim$IB,RIB$BigSim$SC,limits = limits,xlab="",ylab="",textoffset = toffset,legendoffset = loffset)
  axis(2,at=tickloc,labels=T,tick=T)
  axis(1,at=tickloc,label=T,tick=T,gap.axis = 0)
  title(xlab=expression("IB"^"WO"),line=2.5)
  title(ylab=expression("IB"^"W"),line=2)
  text(0,0.3,adj=c(0,3),labels="c)",cex=TextCEX)
  
  par(mar = c(3.5,0.5,1,1)+0.1)
  IB_Scatter_clean(AllR$BigSim$IB,AllR$BigSim$SC,limits = limits,xlab="",ylab="",textoffset = toffset,legendoffset = loffset)
  axis(2,at=tickloc,labels=F,tick=T,tcl=-0.25)
  axis(1,at=tickloc,label=T,tick=T)
  title(xlab=expression("IB"^"WO"),line=2.5)
  text(0,0.3,adj=c(0,3),labels="d)",cex=TextCEX)
  dev.off()
}

SystemHist <- function(){
  Risks <- ScenariowiseLosses()
  Folder <- "./data/"
  Trial <- "CovidShock"
  Equities <- readRDS(paste0(Folder, "Equities_",Trial,".rds"))
  DIRisks <- colSums(t(Risks$Direct_risks)*Equities/sum(Equities))
  SCRisks <- colSums(t(Risks$Direct_risks+Risks$MSC_W_risks)*Equities/sum(Equities))
  IBRisks <- colSums(t(Risks$Direct_risks+Risks$MSC_W_risks+Risks$MIB_W_risks)*Equities/sum(Equities))
  
  PPar <- GetPlotParams()
  TargetFolder <- "./pics/"
  Target <- "SystemLossHistogram"
  pdf(file = paste0(TargetFolder, Target,".pdf"),  width = 8, height = 4.5,# horizontal = FALSE,
      paper = "special", bg="white")
  par(mar = c(3.6,3.9,1,1)+0.1) #c(bottom, left, top, right)
  par(mgp=c(2.5,1,0))
  
  PPar$cex_ax <- 1.2
  PPar$cex.names <- 1.2
  PPar$cex_points <- 1.2
  PPar$cex_xylab <- 1.2
  cex_ax <- 1.2
  
  hist(0.0,xlim=c(0,0.12),ylim=c(0,100),lty="blank", xlab = "Financial system losses",ylab="Frequency",main="",col="white",
       cex.lab = PPar$cex_ax, cex.axis = PPar$cex_ax)
  colors = viridis(3,begin=0.1,end=0.9)
  lines(x = c(quantile(DIRisks,0.95),quantile(DIRisks,0.95)), y = c(0,100), lwd=2 , col = colors[3])
  lines(x = c(quantile(SCRisks,0.95),quantile(SCRisks,0.95)), y = c(0,100), lwd=2 , col = colors[2])
  lines(x = c(quantile(IBRisks,0.95),quantile(IBRisks,0.95)), y = c(0,100), lwd=2 , col = colors[1])
  lines(x = c(mean(DIRisks),mean(DIRisks)), y = c(0,100), lwd=2 , col = colors[3], lty = 2)
  lines(x = c(mean(SCRisks),mean(SCRisks)), y = c(0,100), lwd=2 , col = colors[2], lty = 2)
  lines(x = c(mean(IBRisks),mean(IBRisks)), y = c(0,100), lwd=2 , col = colors[1], lty = 2)
  lines(x = c(mean(DIRisks[DIRisks > quantile(DIRisks,0.95)]),mean(DIRisks[DIRisks > quantile(DIRisks,0.95)])), y = c(0,100), lwd=2 , col = colors[3], lty = 3)
  lines(x = c(mean(SCRisks[SCRisks > quantile(SCRisks,0.95)]),mean(SCRisks[SCRisks > quantile(SCRisks,0.95)])), y = c(0,100), lwd=2 , col = colors[2], lty = 3)
  lines(x = c(mean(IBRisks[IBRisks > quantile(IBRisks,0.95)]),mean(IBRisks[IBRisks > quantile(IBRisks,0.95)])), y = c(0,100), lwd=2 , col = colors[1], lty = 3)
  colors <- viridis(3,begin=0.1,end=0.9,alpha=0.5)
  hist(as.numeric(DIRisks),  xlim = c(0,0.12), ylim = c(0,100), add=T,
       breaks = seq(floor(min(DIRisks)*1000)/1000, max(DIRisks)+0.001, by = 0.001), col=colors[3],
       xlab = "financial system losses ", ylab = "frequency", cex.lab = cex_ax, cex.axis = cex_ax, main = " ",
  )
  hist(as.numeric(SCRisks),  xlim = c(0,0.12), ylim = c(0,100), add=TRUE,
       breaks = seq(floor(min(SCRisks)*1000)/1000, max(SCRisks)+0.001, by = 0.001), col=colors[2],
       xlab = "financial system losses ", ylab = "frequency", cex.lab = cex_ax, cex.axis = cex_ax, main = " ",
  )
  hist(as.numeric(IBRisks),  xlim = c(0,0.12), ylim = c(0,100), add=TRUE,
       breaks = seq(floor(min(IBRisks)*1000)/1000, max(IBRisks)+0.001, by = 0.001), col=colors[1],
       xlab = "financial system losses ", ylab = "frequency", cex.lab = cex_ax, cex.axis = cex_ax, main = " ",
  )
  colors = viridis(3,begin=0.1,end=0.9)
  legend("topright",c("DI","DI+SC","DI+SC+IB","EL","VaR","ES"),col=c(colors[c(3,2,1)],rep("black",3)),pch=c(15,15,15,NA,NA,NA),lty=c(NA,NA,NA,'dashed','solid','dotted'),lwd=c(NA,NA,NA,2,2,2),cex=PPar$cex_ax,bty="n")
  dev.off()
}