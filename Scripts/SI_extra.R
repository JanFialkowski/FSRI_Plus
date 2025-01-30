source("./Scripts/DebtRank.R",local=T)
source("./Scripts/FSRIPlusPipeline.R",local=T)

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

OldDebtRanks <- function(IBmat,Equities,ShockIB,ShockSC){
  IB <- CollectiveDebtRanking(t(IBmat*Equities),Equities,ShockIB)
  SC <- CollectiveDebtRanking(t(IBmat*Equities),Equities,ShockSC)
  IB$DebtRank1 <- IB$DebtRank1 - IB$Original
  IB$DebtRank2 <- IB$DebtRank2 - IB$Original
  SC$DebtRank1 <- SC$DebtRank1 - SC$Original
  SC$DebtRank2 <- SC$DebtRank2 - SC$Original
  return(list(IB=IB,SC=SC))
}

SolvingIB <- function(IBmat,Equities,ShockIB,ShockSC){
  IB_solved <- t(solve(dim(IBmat)[1]-IBmat)%*%t(ShockIB))
  SC_solved <- t(solve(dim(IBmat)[1]-IBmat)%*%t(ShockSC))
  return(list(IB=IB_solved-ShockIB,SC=SC_solved-ShockSC))
}

BigSimulation <- function(IBmat,Equities,ShockIB,ShockSC,Namemod=""){
  Data <- LoadRandomData(Banks=4,Assets=15)
  Data$State$InterbankMarkets[[1]] <- IBmat*Equities
  Data$State$Equity <- Equities
  Data$State$StartingEquity <- Equities
  Data$State$Assets[15,] <- rowSums(IBmat*Equities)
  Data$State$StartingAssets <- Data$State$Data$Assets
  path <- paste0(getwd(),"/")
  Model = "DebtRank" #Andras or DebtRank, chooses the Contagion model
  type = paste0("FactorBreaking_IB",Namemod) #Unique identifier for the run. Appears in the filenames and can potentially overwrite old stuff
  ShockedProperties = c("Equity") # Options: Equity, Liquidity, Inflow
  RelativeShocks = c(T) # Vector of bools, true denotes that the shock is relative to the current value, false is an absolute shock.
  Scenarios <- ShockIB
  
  FSRIs <- FSRIPlus(type, Scenarios, Model = Model, ShockedProperties = ShockedProperties, 
                    RelativeShocks = RelativeShocks,
                    path = path, overwrite = T,InputData=Data)
  IB_Losses <- readRDS(paste0(path,"data/TotalLosses_With_Cascade_",type,".rds"))
  IB_Losses <- t(t(IB_Losses)/Equities)
  
  type=paste0("FactorBreaking_SC",Namemod)
  Scenarios <- ShockSC
  
  FSRIs <- FSRIPlus(type, Scenarios, Model = Model, ShockedProperties = ShockedProperties, 
                    RelativeShocks = RelativeShocks,
                    path = path, overwrite = T,InputData=Data)
  
  SC_Losses <- readRDS(paste0(path,"data/TotalLosses_With_Cascade_",type,".rds"))
  SC_Losses <- t(t(SC_Losses)/Equities)
  IB_Losses_Clamp <- pmax(pmin(IB_Losses+ShockIB,1)-ShockIB,0)
  SC_Losses_Clamp <- pmax(pmin(SC_Losses+ShockSC,1)-ShockSC,0)
  return(list(IB=IB_Losses,SC=SC_Losses,IB_Clamp = IB_Losses_Clamp, SC_Clamp = SC_Losses_Clamp))
}

SingleStep <- function(IBmat,Equities,ShockIB,ShockSC){
  IB <- t(IBmat %*% t(ShockIB))
  SC <- t(IBmat %*% t(ShockSC))
  return(list(IB=IB,SC=SC))
}

RunAll <- function(IBmat,Equities,ShockIB,ShockSC,Namemod=""){
  output = list()
  output$BigSim = BigSimulation(IBmat,Equities,ShockIB,ShockSC,Namemod)
  output$Solve = SolvingIB(IBmat,Equities,ShockIB,ShockSC)
  output$OldDRs = OldDebtRanks(IBmat,Equities,ShockIB,ShockSC)
  output$SingleStep = SingleStep(IBmat,Equities,ShockIB,ShockSC)
  return(output)
}

Risks <- ScenariowiseLosses()
Equities <- readRDS("./data/Equities_CovidShock.rds")

BankShock_IB <- Risks$Direct_risks
BankShock_SC <- Risks$Direct_risks+Risks$MSC_W_risks

# Random IBmat, but smaller ####
RIBs <- matrix(runif(t*t),nrow=t)/20
diag(RIBs) <- 0
RIBDatas <- RunAll(RIBs,Equities,BankShock_IB,BankShock_SC,Namemod="Random_smallIB")
saveRDS(RIBDatas,"./data/RandomSmallIB.rds")

# Less uniform random bankshocks ####
means_IB <- apply(BankShock_IB,FUN = mean,MARGIN = 2)
sd_IB <- apply(BankShock_IB,FUN = sd,MARGIN = 2)
means_SC <- apply(BankShock_SC,FUN = mean,MARGIN = 2)
sd_SC <- apply(BankShock_SC,FUN = sd,MARGIN = 2)
HeteroRIB <- rnorm(1000,mean=means_IB[1],sd=sd_IB[1])
HeteroRSC <- rnorm(1000,mean=means_SC[1],sd=sd_SC[1])
for(i in 2:4){
  HeteroRIB <- cbind(HeteroRIB,rnorm(1000,mean=means_IB[i],sd=sd_IB[i]))
  HeteroRSC <- cbind(HeteroRSC,rnorm(1000,mean=means_SC[i],sd=sd_SC[i]))
}
HeteroRIB <- pmin(pmax(HeteroRIB,0),1)
HeteroRSC <- pmin(pmax(HeteroRIB,HeteroRSC),1)
HeteroRData <- RunAll(IBmat,Equities,HeteroRIB,HeteroRSC,Namemod="HeteroRandomShocks")
saveRDS(HeteroRData,"./data/HeteroShocks.rds")

# Hetero Schocks, random IB ####
HetIBData <- RunAll(RIBs,Equities,HeteroRIB,HeteroRSC,Namemod="HeteroIB")
saveRDS(HetIBData,"./data/RandomIBHetShocks.rds")