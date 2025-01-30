library("Matrix")

source(paste0(getwd(),"/Scripts/BankContagion.R"),local=T)
source(paste0(getwd(),"/Scripts/DebtRankContagion.R"),local=T)
source(paste0(getwd(),"/Scripts/SimulationSetup.R"),local=T)

ShockToState <- function(RelativeShock,State,Shockfunctions){
  ###Takes the shock, the current state and the Properties to shock as inputs and returns the state after the shock. Shock can be a single vector for each bank or a 2-dim matrix with one dimension as long as the number of Properties to shock.
  ShockedProperties <- names(Shockfunctions)
  i=NULL
  if(length(dim(RelativeShock))>1){i=which(dim(RelativeShock)==length(ShockedProperties))}
  if(is.null(i)){
    for(j in seq_along(ShockedProperties)){
      State <- Shockfunctions[[ShockedProperties[j]]](RelativeShock,State)
    }
  } else{
    IndexMatrix <- matrix(rep(seq_len(length(RelativeShock)/2)),ncol=2)
    for(j in seq_along(ShockedProperties)){
      IndexMatrix[,i] <- j
      Localshock <- RelativeShock[IndexMatrix]
      State <- Shockfunctions[[ShockedProperties[j]]](Localshock,State)
    }
  }
  return(State)
}

CreateSubsetArgumentList <- function(ArrayToSubset,DimensionToIterate=1,IterationID=1){
  output <- list(ArrayToSubset)
  Dims <- dim(ArrayToSubset)
  for(i in seq_along(Dims)){
    output[[i+1]] <- if(i==DimensionToIterate){IterationID}else{seq_len(Dims[i])}
  }
  return(output)
}

RunMultipleScenarios <- function(InitialState,Scenarios,ShockFunctions,type, Strategies,Ranks, path=paste0(getwd(),"/"),ScenariosIterationIndex = 1){
  NumberScenarios <- dim(Scenarios)[ScenariosIterationIndex]
  NumberBanks <- length(InitialState$Equity)
  FSRI <- Matrix(0,nrow = NumberScenarios, ncol = NumberBanks)
  FSRIPlus <- Matrix(0,nrow = NumberScenarios, ncol = NumberBanks)
  LCRs <- array(dim=c(NumberScenarios,NumberBanks,3))
  CARs <- array(dim=c(NumberScenarios,NumberBanks,3))
  Losses <- array(dim=c(NumberScenarios,NumberBanks,length(Strategies)))
  IterationCounter <- rep(0,NumberScenarios)
  ScenariosIterationArgument <- CreateSubsetArgumentList(Scenarios,DimensionToIterate = ScenariosIterationIndex)
  for(i in 1:NumberScenarios){
    ScenariosIterationArgument[[ScenariosIterationIndex+1]] <- i
    LocalScenario <- do.call("[",ScenariosIterationArgument)
    ShockedState <- ShockToState(LocalScenario,InitialState,ShockFunctions)
    ShockedStatePlus <- BankingContagionCleanish(ShockedState,Strategies,Ranks,Debugging=F)
    ShockedResultsPlus <- AnalyseSingleFSRICascade(ShockedStatePlus,ShockedState)
    IterationCounter[i] <- ShockedResultsPlus$Iters
    LCRs[i,,1] <- InitialState$Liquidity/InitialState$NetOutflow
    LCRs[i,,2] <- ShockedState$Liquidity/ShockedState$NetOutflow
    LCRs[i,,3] <- ShockedResultsPlus$LCR
    CARs[i,,1] <- InitialState$Equity/InitialState$RiskWeightedAssets
    CARs[i,,2] <- ShockedState$Equity/InitialState$RiskWeightedAssets
    CARs[i,,3] <- ShockedResultsPlus$Equity/(InitialState$RiskWeightedAssets+colSums(ShockedResultsPlus$Assets*ShockedResultsPlus$AssetPrices-InitialState$Assets*InitialState$AssetPrices))
    FSRI[i,] <- InitialState$Equity - ShockedState$Equity
    FSRIPlus[i,] <- ShockedResultsPlus$Losses
    Losses[i,,] <- as.matrix(ShockedResultsPlus$StrategyLosses)
  }
  saveRDS(IterationCounter,paste0(path,"data/IterationCounts",type,".rds"))
  saveRDS(CARs,paste0(path,"data/CARs_",type,".rds"))
  saveRDS(LCRs,paste0(path,"data/LCRs_",type,".rds"))
  saveRDS(InitialState$Equity,paste0(path,"data/Equities_",type,".rds"))
  saveRDS(FSRI,paste0(path,"data/TotalLosses_No_Cascade_",type,".rds"))
  saveRDS(FSRIPlus,paste0(path,"data/TotalLosses_With_Cascade_",type,".rds"))
  saveRDS(Losses,paste0(path,"data/StrategyLosses_",type,".rds"))
  FSRI <- apply(FSRI,FUN="sum",MARGIN=1)
  FSRIPlus <- apply(FSRIPlus,FUN="sum",MARGIN=1)
  TotalEquity <- sum(InitialState$Equity)
  FSRI <- FSRI/TotalEquity
  FSRIPlus <- FSRIPlus/TotalEquity
  output <- cbind(FSRI,FSRIPlus)
  saveRDS(output, paste0(path,"data/FSRIs_",type,".rds"))
  return(cbind(FSRI,FSRIPlus))
}

AnalyseSingleFSRICascade <- function(ResultState,InitialState){
  ### This analyzes the difference between the IntialState and the ResultsState
  output = list()
  output$Losses <- colSums(InitialState$AssetPrices * InitialState$Assets - ResultState$AssetPrices * ResultState$Assets)
  output$FireSaleLosses <- colSums((InitialState$AssetPrices-ResultState$AssetPrices)*InitialState$Assets)
  output$StartingAssets <- InitialState$Assets
  output$Assets <- ResultState$Assets
  output$AssetPrices <- ResultState$AssetPrices
  output$StrategyLosses <- ResultState$Losses
  output$LGD <- apply(ResultState$Losses, MARGIN = 1, FUN = "max")
  output$Defaults <- ResultState$Defaults
  LastAssetIndex <- length(InitialState$AssetPrices)-length(InitialState$InterbankMarkets)
  LiquidChange <- colSums(ResultState$AssetPrices[-1]*ResultState$Assets[-1,]*(1-ResultState$Haircut[-1])) 
  + colSums(ResultState$AssetPrices[2:LastAssetIndex]*(InitialState$Assets[2:LastAssetIndex,] - ResultState$Assets[2:LastAssetIndex,])) 
  - colSums(InitialState$Assets[-1,]*(1-ResultState$Haircut[-1]))
  output$Liquidity <- ResultState$Liquidity+LiquidChange
  output$NetOutflow <- ResultState$NetOutflow
  output$LCR <- output$Liquidity/output$NetOutflow
  InterbankIndex = length(ResultState$AssetPrices)-length(ResultState$InterbankMarkets)+1
  InterbankLiquidityChanges <- ResultState$Assets[InterbankIndex,]-InitialState$Assets[InterbankIndex,]
  NostroChanges <- ResultState$Assets[11,]-InitialState$Assets[11,]
  output$Inflow <- pmax(ResultState$Inflow+NostroChanges+InterbankLiquidityChanges,0)
  output$Outflow <- ResultState$Outflow
  output$StartingEquity <- InitialState$Equity
  output$Equity <- InitialState$Equity - output$Losses
  output$Interbank1Change <- (InitialState$Assets[LastAssetIndex+1,]-ResultState$Assets[LastAssetIndex+1,])
  output$Iters <- ResultState$Iters
  return(output)
}

FSRIPlus <- function(Name,InputData,Scenarios,Model="DebtRank",
                     ShockedProperties = c("Equity"), RelativeShocks=c(T),
                     path=paste0(getwd(),"/"),
                     overwrite = F){
  
  EqPath <- paste0(path,"data/Equities_",Name,".rds")
  EqExistence <- file.exists(EqPath)&&!dir.exists(EqPath)
  if(!overwrite && EqExistence){
    stop("The identifier ",Name," seems to exist already. Since overwrite is set to FALSE execution is terminated to avoid unexpected loss of data")
  }
  
  InitialState <- InputData$State
  Strategies <- StrategySetup(Model)
  Ranks <- RankSetup(InputData,Model)
  Shockfunctions <- ShockFunctionSetup(ShockedProperties,RelativeShocks)

  FSRIs <- RunMultipleScenarios(InitialState,Scenarios,ShockFunctions=Shockfunctions,Name,Strategies=Strategies,Ranks=Ranks,path = path)
  return(FSRIs)
}


#Scenarios <- FSRI_loss_mat
#Scenarios <- fsri_loss_multiple_weighted # Should be a matrix of size Scenarios * banks (19)
#Scenarios <- readRDS("./data/IntermediateData/FSRI_250923.rds")
#Scenarios <- Scenarios$results_either$losses_over_equity_fb_mat
#Scenarios <- Scenarios[,MappingOrder]
#Scenarios <- pmin(Scenarios,1)
#Scenarios <- TopXESRIFirms(GetBankMapping(),top=2500,All=T)
#Scenarios <- DirectCovidShockScenario()
#Scenarios <- CovidShockScenario()

Model = "DebtRank" #Andras or DebtRank, chooses the Contagion model
Name = "NewCovidShock_EquityTermination_11012024" #Unique identifier for the run. Appears in the filenames and can potentially overwrite old stuff
ShockedProperties = c("Equity") # Options: Equity, Liquidity, Inflow
RelativeShocks = c(F) # Vector of bools, true denotes that the shock is relative to the current value, false is an absolute shock.
Path <- paste0(getwd(),"/")

#FSRIs <- FSRIPlus(Name, Scenarios, Model = Model, ShockedProperties = ShockedProperties, 
#                 RelativeShocks = RelativeShocks, 
#                 path = Path, overwrite = F)
