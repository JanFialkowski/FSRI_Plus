library(Matrix)
# This code is a framework around the model from this paper: Borsos, Andras, and Bence Mero. Shock propagation in the banking system with real economy feedback. No. 2020/6. MNB Working Papers, 2020.
# By changing the adjustment strategies this can be used to also model DebtRank

# Data needed:
# Assets of each bank at current time, matrix "Assets"
# Prices of each Asset at current time,  vector "AssetPrices"
# Risk weights of each asset, vector "RiskWeights"
# Equity of Banks at start, vector "Equity"
# Risk weighted assets of each bank at start, vector "RiskWeightedAssets"
# CAR regulatory goal, scalar "CAR"
# Haircut parameter for each asset, vector "Haircut"
# Interbank Market matrices, as a list containing matrices, "InterbankMarkets" 
# Liquidity goal from regulators, scalar, "LCR"
# Expected Income of each bank, Vector, "Inflow"
# Expected Expenditures of each bank, Vector, "Outflow"
# Minimum Prices for assets if the entire banking sector liquidates it, "MinimumPrices"
# Number of Assets held outside of the banking sector, "Outside Assets"
# Liquidity of a bank "Liquidity"

DifferenceToInitialStateDirectly <- function(Assets,Values,StartingAssets,StartingValues){
  return(colSums(Assets*Values-StartingAssets*StartingValues)) #Observed quantities at time t=t - Observed Quantities at time t=0
}

DifferenceToInitialState <- function(State,InitialState,PriceModifier){
  return(DifferenceToInitialStateDirectly(State$Assets,State$AssetPrices*PriceModifier,InitialState$Assets,InitialState$AssetPrices*PriceModifier))
}

SellAssetsOfGivenRank <- function(Assets, Values, Order, CurrentRank, Targets) {
  #For a given Rank of Asset Liquidation this function reduces the assets of banks until Assets*Price=Targets, or Assets are 0. Returns a list of the gained Value and the Change in Assets
  ValueInAssetsAtCurrentRank <- colSums(Assets*Values*(Order==CurrentRank)) #Vector of length = banks, Amount of money for the banks that can be gained by selling all Assets in the current rank
  FractionOfAssetsToSell <- pmax(pmin((Targets)/ValueInAssetsAtCurrentRank,1,na.rm=TRUE),0) #Vector of Length = banks, Clamp the fraction of assets that are sold to something sensible. Can't sell more than you have, and you don't buy anything in our model
  ChangeInAssetsThisRank <- t(t(Assets)*FractionOfAssetsToSell)*(Order == CurrentRank)# R does not Broadcast, just recycles, so the t(t()) is needed to fix that, in other words: FractionToSell has to be applied rowwise, thus the t(Assets), the second one lets us go back to the original convention 
  CashGainedBySellingAssetsThisRank <- colSums(ChangeInAssetsThisRank*Values)
  output <- list()
  output$ValueGained <- CashGainedBySellingAssetsThisRank
  output$ChangeInAssets <- ChangeInAssetsThisRank
  return(output)
}

SellAssetsUntilTargetReached <- function(Assets,Values,Targets,Order){
  # This function sells of Assets of each bank and gains its Values in the order given by Ranks until the Targets are reached. Assets at the same rank are sold in the same proportion
  RanksPossiblyUsed <- sort(unique(Order)) #We iterate over each rank until either all banks have reached their targets or we run out of assets
  CashGainedBySellingAssets <- rep(0,length(Targets)) #Keep track of the amount of money each bank gains by selling their assets
  if(RanksPossiblyUsed[1]==0){RanksPossiblyUsed <- RanksPossiblyUsed[2:length(RanksPossiblyUsed)]} #Assets with rank 0 are not sold (e.g. risk free accounts for CAR)
  ChangeInAssets <- matrix(0,nrow=dim(Assets)[1],ncol=dim(Assets)[2]) #Keeps Track of the sold Assets
  
  for (CurrentRank in RanksPossiblyUsed){
    if (sum(CashGainedBySellingAssets<Targets)!=0){#Check if any bank with unresolved targets still exists
      dummy <- SellAssetsOfGivenRank(Assets,Values,Order,CurrentRank,Targets - CashGainedBySellingAssets)
      ChangeInAssets <- ChangeInAssets + dummy$ChangeInAssets
      CashGainedBySellingAssets <- CashGainedBySellingAssets + dummy$ValueGained
    }
  }
  output <- list()
  output$ChangeInAssets <- ChangeInAssets
  output$CashGainedBySellingAssets <- CashGainedBySellingAssets
  return(output)
}

CalculateCARTargets <- function(State,InitialState){
  RWA <- InitialState$RiskWeightedAssets + DifferenceToInitialState(State,InitialState,State$RiskWeights)
  return(State$CAR*RWA - State$Equity)
}

CARAdjustment <- function(State,InitialState,Order){
  SuperfluousRisk <- 1/State$CAR * CalculateCARTargets(State,InitialState)
  Results <- SellAssetsUntilTargetReached(State$Assets,State$AssetPrices*State$RiskWeights,SuperfluousRisk,Order)
  return(-Results$ChangeInAssets)# The minus keeps the convention intact where A(t)=A(t-1)+Change. We only sell, so changes are negative
}

CARDefaults <- function(State,InitialState,Results){
  RWA <- InitialState$RiskWeightedAssets + DifferenceToInitialStateDirectly(State$Assets+Results,State$AssetPrices*State$RiskWeights, InitialState$Assets,InitialState$AssetPrices*InitialState$RiskWeights)
  Cars <- State$Equity/RWA
  if(NA %in% Cars){stop("One or more CAR values are NA, usually doesn't happen, so check if something broke")}
  Defaults <- Cars < State$CAR*0.999999
  Losses <- Cars < State$CAR/2*0.999999
  output <- list(Defaults = Defaults, Losses = Losses)
  return(output)
}

CalculateCurrentLiquidity <- function(State,InitialState){
  # Initial Liquidity - observed Liquidity at time=0 + observed Liquidity at time=t
  return(InitialState$Liquidity + DifferenceToInitialState(State,InitialState,1-State$Haircut))
}

CalculateLiquidityTargets <- function(LCR,NetOutflow,Liquidity){
  return(LCR*NetOutflow-Liquidity)
}

LiquidityAdjustment <- function(State,InitialState,Ranks){
  CurrentLiquidity <- CalculateCurrentLiquidity(State,InitialState)
  MissingLiquidity <- CalculateLiquidityTargets(State$LCR,State$NetOutflow,CurrentLiquidity)
  Results <- SellAssetsUntilTargetReached(State$Assets,State$AssetPrices*State$Haircut,MissingLiquidity,Ranks)
  return(-Results$ChangeInAssets)
}

LiquidityDefaults <- function(State,InitialState,Results){
  CurrentLiquidity <- CalculateCurrentLiquidity(State,InitialState)
  Defaults <- CurrentLiquidity/State$NetOutflow < State$LCR
  Losses <- CurrentLiquidity/State$NetOutflow < State$LCR/2
  output <- list(Defaults = Defaults, Losses = Losses)
}

CalculateRelativeLGD <- function(MatrixOfLosses){
  MaximumOfLosses <- apply(MatrixOfLosses,1,max) #Columns have to be the individual banks and rows the different scenarios
  Remaining <- MaximumOfLosses
}

BuildNetOutflowCalculation <- function(NostroIndex){
  function(State,InitialState){
    InterbankIndex = length(State$AssetPrices)-length(State$InterbankMarkets)+1
    InterbankLiquidityChanges <- State$Assets[InterbankIndex,]-InitialState$Assets[InterbankIndex,]
    NostroChanges <- State$Assets[NostroIndex,]-InitialState$Assets[NostroIndex,]
    NewNetOutflow <- State$Outflow - pmin(pmax(State$Inflow+NostroChanges+InterbankLiquidityChanges,0),0.75*State$Outflow)
    return(NewNetOutflow)
  }
}

CalculateNewNetOutflow <- function(State,InitialState, NostroIndex=11){
  InterbankIndex = length(State$AssetPrices)-length(State$InterbankMarkets)+1
  InterbankLiquidityChanges <- State$Assets[InterbankIndex,]-InitialState$Assets[InterbankIndex,]
  NostroChanges <- State$Assets[NostroIndex,]-InitialState$Assets[NostroIndex,]
  NewNetOutflow <- State$Outflow - pmin(pmax(State$Inflow+NostroChanges+InterbankLiquidityChanges,0),0.75*State$Outflow)
  return(NewNetOutflow)
}

CalculateNewPrices <- function(State,InitialState){
  Alphas <- log(InitialState$MinimumPrices)/rowSums(InitialState$Assets + InitialState$OutsideAssets)
  Alphas[is.na(Alphas)] <- 0
  return(exp(Alphas*(State$OutsideAssets + rowSums(InitialState$Assets-State$Assets))))
}

CalculateCashGain <- function(State,InitialState,AssetChange){
  N_Interbank <- length(State$InterbankMarkets)
  N_Assets <- length(State$AssetPrices)
  LastAsset <- N_Assets - N_Interbank
  return(colSums(State$AssetPrices[2:LastAsset]*(InitialState$Assets[2:LastAsset,]-State$Assets[2:LastAsset,])))
}

CheckIfAssetChangeSmall <- function(State,AssetChange){
  sum(colSums(AssetChange))/sum(State$Equity)<=0.01 #relative change of equity
}

CheckLiquiditySituation <- function(State,InitialState,Ranks){
  Ranks <- Ranks[[length(Ranks)]]
  LastAssetIndex <- length(State$AssetPrices)-length(State$InterbankMarkets)
  LiquidChange <- colSums(State$AssetPrices[-1]*State$Assets[-1,]*(1-State$Haircut[-1])) #Liquidity still in Assets
     + colSums(State$AssetPrices[2:LastAssetIndex]*(InitialState$Assets[2:LastAssetIndex,] - State$Assets[2:LastAssetIndex,])) #Cash gained by selling
     - colSums(InitialState$Assets[-1,]*(1-State$Haircut[-1])) #Liquidity that was there beforehand
  testLiquidity <- State$Liquidity + LiquidChange
  testLCR <- testLiquidity/State$NetOutflow
  LCRFailure <- testLCR < State$LCR
  LiquidityNotEmpty <- colSums(State$Assets*State$Haircut*Ranks) > 0
  BankFailsLiquidity <- sum(LCRFailure & LiquidityNotEmpty)
  return(BankFailsLiquidity)
}

BankingContagionCleanish <- function(Input,Strategies,Ranks,Debugging = F){
  InitialState <- Input
  State <- Input
  State$NetOutflow <- State$Outflow - pmin(0.75*State$Outflow,State$Inflow)
  Looping = TRUE
  Iter <- 0
  Defaults <- Matrix(0,ncol=length(Strategies),nrow=length(State$Equity)) # Matrix, Banks X Strategies, indicating default events per Bank for each separate issue
  Losses <- Matrix(0,ncol=length(Strategies),nrow=length(State$Equity)) # Matrix, Banks X Strategies, indicating relative losses for each Bank for each separate issue
  while(Looping){
    Iter <- Iter + 1
    ActualizedAssetChangeThisRound <- 0*State$Assets
    State$Equity <- InitialState$Equity + colSums(State$AssetPrices*State$Assets - InitialState$AssetPrices*InitialState$Assets)
    for(i in 1:length(Strategies)){
      StrategyAssetChanges <- Strategies[[i]][[1]](State,InitialState,Ranks[[i]])
      if(Debugging){print(StrategyAssetChanges)}
      DefaultResults <- Strategies[[i]][[2]](State,InitialState,StrategyAssetChanges)
      if(Debugging){print(DefaultResults)
      print(dim(Defaults))
      print(dim(DefaultResults$Defaults))}
      Defaults[,i] <- DefaultResults$Defaults
      if(Debugging){print("Calculation of Defaults happened")
      print(dim(Losses))
      print(dim(DefaultResults$Losses))}
      Losses [,i] <- DefaultResults$Losses
      if(Debugging){print("Calculation of Losses happened")}
      ActualizedAssetChangeThisRound <- pmin(ActualizedAssetChangeThisRound,StrategyAssetChanges)
    }
    
    State$Assets <- State$Assets + ActualizedAssetChangeThisRound
    
    for(i in 1:length(State$InterbankMarkets)){
      Remaining <- 1-CalculateRelativeLGD(Losses)
      NewInterbank <- State$InterbankMarkets[[i]]%*%diag(Remaining)
      #Adjust the assets corresponding to the interbank exposure network. These loans are part of the assets. Can probably be more beautiful.
      ActualizedAssetChangeThisRound[length(State$AssetPrices)-length(State$InterbankMarkets)+i,] <- rowSums(NewInterbank)-State$Assets[length(State$AssetPrices)-length(State$InterbankMarkets)+i,]
      State$Assets[length(State$AssetPrices)-length(State$InterbankMarkets)+i,]<-rowSums(NewInterbank)
    }
    State$NetOutflow <- CalculateNewNetOutflow(State,InitialState)
    State$AssetPrices <- CalculateNewPrices(State,InitialState)
    #Assets are sold down here, all of them at the newest price.
    NewCash <- CalculateCashGain(State,InitialState,ActualizedAssetChangeThisRound)
    NewAsset1 <- InitialState$Assets[1,] + NewCash
    ActualizedAssetChangeThisRound[1,] <- NewAsset1 - State$Assets[1,]
    State$Assets[1,] <- NewAsset1
    if(CheckIfAssetChangeSmall(State,ActualizedAssetChangeThisRound)){Looping=FALSE}
    if(CheckLiquiditySituation(State,InitialState,Ranks)){Looping = TRUE}
  }
  output <- c(State,Defaults = Defaults, Losses = Losses,Iters = Iter)
  return(output)
}
