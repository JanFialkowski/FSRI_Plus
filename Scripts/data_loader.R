library(data.table)

DataConsistencyChecks <- function(Data,numIBMarkets = 1){
  print("Checking MinimumPrices")
  Ib <- numIBMarkets
  TotalAss <- length(Data$CARRanks)
  Banks <- dim(Data$State$Assets)[2]
  print("Cash does not lose value")
  print(Data$State$MinimumPrices[1]==1)
  print("Interbankmarkets do not lose value")
  print(all.equal(Data$State$MinimumPrices[(TotalAss-Ib+1):TotalAss],rep(1,length.out=Ib)))
  print(paste0("Checking Interbank, number of markets currently expected: ",numIBMarkets))
  print(length(Data$State$InterbankMarkets)==numIBMarkets)
  print("Checking consistency of Assets with Interbank matrices")
  for(MarketID in seq_along(Data$State$InterbankMarkets)){
    print(all.equal(Data$State$Assets[TotalAss-Ib+MarketID,],rowSums(Data$State$InterbankMarkets[[MarketID]])))
  }
  print("Checking Diagonals of Interbankmarkets (you don't owe yourself money")
  for(MarketID in seq_along(Data$State$InterbankMarkets)){
    print(all.equal(diag(Data$State$InterbankMarkets[[MarketID]]),rep(0,Banks)))
  }
  print("Making sure the Interbanks are not filled with the same value")
  for(Market in Data$State$InterbankMarkets){
    print(length(table(Market))!=1)
  }
  print("Making sure cash and interbank assets are not sold")
  print("CARRanks for Cash and Interbanks")
  print(Data$CARRanks[1]==0)
  print(all.equal(Data$CARRanks[(TotalAss-Ib+1):TotalAss],rep(0,length.out=Ib)))
  print("LiquidityRanks for Cash and Interbanks")
  print(Data$LiquidityRanks[1]==0)
  print(all.equal(Data$LiquidityRanks[(TotalAss-Ib+1):TotalAss],rep(0,length.out=Ib)))
  print("Sanity checking Liquidity")
  liq <- Data$State$Liquidity
  print("Number of Liquidites fits to Assets dimensions")
  print(length(liq)==Banks)
  print("Positive Liquidites")
  print(all.equal(liq>0,rep(T,Banks)))
  print("Liquidites not unreasonably big")
  print(all.equal(liq/Data$State$NetOutflow<50,rep(T,Banks)))
  print("Sanity checking CAR")
  equ <- Data$State$Equity
  print("Number of Equities fits the Asset structure")
  print(length(equ)==Banks)
  print("Positive Equities")
  print(all.equal(equ>0,rep(T,Banks)))
  print("Equities not unreasonably big")
  print(all.equal(equ/Data$State$RiskWeightedAssets<10,rep(T,Banks)))
  print("Testing for Haircuts to be constrained to 0 and 1")
  print(all.equal(Data$State$Haircut>=0,Data$State$Haircut<=1))
  print("Testing for RiskWeights to be constrained to 0 and 1")
  print(all.equal(Data$State$RiskWeights>=0,Data$State$RiskWeights<=1))
  print("Testing for Cash to be Haircut- and Riskfree")
  print(all.equal(Data$State$Haircut[1]==0,Data$State$RiskWeights[1]==0))
  print("Checking if RiskWeightedAssets is not smaller than what we observe, might not matter")
  print(all.equal(Data$State$RiskWeightedAssets>=colSums(Data$State$Assets*Data$State$RiskWeights),rep(T,length.out=Banks)))
  print("Checking if Liquidity is not smaller than what we observe, might not matter")
  print(all.equal(Data$State$RiskWeightedAssets>=colSums(Data$State$Assets*(1-Data$State$Haircut)),rep(T,length.out=Banks)))
}

LoadRandomData <- function(Banks = 4,Assets = 4){
  State <- list()
  State$Assets <- runif(Banks*Assets)*4
  dim(State$Assets) <- c(Assets,Banks)
  State$CAR <- 0.08
  State$AssetPrices <- rep(1,Assets)
  State$StartingAssetPrices <- State$AssetPrices
  State$Equity <- runif(Banks)
  State$StartingEquity <- State$Equity
  State$RiskWeights <- sample(c(0,0.05,0.5,1),size=Assets,replace=T)
  State$RiskWeightedAssets <- colSums(State$Assets * State$RiskWeights * 1.1)
  State$Haircut <- sample(c(1,0,0.5),size=Assets,replace=T)
  State$MinimumPrices <- rep(0.75,Assets)
  State$OutsideAssets <- runif(Assets)*0
  InterbankMarkets <- runif(Banks^2)
  dim(InterbankMarkets) <- c(Banks,Banks)
  diag(InterbankMarkets) <- 0
  State$InterbankMarkets <- list(InterbankMarkets)
  State$LCR <- 1
  State$Liquidity <- runif(Banks)
  State$Inflow <- runif(Banks)
  State$Outflow <- runif(Banks)
  State$NetOutflow <- State$Outflow - pmin(0.75*State$Outflow,State$Inflow)
  
  State$Assets[Assets,] <- rowSums(State$InterbankMarkets[[1]])
  State$StartingAssets <- State$Assets
  State$MinimumPrices[1] <- 1
  State$MinimumPrices[Assets] <- 1
  
  LiquidityRanks <- sample(c(1,2,3),size=Assets,replace=T)
  CARRanks <- sample(c(1,2,3),size=Assets,replace=T)
  CARRanks[1] <- LiquidityRanks[1] <- 0
  CARRanks[Assets] <- LiquidityRanks[Assets] <- 0
  return(list(State=State,LiquidityRanks=LiquidityRanks,CARRanks=CARRanks))
}
