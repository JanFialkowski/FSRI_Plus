StrategySetup <- function(Type = c("Andras","DebtRank")){
  Type <- match.arg(Type)
  Strategies <- switch(Type,
                       Andras = list(list(CARAdjustment,CARDefaults),list(LiquidityAdjustment,LiquidityDefaults)),
                       DebtRank = list(c(DebtRankAdjustment,DebtRankLosses)))
  return(Strategies)
}

RankSetup <- function(Data, Type = c("Andras","DebtRank")){
  Type <- match.arg(Type)
  Ranks <- switch(Type,
                  Andras = list(Data$CARRanks,Data$LiquidityRanks),
                  DebtRank = list(rep(1,length(Data$State$AssetPrices))))
  return(Ranks)
}

BankingShock <- function(State,Shock,ShockedProperty){
  State[[ShockedProperty]] <- State[[ShockedProperty]] - Shock
  return(State)
}

ShockFactory <- function(Property=c("Equity","Liquidity","Inflow"),RelativeShock=NULL){
  ### Takes the name of Properties of banks e.g. Liquidity and Equity as a character vector and returns a named list with names equal to that vector and values being function that apply the shock to the named property.
  Property <- match.arg(Property)
  SimpleShockFactory <- function(Property,RelativeShocks=NULL){
    if(is.null(RelativeShocks)){
      return(function(Shock,State){
        if(max(Shock)<=1){print("No shock value greater than 1 detected and no flag set, treating the shock as a relative one"); Shock <- State[[Property]]*Shock}
        return(BankingShock(State,Shock,Property))
      })
    }else if(RelativeShocks){
      return(function(Shock,State){
        Shock <- State[[Property]]*Shock
        return(BankingShock(State,Shock,Property))
      })
    }else{
      return(function(Shock,State){
        return(BankingShock(State,Shock,Property))
      })
    }
    
  }
  Equity <- SimpleShockFactory("Equity",RelativeShock)
  Liquidity <- SimpleShockFactory("Liquidity",RelativeShock)
  Inflow <- function(Shock,State){
    Shockfunction <- SimpleShockFactory("Inflow",RelativeShock)
    State <- Shockfunction(State,Shock)
    State$NetOutflow <- State$Outflow - pmin(0.75*State$Outflow,State$Inflow)
    return(State)
  }
  return(get(Property))
}

ShockFunctionSetup <- function(ShockedProperties,RelativeShocks = rep(T,length(ShockedProperties))){
  if(length(RelativeShocks==1)){RelativeShocks <- rep(RelativeShocks,length(ShockedProperties))}
  
  output <- list()
  for(i in seq_along(ShockedProperties)){output[[ShockedProperties[[i]]]]<-ShockFactory(ShockedProperties[[i]],RelativeShocks[i])}
  return(output)
}
