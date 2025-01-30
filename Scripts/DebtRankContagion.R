DebtRankAdjustment <- function(State,InitialState,Ranks){
  return(0)
}

DebtRankLosses <- function(State,InitialState,Results){
  Defaults <- State$Equity < 10^-12
  Losses <- pmin((State$StartingEquity-State$Equity)/State$StartingEquity,1)
  return(list(Defaults=Defaults,Losses=Losses))
}
