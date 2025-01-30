library(RSpectra)

DebtRank2_2 <- function(h_0,leverages){
  L <- leverages
  h_out <- h_0
  h_p <- 0
  h_pp <- 0
  failure <- FALSE
  UpdateMatrix <- function(Matrix, h){
    Matrix[,h==1]<-0
    Matrix[h==1,]<-0
    return(Matrix)
  }
  StepOneTime <- function(h,ht1,Matrix){
    h<-h+Matrix %*% (h-ht1)
    h<-pmin(1,h)
  }
  CalcConvergence <- function(h, L){
    solve(diag(dim(L)[1])-L)%*%h
  }
  while(TRUE %in% ((h_out-h_p)!=0)){
    if(failure){
      L <- UpdateMatrix(L,h_pp)
      failure <- FALSE
      if(Mod(eigs(L,k=1,which="LM")$values)<1){
        dummy <- CalcConvergence(h_p-h_pp,L)+h_p
        dummy <- pmin(1,dummy)
        if(sum(dummy==1)==sum(h_p==1)){
          return(dummy)
        }
      }
    }
    if(sum(h_p==1)>sum(h_pp==1)){
      failure <- TRUE
    }
    h_pp <- h_p
    h_p <- h_out
    h_out <- StepOneTime(h_p,h_pp,L)
  }
  return(h_out)
}

DebtRank1 <- function(h_0,leverages){
  L <- pmin(1,leverages)
  dim(L)<-dim(leverages)
  h_out <- h_0
  h_p<-0
  UpdateMatrix <- function(Matrix, h){
    Matrix[,h>0]<-0
    return(Matrix)
  }
  StepOneTime <- function(h,Matrix){
    h<-h+Matrix %*% h
    h<-pmin(1,h)
  }
  while(TRUE %in% ((h_out-h_p)!=0)){
    h_p <- h_out
    h_out<-StepOneTime(h_p,L)
    L<-UpdateMatrix(L,h_p)
  }
  return(h_out)
}
TurnLiabilitesIntoLeverages <- function(L,E){
  return(as.matrix(t(L)/E))#L_ij is the money owed from i to j, t(L) is "asset space"
}
CollectiveDebtRanking <- function(Liabilities,Equities,H){
  Leverages <- TurnLiabilitesIntoLeverages(Liabilities,Equities)
  output <- list(Original = H,DebtRank1 = H, DebtRank2 = H)
  for (i in 1:dim(H)[1]){
    output$DebtRank1[i,] <- DebtRank1(H[i,],Leverages)
    output$DebtRank2[i,] <- DebtRank2_2(H[i,],Leverages)
  }
  pmin(1,output$DebtRank2)
  dim(output$DebtRank2)<-dim(H)
  return(output)
}