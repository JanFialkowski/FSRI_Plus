library(Matrix)
library("data.table")

#zlatas vetka
#0#
H_pmax <- function(results_cascade = results_h){
  
  n<-dim(results_cascade$hd_T_mat)[2]
  m<-dim(results_cascade$hd_T_mat)[1]
  
  blocks<-c(seq(1,n,1000),n+1)
  
  H<-Matrix(0, ncol=n, nrow = m)
  
  for (i in 1:length(blocks[-1])) {
    cat(i, " ")
    cat(as.character(Sys.time()), " ", "\n")
    H[ , blocks[i]:(blocks[i+1]-1)] <- as.matrix( pmax( results_cascade$hd_T_mat[ ,blocks[i]:(blocks[i+1]-1)],
                                                        results_cascade$hu_T_mat[ ,blocks[i]:(blocks[i+1]-1)]))
  }
  
  return(H)
}

#1#
financial_shock<-function(max_hd_hu = H,
                          d = firm_list,
                          param_eq=1,
                          param_cash=1,
                          param_liq=1,
                          param_short=1){
  
  rev<-d$revenue
  mat_costs<-d$material_costs
  equity<-d$equity
  liquidity<-d$short_term_assets - d$short_term_liabilities
  
  N <- dim(max_hd_hu)[1]
  bad_either <- Matrix(0, nrow = N , ncol = N)
  bad_liq <- Matrix(0, nrow = N , ncol = N)
  bad_eq <- Matrix(0, nrow = N , ncol = N)
  
  eq_inv <- as.numeric( ((equity== 0)*1 + (equity != 0)*equity)^(-1))*(equity != 0)
  liq_inv <- as.numeric(( (liquidity== 0)*1 + (liquidity != 0)*liquidity )^(-1))*(liquidity!= 0)
  
  re_minus_co <- as.numeric(rev-mat_costs) # profit
  after <- max_hd_hu * re_minus_co # Delta profit
  
  bad_eq <- ((after * eq_inv) >= param_eq)
  bad_liq <- ((after * liq_inv) >= param_liq) 
  bad_either <- (bad_eq + bad_liq)>0
  
  return(list(bad_either=bad_either,
              bad_liq=bad_liq,
              bad_eq=bad_eq
  ))
  
}


#2#
bank_shock<-function(bad,
                     B_outst = bank_firm_matrix_out){
  
  ### Calculates the total loss as a results of firms defaulting on their outstanding loans
  ### Outputs total loss per scenario, with every row being a scenario
  bad_T_eq<-t(bad$bad_eq) * 1
  shock_outst_eq<-bad_T_eq %*% B_outst#column is a bank, row is a firm
  
  bad_T_either<-t(bad$bad_either) * 1
  shock_outst_either<-bad_T_either %*% B_outst#column is a bank, row is a firm
  
  bad_T_liq<-t(bad$bad_liq) * 1
  shock_outst_liq<-bad_T_liq %*% B_outst#column is a bank, row is a firm
  
  return(list(loss_out_eq=shock_outst_eq,
              loss_out_either=shock_outst_either,
              loss_out_liq=shock_outst_liq
              ))
  
}


#3#
losses_scenarios <- function(loss_out,
                             bad,
                             B = bank_firm_matrix_out,
                             equity_banks = equity_banks){
  ###
  
  loans_total <- sum(B)
  equity_total <- sum(equity_banks)
  loans_per_bank <- colSums(B)
  
  equity_inv <- (equity_banks != 0)*as.numeric( ((equity_banks == 0)*1 + 
                                                 (equity_banks != 0)*equity_banks)^(-1) )
  
  losses_over_total_loans_fb_mat <- loss_out * loans_total^(-1)

  losses_over_equity_bf_mat <- t(loss_out) * equity_inv
  
  losses_over_equity_bf_mat2 = pmin(losses_over_equity_bf_mat, 1)
  
  equity_loss_f_vec <- crossprod( losses_over_equity_bf_mat2, (equity_banks / equity_total) )
  loan_weighted_equity_loss_f_vec <- crossprod( losses_over_equity_bf_mat2, (loans_per_bank / loans_total) )
  
  cat(as.character(Sys.time()), "\n")
  
  return(list(losses_over_total_loans_fb_mat = losses_over_total_loans_fb_mat,
              losses_over_equity_fb_mat = t(losses_over_equity_bf_mat),
              equity_loss_f_vec = equity_loss_f_vec,
              loan_weighted_equity_loss_f_vec = loan_weighted_equity_loss_f_vec,
              bad = bad))
}

#4#
general_outputs <- function(B = bank_firm_matrix_out,
                            equity_banks = equity_banks, 
                            H = H){
  
  loans_over_total_loans_fb_mat <- B / sum(B)
  
  loans_total <- sum(B)
  equity_total <- sum(equity_banks)
  loans_per_bank <- colSums(B)
  
  loans_over_total_equity_fb_mat <- B / equity_total
  
  equity_inv <- (equity_banks != 0)*as.numeric( ((equity_banks == 0)*1 + 
                                                  (equity_banks != 0)*equity_banks)^(-1) )
  
  loans_over_bank_equities_bf_mat <- t(B) * equity_inv
  
  h_times_bank_losses_fb_mat <- crossprod( H , B )
  
  h_losses_over_equity_bf_mat <- t(h_times_bank_losses_fb_mat) * equity_inv
  
  h_losses_over_equity_bf_mat2 = pmin(h_losses_over_equity_bf_mat, 1)
  
  equity_h_loss_f_vec <- crossprod( h_losses_over_equity_bf_mat2, (equity_banks / equity_total) )
  loan_weighted_equity_h_loss_f_vec <- crossprod( h_losses_over_equity_bf_mat2, (loans_per_bank / loans_total) )
  
  
  equity_ratio_b_vec <- equity_banks / equity_total
  
  return(list(loans_over_total_loans_fb_mat = loans_over_total_loans_fb_mat,
         loans_over_bank_equities_fb_mat = t(loans_over_bank_equities_bf_mat),
         loans_over_total_equity_fb_mat = loans_over_total_equity_fb_mat,
         h_times_bank_losses_fb_mat = h_times_bank_losses_fb_mat,
         h_losses_over_equity_fb_mat = t(h_losses_over_equity_bf_mat),
         equity_h_loss_f_vec = equity_h_loss_f_vec,
         loan_weighted_equity_h_loss_f_vec = loan_weighted_equity_h_loss_f_vec,
         equity_ratio_b_vec = equity_ratio_b_vec))
}
