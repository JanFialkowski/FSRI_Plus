# Set working directory
setwd("/Users/janfialkowski/Documents/Work/FSRI_Plus/")

# Load required libraries
library(Matrix)
library(data.table)
library(GLcascade)
library(fastcascade)
library(igraph)
library(colorspace)
library(parallel)

# Data preparation ####

source("./Scripts/data_loader.R")

# Firm network: edge list with supply chain links
# Firm list: data table with firm characteristics
# - id
# - out strength
# - nace4 (industry code)
# - financial variables: net profit, operational profit, retained earnings, equity, liquidity

# Banks equity: vector of bank equities
# Firm-bank loans: matrix of firms x banks with loan amounts

# Firms
n <- 15  # Number of firms

# Create firm list
firm_list <- data.table(1:n)
colnames(firm_list) <- c("id")
firm_list$nace4 <- c(1920, 3523, 4221, 710, 4334, 2410, 2410, 5221, 5221, 4531, 2410, 2410, 4612, 2410, 4612)

# Financial variables (values in 1000 HUF)
#(revenue,material_costs,cash,net_short_term_liquidity,short_term_assets)
firm_list$material_costs <- c(500, 400, 500, 200, 150, 40, 100, 350, 100, 300, 200, 150, 100, 50, 5)
firm_list$revenue <- c(600, 450, 520, 300, 160, 80, 600, 400, 110, 190, 280, 80, 20, 55, 10)
firm_list$short_term_liabilities <- c(175, 120, 70, 40, 10, 50, 40, 20, 20, 20, 16, 16, 15, 15, 8)
firm_list$equity <- c(500, 200, 300, 200, 100, 100, 50, 100, 50, -100, 200, 80, 60, 50, 50)
firm_list$short_term_assets <- c(300, 100, 300, 100, 80, 150, 80, 50, 60, 200, 150, 100, 50, 60, 10)

# Edge list: (supplier id, buyer id, link weight)
el <- rbind(
  c(1, 8, 300),
  c(1, 9, 500),
  c(1, 3, 500), 
  c(2, 10, 1000),
  c(2, 7, 800),
  c(2, 3, 500),
  c(3, 7, 2000),
  c(3, 8, 1500), 
  c(4, 2, 800),
  c(4, 3, 1000),
  c(4, 6, 600),
  c(5, 3, 1000),
  c(5, 7, 500),
  c(6, 5, 400), 
  c(7, 8, 200),
  c(7, 9, 300),
  c(7, 11, 200),
  c(8, 10, 2000),
  c(9, 10, 600),
  c(11, 15, 30),
  c(12, 8, 100),
  c(12, 11, 100),
  c(13, 2, 150),
  c(14, 11, 20),
  c(14, 12, 100),
  c(14, 13, 150),
  c(15, 1, 100)
)

firm_network <- data.table(el)
colnames(firm_network) <- c("supplier", "buyer", "weight")

# Assign NACE4 codes from firm_list to buyers and suppliers
firm_network <- merge(firm_network, firm_list[, .(id, buyer_nace4 = nace4)], by.x = "buyer", by.y = "id")
firm_network <- merge(firm_network, firm_list[, .(id, supplier_nace4 = nace4)], by.x = "supplier", by.y = "id")

# Adjacency matrix of the supply chain network (SCN)
W <- Matrix::sparseMatrix(i = el[, 1], j = el[, 2], x = el[, 3], dims = c(n, n))

# Calculate out-strength of firms
firm_list$out_strength <- rowSums(W)

# Banks
# Number of banks
t <- 4
equity_banks <- c(1000, 100, 100, 1000)  # Equity in 1000 HUF

# Firm-bank links with outstanding principal
bl <- rbind(
  c(1, 2, 100),
  c(1, 3, 50),
  c(1, 5, 50), 
  c(2, 3, 60),
  c(2, 7, 40),
  c(2, 9, 10),
  c(3, 10, 20),
  c(3, 8, 30), 
  c(4, 10, 10),
  c(4, 8, 30),
  c(4, 11, 10),
  c(4, 12, 30)
)

# Firm x bank matrix with loans
bank_firm_matrix_out <- Matrix::sparseMatrix(i = bl[, 2], j = bl[, 1], x = bl[, 3], dims = c(n, t))

# Filter out firms with bad financial data
bank_firm_matrix_out[which(firm_list$equity<0),] <- 0 # negative Equity
bank_firm_matrix_out[which((firm_list$revenue-firm_list$material_costs)<0),] <- 0 # negative profit
bank_firm_matrix_out[which((firm_list$short_term_assets-firm_list$short_term_liabilities)<0),] <- 0 # negative liquidity

# Firm x bank matrix weighted by bank equities
loan_fb_matrix_bank_equities <- t(t(bank_firm_matrix_out) / equity_banks)

# Bank to Bank edgelist
L_el <- rbind(
  c(2, 1, 100),
  c(1, 3, 50),
  c(2, 4, 60),
  c(4, 1, 100)
)

bank_bank_matrix <- Matrix::sparseMatrix(i = L_el[, 2], j = L_el[, 1], x = L_el[, 3], dims = c(t, t))

# Create random banking data and replace the relevant fields.
InputData <- LoadRandomData(Banks=4, Assets=15)
InputData$State$Equity <- equity_banks
InputData$State$StartingEquity <- equity_banks
InputData$State$InterbankMarkets[[1]] <- bank_bank_matrix
InputData$State$Assets[15,] <- rowSums(bank_bank_matrix)
InputData$State$StartingAssets[15,] <- rowSums(bank_bank_matrix)

# FSRI ####

source("./Scripts/functions_supply_chain_contagion.R")
source("./Scripts/financial_systemic_risk_functions_II.R")
source("./Scripts/FSRIPlusPipeline.R")

ESRI <- GL_cascade(W,
                   firm_list$nace4,
                   ess_mat_sec=ess_mat_sec_gl,
                   use_rcpp = F,
                   ncores = 0,
                   track_h = T)

H <- H_pmax(ESRI)
Firm_Shock <- financial_shock(H,firm_list)
FSRI_Bank_Shock <- bank_shock(Firm_Shock,loan_fb_matrix_bank_equities)

Model <- "DebtRank" #Andras or DebtRank, chooses the Contagion model
Name <- "SingleFirm_FSRI" #Unique identifier for the run. Appears in the filenames and can potentially overwrite old stuff
ShockedProperties <- c("Equity") # Options: Equity, Liquidity, Inflow
RelativeShocks <- c(T) # Vector of bools, true denotes that the shock is relative to the current value, false is an absolute shock.
Path <- paste0(getwd(),"/")

# Runs the DebtRank cascade and does the analysis. Results are saved in the data folder
FSRIPlus(Name, InputData = InputData, Scenarios=FSRI_Bank_Shock$loss_out_either,
        Model = Model, ShockedProperties = ShockedProperties, 
        RelativeShocks = RelativeShocks,
        path = Path, overwrite = F)

# Sampling Covid shocks ####

source("./Scripts/sample_synthetic_firm_level_shocks.R")

# 100% shock to firm 7, 0% to all other firms
shocked_firm <- 7
xi_0 <- c(rep(0, shocked_firm-1), 1, rep(0, n - shocked_firm))

# in- and out-strength of firms
s_in <- colSums(W)
s_out <- rowSums(W)

set.seed(100)
xi_synth <- sample_firm_lev_shocks(psi_k_mat = NULL, # named with nace category, percentage shock to sector k instrength in the first row, and out strength in the second row
                                   firm_lev_shock = xi_0, # n dim vector, elements \in [0,1], empirical shock for each firm
                                   s_in = s_in, # instrengths of firms of sector k 
                                   s_out = s_out, # outstrengths of firms of sector k
                                   n_scen = 1000, # number of shocks for the sector,
                                   #m_secs = m_secs, # number of firms within the sector
                                   nace_cat = substring(as.character(firm_list$nace4),1,2), # vector with firms nace2 categories
                                   tracker = TRUE, 
                                   sample_mode = "empirical",
                                   silent = FALSE)

Covid_Direct_Shock <- financial_shock(xi_synth$psi_mat,firm_list)
Covid_Direct_Bank_Shock <- bank_shock(Covid_Direct_Shock,loan_fb_matrix_bank_equities)

Covid_Cascade_Shock <- GL_cascade(W,
                   firm_list$nace4,
                   psi_mat = xi_synth$psi_mat,
                   ess_mat_sec=ess_mat_sec_gl,
                   use_rcpp = F,
                   ncores = 0,
                   track_h = T)

H_Covid <- H_pmax(Covid_Cascade_Shock)
Covid_Cascade_Shock <- financial_shock(H_Covid,firm_list)
Covid_Cascade_Bank_Shock <- bank_shock(Covid_Cascade_Shock,loan_fb_matrix_bank_equities)

Name_CSD <- "CovidShock_Direct"
# Runs the DebtRank cascade and does the analysis. Results are saved in the data folder
FSRIPlus(Name_CSD, InputData = InputData,
         Scenarios=Covid_Direct_Bank_Shock$loss_out_either,
         Model = Model, ShockedProperties = ShockedProperties, 
         RelativeShocks = RelativeShocks,
         path = Path, overwrite = T)


Name_CS <- "CovidShock"
# Runs the DebtRank cascade and does the analysis. Results are saved in the data folder
FSRIPlus(Name_CS, InputData = InputData,
         Scenarios=Covid_Cascade_Bank_Shock$loss_out_either,
         Model = Model, ShockedProperties = ShockedProperties, 
         RelativeShocks = RelativeShocks,
         path = Path, overwrite = T)

# SI Analysis ####

# DebtRanks for the individual banks
source("./Scripts/DebtRank.R")
DR <- matrix(0,t,t)
for(i in 1:t){
  h0 <- c(rep(0,i-1),1,rep(0,t-i))
  DR[i,] <- DebtRank1(h0,t(bank_bank_matrix)/equity_banks)
}
saveRDS(DR,"./data/DebtRank.rds")

# This repeats the Cascades with a combination of random shocks and random IB networks
IBmat <- bank_bank_matrix/equity_banks
source("./Scripts/SI_extra.R")



# Redoing the figures ####

source("./Scripts/ProperPicsing.R")
Figure2FSRI()
Figure2_FSRI_CCDF()
Figure2_SI()
Figure2_DebtRanks()
Figure3()
Figure3_SI_OtherRisks()
Figure3_SystemSI()
Figure4()
Figure4aAggregated()
Figure4aLogged()
Figure4_AmpCCDF()
Figure4_FactorBreaking()
Figure4bc_SI()
SI_Figure_BankingShock()
SI_Figure_BankingShock_Logged()
SystemHist()
