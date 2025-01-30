#==============================================================================#
####################### sample synthetic shocks ################################
#==============================================================================#

# Code corresponding to algorithm in 
# SI Section 8. Constructing synthetic firm-level shocks with same sector level impact

# 1. wrapper function that performs function 2. for each sector supplied
# 2. sample synthetic firm-level shock corresponding to the sector shock
# 3. impute shocks for firms where no shock data was available with shocks other firms in the sector received

# provides the generalized inverse function
library(MASS)

#==============================================================================#
######################## 1. wrapper for all sectors ############################
#==============================================================================#

# takes an empirical firm-level, or an aggregated sector level shock and samples
# synthetic shocks affecting firms within sectors differently, but aggregated sectors
# the same for each sampled shock vector

sample_firm_lev_shocks <- function(psi_k_mat = NULL, # named with nace category, percentage shock to sector k instrength in the first row, and out strength in the second row
                                   firm_lev_shock = firm_lev_shock, # n dim vector, elements \in [0,1], empirical shock for each firm
                                   s_in = s_in, # instrengths of firms of sector k 
                                   s_out = s_out, # outstrengths of firms of sector k
                                   n_scen = n_scen, # number of shocks for the sector,
                                   #m_secs = m_secs, # number of firms within the sector
                                   nace_cat = nace_cat, # vector with firms nace categories
                                   tracker = FALSE, 
                                   sample_mode = "uniform",
                                   silent = TRUE
){
  
  # general numbers
  n <- length(nace_cat)
  
  if((length(psi_k_mat)>0) & (length(firm_lev_shock) > 1)){
    cat("SOFT WARNING: \n
        Empirical firm level shocks and sector level perecentage shocks have been supplied: \n 
        Using empirical firm level shocks as default!")
    psi_k_mat <- NULL
  }
  
  # if firm shock creation is initialized with percentage shocks from each sector
  if(length(names(psi_k_mat)) > 0){
    
    shocked_naces <- names(psi_k_mat) 
    
    sectors <- sort(unique(nace_cat))
    m_secs <- length(sectors)
    cat(m_secs, "industry sectors contained in p", "\n")
    
  }
  
  # if firm shock creation is initialized with an empirical firm level shocks, e.g. covid-19 employment shock
  if(length(firm_lev_shock) > 1){
    
    # identify the sectors that are in the data and will be shocked
    sectors <- sort(unique(nace_cat))
    m_secs <- length(sectors)
    cat(m_secs, "industry sectors contained in nace_cat", "\n")
    
    
    # calculate the s_in_k size for each sector and the s_out_k size for each sector 
    s_in_sec_vec  <- sapply(sectors, function(x) s_in  %*% (nace_cat == x))
    s_out_sec_vec <- sapply(sectors, function(x) s_out %*% (nace_cat == x))
    
    # calculate the absolute shock s_in_k size for each sector and the absolute shock s_out_k size for each sector 
    s_in_shock_abs_sec_vec  <- sapply(sectors, function(x) s_in  %*% (firm_lev_shock * (nace_cat == x)))
    s_out_shock_abs_sec_vec <- sapply(sectors, function(x) s_out %*% (firm_lev_shock * (nace_cat == x)))
    
    # calculate the relative shock s_in_k size for each sector and the relative shock s_out_k size for each sector 
    s_in_shock_rel_sec_vec   <- ifelse(s_in_sec_vec  < (10^-2), 0,   s_in_shock_abs_sec_vec  / s_in_sec_vec)
    s_out_shock_rel_sec_vec  <- ifelse(s_out_sec_vec < (10^-2), 0,   s_out_shock_abs_sec_vec / s_out_sec_vec)
    
    
    # list for the shock distribution for each sector
    firm_lev_shock_sector_list <- lapply(sectors, function(x) firm_lev_shock[nace_cat == x])
    names(firm_lev_shock_sector_list) <- sectors
    
    psi_k_mat <- rbind(s_in_shock_rel_sec_vec,
                       s_out_shock_rel_sec_vec)
  }
  
  # sector shocks from firm level shock
  psi_inds <- matrix(0, nrow=0, ncol=3)
  track_mat <-  matrix(0, nrow=0, ncol=4)
  
  for(i in 1:m_secs){
    
    target_ids <- which(nace_cat == sectors[i] )
    # number of firms in the sector
    k_sec <- length(target_ids)
    
    names(target_ids) <- as.character(1:k_sec)
    
    cat("SECTOR: ", i, " out of ", m_secs, "| sector name:", sectors[i] ,"| sector size:", k_sec , "\n")
    
    
    # distribution of positive shocks to sample from
    emp_firm_shock_dist <- firm_lev_shock_sector_list[[i]][firm_lev_shock_sector_list[[i]] > 0]
    
    
    
    # create shocks for one sector
    psi_k_inds <- sample_firm_lev_shocks_of_sector(psi_k_in  = psi_k_mat[1, i], # percentage shock to sector k
                                                   psi_k_out = psi_k_mat[2, i],
                                                   s_in  = s_in[target_ids], # instrengths of firms of sector k 
                                                   s_out = s_out[target_ids], # outstrengths of firms of sector k
                                                   n_scen = n_scen, # number of shocks for the sector,
                                                   k_sec = k_sec, # number of firms within the sector,
                                                   target_ids = target_ids,
                                                   tracker = tracker,
                                                   emp_firm_shock_dist = emp_firm_shock_dist,
                                                   sample_mode = sample_mode,
                                                   silent = silent
    )
    
    # bind output together
    psi_inds <- rbind(psi_inds, psi_k_inds$arr_inds)
    
    if(tracker == TRUE){
      track_mat <- rbind(track_mat, cbind(psi_k_inds$track_mat, rep(i, n_scen)) )
    }
    
  }
  
  # make psi_mat out of the row and column indices
  
  psi_mat <- sparseMatrix(i = psi_inds[,1],
                          j = psi_inds[,2],
                          x = psi_inds[,3],
                          dims = c(n, n_scen))
  
  outlist <- list(psi_mat = psi_mat)
  
  if(tracker == TRUE){
    colnames(track_mat) <- c("nr_iterations", "in_str_deviation_perc", "out_str_deviation_perc", "sector")
    outlist <- c(outlist, 
                 list(track_mat = track_mat)
    )
  }
  
  outlist
}


#==============================================================================#
####################### 2. sampler for a single sector #########################
#==============================================================================#



sample_firm_lev_shocks_of_sector <- function(psi_k_in = psi_k_in, # percentage shock to sector k
                                             psi_k_out = psi_k_out, # percentage shock to sector k
                                             s_in = s_in, # instrengths of firms of sector k 
                                             s_out = s_out, # outstrengths of firms of sector k
                                             n_scen = n_scen, # number of shocks for the sector,
                                             k_sec = k_sec, # number of firms within the sector,
                                             target_ids = target_ids, # ids of firms in sector k
                                             tracker = tracker, # should calculations be tracked
                                             emp_firm_shock_dist = emp_firm_shock_dist, # empirical distribution of firm level shocks
                                             sample_mode = sample_mode, # should empirical shocks be used or shocks be drawn from uniform distributions
                                             silent = silent){
  
  
  
  # track outputs
  if(tracker == TRUE){
    track_mat <- matrix(0, nrow = n_scen, ncol = 3)
  }
  
  ### special cases: 
  
  # if there is only one firm in the sector 
  if(k_sec == 1){
    # track outputs
    arr_inds <- cbind(rep(target_ids, n_scen), 1:n_scen, mean(psi_k_in, psi_k_out)) # should be that psi_k_in = psi_k_out if there is only one firm
    
    outlist <- list(arr_inds = arr_inds)
    
    if(tracker == TRUE){
      outlist <- c(outlist, 
                   list(track_mat = track_mat)
      )
    }
    
    cat("only one firm per sector", "\n")
    return(outlist)
    
  }
  
  if( mean(c(psi_k_in, psi_k_out)) <= 10^-15){
    # track output
    arr_inds <- cbind(rep(target_ids, n_scen), rep(1:n_scen, each = length(target_ids)), 0)
    
    
    outlist <- list(arr_inds = arr_inds)
    
    if(tracker == TRUE){
      outlist <- c(outlist, 
                   list(track_mat = track_mat)
      )
    }
    
    cat("sector shock was zero", "\n")
    
    return(outlist)
    
  }
  
  ### pre calculations 
  
  
  
  # absolute shock to in strength and out strength
  s_in_k_shock_abs <- sum(s_in) * psi_k_in
  s_in_k_shock_abs
  s_out_k_shock_abs <- sum(s_out) * psi_k_out
  s_out_k_shock_abs
  
  # which firms are having more in strength than out strength and vice versa
  # in_heavy  <- (s_in / (s_out + (s_out <= 0)) ) >  ( sum(s_in) / sum(s_out))
  # out_heavy <- (s_in / (s_out + (s_out <= 0)) ) <= ( sum(s_in) / sum(s_out))
  
  # which firms are having more in strength than out strength than the sector level shock and vice versa 
  in_heavy  <- (s_in / (s_out + (s_out <= 0)) ) >  ( s_in_k_shock_abs / (s_out_k_shock_abs + (s_out_k_shock_abs <=0)) )
  out_heavy <- (s_in / (s_out + (s_out <= 0)) ) <= ( s_in_k_shock_abs / (s_out_k_shock_abs + (s_out_k_shock_abs <=0)) )
  
  if(silent == FALSE){
    cat("number of in_heavy: ", sum(in_heavy), "\n")
    cat("number of out_heavy: ", sum(out_heavy), "\n") 
  }
  
  
  if( (sum(in_heavy)==0) | (sum(out_heavy) == 0)){
    
    cat("either no in_heavy or no out_heavy firms  \n")
    
    
    if((sum(in_heavy) == 0)){
      in_heavy <- abs( (s_in / (s_out + (s_out <= 0)) ) -  ( s_in_k_shock_abs / (s_out_k_shock_abs + (s_out_k_shock_abs <=0)) )) < 10^-10
      
      out_heavy <- (! in_heavy) & out_heavy
      
      if(sum(in_heavy) > 0){
        cat("added edge case to in_heavy and deleted from out_heavy \n")
      }
    }
    
    if((sum(out_heavy) == 0)){
      out_heavy <- abs( ( s_in / (s_out + (s_out <= 0)) ) -  
                          ( s_in_k_shock_abs / (s_out_k_shock_abs + (s_out_k_shock_abs <=0)) )
      ) < 10^-10
      
      in_heavy <- (! out_heavy) & in_heavy
      
      if(sum(out_heavy) > 0){
        cat("added edge case to out_heavy and deleted from in_heavy  \n")
      }
      
    }
    
  }
  
  ### sample shocks for each scenario such that the aggregated shocks are larger than s_in_k_shock_abs and s_out_k_shock_abs
  
  # permute the sequence of firms for each scenario to calculate (properly randomize who gets a shock)
  n_sec_permutations <- matrix(sapply(1:n_scen, function(x) sample(1:k_sec )),
                               nrow = k_sec, 
                               ncol = n_scen)
  
  # order to permute the firms back 
  n_sec_inv_permutations <- sapply(1:n_scen, function(x) order(n_sec_permutations[ ,x]))
  # sapply(1:n_scen, function(x) n_sec_permutations[n_sec_inv_permutations[,x] ,x])
  
  # k_sec x n_scen matrices with permuted in and out strength vectors for each scenario
  s_in_perm_mat  <- sapply(1:n_scen, function(x) s_in[n_sec_permutations[,x]])
  s_out_perm_mat <- sapply(1:n_scen, function(x) s_out[n_sec_permutations[,x]])
  
  
  
  # initialization for adding additional shocks until size is large enough
  psi_mat_k <- matrix(0, nrow = k_sec, ncol=n_scen)
  s_in_psi_cutoff <- s_out_psi_cutoff <- rep(NA, n_scen)
  t = 0
  
  while( sum( is.na(s_in_psi_cutoff) | is.na(s_out_psi_cutoff)) > 0 ){
    t <-  t + 1 
    
    if(silent==FALSE){
      cat(t, "   ")
    }
    # which scenarios need additional shocks
    scen_to_add <- is.na(s_in_psi_cutoff) | is.na(s_out_psi_cutoff)
    
    
    # sample shocks for scenarios where aggegated shocks are still smaller than shock target; alternative sample from empirical distribution
    if(sample_mode == "uniform"){
      psi_mat_k[, scen_to_add] <-  pmin(1, psi_mat_k[, scen_to_add] +  matrix(runif(sum(scen_to_add) * k_sec), 
                                                                              nrow = k_sec, 
                                                                              ncol=sum(scen_to_add)))
    }
    
    
    if(sample_mode == "empirical"){
      psi_mat_k[, scen_to_add] <-  pmin(1, psi_mat_k[, scen_to_add] +  matrix( emp_firm_shock_dist[sample(1:length(emp_firm_shock_dist), sum(scen_to_add) * k_sec, replace=TRUE)], 
                                                                               nrow = k_sec, 
                                                                               ncol=sum(scen_to_add)))
    }
    
    # absolute shock sizes for each scenario
    s_in_shock_abs  <- psi_mat_k * s_in_perm_mat
    s_out_shock_abs <- psi_mat_k * s_out_perm_mat
    
    # cumulative absolute shock sizes for each scenario
    s_in_shock_abs_cum  <- apply(s_in_shock_abs,  MARGIN=2, FUN=cumsum)
    s_out_shock_abs_cum <- apply(s_out_shock_abs, MARGIN=2, FUN=cumsum)
    
    
    # vector that for each scenario (column) determine at which row index the sector shock size is reached
    s_in_psi_cutoff  <- sapply(1:n_scen, function(x) match(TRUE, s_in_shock_abs_cum[,x]  >= s_in_k_shock_abs) )
    s_out_psi_cutoff <- sapply(1:n_scen, function(x) match(TRUE, s_out_shock_abs_cum[,x] >= s_out_k_shock_abs) )
    
  }
  
  
  cutoff <- pmax(s_in_psi_cutoff, s_out_psi_cutoff)
  
  # set to zero the shock vec for firms after reaching the shock threshold
  # set_to_zero <- pmin(k_sec, cutoff+1) # largest row index to keep the shock
  # set_to_zero <- sapply(1:n_scen, function(x) (1:k_sec) < set_to_zero[x] )
  
  # set to zero the shock vec for firms after reaching the shock threshold
  set_to_zero <- sapply(1:n_scen, function(x) (1:k_sec) <= cutoff[x] )
  
  psi_mat_k <- psi_mat_k * set_to_zero
  psi_mat_k <- sapply(1:n_scen, function(x) psi_mat_k[n_sec_inv_permutations[,x] ,x])
  
  # check if criterion worked
  # s_in %*% psi_mat_k >= s_in_k_shock_abs
  # s_out %*% psi_mat_k >= s_out_k_shock_abs
  
  
  psi_k_out_mat <- matrix(0, nrow = k_sec, ncol = n_scen)
  
  # do the rescaling of shocks to firms such that affected in and out strength have the right size as in the original sector level shock
  for(k_it in 1:n_scen){
    
    
    if(silent == FALSE){
      cat("\n","\n",  "Scenario:", k_it, "\n", "\n")
    }
    
    psi <- psi_mat_k[, k_it]
    psi_full <- rep(FALSE, k_sec)
    
    crit <- TRUE
    t <- 0
    while(crit){ # rescale and add shocks until aggregated firm level shock matches sector shock
      t <- t+1
      #cat(t, " ")
      
      # include only psi elements that are smaller than 1, and also only rescale them. 
      # need to subtract the ones with psi=1 on the right hand side too
      
      
      # the remaining absolute shock size to distribute after deducting 100% shocked firms
      b_ls <- c(s_in_k_shock_abs - psi_full %*% s_in, 
                s_out_k_shock_abs - psi_full %*% s_out)
      
      # update which firms are in_heavy and which are out_heavy when remaining shock size changes
      in_heavy  <- (s_in / (s_out + (s_out <= 0)) ) >  ( b_ls[1] / (b_ls[2] + (b_ls[2] <=0)) )
      out_heavy <- (s_in / (s_out + (s_out <= 0)) ) <= ( b_ls[1] / (b_ls[2] + (b_ls[2] <=0)) )
      
      
      # check for the edge case of in and out heavy:
      if( (sum(in_heavy)==0) | (sum(out_heavy) == 0)){
        cat("no in_heavy or out_heavy firms  \n")
        if((sum(in_heavy) == 0)){
          in_heavy <- abs( (s_in / (s_out + (s_out <= 0)) ) -  ( s_in_k_shock_abs / (s_out_k_shock_abs + (s_out_k_shock_abs <=0)) )) < 10^-10
          out_heavy <- (! in_heavy) & out_heavy
          if(sum(in_heavy) > 0){
            cat("added edge case to in_heavy and deleted from out_heavy \n")
          }
        }
        if((sum(out_heavy) == 0)){
          out_heavy <- abs((s_in / (s_out + (s_out <= 0)) ) -  ( s_in_k_shock_abs / (s_out_k_shock_abs + (s_out_k_shock_abs <=0)) )) < 10^-10
          in_heavy <- (! out_heavy) & in_heavy
          if(sum(out_heavy) > 0){
            cat("added edge case to out_heavy and deleted from in_heavy  \n")
          }
        }
      }
      
      #barplot(rbind(in_heavy, out_heavy), col=c("red", "blue"))
      
      
      # first row contains in strength shock sizes, second row out strength shock sizes
      A_ls <- rbind((s_in  * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy),
                    (s_out * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy) )
      
      
      # check if it is solvable
      #qr(A_ls1)$rank == qr(cbind(A_ls1, c(s_in_k_shock_abs ,s_out_k_shock_abs )))$rank
      
      solvable <- qr(A_ls)$rank == qr(cbind(A_ls, b_ls))$rank
      
      if(silent == FALSE){
        if(solvable){
          cat("Linear systems has a solution", "\n")
        }else{
          
          cat("Linear system does not have a solution: add shocks", "\n")
          
          if(sum(A_ls[,'in_heavy'])<10^-15){
            cat("in_heavy firms did not receive shocks \n")
          }
          
          if(sum(A_ls[,'out_heavy'])<10^-15){
            cat("out_heavy firms did not receive shocks \n")
          }
          
          if(sum(A_ls[1,])<10^-15){
            cat("sum of in_str shocks is zero \n")
          }
          
          if(sum(A_ls[2,])<10^-15){
            cat("sum of out_str shocks is zero \n")
          }
        }      
      }
      
      
      # if there is no in heavy firm affected shock an additional in heavy firm
      if( sum(A_ls[,'in_heavy']) <=0){ 
        # shock by how much or how many?}
        if(silent == FALSE){
          cat("added shock to in_heavy firms", "\n")
        }
        # firms i can shock additionally
        firms_to_shock <- which( (in_heavy * ( psi <= 0)) > 0)
        add_shock <- sample(1:length(firms_to_shock), 1)
        psi[firms_to_shock[add_shock]] <- runif(1)
        
        A_ls <- rbind((s_in  * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy),
                      (s_out * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy) )
      }
      
      
      # if there is no out heavy firm affected shock an additional in heavy firm
      if( sum(A_ls[,'out_heavy']) <=0){ 
        if(silent == FALSE){
          cat("added shock to out_heavy firms", "\n")
        }
        # firms i can shock additionally
        firms_to_shock <- which( (out_heavy * ( psi <= 0)) > 0)
        add_shock <- sample(1:length(firms_to_shock), 1)
        psi[firms_to_shock[add_shock]] <- runif(1)      
        
        A_ls <- rbind((s_in  * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy),
                      (s_out * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy) )
      }
      
      
      # if there is no in heavy firm affected shock an additional in heavy firm
      if( sum(A_ls[1,]) <=0){ 
        # shock by how much or how many?}
        if(silent == FALSE){
          cat("added shock to in_str firm", "\n")
        }
        # firms i can shock additionally
        firms_to_shock <- which( (s_in * ( psi <= 0)) > 0)
        add_shock <- sample(1:length(firms_to_shock), 1)
        psi[firms_to_shock[add_shock]] <- runif(1)
        
        A_ls <- rbind((s_in  * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy),
                      (s_out * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy) )
      }
      
      
      # if there is no out heavy firm affected shock an additional in heavy firm
      if( sum(A_ls[2,]) <=0){ 
        if(silent == FALSE){
          cat("added shock to out_str firm", "\n")
        }
        # firms i can shock additionally
        firms_to_shock <- which( (s_out * ( psi <= 0)) > 0)
        add_shock <- sample(1:length(firms_to_shock), 1)
        psi[firms_to_shock[add_shock]] <- runif(1)      
        
        A_ls <- rbind((s_in  * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy),
                      (s_out * psi * (!psi_full)) %*% cbind(in_heavy, out_heavy) )
      }
      
      
      
      # constraint matrix for solving undetermined system 
      # A_ls2 <- rbind(s_in  * psi * (!psi_full) ,
      #               s_out * psi * (!psi_full)  )
      # 
      # 
      # rescale_factors <- qr.solve(A_ls2, c(s_in_k_shock_abs - psi_full %*% s_in,
      #                                  s_out_k_shock_abs - psi_full %*% s_out))
      # 
      # psi * rescale_factors
      
      # rescale_factors1 <- solve(A_ls1, c(s_in_k_shock_abs, s_out_k_shock_abs))
      
      
      A_ls <- round(A_ls, 10)
      b_ls <- round(b_ls, 10)
      
      
      #xxx <- try(rescale_factors <- solve(a = A_ls, b = b_ls))
      #rescale_factors <- solve(a = A_ls, b = b_ls)
      # 
      rescale_factors <- as.vector(MASS::ginv(A_ls) %*% b_ls)
      
      # is this a good choice?
      rescale_factors <- pmax(0, round(rescale_factors, 15))
      
      if(silent == FALSE){
        cat("rescaling factors: \n")
        print(rescale_factors)
      }
      
      # rescale shocks that were smaller than 1 and leave the ones that are 1 the same as 1.
      # psi1 <- ((psi * rescale_factors1['in_heavy'])  * in_heavy + 
      #           (psi * rescale_factors1['out_heavy']) * out_heavy) 
      
      psi <- (!psi_full) * ((psi * rescale_factors[1])  * in_heavy + 
                              (psi * rescale_factors[2]) * out_heavy) +
        psi_full
      
      if(silent == FALSE){
        cat("psi_i are larger than one: ", sum(psi > 1), "\n")
      }
      # check:
      # out and in shock too large
      #print( c(sum(s_in)  * psi_k - s_in  %*% psi1, sum(s_out) * psi_k - s_out %*% psi1))
      #print( c(sum(s_in)  * psi_k - s_in  %*% psi, sum(s_out) * psi_k - s_out %*% psi))
      
      
      psi_full <- (psi >= 1) | (psi < 0)
      psi <- pmax(0, pmin(1, psi))
      
      s_in_objective  <- s_in_k_shock_abs  - s_in  %*% psi
      s_out_objective <- s_out_k_shock_abs - s_out %*% psi
      
      if( (sum(psi > 1) == 0) & 
          (abs(s_in_objective) < 10^-2) & 
          (abs(s_out_objective) < 10^-2) ){
        # done and break the loop
        crit = FALSE
      }
      # else{
      #   # repeat with no elements larger than 1 and smaller than 0
      #   psi <- pmax(0, pmin(1, psi))
      #   
      #   # stopping criterion is that difference between firm shock and sector shock is zero and all psi_i \in [0,1]
      # }
      
      
      
      if(silent == FALSE){
        cat("s_in_objective: ", s_in_objective, ",  s_out_objective: ", s_out_objective, "\n")
      }
      
      if(t==max(1000,k_sec)){break}
    } # end while loop
    
    s_in_objective  <- s_in_k_shock_abs  - s_in  %*% psi
    s_out_objective <- s_out_k_shock_abs - s_out %*% psi
    
    if(tracker == TRUE){
      track_mat[k_it, ] <- c(t, 
                             s_in_objective ,#/ (s_in_k_shock_abs + (s_in_k_shock_abs<= 0)) , 
                             s_out_objective# / (s_out_k_shock_abs + (s_out_k_shock_abs <= 0)) 
      )
    }
    psi_k_out_mat[ , k_it] <- psi 
    
    #plot(track_mat[,1]); plot(track_mat[,2]); plot(track_mat[,3]);
    
  }# end for loop
  
  if(silent == FALSE){
    cat("check if all psi are smaller than one: ", sum(psi_mat_k > 1), "\n")
  }
  
  arr_inds <- cbind(which(psi_k_out_mat > 0 , arr.ind = TRUE), psi_k_out_mat[which(psi_k_out_mat > 0)])
  arr_inds[,1] <- target_ids[as.character(arr_inds[,1])] # use original indices of firms 
  
  
  outlist <- list(arr_inds = arr_inds)
  
  if(tracker == TRUE){
    outlist <- c(outlist, 
                 list(track_mat = track_mat)
    )
  }
  
  outlist
}













#==============================================================================#
####################### 3. impute shocks for unshocked firms ###################
#==============================================================================#

# create for an empirical shock distribution across sectors new shock distributions that
# shock also unmatched firms, i.e. sampling up the empirical scenario

shock_non_matched_firms <- function(firm_sec_memb, # sector affiliation of firms
                                    firm_matched,  # vector that is 1 if firm is matched 0 if not matched
                                    firm_lev_shock,# the empirical shock vector 
                                    n_scen = 10,   # number of extended scenarios to calculate 
                                    rsl_wd = "",
                                    save_psi_mat = TRUE){ # should psi_mat be saved directly to the path
  
  cat( as.character(Sys.time()), "\n")
  
  # identify the sectors that are present 
  sectors <- sort(unique(firm_sec_memb))
  
  # count firms per sector
  sectors_firm_count <- sapply(sectors, function(x) sum(firm_sec_memb==x))
  
  # count matched firms per sector
  sectors_firm_matched_count <- sapply(sectors, function(x) firm_matched %*% (firm_sec_memb == x))
  
  # split up the overall shock distribution to sector level shock distributions (to draw from them)
  matched_firm_lev_shock_sector_list <- lapply(sectors, function(x) firm_lev_shock[ (firm_sec_memb == x) & firm_matched ] )
  names(matched_firm_lev_shock_sector_list) <- sectors
  
  # firm_lev_shock_sector_list <- lapply(sectors, function(x) firm_lev_shock[ (firm_sec_memb == x) ] )
  #
  # # for each sector the number of firms that could be matched and received a shock 
  # sectors_firm_matched_shocked_count <- sapply(firm_lev_shock_sector_list, function(x) sum(x > 0))
  # 
  # # for each sector the fraction of firms that could be matched and received a shock out of the ones that could be matched
  # sectors_firm_matched_shocked_fraction <- sectors_firm_matched_shocked_count /
  #                                           ifelse(sectors_firm_matched_count > 0, sectors_firm_matched_count, 1)
  # 
  # # the number of additional firms to shock in each sector ### maybe round instead of floor is also ok?
  # firm_unmatched_tobe_shocked_count <- round(sectors_firm_matched_shocked_fraction * (sectors_firm_count - sectors_firm_matched_count))
  # 
  
  
  
  # generate the size empty psi_mat (each column contains the shock vector of one scenario)
  n <- length(firm_sec_memb)
  psi_mat <- Matrix(rep(firm_lev_shock, n_scen), 
                    ncol = n_scen, 
                    nrow = n)
  
  
  
  for(k_sec in 1:length(sectors)){
    
    psi_mat_k_sec <- Matrix(0, nrow=n, ncol=n_scen)
    
    
    cat("sector", sectors[k_sec], k_sec, "  sector size",  sectors_firm_count[k_sec], "\n")
    
    # if for a sector there is no unmatched firm or non of the matched firm got a shock skip
    if(length(matched_firm_lev_shock_sector_list[[k_sec]])==0){ next }
    
    #ids of firms within a sector to insert them i.e. unmatched firms
    target_ids <- (1:n)[(firm_sec_memb == sectors[k_sec]) & !firm_matched]
    
    
    # firm_ids_to_shock <- matrix(sapply(1:n_scen, function(x) target_ids[ sample(1:length(target_ids), 
    #                                                        firm_unmatched_tobe_shocked_count[k_sec] )]),
    #        nrow = firm_unmatched_tobe_shocked_count[k_sec],
    #        ncol = n_scen)
    
    # sample a shock from the empirical distribution of a sector for each unmatched firm (can be positive or zero)
    firm_shocks_drawn <- Matrix( sample(matched_firm_lev_shock_sector_list[[k_sec]],
                                        size = length(target_ids)*n_scen,
                                        replace = TRUE), 
                                 nrow = length(target_ids),
                                 ncol = n_scen, sparse = TRUE)
    
    # assign the psi_mat
    psi_mat_k_sec[target_ids, ] <- firm_shocks_drawn
    psi_mat <- psi_mat + psi_mat_k_sec
  }
  
  if(save_psi_mat){
    
    cat("\n", "save psi_mat of size", object.size(psi_mat)/10^6, "MB", '\n')
    cat("start save", as.character(Sys.time()), "\n")
    
    today <- gsub(x = substr(date(), 1,10), pattern = " ", replacement = "")
    
    saveRDS(psi_mat, file = paste0(rsl_wd, 
                                   paste0("psi_mat_for_empirical_shock_extension_", today  ,".RDS")) )
    cat("end save", as.character(Sys.time()), "\n")
    
  }
  
  psi_mat
}