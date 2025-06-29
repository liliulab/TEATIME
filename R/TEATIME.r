###TEATIME####
library(dplyr)
library(RBesT)
#library(MAGOS)
library(stats)
library(strucchange)
library(rcompanion)
library(likelihoodExplore)
library(tidyr)

source(file.path("R", "MAGOS.R"))



#' Rbest to find peak in main cluster
#' @return save Rbest info in files
Rbest_classify<-function(){
  depth=TEATIME$depth
  beta=TEATIME$beta
  output.folder=TEATIME$outputfolder
  output.prefix=TEATIME$outputprefix
  sample_name= TEATIME$id
  main.vaf=TEATIME$main_cluster_vaf

  working.dir <- paste(output.folder, '/', sep='');
  Rbest.file.name <- paste(working.dir, output.prefix, '.rbest.info.txt', sep='');

  possibleError <- tryCatch({
    m<-automixfit(main.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf)/5)), Niter.max=10000)

  },error=function(e){
    e
  })

  if(!inherits(possibleError, "error")){
    a=m["a",]
    b=m["b",]
    mean.a.b=a/(a+b)
    minvaf=min(mean.a.b)
    maxvaf=max(mean.a.b)
    #must_have_one
    mean.a.b.sub<-mean.a.b[abs(mean.a.b-0.5)>0.005]

    len_2=length(mean.a.b)
    len=length(mean.a.b.sub)+1
  }else{

    minvaf=mean(main.vaf)
    maxvaf=mean(main.vaf)
    #
    len_2=1
    len=1
  }

  ##adjust to 0.5
  adjustment <- 0.5 - mean(main.vaf)
  main.vaf <- main.vaf + adjustment
  possibleError <- tryCatch({
    m<-automixfit(main.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf)/5)), Niter.max=10000)


  },error=function(e){
    e
  })
  if(!inherits(possibleError, "error")){
    a=m["a",]
    b=m["b",]
    mean.a.b=a/(a+b)
    minvaf_adj=min(mean.a.b)
    maxvaf_adj=max(mean.a.b)
    len_adj=length(mean.a.b)
    mean.a.b.sub<-mean.a.b[abs(mean.a.b-0.5)>0.005]
    len_adj_2=length(mean.a.b)
  }else{
    minvaf_adj=mean(main.vaf)
    maxvaf_adj=mean(main.vaf)
    len_adj=1
    len_adj_2=1
  }

  Rbest_info<-data.frame(
    samplename=sample_name,
    minvaf=minvaf,
    maxvaf=maxvaf,
    len=len,
    len_2=len_2,
    minvaf_adj=minvaf_adj,
    maxvaf_adj=maxvaf_adj,
    len_adj=len_adj,
    len_adj_2=len_adj_2
  )
  write.table(Rbest_info, Rbest.file.name, sep='\t', row.names=F, quote=F);
  #saveRDS(data,paste0(savepath,'/Rbest.classify.depth',depth,'.rds'))
}

#' Mu Estimation from Peak
#'
#' This function estimates the mu values from the peak of the VAF distribution.
#'
#' @param vaf_set A numeric vector of VAF values.
#' @param p A numeric value representing the theoretical peak.
#' @param num_decimal An integer representing the number of decimal places.
#' @param p_thre A numeric value representing the threshold for p-values.
#' @param celldivlist An optional numeric vector of cell division values. Default is NA, Bac case.
#' @param collect_data_check An optional Dataframe, we prepare it for normal case.
#' @return A data frame containing the estimated mu values and related metrics.
mu_estimation_from_peak<-function(vaf_set,p,num_decimal,p_thre,celldivlist=NULL,collect_data_check=NULL){
  depth=TEATIME$depth
  beta=TEATIME$beta
  Min_Sample_size=6


  collect.data<-data.frame()
  right_most.vaf <- p/2 + (1-p)/(2 * exp(log(2) * beta ))
  abs_diff <- abs(vaf_set - right_most.vaf)
  right_df<-data.frame(vaf=vaf_set,abs_diff=abs_diff)

  right_df$prob<-pbeta(vaf_set,depth*right_most.vaf,depth-depth*right_most.vaf)

  right_df$score<-right_df$prob/max(right_df$prob)-right_df$abs_diff/max(right_df$abs_diff)
  right_df <- right_df[order(-right_df$score), ]
  collect=F
  if(is.null(collect_data_check)){
    collect=T
  }else{
    collect.data<-collect_data_check
  }
  ##is.na means Bac
  ##Fit, Bac should process collect.data
  if(collect){
    collect.data=data.frame()
    if(is.null(celldivlist)){

      for (mu in 3:length(vaf_set)){

        k<-mu
        # Monte Carlo simulation to generate VAFs using the binomial distribution
        right_result<-peak_test(right_df,Min_Sample_size,right_most.vaf,k,depth,num_decimal)
        collect.data <- rbind(collect.data, data.frame(
          cell.div=NA,
          mu_est = mu,
          simu.left.mean = NA,
          simu.right.mean = right_result[1],
          left.mean = NA,
          right.mean = right_result[2],
          left.cd = NA,
          right.cd = right_result[3],
          left.p = NA,
          right.p = right_result[4]
        ))



      }
    }else{
      right_most.vaf=p/2 + (1-p)/(2 * exp(log(2) * beta * 1))
      for (cell.div in celldivlist){

        i_values <- 1:cell.div
        result_vector <- sapply(i_values, function(i) {
          calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
          return(calculated_value)
        })
        result_vector<-c(0.5,p/2,result_vector)
        a<-depth*result_vector
        b<-depth-a
        probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
        df<-data.frame(prob=probs,vaf=vaf_set)
        df<-beta_reassign(df)
        mu1<-max(nrow(df[df$cluster>=3,])/cell.div,1)
        left_df <- df[df$cluster==2+cell.div,]
        right_df <-df[df$cluster==3,]
        if(nrow(right_df)<nrow(left_df)){right_df<-df[df$cluster==2,]}
        if(nrow(right_df)>mu1 & mu1>3){
          left_most_vaf=p/2 + (1-p)/(2 * exp(log(2) * beta * cell.div))

          left_df$abs_diff <- abs(left_df$vaf - left_most_vaf)
          left_df$prob<-pbeta(left_df$vaf,depth*left_most_vaf,depth-depth*left_most_vaf)
          left_df$score<-left_df$prob/max(left_df$prob)-left_df$abs_diff/max(left_df$abs_diff)
          left_df <- left_df[order(-left_df$score), ]

          right_df$abs_diff <- abs(right_df$vaf - right_most.vaf)
          right_df$prob<-pbeta(right_df$vaf,depth*right_most.vaf,depth-depth*right_most.vaf)
          right_df$score<-right_df$prob/max(right_df$prob)-right_df$abs_diff/max(right_df$abs_diff)
          right_df <- right_df[order(-right_df$score), ]


          k<-round(mu1)
          # Monte Carlo simulation to generate VAFs using the binomial distribution
          #df,Min_Sample_size,vaf,k,depth,num_decimal
          left_result<-peak_test(left_df,Min_Sample_size,left_most_vaf,k,depth,num_decimal)
          right_result<-peak_test(right_df,Min_Sample_size,right_most.vaf,k,depth,num_decimal)
          # Print p-value
          collect.data <- rbind(collect.data, data.frame(
            cell.div=cell.div,
            mu_est = mu1,
            simu.left.mean = left_result[1],
            simu.right.mean = right_result[1],
            left.mean = left_result[2],
            right.mean = right_result[2],
            left.cd = left_result[3],
            right.cd = right_result[3],
            left.p = left_result[4],
            right.p = right_result[4]
          ))


        }


      }
    }
  }




  Temp_keep<-collect.data
  if (all(is.na(Temp_keep$left.p))) {
    if (all(is.na(Temp_keep$right.p))) {
      collect.data<-collect.data
    }else{
      collect.data<-Temp_keep[Temp_keep$right.p>p_thre,]
      }

    if(nrow(collect.data)==0){collect.data=Temp_keep}


    collect.data$right_rank<-rank(collect.data$right.cd)

    n_rows <- nrow(collect.data)
    if(collect){
      top_percentage <- if (n_rows > 4) 0.25 else if (n_rows > 2) 0.5 else 1
    }else{
      top_percentage <- if (n_rows > 2) 0.5 else 1
    }
    top_n <- floor(top_percentage * n_rows)
    sorted_df <- collect.data[order(collect.data$right_rank), ]
    top_pick <- head(sorted_df, top_n)

  }else{
    collect.data<-Temp_keep[Temp_keep$left.p>p_thre | Temp_keep$right.p>p_thre,]
    collect.data<-collect.data[complete.cases(collect.data[ , 3]),]
    collect.data$left.cd<-ifelse(is.na(collect.data$left.cd),10,collect.data$left.cd)
    n_rows <- nrow(collect.data)
    top_percentage <- if (n_rows > 2) 0.5 else 1
    top_n <- floor(top_percentage * n_rows)

    if (all(collect.data$left.cd == 10)){
      collect.data$right_rank<-rank(collect.data$right.cd)
      sorted_df <- collect.data[order(collect.data$right_rank), ]
      top_pick <- head(sorted_df, top_n)

    }else{
      collect.data$left_rank<-rank(collect.data$left.cd)
      collect.data$right_rank<-rank(collect.data$right.cd)
      sorted_df <- collect.data[order(collect.data$left_rank), ]
      top_left <- head(sorted_df, top_n)

      sorted_df <- collect.data[order(collect.data$right_rank), ]
      top_right <- head(sorted_df, top_n)

      top_pick<-rbind(top_left,top_right)
      top_pick<-top_pick[!duplicated(top_pick),]

    }


  }


  top_pick_p<-Compare_real_simu_peak(top_pick,p,vaf_set,num_decimal,right_df=right_df,right_most.vaf=right_most.vaf)

  top_pick_p <- as.data.frame(t(top_pick_p))
  top_pick_p$V1<-ifelse(top_pick_p$V1<0.05,0,top_pick_p$V1)
  top_pick$r1<-top_pick_p$V1

  top_pick$loglike<-top_pick_p$V2
  top_pick$aic<-top_pick_p$V3
  top_pick$bic<-top_pick_p$V4

  top_pick$top_pick_score<-top_pick$r1+1/top_pick$right_rank
  if(is.null(celldivlist)){
    top_pick<- top_pick[order(-top_pick$top_pick_score), ]
  }else{
    if(nrow(top_pick)>1){
      top_pick<- top_pick[top_pick$top_pick_score>median(top_pick$top_pick_score), ]
    }
    top_pick <- top_pick[order(top_pick$right_rank), ]
  }


  pick_mu_cell.div<-top_pick[1:min(10,nrow(top_pick)),] ### First 10 guess that is good

  return(pick_mu_cell.div)

}


#' Peak compare, real vs simu
#' @param peakdata A dataframe contain cell.div, guess mu, etc. If have cell.div Fit case. Else, Bac case.
#' @param p proportion.
#' @param vaf_set Vaf seq used for analysis.
#' @param num_decimal An integer representing the number of decimal places.
#' @param right_df Used in Bac case.
#' @param right_most.vaf Used in Bac case.
#' @param mu_small Used in if mu should be very small (<3).
#' @param bac Used in Bac case.
#' @return A list of comparison.
Compare_real_simu_peak<-function(peakdata,p,vaf_set,num_decimal,right_df=NULL,right_most.vaf=NA,mu_small=F,bac=F){
  depth=TEATIME$depth
  beta=TEATIME$beta
  top_pick_p <- apply(peakdata, 1, function(row) {
    cell_div_value <- row['cell.div']
    mu_value <- row['mu_est']
    if(is.na(cell_div_value)){
      largest_right_values <- right_df[1:mu_value,]$vaf
      p_v_1_values = c()
      for (try in 1:50){
        # Simulate VAFs for each element in vaf_list and round them
        all_simulated_vafs_list <- round(rbinom(round(mu_value), depth, right_most.vaf) / depth, num_decimal)
        # Concatenate all vectors into a single vector
        all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
        p_v_1=wilcox.test(largest_right_values, all_simulated_vafs)$p.value
        p_v_1_values = c(p_v_1_values, p_v_1)
      }
      mean_p_v_1 = mean(p_v_1_values)

      log_likelihood = log_likelihood_mixture(largest_right_values, right_most.vaf,depth) # example
      num_params = 1
      # Sample size
      sample_size = length(largest_right_values)
      AIC_value = compute_AIC(log_likelihood, num_params)
      BIC_value = compute_BIC(log_likelihood, num_params, sample_size)
      return(c(mean_p_v_1,log_likelihood,AIC_value,BIC_value))

    }else{

      if (!bac) {
        i_values <- if (mu_small) 1:20 else 1:cell_div_value
        result_vector <- sapply(i_values, function(i) {
          calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
          return(calculated_value)
        })
        result_vector<-c(0.5,p/2,result_vector)
        a<-depth*result_vector
        b<-depth-a
        probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
        df<-data.frame(prob=probs,vaf=vaf_set)
        df<-beta_reassign(df)
        data<-df[df$cluster>=2,]$vaf
      }else{
        data<-vaf_set
      }

      vaf_list<-c(p/2)
      for(i in 1:cell_div_value){
        cell.num<-exp(log(2)*beta*i)
        vaf_list<-c(vaf_list,p/2 + (1-p)/(2*cell.num))
      }
      p_v_1_values = c()
      for (try in 1:50){
        # Simulate VAFs for each element in vaf_list and round them
        all_simulated_vafs_list <- sapply(vaf_list, function(vaf) {
          round(rbinom(round(mu_value), depth, vaf) / depth, num_decimal)
        })

        # Concatenate all vectors into a single vector
        all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
        p_v_1=wilcox.test(data, all_simulated_vafs)$p.value
        p_v_1_values = c(p_v_1_values, p_v_1)
      }
      mean_p_v_1 = mean(p_v_1_values)
      log_likelihood = log_likelihood_mixture(data, vaf_list,depth) # example
      num_params = length(vaf_list)
      sample_size = length(data)
      AIC_value = compute_AIC(log_likelihood, num_params)
      BIC_value = compute_BIC(log_likelihood, num_params, sample_size)
      return(c(mean_p_v_1,log_likelihood,AIC_value,BIC_value))

    }

  })
  return(top_pick_p)
}



#' Peak compare, real vs simu
#' @param df A numeric vector of VAF values.
#' @param Min_Sample_size Check if we need bootstrap.
#' @param vaf right most vaf.
#' @param k mu values.
#' @param depth depth sequence.
#' @param num_decimal An integer representing the number of decimal places.
#' @param temp_check if NA is required
#' @return A list of comparison.
peak_test<-function(df,Min_Sample_size,vaf,k,depth,num_decimal,temp_check=F){
  depth=TEATIME$depth
  beta=TEATIME$beta
  if(temp_check){
    temp_k<-0
  }else{
    temp_k<-k
  }
  if(nrow(df)<temp_k){
    return(c(NA,NA,NA,NA))
  }else{
    #Define start.div and end.div
    n_simulations <- 100000
    simulated_vafs.right <- round(rbinom(n_simulations, depth, vaf) / depth,num_decimal)

    # For left_df
    if(k<Min_Sample_size){
      largest_right_values <- generate_bootstrap_samples(df[1:k,]$vaf,Min_Sample_size,num_decimal)
    }else{

      largest_right_values <- df[1:k,]$vaf
    }
    largest_right_values <- na.omit(largest_right_values)
    largest_right_values <- as.vector(largest_right_values)

    w.right.result <- wilcox.test(simulated_vafs.right, largest_right_values)

    data1_long <- data.frame(value = simulated_vafs.right, group = "data1")
    data2_long <- data.frame(value = largest_right_values, group = "data2")
    combined_data <- rbind(data1_long, data2_long)

    cd.right=cliffDelta(value ~ group, data = combined_data)
    return(c(mean(simulated_vafs.right),mean(largest_right_values),abs(cd.right),w.right.result$p.value))
  }



}


#' Generate Bootstrap Samples
#'
#' This function generates bootstrap samples from the original data.
#'
#' @param original_data A numeric vector of the original data.
#' @param n_samples An integer representing the number of bootstrap samples to generate.
#' @param num_decimal An integer representing the number of decimal places to round the results.
#' @return A numeric vector containing the bootstrap samples.
generate_bootstrap_samples <- function(original_data, n_samples,num_decimal) {
  # Initialize an empty vector to store bootstrap samples
  depth=TEATIME$depth
  beta=TEATIME$beta
  bootstrap_samples <- numeric(n_samples)
  # Generate n_samples bootstrap samples
  for (i in 1:n_samples) {
    resample <- sample(original_data, size = length(original_data), replace = TRUE)
    bootstrap_samples[i] <- mean(resample)
  }

  return(bootstrap_samples)
}

#' Log Likelihood Mixture
#'
#' This function calculates the log likelihood of a mixture model for VAF data.
#'
#' @param data A numeric vector of VAF values.
#' @param p_vec A numeric vector of probabilities.
#' @param depth A numeric value representing the depth.
#' @return A numeric value representing the total log likelihood.
log_likelihood_mixture <- function(data, p_vec,depth) {

  depth=TEATIME$depth
  beta=TEATIME$beta
  a<-depth*p_vec
  b<-depth-a
  probs <- sapply(1:length(a), function(i) dbeta(data, a[i], b[i]))
  df<-data.frame(prob=probs,vaf=data)

  if(length(p_vec)>1){
    # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
    df<-beta_reassign(df)
    df <- df %>%
      mutate(p = p_vec[cluster])
  }else{ df$p<-p_vec}

  df <- df %>%
    mutate(log_likelihood = likbeta(x = vaf,
                                    shape1 = p*depth, shape2 = depth-p))

  total_log_likelihood <- sum(df$log_likelihood)

  return(total_log_likelihood)
}

#' Compute AIC
#'
#' This function computes the Akaike Information Criterion (AIC) for a model.
#'
#' @param log_likelihood A numeric value representing the log likelihood of the model.
#' @param num_params An integer representing the number of parameters in the model.
#' @return A numeric value representing the AIC.
compute_AIC <- function(log_likelihood, num_params) {
  return(-2 * log_likelihood + 2 * num_params)
}

#' Compute BIC
#'
#' This function computes the Bayesian Information Criterion (BIC) for a model.
#'
#' @param log_likelihood A numeric value representing the log likelihood of the model.
#' @param num_params An integer representing the number of parameters in the model.
#' @param sample_size An integer representing the sample size.
#' @return A numeric value representing the BIC.
compute_BIC <- function(log_likelihood, num_params, sample_size) {
  return(-2 * log_likelihood + num_params * log(sample_size))
}

#' This function estimates the p values, normal case.
#'
#' @param second.vaf cluster near clonal vafs.
#' @param mean.a.b peak contained in clonal vafs.
#' @param upper_clonal_vaf peak that is near the second cluster(mixed maybe).
#' @param clonal.vaf.left all vafs that belongs to upper_clonal_vaf.
#' @return Guess about p with clear, meaning if clonal.vaf>1 peak, take care.
Find_p_process<-function(second.vaf,mean.a.b,upper_clonal_vaf,clonal.vaf.left){
  depth=TEATIME$depth
  beta=TEATIME$beta
  m<-automixfit(second.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(second.vaf)/3)), Niter.max=10000)
  a1=m["a",]
  b1=m["b",]
  mean.a.b1=a1/(a1+b1)
  second_clone_left_vaf=min(mean.a.b1)
  clear=ifelse(length(mean.a.b)>1,FALSE,TRUE)

  left_limit <- (min(second.vaf) + second_clone_left_vaf) / 2



  if(length(mean.a.b)>1){
    right_limit=max(((upper_clonal_vaf+max(clonal.vaf.left))/2*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)/2,(max(second.vaf)+second_clone_left_vaf)/2)
  }else{
    #Try from min to a left_limit to right_limit
    right_limit=((max(second.vaf)+upper_clonal_vaf)/2*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)/2
    if(right_limit>0.45 | right_limit<left_limit){
      right_limit=(max(second.vaf)+second_clone_left_vaf)/2
    }

  }
  ##if too few point
  if(length(second.vaf[second.vaf>right_limit])<=6){
    sorted_vaf <- sort(second.vaf, decreasing = TRUE)
    right_limit <- sorted_vaf[max(round(0.25*length(sorted_vaf)),6)]
  }

  left <- left_limit
  middle <- (left_limit + right_limit) / 2
  right <- right_limit
  # Quantile points sequence
  my_seq <- c(left,middle,right) ##3 guess
  my_seq<-my_seq[my_seq>0.5/(1+exp(log(2)*beta))]
  return(list(my_seq = my_seq, clear = clear))
}

#' Reassign VAF Groups Based on Beta Distribution Probabilities
#'
#' This function processes a data frame of VAF values and their associated probabilities,
#' reassigning each VAF to a cluster based on the maximum probability. It calculates the
#' frequency of each VAF and distributes the VAFs into clusters accordingly.
#'
#' @param df A data frame containing VAF values and their associated probabilities.
#'        The data frame must include columns named `vaf` and `prob.<n>` where `<n>`
#'        represents different clusters.
#' @return A data frame with columns `vaf`, `cluster`, and `freq`, where `cluster` is the
#'         reassigned cluster number and `freq` is the frequency count of VAF values in each cluster.
beta_reassign<-function(df){
  depth=TEATIME$depth
  beta=TEATIME$beta
  # Process each VAF group
  # Grouping and summarizing
  freq_df <- df %>%
    group_by(vaf) %>%
    summarise(n = n(), .groups = 'drop')
  # Calculate proportions and join back with frequency data
  proportion_df <- df %>%
    group_by(vaf) %>%
    summarise(across(starts_with("prob"), sum), .groups = 'drop') %>%
    rowwise() %>%
    mutate(
      total = sum(c_across(starts_with("prob"))),
      across(starts_with("prob"), ~round(. / total * freq_df[freq_df$vaf == vaf, ]$n, 0))
    ) %>%
    left_join(freq_df, by = "vaf")
  final_df <- proportion_df %>%
    pivot_longer(cols = starts_with("prob"), names_to = "cluster", values_to = "freq") %>%
    mutate(cluster = as.integer(gsub("prob\\.", "", cluster))) %>% # Convert cluster to integer
    #{print(head(.)); .} %>%  # Print the head of the dataframe
    uncount(freq, .remove = FALSE)

  df_not_in_final<-df[!(df$vaf %in% final_df$vaf),]
  if(nrow(df_not_in_final)>0){
    result <- df_not_in_final %>%
      rowwise() %>%
      mutate(
        cluster = which.max(c_across(starts_with("prob")))
      ) %>%
      ungroup()
    result$freq<-1
    final_df<-rbind(final_df[,c("vaf","cluster","freq")],result[,c("vaf","cluster","freq")])
  }else{
    final_df<-final_df[,c("vaf","cluster","freq")]
  }

  return(final_df)
}

#' Given a guessed P, find mu from normal case
#'

#'
#' @param clear if we need further process
#' @param give.vaf Given P, coresponding Vafs
#' @param upper_clonal_vaf peak that is near the second cluster(mixed maybe).
#' @param clonal.vaf.left all vafs that belongs to upper_clonal_vaf.
#' @param second.vaf cluster near clonal vafs.
#' @param num_decimal An integer representing the number of decimal places to round the results.
#' @param p_thre cutoff
#' @return A data frame calculated mu.
Iterate_P_optimize<-function(clear,give.vaf,upper_clonal_vaf,clonal.vaf.left,second.vaf,num_decimal,p_thre=0.05){
  depth=TEATIME$depth
  beta=TEATIME$beta
  Min_Sample_size=6
  i_values <- 1:20

  p=give.vaf*2
  reliable=1

  suppose_right_vaf=p/2+(1-p)/(2*exp(log(2)*beta))
  first_vaf_list<-c()
  if(!clear){

    a <- if (suppose_right_vaf > upper_clonal_vaf) {
      c(depth * 0.5, depth * suppose_right_vaf)
    } else {
      c(depth * upper_clonal_vaf, depth * suppose_right_vaf)
    }
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(clonal.vaf.left, a[i], b[i]))
    df<-data.frame(prob=probs,vaf=clonal.vaf.left)
    df<-beta_reassign(df)
    first_vaf_list<-df[df$cluster==2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
  }

  vaf_set<-c(first_vaf_list,second.vaf[second.vaf>give.vaf])

  result_vector <- sapply(i_values, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  a<-depth*result_vector
  b<-depth-a
  probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
  df<-data.frame(prob=probs,vaf=vaf_set)
  df<-beta_reassign(df)
  df<-df[!duplicated(df),]

  df_count_freq <- df %>%
    group_by(cluster) %>%
    summarise(count = sum(freq), .groups = 'drop')
  if(min(df_count_freq$cluster)>1){
    df_count_freq <- rbind(data.frame(cluster = 1, count = 0), df_count_freq)
    # Sort the dataframe by cluster
    df_count_freq <- df_count_freq[order(df_count_freq$cluster), ]
  }

  max_mu=df_count_freq[1,]$count+df_count_freq[2,]$count

  if(nrow(df_count_freq)>6 & df_count_freq[2,"count"]/df_count_freq[3,"count"]<=4){
    df_count_freq<-df_count_freq[2:3,]
  }else{
    df_count_freq<-df_count_freq[2,]
  }
  min_mu=max(3,min(df_count_freq$count/4))
  total_count_temp=length(vaf_set)
  start.div=round(total_count_temp/max_mu)
  end.div<-round(total_count_temp/min_mu)

  cat('cell.div for normal case...\n');
  cat('   ', start.div,end.div,'...'); flush.console();
  #fit_slope_check<-function(p,start.div,num_decimal,end.div=NA,vaf_set=NA)
  collect.data<-fit_slope_check(p,start.div,num_decimal,end.div=end.div,vaf_set_check=vaf_set)
  collect.data<-collect.data[collect.data$mu_real>3,]
  #print('first done')
  cat('Slope method Done for normal case...\n');
  cat('   ', nrow(collect.data),'...'); flush.console();


  if(nrow(collect.data[abs(collect.data$z_score)<=1.96,])>0){
    ##if over 20 ( very hard to distinguish)
    collect.data<-collect.data[abs(collect.data$z_score)<=1.96,]

    first_div_vaf_data<-df[df$cluster==df_count_freq[1,]$cluster,]
    sec.div.vaf<- p/2 + (1-p)/(2 * exp(log(2) * beta * 2))
    first_div_vaf_data$abs_diff <- abs(first_div_vaf_data$vaf - sec.div.vaf)
    first_div_vaf_data$prob<-pbeta(first_div_vaf_data$vaf,depth*sec.div.vaf,depth-depth*sec.div.vaf)
    first_div_vaf_data$score<-first_div_vaf_data$prob/max(first_div_vaf_data$prob)-first_div_vaf_data$abs_diff/max(first_div_vaf_data$abs_diff)
    right_df <- first_div_vaf_data[order(-first_div_vaf_data$score), ]
    celldivlist<-collect.data$cell.div
    # Sort the data frame by absolute differences
    collect.data1<-data.frame()
    for (cell.div in celldivlist){
      mu_suppose=length(vaf_set)/cell.div
      k<-round(mu_suppose)
      #peak_test<-function(df,Min_Sample_size,vaf,k,depth,num_decimal)
      right_result<-peak_test(df=right_df,Min_Sample_size=Min_Sample_size,vaf=sec.div.vaf,k=k,depth=depth,num_decimal=num_decimal,temp_check=F)
      collect.data1 <- rbind(collect.data1, data.frame(
        cell.div=cell.div,
        mu_est = mu_suppose,
        simu.right.mean = right_result[1],
        right.mean = right_result[2],
        right.cd = right_result[3],
        right.p = right_result[4]
      ))

    }
    #mu_estimation_from_peak<-function(vaf_set,p,num_decimal,p_thre,celldivlist=NA,collect.data=NULL)
    #collect.data2<-Second_run_test(first_vaf_list,second.vaf,give.vaf,sec.div.vaf,right_df,dynamic,celldivlist,p,beta,depth,num_decimal,Min_Sample_size=6,p_thre)
    collect.data2<-mu_estimation_from_peak(vaf_set=vaf_set,p=p,num_decimal=num_decimal,p_thre=p_thre,celldivlist=celldivlist,collect_data_check=collect.data1)
    cat('Peak method Done for normal case...\n');
    cat('   ', nrow(collect.data2),'...'); flush.console();

    flush.console();

    range1 <- ifelse(nrow(collect.data2) == 1, collect.data2$mu_est - 1, min(collect.data2$mu_est))
    range2 <- ifelse(nrow(collect.data2) == 1, collect.data2$mu_est + 1, max(collect.data2$mu_est))

    over12<-collect.data[collect.data$mu_real>=range1 & collect.data$mu_real<=range2,]

    if(nrow(over12)>0){
      overlap.pick<-over12
    }else{
      overlap.pick<-NULL
    }
    #overlap.pick<-ifelse(nrow(over12)>0,over12,NULL)


    if(nrow(df_count_freq)>1){
      max_mu=df_count_freq[df_count_freq$cluster==df_count_freq[2,]$cluster,]$count
      min_mu=max(1,df_count_freq[df_count_freq$cluster==3,]$count/2)

      start.div=round(total_count_temp/max_mu)
      end.div<-round(total_count_temp/min_mu)
      first_div_vaf_data<-df[df$cluster==3,]
      sec.div.vaf<- p/2 + (1-p)/(2 * exp(log(2) * beta * 3))
      first_div_vaf_data$abs_diff <- abs(first_div_vaf_data$vaf - sec.div.vaf)
      first_div_vaf_data$prob<-pbeta(first_div_vaf_data$vaf,depth*sec.div.vaf,depth-depth*sec.div.vaf)
      first_div_vaf_data$score<-first_div_vaf_data$prob/max(first_div_vaf_data$prob)-first_div_vaf_data$abs_diff/max(first_div_vaf_data$abs_diff)
      right_df <- first_div_vaf_data[order(-first_div_vaf_data$score), ]


      collect.data1<-data.frame()
      for (cell.div in celldivlist){
        mu_suppose=length(vaf_set)/cell.div
        k<-round(mu_suppose)
        #peak_test<-function(df,Min_Sample_size,vaf,k,depth,num_decimal)
        right_result<-peak_test(df=right_df,Min_Sample_size=Min_Sample_size,vaf=sec.div.vaf,k=k,depth=depth,num_decimal=num_decimal)
        collect.data1 <- rbind(collect.data1, data.frame(
          cell.div=cell.div,
          mu_est = mu_suppose,
          simu.right.mean = right_result[1],
          right.mean = right_result[2],
          right.cd = right_result[3],
          right.p = right_result[4]
        ))

      }
      collect.data3<-mu_estimation_from_peak(vaf_set=vaf_set,p=p,num_decimal=num_decimal,p_thre=p_thre,celldivlist=celldivlist,collect_data_check=collect.data1)
      #collect.data3<-Second_run_test(first_vaf_list,second.vaf,give.vaf,sec.div.vaf,right_df,dynamic,celldivlist,p,beta,depth,num_decimal,Min_Sample_size=6,p_thre)

      range1 <- ifelse(nrow(collect.data3) == 1, collect.data3$mu_est - 1, min(collect.data3$mu_est))
      range2 <- ifelse(nrow(collect.data3) == 1, collect.data3$mu_est + 1, max(collect.data3$mu_est))

      over23<-collect.data[collect.data$mu_real>=range1 & collect.data$mu_real<=range2,]
      overlap_values <- intersect(over12$mu_real, over23$mu_real)

      if(nrow(over23)>0 | nrow(over12)>0){
        if(length(overlap_values) > 0){
          overlap.pick<-over23[over23$mu_real %in% overlap_values, ]
        }else{
          overlap.pick<-rbind(over12,over23)
        }
      }
    }
    ##check three overlap
    if(is.null(overlap.pick)){
      overlap.pick<-collect.data
    }
    overlap.pick<-na.omit(overlap.pick)

  }else{
    overlap.pick<-collect.data
    reliable=0
  }

  #vaf_set=c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  ##if mu too small

  collect.data4<-mu_find_small(vaf_set,p)
  mu_small_selection=collect.data4[collect.data4$p_value>0.05 & collect.data4$wx_p>0.05 ,]

  pick.cell.div=overlap.pick$cell.div
  pick.mu=overlap.pick$mu_real
  z_score=overlap.pick$z_score

  overlap.pick$mu_est=overlap.pick$mu_real
  top_pick_p<-Compare_real_simu_peak(overlap.pick,p,vaf_set,num_decimal)
  top_pick <- as.data.frame(t(top_pick_p))
  #overlap.pick$loglike<-
  overlap.pick$loglike<-top_pick$V2
  overlap.pick$aic<-top_pick$V3
  overlap.pick$bic<-top_pick$V4
  pick.log<-overlap.pick$loglike
  pick.bic<-overlap.pick$bic
  pick.aic<-overlap.pick$aic


  if(nrow(mu_small_selection)>0){
    num_df<-length(vaf_set)
    mu_small_selection$cell.div<-round(num_df/mu_small_selection$mu_est)
    small_mu_pick<-Compare_real_simu_peak(peakdata =mu_small_selection,p=p,vaf_set=vaf_set,num_decimal=num_decimal,bac = T)
    small_mu_pick <- as.data.frame(t(small_mu_pick))
    small_mu_pick$V1<-ifelse(small_mu_pick$V1<0.05,0,small_mu_pick$V1)
    mu_small_selection$r1<-small_mu_pick$V1
    mu_small_selection$loglike<-small_mu_pick$V2
    mu_small_selection$aic<-small_mu_pick$V3
    mu_small_selection$bic<-small_mu_pick$V4

    mu_small_selection <- mu_small_selection[which(mu_small_selection$bic == min(mu_small_selection$bic)), ]

    if(mu_small_selection$bic<min(data$bic)*0.5){
      pick.mu=mu_small_selection$mu_est
      pick.bic=mu_small_selection$bic
      pick.aic=mu_small_selection$aic
      pick.log=mu_small_selection$loglike
      z_score<-rep(1,length(mu_small_selection$aic))
    }

  }

  normal_pick<-data.frame(
    cell.div=pick.cell.div,
    mu=pick.mu,
    loglike=pick.log,
    bic=pick.bic,
    aic=pick.aic,
    p=p,
    z_score=z_score,
    reliable=reliable
  )
  return(normal_pick)

}


#' Estimate mu from main cluster, from eq in paper, normal case
#' @param p_thre Check if input is from magos analysis
#' @return A list containing the analysis results.
Run.normal.maincluster = function(p_thre=0.01) {
  num_decimal <- nchar(as.character(TEATIME$depth))
  depth=TEATIME$depth
  beta=TEATIME$beta
  vafdata.summary<-TEATIME$vafdata_summary
  vafdata.summary.filter<-TEATIME$vafdata_summary_filter
  vafdata<-TEATIME$vafdata
  second.vaf<-TEATIME$second_cluster_vaf$normal
  color <- vafdata.summary$colors[vafdata.summary$max %in% vafdata.summary.filter[1, "max"] | vafdata.summary$min %in% vafdata.summary.filter[1, "min"]]

  clonal.vaf<-vafdata[vafdata$colors %in% color,]$vaf.1 ##this should be same as TEATIME$main_cluster_vaf
  m<-automixfit(clonal.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(clonal.vaf)/3)), Niter.max=10000)
  a=m["a",]
  b=m["b",]
  mean.a.b=a/(a+b)
  upper_clonal_vaf=min(mean.a.b)
  approx_vaf_index=which.min(mean.a.b)
  probs <- sapply(1:length(a), function(i) dbeta(clonal.vaf, a[i], b[i]))
  df<-data.frame(prob=probs,vaf=clonal.vaf)
  df<-beta_reassign(df)

  clonal.vaf.left<-df[df$cluster==approx_vaf_index,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
  second.try <- nrow(vafdata.summary.filter) > 1
  cat('proportion needs run twice?...\n');
  cat('   ', second.try ,'...'); flush.console();
  if(nrow(vafdata.summary.filter)>1){
    second.vaf=second.vaf
  }else{
    temp_vaf<-(upper_clonal_vaf+0.5)/2
    second.vaf<-clonal.vaf[clonal.vaf<temp_vaf]
  }
  #Find_p_process<-function(second.vaf,mean.a.b,upper_clonal_vaf,clonal.vaf.left)
  seq_data<-Find_p_process(second.vaf,mean.a.b,upper_clonal_vaf,clonal.vaf.left)
  my_seq<-seq_data$my_seq
  clear<-seq_data$clear
  if(length(my_seq)<=1 & second.try==T ){
    temp_vaf<-(upper_clonal_vaf+0.5)/2
    second.vaf<-clonal.vaf[clonal.vaf<temp_vaf]
    seq_data<-Find_p_process(second.vaf,mean.a.b,upper_clonal_vaf,clonal.vaf.left)
    my_seq<-seq_data$my_seq
    clear<-seq_data$clear
  }

  all.p.data.final <- data.frame(
    cell.div = NA, mu = NA, loglike = NA, bic = NA, aic = NA, p = NA,
    z_score = NA, reliable = NA, z_score1 = NA, bicrank = NA, score = 0, clear = NA
  )
  if(length(my_seq)<=0){
    cat('Unable to find point to guess P for normal case...\n');flush.console();

  }else{
    cat('P for normal case...\n');
    cat('   ', my_seq,'...'); flush.console();
    all.p.data<-data.frame()
    start_time <- Sys.time()
    for(give.vaf in my_seq){
      tryCatch({
        opt.data<-Iterate_P_optimize(clear,give.vaf,upper_clonal_vaf,clonal.vaf.left,second.vaf,num_decimal,p_thre)
        if (!is.null(opt.data)) {
          all.p.data <- rbind(all.p.data, opt.data)
        }
      }, error = function(e) {
        e
      })

    }
    end_time <- Sys.time()
    time_durations<-end_time-start_time
    cat('Fast run time duration...\n');
    cat('   ', time_durations,'...'); flush.console();

    if(nrow(all.p.data)>0){
      all.p.data$mu<-round(all.p.data$mu,num_decimal)
      all.p.data <- if (nrow(all.p.data[abs(all.p.data$z_score) <= 1.96, ]) >= 2) {
        all.p.data[abs(all.p.data$z_score) <= 1.96, ]
      } else {
        all.p.data
      }

      if(nrow(all.p.data)>=3){
        # Calculate the mean and standard deviation without the maximum value
        mean_filtered <- median(all.p.data$mu)
        sd_filtered <-mad(all.p.data$mu)
        if(sd_filtered>0){
          all.p.data$z_score1 <- abs((all.p.data$mu - mean_filtered) / sd_filtered)
          all.p.data <- all.p.data[all.p.data$z_score1 < 2,]
        }
      }
      all.p.data$z_score1<-1
      all.p.data$z_score<-abs(all.p.data$z_score)
      all.p.data$bicrank<-rank(all.p.data$bic)
      all.p.data$score<-1/all.p.data$bicrank+(1-all.p.data$z_score)*0.25
      all.p.data$clear<-clear
      all.p.data.final<-all.p.data
      #print(all.p.data)
    }
  }

  return(all.p.data.final)



}


#' Estimate mu from main cluster, from eq in paper, Bac case
#' @param p_thre Check if input is from magos analysis
#' @return A list containing the analysis results.
Run.bac.maincluster = function(p_thre=0.01) {

  num_decimal <- nchar(as.character(TEATIME$depth))
  second.vaf<-TEATIME$second_cluster_vaf$bac
  main.vaf<-TEATIME$main_cluster_vaf
  depth=TEATIME$depth
  beta=TEATIME$beta
  m<-automixfit(second.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(second.vaf)/3)), Niter.max=10000)
  a=m["a",]
  b=m["b",]
  mean.a.b=a/(a+b)
  left_most_vaf=min(mean.a.b)
  right_most.vaf=max(mean.a.b)
  right_vaf_index=which.max(mean.a.b)
  keep.right.vaf<-right_most.vaf
  probs <- sapply(1:length(a), function(i) dbeta(second.vaf, a[i], b[i]))
  df<-data.frame(prob=probs,vaf=second.vaf)
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  if(length(mean.a.b)>1){
    df<-beta_reassign(df)
    second.update.vaf<-df[df$cluster==right_vaf_index,]$vaf
  }else{
    second.update.vaf<-df$vaf
  }
  cat('Find right most VAF...\n');
  cat('   ', mean.a.b, '...'); flush.console();
  while (length(mean.a.b)>1){
    possibleError <- tryCatch({
      m<-automixfit(second.update.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(second.update.vaf)/3)), Niter.max=10000)

    },error=function(e){
      e
    })
    if(!inherits(possibleError, "error")){
      #m<-automixfit(second.update.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(second.update.vaf)/3)), Niter.max=10000)
      a=m["a",]
      b=m["b",]
      mean.a.b=a/(a+b)
      left_most_vaf=min(mean.a.b)
      right_most.vaf=max(mean.a.b)
      if(length(mean.a.b)<2){break}
      right_vaf_index=which.max(mean.a.b)
      probs <- sapply(1:length(a), function(i) dbeta(second.update.vaf, a[i], b[i]))
      df<-data.frame(prob=probs,vaf=second.update.vaf)
      df<-beta_reassign(df)
      second.update.vaf<-df[df$cluster==right_vaf_index,]$vaf
      if (length(second.update.vaf) == 0 || (length(unique(second.update.vaf)) == 1)) {break}
    }else{break}
  }
  right_most.vaf<-max(right_most.vaf,keep.right.vaf)
  p=(right_most.vaf*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)
  p <- ifelse(p < 0, 1 - 2 * exp(log(2) * beta) * right_most_vaf, p)
  p <- ifelse(p<0, 0.01,p)


  righta=c(0.5*depth,right_most.vaf*depth)
  rightb=depth-righta
  probs <- sapply(1:length(righta), function(i) dbeta(main.vaf, righta[i], rightb[i]))
  df<-data.frame(prob=probs,vaf=main.vaf)
  df<-beta_reassign(df)

  vaf_set<-c(df[df$cluster==2,]$vaf,second.update.vaf)
  start.div<-ifelse(left_most_vaf == right_most.vaf,1,2)
  cat('Start Running Peak Method For Bac Case...\n');
  cat('   ', length(vaf_set), '...'); flush.console();

  pick_mu_cell.div=mu_estimation_from_peak(vaf_set,p,num_decimal,p_thre,celldivlist=NULL)



  vaf_set=c(df[df$cluster==2,]$vaf,second.vaf)
  cat('Checking if mu is small...\n');flush.console();
  result_df<-mu_find_small(vaf_set,p)
  mu_small_selection=result_df[result_df$p_value>0.05 & result_df$wx_p>0.05 ,]
  pick.mu=pick_mu_cell.div$mu_est
  pick.bic=pick_mu_cell.div$bic
  pick.aic=pick_mu_cell.div$aic
  pick.log=pick_mu_cell.div$loglike


  if(nrow(mu_small_selection)>0){
    num_df<-length(vaf_set)
    mu_small_selection$cell.div<-round(num_df/mu_small_selection$mu_est)
    small_mu_pick<-Compare_real_simu_peak(peakdata =mu_small_selection,p=p,vaf_set=vaf_set,num_decimal=num_decimal,bac = T)

    small_mu_pick <- as.data.frame(t(small_mu_pick))
    small_mu_pick$V1<-ifelse(small_mu_pick$V1<0.05,0,small_mu_pick$V1)
    mu_small_selection$r1<-small_mu_pick$V1
    mu_small_selection$loglike<-small_mu_pick$V2
    mu_small_selection$aic<-small_mu_pick$V3
    mu_small_selection$bic<-small_mu_pick$V4


    mu_small_selection <- mu_small_selection[which(mu_small_selection$bic == min(mu_small_selection$bic)), ]

    if(mu_small_selection$bic<min(data$bic)*0.5){
      pick.mu=mu_small_selection$mu_est
      pick.bic=mu_small_selection$bic
      pick.aic=mu_small_selection$aic
      pick.log=mu_small_selection$loglike
    }
  }


  bac_pick<-data.frame(
    cell.div=pick.mu,
    mu=pick.mu,
    loglike=pick.log,
    bic=pick.bic,
    aic=pick.aic,
    p=p
  )

  return(bac_pick)


}
#' Estimate mu from main cluster, from eq in paper, Fit case
#' @param p_thre Check if input is from magos analysis
#' @return A list containing the analysis results.
Run.fit.maincluster = function(p_thre=0.01) {
  num_decimal <- nchar(as.character(TEATIME$depth))
  border_vaf<-calculate_border_vaf()
  tea.result<-calculate_left_right_most_vaf(border_vaf=border_vaf,num_decimal=num_decimal)
  cat('Start Running Slope Method For Fit Case...\n');flush.console();
  #fit_slope_check<-function(p,start.div,num_decimal,end.div=NA,vaf_set=NA)
  p=tea.result$p
  start.div=tea.result$start.div
  first_check<-fit_slope_check(p=p,start.div=start.div,num_decimal,end.div=NA,vaf_set_check=NULL)
  first_check<-first_check[first_check$mu_real>3,]
  if(nrow(first_check[abs(first_check$z_score) <= 1.96, ]) > 0){
    collect.data<-first_check[abs(first_check$z_score) <= 1.96, ]
  }else{
    collect.data<-first_check
  }
  #collect.data <- ifelse(nrow(first_check[abs(first_check$z_score) <= 1.96, ]) > 0,first_check[abs(first_check$z_score) <= 1.96, ],first_check)
  celldivlist<-collect.data$cell.div
  celldivlist <- celldivlist[!is.na(celldivlist) & !is.nan(celldivlist)]

  # Sort the data frame by absolute differences
  vaf_set=TEATIME$main_cluster_vaf
  cat('Start Running Peak Method For Fit Case...\n');
  cat('   ', length(celldivlist), '...'); flush.console();
  if(debug_mode){
  cat("   ", celldivlist, "\n");flush.console();
  }
  collect.data2<- mu_estimation_from_peak(vaf_set,p=tea.result$p,num_decimal,p_thre,celldivlist=celldivlist)
  range1 <- ifelse(nrow(collect.data2) == 1, collect.data2$mu_est - 1, min(collect.data2$mu_est))
  range2 <- ifelse(nrow(collect.data2) == 1, collect.data2$mu_est + 1, max(collect.data2$mu_est))

  if(range2>median(collect.data$mu)){
    over12<-collect.data[collect.data$mu_real>=range1 & collect.data$mu_real<=range2,]
  }else{
    over12<-collect.data[collect.data$mu_real>=range1,]
  }

  if(nrow(over12)>0){
    overlap.pick<-over12
  }else{
    overlap.pick<-collect.data
  }
  #overlap.pick<-ifelse(nrow(over12)>0,over12,collect.data)

  overlap.pick<-na.omit(overlap.pick)

  ##if mu is very small
  cat('Checking if mu is small...\n');flush.console();
  all.p.data<-mu_estimation_small(overlap.pick,p,num_decimal,p_thre)
  if(nrow(all.p.data)>=3){
    # Calculate the mean and standard deviation without the maximum value
    mean_filtered <- median(all.p.data$mu)
    sd_filtered <-mad(all.p.data$mu)
    if(sd_filtered>0){
      all.p.data <- all.p.data %>%
        mutate(z_score1 = abs((mu - mean_filtered) / sd_filtered)) %>%
        filter(z_score1 < 2)
    }
  }

  all.p.data <- all.p.data %>%
    mutate(
      z_score1 = 1,
      z_score = abs(z_score),
      bicrank = rank(bic),
      score = 1 / bicrank + (1 - z_score) * 0.25,
      lowerbound1 = max(collect.data2$mu_est),
      lowerbound2 = min(collect.data2$mu_est)
    )

  all.p.data <- all.p.data[order(-all.p.data$score),]
  return(all.p.data)

}


#' Calculate Border VAF
#' @return A numeric value representing the border VAF.
calculate_border_vaf <- function() {
  main_vaf=TEATIME$main_cluster_vaf
  depth=TEATIME$depth
  a <- depth * 0.5
  b <- depth - a
  probs <- pbeta(main_vaf, a, b)
  combined_df <- data.frame(main_vaf = main_vaf, probs = probs)
  if (nrow(combined_df[combined_df$probs <= 0.05,]) > 0) {
    border_vaf <- max(combined_df[combined_df$probs <= 0.05, "main_vaf"])
  } else {
    border_vaf <- min(combined_df$main_vaf)
  }
  return(border_vaf)
}

#' Calculate Proportion from left, right most vaf
#' @param border_vaf calulated from border_vaf
#' @param num_decimal An integer representing the number of decimal places.
#' @return A numeric value representing the Proportion.
calculate_left_right_most_vaf <- function(border_vaf,num_decimal) {
  main_vaf=TEATIME$main_cluster_vaf
  depth=TEATIME$depth
  beta=TEATIME$beta
  m <- automixfit(main_vaf, type = "beta", Nc = 2:10, thresh = 0, k = 6, Ninit = min(50, round(length(main_vaf) / 5)), Niter.max = 10000)
  a=m["a",]
  b=m["b",]
  mean.a.b=a/(a+b)
  approx_vaf=min(mean.a.b)
  approx_vaf_index=which.min(mean.a.b)
  left_most_vaf <- round(approx_vaf, num_decimal)

  inita=c(0.5*depth,a[approx_vaf_index])
  initb=c(0.5*depth,b[approx_vaf_index])
  probs <- sapply(1:length(inita), function(i) dbeta(main_vaf, inita[i], initb[i]))
  df<-data.frame(prob=probs,vaf=main_vaf)
  df<-beta_reassign(df)
  main_vaf_update<-df[df$cluster<2,]$vaf

  right_most.vaf<-approx_vaf
  right_save<-right_most.vaf
  iteration=1
  while (iteration<=10000){
    m<-automixfit(main_vaf_update, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main_vaf_update)/5)), Niter.max=10000)
    updatea=m["a",]
    updateb=m["b",]
    mean.a.b=updatea/(updatea+updateb)
    approx_vaf=min(mean.a.b)
    if(approx_vaf>=0.5 | length(mean.a.b)==1){break}
    approx_vaf_index=which.min(mean.a.b)
    currenta=c(0.5*depth,updatea[approx_vaf_index])
    currentb=c(0.5*depth,updateb[approx_vaf_index])
    probs <- sapply(1:length(currenta), function(i) dbeta(main_vaf_update, currenta[i], currentb[i]))
    df<-data.frame(prob=probs,vaf=main_vaf_update)
    df<-beta_reassign(df)
    main_vaf_update<-df[df$cluster<2,]$vaf
    #print(diff)
    if(approx_vaf<0.5){
      right_most.vaf<-approx_vaf
      if(approx_vaf>right_save){right_save<-approx_vaf}
    }
    iteration=iteration+1
  }

  right.th.vaf<-left_most_vaf+(0.5-left_most_vaf)/exp(log(2)*beta)
  right_most.vaf=round(max(right_most.vaf,border_vaf,right.th.vaf,right_save),num_decimal)

  p<-(right_most.vaf*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)
  if(right_most.vaf == left_most_vaf){
    start.div<-1
  }else{
    start.div<-2
  }

  tea.result <- list(
    p = p,
    start.div = start.div
  )

  return(tea.result)
}

#' Get mu guess from slope method
#' @param p Proportion
#' @param start.div minimum cell.div
#' @param num_decimal An integer representing the number of decimal places.
#' @param end.div Set upper limit for cell.div
#' @param vaf_set_check provide vafs for main cluster
#' @return Guess from slope method.
fit_slope_check<-function(p,start.div,num_decimal,end.div=NA,vaf_set_check=NULL){

  collect.data=data.frame()
  min_mu=3
  if(is.null(vaf_set_check)){
    vaf_set=TEATIME$main_cluster_vaf
  }else{
    vaf_set=vaf_set_check
  }
  #vaf_set=ifelse(is.na(vaf_set_check),TEATIME$main_cluster_vaf,vaf_set_check)
  #print(length(vaf_set))
  depth=TEATIME$depth
  beta=TEATIME$beta
  #cell.div<-ifelse(is.na(end.div),start.div:round(length(vaf_set)/min_mu),start.div:end.div)
  if (is.na(end.div)) {
    cell.div <- start.div:round(length(vaf_set) / min_mu)
  } else {
    cell.div <- start.div:end.div
  }
  #cat('   ', length(cell.div), '...'); flush.console();
  #cell.div<-start.div:round(length(vaf_set)/min_mu)
  mu_list<-round(length(vaf_set)/cell.div)
  mu_data <- data.frame(mu_list, cell.div)

  # randomly select one corresponding cell.div
  unique_mu_list <- unique(mu_list)
  selected_cell_div <- sapply(unique_mu_list, function(x) {
    candidates <- mu_data$cell.div[mu_data$mu_list == x]
    if (length(candidates) > 1) {
      return(sample(candidates, 1))
    } else {
      return(candidates)
    }
  })

  # Create a new data frame with unique mu_list and corresponding cell.div
  unique_df <- data.frame(mu_list = unique_mu_list, cell.div = selected_cell_div)
  div.list <- unique_df$cell.div[
    unique_df$cell.div <= quantile(unique_df$cell.div, 0.75) &
      unique_df$cell.div >= quantile(unique_df$cell.div, 0.25)]

  if(is.na(end.div)){
    vaf_set=vaf_set[vaf_set>p/2]
    cat('Fit case slope running, Total cell.div count ...\n');
    cat('   ', length(div.list), '...'); flush.console();}else{
      cat('normal case slope running, Total cell.div count ...\n');
      cat('   ', length(div.list), '...'); flush.console();
    }

  for (cell.div in div.list){
    possibleError <- tryCatch({
      mu_est=length(vaf_set)/cell.div
      #print(mu_est)
      #print(length(vaf_set))
      if (mu_est <= 3) next

      i_values <- 1:cell.div
      result_vector <- sapply(i_values, function(i) {
        calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
        return(calculated_value)
      })

      a<-depth*result_vector
      b<-depth-a
      probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
      df<-data.frame(prob=probs,vaf=vaf_set)
      df<-beta_reassign(df)
      df<-df[!duplicated(df),]

      mulist<-c()
      plist<-c()

      for(try in 1:3){
        mu_from_real<-get_slope(df,p,result_vector,num_decimal)
        mu_from_simu<-slope_simu(cell.div,mu_est,p,num_decimal)
        mean_list <- mean(mu_from_simu)
        sd_list <- max(sd(mu_from_simu),1.5)
        z_score <- (mu_from_real - mean_list) / sd_list
        mulist<-c(mulist,mu_from_real)
        plist<-c(plist,abs(z_score))
      }

      closest_index <- which.min(plist)
      mu_from_real <- mulist[closest_index]
      z_score<-plist[closest_index]
    },error=function(e){
      e
    })
    if(!inherits(possibleError, "error")){
      collect.data <- rbind(collect.data, data.frame(
        mu = mu_est,
        mu_real = mu_from_real,
        mu_simu = mean_list,
        cell.div = cell.div,
        z_score = z_score
      ))
    }else{
      collect.data <- rbind(collect.data, data.frame(
        mu = NA,
        mu_real = NA,
        mu_simu = NA,
        cell.div = NA,
        z_score = NA
      ))
    }

  }
  if(nrow(na.omit(collect.data))>0){
    collect.data<-na.omit(collect.data)
  }
  collect.data$p<-p


  return(collect.data)

}

#' Prepare data for slope method
#' @param df Vaf and its cluster
#' @param p proportion
#' @param result_vector Theoretical VAF from cell.div
#' @param num_decimal An integer representing the number of decimal places.
#' @return mu Guess.
get_slope<-function(df,p,result_vector,num_decimal){
  depth=TEATIME$depth
  beta=TEATIME$beta
  G1<-df%>%
    group_by(cluster) %>%
    summarise(count = sum(freq))
  G1<-G1[order(G1$cluster),]

  G1 <- G1 %>%
    arrange(cluster) %>% # Sort by 'cluster' if it's not already sorted
    mutate(vaf = result_vector[cluster]) # Assuming 'result_vector' and 'cluster' align in terms of index.


  G1.sub <- G1 %>%
    mutate(vaf = format(vaf, nsmall = num_decimal)) %>%
    arrange(desc(vaf), desc(count)) %>%
    mutate(cumsum = cumsum(count)) %>%
    group_by(vaf) %>%
    slice_max(cumsum) %>%
    ungroup()

  G1.sub$vaf<-as.numeric(G1.sub$vaf)

  G1.sub$lnf<-2*G1.sub$vaf-p
  G1<-G1.sub[G1.sub$lnf>0,]
  G1$x<-log(G1$lnf)
  G1<-G1[order(G1$cluster),]
  if(nrow(G1)>=6){
    G1<-G1[2:floor(nrow(G1) * 0.5),]
  }
  mu_1=calculate_mu(G1)
  return(mu_1)

}

#' Calculate mu
#' @param data processed vaf data
#' @return mu
calculate_mu <- function(data) {
  depth=TEATIME$depth
  beta=TEATIME$beta
  possibleError <- tryCatch({
    mu.turn <-breakpoints(cumsum~x,data=data,h=3/nrow(data))

  },error=function(e){
    e
  })
  if(!inherits(possibleError, "error")){
    #REAL WORK
    mu.turn <-breakpoints(cumsum~x,data=data,h=3/nrow(data))
    bf <- breakfactor(mu.turn)
    # Initialize a vector to store AIC values
    AIC_values <- c()
    BIC_values<- c()
    slopes<-c()
    # Loop over each segment and calculate the AIC
    for(i in unique(bf)) {
      segment_data <- subset(data, bf == i)
      fm_segment <- lm(cumsum ~ x, data = segment_data)
      coeffs <- coef(fm_segment)
      slope <- coeffs[-1]*(-1)*beta*log(2)
      slopes<-c(slopes,slope)
      AIC_values <-c(AIC_values, AIC(fm_segment))
      BIC_values<-c(BIC_values, BIC(fm_segment))
    }

    best_segment_AIC <- which.min(AIC_values)
    best_segment_BIC <- which.min(BIC_values)
    # Calculate the average slope of the best segments based on AIC and BIC
    mu <- (slopes[best_segment_AIC] + slopes[best_segment_BIC]) / 2
  }else{
    fm<-lm(cumsum~x,data=data)
    mu=summary(fm)$coefficients[2,1]*(-1)*beta*log(2)
  }
  return(mu)  # Return mu_high
}

#' Calculate mu
#' @param cell.div processed vaf data
#' @param mu iterated mu_est from cell.div
#' @param p proportion
#' @param num_decimal An integer representing the number of decimal places.
#' @return mu list
slope_simu<-function(cell.div,mu,p,num_decimal){
  depth=TEATIME$depth
  beta=TEATIME$beta
  i_values <- 1:cell.div
  result_vector <- sapply(i_values, function(i) {
    calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  a<-depth*result_vector
  b<-depth-a
  # Pre-allocate storage
  mu_est_list <- numeric(3)

  # This is constant within the loop, so compute once
  result_matrix <- matrix(result_vector, nrow=length(result_vector), ncol=round(mu), byrow=TRUE)

  # Loop to generate simulated values
  for (try in 1:3) {
    # Use matrix operations for speed
    all_simulated_vafs_list <- sapply(result_vector, function(vaf) {
      round(rbinom(round(mu), depth, vaf) / depth, num_decimal)
    })
    all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
    probs <- sapply(1:length(a), function(i) dbeta(all_simulated_vafs, a[i], b[i]))
    df <- data.frame(prob = probs, vaf = all_simulated_vafs)
    df <- beta_reassign(df)
    df<-df[!duplicated(df),]
    mu_est_list[try] <- get_slope(df,p,result_vector,num_decimal)
  }


  return(mu_est_list)
}

#' This function estimates the small mu values
#'
#' @param data A data frame containing VAF data and associated metrics.
#' @param p proportion.
#' @param num_decimal An integer representing the number of decimal places.
#' @param p_thre A numeric value representing the threshold for p-values.
#' @return A data frame containing the estimated mu values and related metrics.
mu_estimation_small<-function(data,p,num_decimal,p_thre){
  depth=TEATIME$depth
  beta=TEATIME$beta
  vaf_set=TEATIME$main_cluster_vaf
  i_values <- 1:20
  result_vector <- sapply(i_values, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  result_vector<-c(0.5,p/2,result_vector)

  a<-depth*result_vector
  b<-depth-a
  probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
  df<-data.frame(prob=probs,vaf=vaf_set)
  df<-beta_reassign(df)
  result_df<-mu_find_small(df[df$cluster>=2,]$vaf,p)

  mu_small_selection=result_df[result_df$p_value>0.05 & result_df$wx_p>0.05 ,]

  data$mu_est<-data$mu_real
  top_pick_p<-Compare_real_simu_peak(peakdata =data,p=p,vaf_set=vaf_set,num_decimal=num_decimal)

  top_pick <- as.data.frame(t(top_pick_p))
  data$loglike<-top_pick$V2
  data$aic<-top_pick$V3
  data$bic<-top_pick$V4

  pick.cell.div=data$cell.div
  pick.mu=data$mu_real
  z_score=data$z_score
  pick.log<-data$loglike
  pick.bic<-data$bic
  pick.aic<-data$aic

  if(nrow(mu_small_selection)>0){
    num_df<-length(df[df$cluster>=2,]$vaf)
    mu_small_selection$cell.div<-round(num_df/mu_small_selection$mu_est)
    small_mu_pick<-Compare_real_simu_peak(peakdata =mu_small_selection,p=p,vaf_set=vaf_set,num_decimal=num_decimal,mu_small=T)

    small_mu_pick <- as.data.frame(t(small_mu_pick))
    small_mu_pick$V1<-ifelse(small_mu_pick$V1<0.05,0,small_mu_pick$V1)
    mu_small_selection$r1<-small_mu_pick$V1
    mu_small_selection$loglike<-small_mu_pick$V2
    mu_small_selection$aic<-small_mu_pick$V3
    mu_small_selection$bic<-small_mu_pick$V4
    mu_small_selection <- mu_small_selection[which(mu_small_selection$bic == min(mu_small_selection$bic)), ]

    if(mu_small_selection$bic<min(data$bic)*0.5){
      pick.cell.div=round(length(df[df$cluster>=2,]$vaf)/mu_small_selection$mu_est)
      pick.mu=mu_small_selection$mu_est
      pick.bic=mu_small_selection$bic
      pick.aic=mu_small_selection$aic
      pick.log=mu_small_selection$loglike
      z_score<-rep(1,length(mu_small_selection$aic))
    }

  }
  all.p.data<-data.frame(
    cell.div=pick.cell.div,
    mu=pick.mu,
    loglike=pick.log,
    bic=pick.bic,
    aic=pick.aic,
    z_score=z_score,
    p=p
  )

  return(all.p.data)

}


#' Find Small Mu Estimates
#'
#' @param vaf_set A numeric vector of VAF values.
#' @param p A numeric value representing the theoretical peak.
#' @return A data frame containing the estimated mu values and related p-values.
mu_find_small<-function(vaf_set,p){
  beta=TEATIME$beta
  depth=TEATIME$depth
  ##Now if mu is very small, then it is a line
  slope1 <- simu_slope(vaf_set, p)
  slope2 <- simu_slope(vaf_set, p + 0.01)
  slope3 <- simu_slope(vaf_set, p - 0.01)
  # Add small random noise in case the slopes are the same
  if (slope1 == slope2 || slope1 == slope3 || slope2 == slope3) {
    noise <- runif(3, min = -1e-5, max = 1e-5) # Change the range of noise as needed
    slope1 <- slope1 + noise[1]
    slope2 <- slope2 + noise[2]
    slope3 <- slope3 + noise[3]
  }

  # Add to mu_collect
  mu_collect <- c(slope1, slope2, slope3)
  # Means to test against
  means_to_test <- 1:3
  # Significance level (e.g., 0.05 for 95% confidence)
  alpha <- 0.05
  # Run the t-tests for each mean using sapply
  test_results <- sapply(means_to_test, function(x) t.test(mu_collect, mu = x)$p.value)
  bac_p=1-p
  sorted_vector <- sort(vaf_set)

  test_results2<-sapply(means_to_test, function(x) {
    i <- 1:round(length(vaf_set)/x)
    # Sort the vector
    # Compute differences between adjacent elements
    diff1<-sorted_vector[(1+x):length(sorted_vector)] - sorted_vector[1:(length(sorted_vector)-x)]
    diff2 <- (bac_p / 2*exp(log(2)*beta*i)) - (bac_p / 2*exp(log(2)*beta*(i+1)))
    wilcox.test(diff2, diff1)$p.value}
  )

  result_df <- data.frame(mu_est = means_to_test, p_value = test_results, wx_p = test_results2)

  return(result_df)

}


#' Simulate Slope for VAF Set
#'
#' @param vaf_set A numeric vector of VAF values.
#' @param p A numeric value representing the theoretical peak.
#' @return A numeric value representing the simulated slope (mu_high).
simu_slope<-function(vaf_set,p){
  depth=TEATIME$depth
  beta=TEATIME$beta
  G1.sub<-data.frame(vaf_set)
  colnames(G1.sub)<-"vaf"

  G1.sub$count<-1

  G1.sub <- G1.sub %>%
    mutate(vaf = format(vaf, nsmall = 3)) %>%
    arrange(desc(vaf), desc(count)) %>%
    mutate(cumsum = cumsum(count)) %>%
    group_by(vaf) %>%
    slice_max(cumsum) %>%
    ungroup()

  G1.sub$vaf<-as.numeric(G1.sub$vaf)

  G1.sub$lnf<-2*G1.sub$vaf-p
  G1<-G1.sub[G1.sub$lnf>0,]
  G1$x<-log(G1$lnf)


  num_rows_to_remove <- round(nrow(G1) * 0.05)
  # Calculate the index of the rows to remove
  rows_to_remove <- c(1:num_rows_to_remove, (nrow(G1) - num_rows_to_remove + 1):nrow(G1))


  # Subset the data frame to remove the rows
  if(nrow(G1[-rows_to_remove,]) > 5){
    G1<-G1[-rows_to_remove,]
  }else{
    G1<-G1
  }


  mu_high=calculate_mu(G1)
  return(mu_high)
}


extract_subp<-function(result){

  result$count<-1
  #result$vaf.1<-result$vaf.1/purity
  meanvaf<-aggregate(result$vaf.1, list(result$colors), mean)
  maxvaf<-aggregate(result$vaf.1, list(result$colors), max)
  minvaf<-aggregate(result$vaf.1, list(result$colors), min)

  sumvaf<-aggregate(result$count, list(result$colors), sum)
  mago.result<-data.frame(max=maxvaf$x,min=minvaf$x,vaf=meanvaf$x,sum=sumvaf$x)
  # Step 1: Drop the cluster with the largest 'vaf'
  max_vaf_index <- which.max(mago.result$vaf)

  mago.result <- mago.result[-max_vaf_index, ]

  # Step 2: Drop the cluster with the smallest 'min'
  min_min_index <- which.min(mago.result$min)
  mago.result <- mago.result[-min_min_index, ]

  min_min_index <- which.max(mago.result$min)
  mago.result <- mago.result[min_min_index, ]
  # Step 3: Filter clusters with 'sum' >= 5
  filtered_result <- mago.result[mago.result$sum >= 0, ]
  if (nrow(filtered_result) > 0) {
    total_sum <- sum(filtered_result$sum)
    weighted_mean_vaf <- sum(filtered_result$vaf * filtered_result$sum) / total_sum
    p=weighted_mean_vaf*2
  }else{
    p=0
  }
  return(p)
}


#' Prepare VAF and Depth Data
#' @param vafdata magos object, OR, manually, it should be dataframe with columns: vaf.1, depth.1,colors
#' colors refers to the cluster corresponding to each mutation
#' @param beta Numeric, a fixed survival rate.
#' @param depth Numeric, a fixed value.
#' @param magos_object Check if input is from magos analysis
#' @param output.folder folder name
#' @param output.prefix prefix name
#' @param id sample name
#' @param write_final Want to write final result?
#' @param debug_mode Want to debug?
#' @return A list containing the preparation results.
prepare.vaf.data = function(vafdata, beta,depth,magos_object,output.folder, output.prefix,id,write_final=write_final,debug_mode=debug_mode,purity_set=purity_set) {

  assign("debug_mode", debug_mode, envir = .GlobalEnv)  # Initialize with TRUE or FALSE

  if(debug_mode){
    cat('Debug mode ON, show error if possible ...\n');flush.console();
  }else{
    cat('Debug mode OFF ...\n');flush.console();
  }
  if(magos_object){
    purity<-min(vafdata$purity,1)
    if(purity_set>0){purity=purity_set}
    vafdata<-vafdata$result
    vafdata$vaf.1<-vafdata$vaf.1*(2-purity)/(2*vafdata$vaf.1*(1-purity)+purity)

    }
  if(is.na(depth)){depth <- round(mean(vafdata$depth.1))}

  ##Extract p value
  magosp<-extract_subp(vafdata)


  #vafdata
  vafdata.summary<- vafdata %>%
    group_by(colors) %>%
    summarise(
      max = max(vaf.1),
      min = min(vaf.1),
      vaf = mean(vaf.1),
      sum = n()
    ) %>%
    mutate(count = 1)

  vafdata.summary<-vafdata.summary[order(vafdata.summary$max,decreasing = T),]

  vafdata.summary.filter<-vafdata.summary[which(vafdata.summary$max>0.25 | vafdata.summary$min >0.25),]
  if (nrow(vafdata.summary.filter) > 1) {
    for (i in 2:nrow(vafdata.summary.filter)) {
      # check if the previous range should be merged
      if (vafdata.summary.filter$sum[i] <= 6) {
        vafdata.summary.filter$min[i-1] <- vafdata.summary.filter$min[i]
        vafdata.summary.filter$count[i] <- 0
      }
    }
  }

  vafdata.summary.filter<-vafdata.summary.filter[vafdata.summary.filter$count>0,]
  vafdata.summary.filter<-vafdata.summary.filter[order(vafdata.summary.filter$max,decreasing = T),]

  max_color <- vafdata.summary$colors[which.max(vafdata.summary$max)]
  main.cluster.vaf <- vafdata$vaf.1[which(vafdata$colors == max_color)] ## Used for F B N
  second.cluster.vaf <- list()

  ##Fit case
  if(nrow(vafdata.summary.filter)>1){
    second_highest_value <- sort(vafdata.summary.filter$max, decreasing = TRUE)[2]
    # Find the color corresponding to the second highest max value
    second_max_color <- vafdata.summary.filter$colors[vafdata.summary.filter$max %in% second_highest_value]
    second.cluster.vaf$fit  <- vafdata$vaf.1[which(vafdata$colors %in% second_max_color)]
    second.cluster.vaf$bac  <- vafdata$vaf.1[which(vafdata$colors %in% second_max_color)]
    second.cluster.vaf$normal  <- vafdata$vaf.1[which(vafdata$colors %in% second_max_color)]

  }else{
    cat('Warning sample, No second cluster found, Fit case failed ...\n');flush.console();
    second.cluster.vaf$fit  <- NULL
    cat('Warning sample, No second cluster found, Bac case recovering ...\n');flush.console();
    vafdata.summary.filter.bac<-rbind(vafdata.summary.filter,vafdata.summary[2,])
    second_highest_value <- sort(vafdata.summary.filter.bac$max, decreasing = TRUE)[2]
    second_max_color <- vafdata.summary.filter.bac$colors[vafdata.summary.filter.bac$max %in% second_highest_value]
    second.cluster.vaf$bac  <- vafdata$vaf.1[which(vafdata$colors %in% second_max_color)]

    cat('Warning sample, No second cluster found, Normal case recovering ...\n');flush.console();
    second_highest_value <-sort(vafdata.summary$max, decreasing = TRUE)[2]
    second_max_color <- vafdata.summary$colors[vafdata.summary$max %in% second_highest_value]
    second.cluster.vaf$normal  <- vafdata$vaf.1[which(vafdata$colors %in% second_max_color)]

  }


  TEATIME <- list(
    vafdata = vafdata,
    vafdata_summary = vafdata.summary,
    vafdata_summary_filter = vafdata.summary.filter,
    main_cluster_vaf = main.cluster.vaf,
    second_cluster_vaf = second.cluster.vaf,
    beta = beta,
    depth = depth,
    outputfolder=output.folder,
    outputprefix=output.prefix,
    id=id,
    magosp=magosp,
    write_final=write_final
  )
  assign("TEATIME", TEATIME, envir = .GlobalEnv)
  return(TEATIME)

}


#' Filter  mu from a range of guess
#' @param data  different cell.div and mu guess
#' @return A list containing the reliable mu.
filter.mu.estimate = function(data) {
  depth=TEATIME$depth
  beta=TEATIME$beta
  ### small mu is very hard to detect,unreliable. Unless bigger mu is not found
  if(nrow(data[data$mu>3,])>1){
    data<-data[data$mu>3,]
  }
  ##remove out
  if(nrow(data)>3){
    # Calculate the mean and standard deviation without the maximum value
    mean_filtered <- median(data$cell.div)
    sd_filtered <-mad(data$cell.div)
    data$z_score <- abs((data$cell.div - mean_filtered) / sd_filtered)
    data <- data[data$z_score < 2,]
  }
  return(data)
}



#' Estimate mu from main cluster, from slope method & Peak method in paper
#' colors refers to the cluster corresponding to each mutation
#' @param p_thre Check if input is from magos analysis
#' @return A list containing the analysis results.
Run.para.estimate.maincluster = function(p_thre=0.01) {
  output.folder=TEATIME$outputfolder
  output.prefix=TEATIME$outputprefix
  sample_name= TEATIME$id
  working.dir <- paste(output.folder, '/', sep='');


  ####Fit case####
  Fit.file.name <- paste(working.dir, output.prefix, '.fit.txt', sep='');
  Fit.file.all.name <- paste(working.dir, output.prefix, '.fit.all.txt', sep='');

  Fit.data=data.frame(
    cell.div = NA, mu = NA, loglike = NA, bic = NA,
    aic = NA, z_score = NA, p = NA, z_score1 = NA,
    bicrank = NA, score = NA, lowerbound1 = NA,
    lowerbound2 = NA, try = NA, name = sample_name, stringsAsFactors = FALSE )

  Fit.select.data=data.frame(
    mu=NA,up=NA,p1=NA,name=sample_name
  )
  Fit.mu.list=Fit.upper.list=Fit.p.list=c()
  for(try in 1:3){
    ##FIT CASE fit.over.check(main.vaf,depth,beta,p_thre)
    possibleError <- tryCatch({
      fit.check<-Run.fit.maincluster(p_thre=p_thre)
    },error=function(e){
      if(debug_mode){
      message('Caught an error in Fit.run !')
      print(e) }else{e}
      })
    #print(possibleError)
    if(!inherits(possibleError, "error")){
      fit.check <- fit.check[order(fit.check$bic),]
      fit.check$try<-try
      fit.check$name<-sample_name
      Fit.data<-rbind(Fit.data,fit.check)
      fit.check<-filter.mu.estimate(fit.check)
      Fit.mu.list=c(Fit.mu.list,fit.check$mu)
      Fit.upper.list=c(Fit.upper.list,fit.check$lowerbound1)
      Fit.p.list=c(Fit.p.list,fit.check$p)

    }

  }
  Fit.select.data$mu=ifelse(length(Fit.mu.list)>0,max(Fit.mu.list),NA)
  Fit.select.data$up=ifelse(length(Fit.mu.list)>0,max(Fit.upper.list),NA)
  Fit.select.data$p1=ifelse(length(Fit.mu.list)>0,Fit.p.list[which.max(Fit.mu.list)],NA)
  if (nrow(Fit.data) > 1) {
    Fit.data <- Fit.data[!is.na(Fit.data$mu),]
  }

  cat('Main cluster: Fit case Finish...\n');flush.console();

  write.table(Fit.select.data, Fit.file.name, sep='\t', row.names=F, quote=F);
  write.table(Fit.data, Fit.file.all.name, sep='\t', row.names=F, quote=F);



  ####Bac case####
  Bac.file.name <- paste(working.dir, output.prefix, '.bac.txt', sep='');
  Bac.file.all.name <- paste(working.dir, output.prefix, '.bac.all.txt', sep='');

  Bac.data=data.frame(
    cell.div = NA, mu = NA, loglike = NA, bic = NA,
    aic = NA, p = NA, try = NA, name = sample_name, stringsAsFactors = FALSE
  )
  Bac.select.data=data.frame(
    mu=NA,up=NA,p1=NA,name=sample_name
  )
  Bac.mu.list=Bac.p.list=c()
  for(try in 1:3){
    ##FIT CASE fit.over.check(main.vaf,depth,beta,p_thre)
    possibleError <- tryCatch({
      bac.check<-Run.bac.maincluster(p_thre=p_thre)
    },error=function(e){
      if(debug_mode){
        message('Caught an error in Bac.run !')
        print(e)
      }else{
        e
      }
      })

    if(!inherits(possibleError, "error")){
      bac.check <- bac.check[order(bac.check$bic),]
      bac.check$try<-try
      bac.check$name<-sample_name
      Bac.data<-rbind(Bac.data,bac.check)
      bac.check<-filter.mu.estimate(bac.check)
      Bac.mu.list=c(Bac.mu.list,bac.check$mu)
      Bac.p.list=c(Bac.p.list,bac.check$p)

    }

  }
  Bac.select.data$mu=ifelse(length(Bac.mu.list)>0,max(Bac.mu.list),NA)
  Bac.select.data$up=length(second.vaf<-TEATIME$second_cluster_vaf$bac)/3
  Bac.select.data$p1=ifelse(length(Bac.mu.list)>0,Bac.p.list[which.max(Bac.mu.list)],NA)
  if (nrow(Bac.data) > 1) {
    Bac.data <- Bac.data[!is.na(Bac.data$mu),]
  }

  cat('Main cluster: Bac case Finish...\n');flush.console();

  write.table(Bac.select.data, Bac.file.name, sep='\t', row.names=F, quote=F);
  write.table(Bac.data, Bac.file.all.name, sep='\t', row.names=F, quote=F);


  ####Normal case####
  Normal.file.name <- paste(working.dir, output.prefix, '.normal.txt', sep='');
  Normal.file.all.name <- paste(working.dir, output.prefix, '.normal.all.txt', sep='');
  ##all point
  Normal.data=data.frame(
    cell.div = NA, mu = NA, loglike = NA, bic = NA,
    aic = NA, p = NA, z_score = NA, reliable = NA,
    z_score1 = NA, bicrank = NA, score = NA, clear = NA,
    try = NA, name = sample_name, stringsAsFactors = FALSE
  )
  Normal.select.data=data.frame(
    mu=NA,up=NA,p1=NA,name=sample_name
  )
  Normal.mu.list=Normal.p.list=c()
  for(try in 1:3){
    ##FIT CASE fit.over.check(main.vaf,depth,beta,p_thre)
    possibleError <- tryCatch({
      normal.check<-Run.normal.maincluster(p_thre=p_thre)
    },error=function(e){
      if(debug_mode){
        message('Caught an error in Normal.run !')
        print(e)
      }else{
        e
      }
      })

    if(!inherits(possibleError, "error")){

      normal.check <- normal.check[order(-normal.check$score),]
      normal.check$try<-try
      normal.check$name<-sample_name
      Normal.data<-rbind(Normal.data,normal.check)

      Normal.mu.list=c(Normal.mu.list,normal.check$mu)
      Normal.p.list=c(Normal.p.list,normal.check$p)

    }

  }
  Normal.select.data$mu=ifelse(length(Normal.mu.list)>0,max(Normal.mu.list),NA)
  Normal.select.data$up=NA
  Normal.select.data$p1=ifelse(length(Normal.mu.list)>0,Normal.p.list[which.max(Normal.mu.list)],NA)
  if (nrow(Normal.data) > 1) {
    Normal.data <- Normal.data[!is.na(Normal.data$mu),]
  }

  cat('Main cluster: Normal case Finish...\n');flush.console();

  write.table(Normal.select.data, Normal.file.name, sep='\t', row.names=F, quote=F);
  write.table(Normal.data, Normal.file.all.name, sep='\t', row.names=F, quote=F);




}



######## S part##########

#' Get S values
#'
#' @param p Numeric, the parameter p
#' @param p_thre Numeric, threshold value for p
#' @return Data frame with estimated S values
Get.S<-function(p,p_thre=1e-6){
  depth=TEATIME$depth
  beta=TEATIME$beta
  result=TEATIME$vafdata
  ##note we have p, mu
  vaf.t1=p/2
  ##s would not be over 2, unrealistic
  min.s.detect<-log(p/(quantile(result$vaf.1,0.01)*2))/(log(2)*beta)-1
  round.cell.max<-exp(log(2)*beta*(1+min.s.detect))
  n <- max(2,floor(1+min.s.detect)) # Calculate the length of the sequence
  # Generate the sequence
  cell.list <-1:n
  vaf_set <- result$vaf.1[which(result$vaf.1 < vaf.t1 & result$vaf.1 > quantile(result$vaf.1, 0.05))]
  #vaf_set<-result[result$vaf.1<vaf.t1 & result$vaf.1 > quantile(result$vaf.1,0.05),]$vaf.1
  if(length(vaf_set) == 0){
    vaf_set<-result$vaf.1[which(result$vaf.1 < vaf.t1)]
  }else{
    vaf_set<-vaf_set
  }
  #vaf_set <- ifelse(length(vaf_set) == 0, result$vaf.1[which(result$vaf.1 < vaf.t1)],vaf_set)
  cat('Start Estimating S...\n');
  cat('   ', length(vaf_set), '...'); flush.console();
  estimate.s.data<-data.frame(s=NA,bic=NA,vaf=NA,start=NA)

  if(length(vaf_set)>0){

    cluster.result=s.dataframe.update(cell.list,p,vaf.t1,vaf_set,min.s.detect)

    if(nrow(cluster.result)==1){
      check_div<-min(cluster.result$cluster)
    }else{
      check_div<-1:nrow(cluster.result)
    }
    #check_div<-ifelse(nrow(cluster.result)==1,min(cluster.result$cluster),1:nrow(cluster.result))
    #print(cluster.result)
    svalue.list<-sapply(check_div, function(give_n) {
      svalue<-s_update_process(give_n,cluster.result,vaf_set,min.s.detect,vaf.t1,p,p_thre)
      return(svalue)
    })
    #print(svalue.list)
    if (any(svalue.list > 0)) {
      estimate.s.data<-evaluate_all_s(svalue.list,vaf.t1,p)
      #print(estimate.s.data)
      estimate.s.data<-estimate.s.data[estimate.s.data$s>0,]
    }else{
      cat('Subclonal cluster unable to estimate s, too few point Or low depth...\n');
      cat('   ', length(vaf_set), '...'); flush.console();
    }
  }else{
    cat('Subclonal cluster unable to estimate s, too few point Or low depth...\n');
    cat('   ', length(vaf_set), '...'); flush.console();
  }



  return(estimate.s.data)


}

#' Keep longest consecutive rows
#'
#' @param df Data frame to filter
#' @return Data frame with the longest consecutive rows
keep_longest_consecutive_rows <- function(df) {
  vec <- df$start
  # Find the difference between consecutive elements
  d <- c(NA, diff(vec))
  # Identify breaks in the sequence
  breaks <- which(d > 1)
  # Get start and end indices for consecutive sequences
  start <- c(1, breaks)
  end <- c(breaks - 1, length(vec))
  # Find the longest sequence
  lengths <- end - start + 1
  longest_seq <- which.max(lengths)

  # Return the rows corresponding to the longest consecutive sequence
  return(df[start[longest_seq]:end[longest_seq], ])
}

#' Pick S values
#'
#' @param df Data frame to pick S values from
#' @return Data frame with picked S values
pick_s<-function(df){
  df$bicrank<-rank(df$bic)
  df$score<-df$vaf
  if(max(df$start)>1){
    df<-keep_longest_consecutive_rows(df)
  }

  pickdata <- df[which(df$bic == min(df$bic)), ]

  if(max(df$start)>=1 ){
    data.pick.sub<-df
    data.pick.sub$truescore <- 0:(nrow(data.pick.sub)-1)
    valid_rows <- which(data.pick.sub$truescore <= data.pick.sub$start & data.pick.sub$start > 0)
    if (length(valid_rows) > 0) {
      data.pick.sub <- data.pick.sub[valid_rows, ]
      pickdata <- df[which(df$bic == min(data.pick.sub$bic)), ]
      #data.pick.sub<-data.pick.sub[data.pick$bicrank <=ceiling(nrow(data.pick) * 0.5),]
      #pickdata<-df[df$bic==min(data.pick.sub$bic),]
    }
  }



  return(pickdata)

}

#' Simulate peak ratio
#'
#' @param p Numeric, the parameter p
#' @param s Numeric, the parameter s
#' @return Numeric vector with simulated peak ratios
Simulate_ratio_peak<-function(p,s){
  depth=TEATIME$depth
  beta=TEATIME$beta
  result=TEATIME$vafdata
  ##note we have p, mu
  vaf.t1=p/2
  n <- floor(1+s)# Calculate the length of the sequence
  # Generate the sequence
  give_n=n
  if(give_n>1){
    cell.list<-generate_sequence(give_n,s,set_max=T,set_max_num=3)
  }else{
    cell.list<-(1:3)*(1+s)
  }
  #cell.list<-ifelse(give_n>1,generate_sequence(give_n,s,set_max=T,set_max_num=3),(1:3)*(1+s))
  vaf_set <- result$vaf.1[which(result$vaf.1 < vaf.t1)]
  result_vector <- sapply(cell.list, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })

  cluster.result<-s.dataframe.update(cell.list=cell.list,p=p,vaf.t1=vaf.t1,vaf_set=vaf_set,min.s.detect=NA,evaluate=F,simulation=T)

  vaf_list<-c(vaf.t1,result_vector)
  cell_list<-c(1,exp(log(2) * beta * cell.list))
  simu_ratio<-simulate_peak(vaf_list,cell_list,vaf_set)

  mu.range<-cluster.result[cluster.result$cluster==2,]$count/simu_ratio/exp(log(2) * beta * cell.list[1])
  return(mu.range)
}

#' Simulate peak values
#'
#' @param vaf_list List of VAF values
#' @param cell_list List of cell values
#' @param vaf_set List of VAF set values
#' @return Numeric vector with simulated peak ratios
simulate_peak<-function(vaf_list,cell_list,vaf_set){
  # Sample data
  depth=TEATIME$depth
  beta=TEATIME$beta
  num_decimal <- nchar(as.character(depth))
  # Assume mu and depth and num_decimal are defined in environment
  iterations <- 10
  ratio <- 0
  lowerratio <- 0
  vaf.div<-c(1:length(vaf_list))
  for(i in 1:iterations){

    all_simulated_vafs_list <- sapply(vaf.div, function(idx) {
      vaf <- vaf_list[idx]
      size <- cell_list[idx]*100
      round(rbinom(size, depth, vaf) / depth, num_decimal)
    })
    all_simulated_vafs <- unlist(all_simulated_vafs_list)
    all_simulated_vafs<-all_simulated_vafs[all_simulated_vafs>min(vaf_set)]
    ratio=ratio+Get_second_peak_ratio(all_simulated_vafs,vaf_list,cell_list,pick.ratio=1)

    sampled_vafs_list <- lapply(vaf.div, function(idx) {
      vaf_vec=all_simulated_vafs_list[idx]
      sample_size <- ceiling(length(vaf_vec)/cell_list[idx])
      sample(vaf_vec, size = sample_size)
    })

    #flattening the list
    sampled_vafs_vector <- unlist(sampled_vafs_list)
    sampled_vafs_vector<-sampled_vafs_vector[sampled_vafs_vector>min(vaf_set)]
    lowerratio=lowerratio+Get_second_peak_ratio(sampled_vafs_vector,vaf_list,cell_list,pick.ratio=1/cell_list[2])

  }
  ratio=ratio/iterations
  lowerratio=lowerratio/iterations
  return(c(ratio,lowerratio))
}

#' Get second peak ratio
#'
#' @param simu_vaf Simulated VAF values
#' @param vaf_list List of VAF values
#' @param cell_list List of cell values
#' @param pick.ratio Numeric, the pick ratio
#' @return Numeric ratio of the second peak
Get_second_peak_ratio<-function(simu_vaf,vaf_list,cell_list,pick.ratio=1){
  depth=TEATIME$depth
  beta=TEATIME$beta

  inita=depth*vaf_list
  initb=depth-inita
  probs <- sapply(1:length(inita), function(i) dbeta(simu_vaf, inita[i], initb[i]))
  df<-data.frame(prob=probs,vaf=simu_vaf)
  df<-beta_reassign(df)
  #if s close to cell 2, 4, 8,(higher less possible)
  df<-df[df$cluster>1,]
  # Sampling function
  cluster.result<- df %>%
    group_by(cluster) %>%
    summarise(
      count = n(),
      mean_vaf = mean(vaf),
      min_vaf=min(vaf)
    )

  ratio=cluster.result[cluster.result$cluster==2,]$count/(cell_list[2]*100*pick.ratio)
  return(ratio)
}


#' Update S values
#'
#' @param give_n Integer, a given cluster
#' @param cluster.result Data frame with cluster results
#' @param vaf_set List of VAF set values
#' @param min.s.detect Numeric, the minimum S detection value
#' @param vaf.t1 Numeric, the VAF at t1 time
#' @param p Numeric, the parameter p
#' @param p_thre Numeric, threshold value for p
#' @return updated S value
s_update_process<-function(give_n,cluster.result,vaf_set,min.s.detect,vaf.t1,p,p_thre=1e-6){
  depth=TEATIME$depth
  beta=TEATIME$beta

  s <- cluster.result$new_s[which(cluster.result$cluster == give_n)]
  i=1
  #print(give_n)
  while(i<=100 ){
    round.cell.max<-exp(log(2)*beta*(1+s))
    max_n=give_n+2

    if(give_n>1){
      cell.list<-generate_sequence(give_n,s)
    }else{
      cell.list<-(1:max_n)*(1+s)
    }
    #cell.list<-ifelse(give_n>1,generate_sequence(give_n,s),(1:max_n)*(1+s))

    cluster.result=s.dataframe.update(cell.list,p,vaf.t1,vaf_set,min.s.detect)
    #print(cluster.result)
    update_s <- cluster.result$new_s[which(cluster.result$cluster == give_n)]
    s=update_s
    if(abs(update_s-s)<p_thre){break}
    i=i+1
  }

  return(update_s)

}



#' Evaluate all S values
#'
#' @param svalue.list List of S values
#' @param vaf.t1 Numeric, the VAF at t1 time
#' @param p Numeric, the parameter p
#' @return Data frame with evaluated S values
evaluate_all_s<-function(svalue.list,vaf.t1,p){
  depth=TEATIME$depth
  beta=TEATIME$beta
  result=TEATIME$vafdata

  n <- length(svalue.list)+1
  end.vaf<-p/(2*exp(log(2)*beta*(1+max(svalue.list))))
  vaf_set <- result$vaf.1[which(result$vaf.1 < vaf.t1 & result$vaf.1 > end.vaf)]
  bic.list<-c()
  vaf.list<-c()
  start.list<-c()
  past_diff<-data.frame()

  s.total.num=length(svalue.list)
  for(i in 1:length(svalue.list)){
    s=0
    if(svalue.list[i]>0){
      s=svalue.list[i]
    }
    #s=ifelse(svalue.list[i]>0,svalue.list[i],0)
    give_n=i
    if(give_n>1){
      cell.list<-generate_sequence(give_n,s,set_max=T,set_max_num=s.total.num)
    }else{
      cell.list<-(1:s.total.num)*(1+s)
    }
    #cell.list<-ifelse(give_n>1,generate_sequence(give_n,s,set_max=T,set_max_num=s.total.num),(1:s.total.num)*(1+s))
    result_list<-s.dataframe.update(cell.list=cell.list,p=p,vaf.t1=vaf.t1,vaf_set=vaf_set,min.s.detect=NA,evaluate=T,simulation=F)
    cluster.result<-result_list$cluster.result
    cluster.result$try<-i
    sampled_df<-result_list$sampled_df

    # get the total log likelihood,sum these up
    total_log_likelihood <- sum(sampled_df$log_likelihood)

    bic<- compute_BIC(total_log_likelihood, 1, length(vaf_set))
    #aic <- compute_AIC(total_log_likelihood, 1)
    bic.list<-c(bic.list,bic)
    vaf.list<-c(vaf.list,sum(abs(cluster.result$realvaf-cluster.result$mean_vaf))/nrow(cluster.result))
    if(i>1){
      end_row=i-1
      down=0
      for(iter in 1:end_row){
        past_diff.pick<-past_diff[past_diff$try==iter,]
        diff_past <- past_diff.pick[1:iter,]$realvaf - past_diff.pick[1:iter,]$mean_vaf
        diff_cluster <- cluster.result[1:iter,]$realvaf - cluster.result[1:iter,]$mean_vaf

        # Compare the differences and get a binary vector (1 if true, 0 otherwise)
        binary_result <- ifelse(diff_past - diff_cluster >= 0, 1, 0)
        #binary_result <- ifelse(diff_past - diff_cluster >= 0, 1, 0)
        down=down+sum(binary_result)
      }
      start.list<-c(start.list,down)

      past_diff<-rbind(past_diff,cluster.result)
    }else{
      past_diff<-rbind(past_diff,cluster.result)
      start.list<-c(start.list,0)
    }

  }

  final.data<-data.frame(
    s=svalue.list,
    bic=bic.list,
    vaf=vaf.list,
    start=start.list
  )
  return(final.data)



}


#' Get a seq from S value
#'
#' @param n current div
#' @param s s value
#' @param set_max if need to set max div
#' @param set_max_num max div
#' @return Data frame with evaluated S values
generate_sequence <- function(n, s,set_max=F,set_max_num=0) {
  depth=TEATIME$depth
  beta=TEATIME$beta
  # Initialize an empty vector to store the result
  cell.list <- c()
  max_n=n+2
  if(set_max){
    max_n=set_max_num
  }
  # Initialize the starting value
  current_val <- 0
  # Initialize a counter to keep track of positions
  # Generate sequence
  for (i in 1:max_n) {  # Generate n * 10 terms as an example
    # Append the current value to the list
    #print(current_val)
    # Check if the counter is divisible by n
    current_val <- if (i %% n == 0) (1 + s) * (i / n) else current_val + 1

    cell.list[i] <- current_val
  }

  return(cell.list)
}


#' Update dataframe for S values
#'
#' @param cell.list List of cell values
#' @param p Numeric, the parameter p
#' @param vaf.t1 Numeric, the VAF at t1 time
#' @param vaf_set List of VAF set values
#' @param min.s.detect Numeric, the minimum S detection value
#' @param evaluate Logical, whether to evaluate
#' @param simulation Logical, whether it is a simulation
#' @return Data frame with updated values
s.dataframe.update<-function(cell.list,p,vaf.t1,vaf_set,min.s.detect,evaluate=F,simulation=F){
  depth=TEATIME$depth
  beta=TEATIME$beta

  result_vector <- sapply(cell.list, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  mix_check<-sapply(result_vector, function(i) {
    # Replace this with your actual calculation
    calculated_value <- ifelse(log((1-p)/2/i)/(log(2) * beta)>0,1,0)
    return(calculated_value)
  })
  cluster_idx <- 1:length(mix_check)

  # Create a named vector for easier matching
  named_mix_vec <- setNames(mix_check, cluster_idx)

  inita=depth*c(vaf.t1,result_vector)
  initb=depth-inita
  probs <- sapply(1:length(inita), function(i) dbeta(vaf_set, inita[i], initb[i]))
  df<-data.frame(prob=probs,vaf=vaf_set)
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  df<-beta_reassign(df)

  if(evaluate){
    p_vec=c(vaf.t1,result_vector)
    df <- df %>%
      mutate(p = p_vec[cluster])
  }
  df<-df[df$cluster>1,]

  # Sampling function
  sample_or_not <- function(data, cluster_id,fraction, mix_vec) {
    if (mix_vec[cluster_id-1] == 1) {
      return(sample_frac(data, size = fraction))
    } else {
      return(data)
    }
  }
  # Apply sampling only to clusters with mix_value of 1
  sampled_df <- df %>%
    group_by(cluster) %>%
    group_modify(~ sample_or_not(.x, .y$cluster,p, named_mix_vec)) %>%
    ungroup()



  if(simulation){
    cluster.result<- sampled_df %>%
      group_by(cluster) %>%
      summarise(
        count = n(),
        mean_vaf = mean(vaf),
        min_vaf=min(vaf)
      )
    return(cluster.result)
  }

  if(evaluate){
    sampled_df <- sampled_df %>%
      mutate(log_likelihood = likbeta(x = vaf,
                                      shape1 = p*depth, shape2 = depth-p))
    cluster.result<- sampled_df %>%
      group_by(cluster) %>%
      summarise(
        count = n(),
        mean_vaf = mean(vaf)
      )
    cluster.result$realvaf<-result_vector[cluster.result$cluster-1]


    return(list(sampled_df = sampled_df, cluster.result = cluster.result))
  }

  cluster.result <- sampled_df %>%
    group_by(cluster) %>%
    summarise(
      count = n(),
      mean_vaf = mean(vaf),
      min_vaf = min(vaf)
    ) %>%
    filter(cluster > 1) %>%
    mutate(
      new_s = log(p / (mean_vaf * 2)) / (log(2) * beta) - 1,
      alt_s = log(p / (min_vaf * 2)) / (log(2) * beta) - 1,
      new_s = ifelse(new_s < 0, alt_s, new_s)
    ) %>%
    filter(new_s <= min.s.detect) %>%
    mutate(cluster = cluster - 1)

  return(cluster.result)



}

#####S part end###

#' All.guess.update: This function processes the given TEATIME object to update guesses
#' @return updated Data
All.guess.update<-function(){
  depth=TEATIME$depth
  beta=TEATIME$beta
  output.folder=TEATIME$outputfolder
  output.prefix=TEATIME$outputprefix
  sample_name= TEATIME$id


  working.dir <- paste(output.folder, '/', sep='');
  Fit.file.name <- paste(working.dir, output.prefix, '.fit.txt', sep='');
  Fit.file.all.name <- paste(working.dir, output.prefix, '.fit.all.txt', sep='');
  Bac.file.name <- paste(working.dir, output.prefix, '.bac.txt', sep='');
  Bac.file.all.name <- paste(working.dir, output.prefix, '.bac.all.txt', sep='');
  Normal.file.name <- paste(working.dir, output.prefix, '.normal.txt', sep='');
  Normal.file.all.name <- paste(working.dir, output.prefix, '.normal.all.txt', sep='');

  output.file.name <- paste(working.dir, output.prefix, '.all.guess.update.txt', sep='');


  fit.all<- read.table(Fit.file.all.name, sep='\t', header=T, stringsAsFactors=F); ##bac.all allpoint
  inter.all<-read.table(Normal.file.all.name, sep='\t', header=T, stringsAsFactors=F);
  #bac.all<-read.table(Bac.file.all.name, sep='\t', header=T, stringsAsFactors=F);

  fit.select.one<-read.table(Fit.file.name, sep='\t', header=T, stringsAsFactors=F);
  inter.select.one<-read.table(Normal.file.name, sep='\t', header=T, stringsAsFactors=F);
  bac.select.one<-read.table(Bac.file.name, sep='\t', header=T, stringsAsFactors=F);##single predict.all


  fitmu=NA
  fitmu_candidate=1 ##In case all case are NA
  fitcell=NA
  fits=NA
  fitp=NA
  cat('Fit case guess...\n');
  cat('   ',fit.select.one$mu,'...'); flush.console();
  if (nrow(fit.all[!is.na(fit.all$cell.div), ]) > 0) {
    fit.all <- fit.all[!is.na(fit.all$cell.div), ]
    fit.data<-Find.s.from.predict(fit.all)

    if(!("mu" %in% colnames(fit.data))) {
      #fit.select.one
      mu.fit.pick=fit.select.one$up
      #fit.all
      fit.select.row<-fit.all[fit.all$lowerbound1==mu.fit.pick,]
      fit.select.row<-fit.select.row[sample(nrow(fit.select.row), 1),]
      p<-fit.select.row$p
      mu.s.update<-fit.data[fit.data$p_value==p,]
      mu.select2<-mu.s.update[sample(nrow(mu.s.update), 1),]$minmu
      df<-Get.S(p, p_thre=1e-6)
      df$p_value <- p
      data<-pick_s(df)

      fitmu=mu.fit.pick
      fitmu_candidate=mu.select2
      fitcell=fit.select.row$cell.div
      fits=data$s
      fitp=p

    }else{
      fitmu=fit.data$mu
      fitmu_candidate=1
      fitcell=fit.data$cell.div
      fits=fit.data$s
      fitp=fit.data$p
    }
  }

  cat('Normal case guess...\n');
  cat('   ',inter.select.one$mu,'...'); flush.console();
  intermu=NA
  intermu_candidate=NA
  intertrust=NA
  intercell=NA
  inters=NA
  interp=NA
  if (nrow(inter.all[!is.na(inter.all$cell.div), ]) > 0) {
    inter.all <- inter.all[!is.na(inter.all$cell.div), ]
    inter.data<-Find.s.from.predict(inter.all)

    if(!("mu" %in% colnames(inter.data))) {
      mu.inter.pick=inter.select.one$mu

      inter.select.row<-inter.all[inter.all$mu==mu.inter.pick,]
      inter.select.row<-inter.select.row[sample(nrow(inter.select.row), 1),]
      p<-inter.select.row$p
      mu.s.update<-inter.data[inter.data$p_value==p,]
      mu.select2<-mu.s.update[sample(nrow(mu.s.update), 1),]$minmu

      df<-Get.S(p, p_thre=1e-6)
      df$p_value <- p
      data<-pick_s(df)

      intermu=mu.select2
      intermu_candidate=mu.inter.pick
      intertrust=(2*abs(mu.select2-mu.inter.pick))/(mu.inter.pick+mu.select2)
      intercell=inter.select.row$cell.div
      inters=data$s
      interp=p


    }else{
      #inter.data$cell.div
      intermu=inter.data$mu
      intermu_candidate=1
      intertrust=0
      intercell=inter.data$cell.div
      inters=inter.data$s
      interp=inter.data$p
    }
  }

  backp=bac.select.one$p1

  data.rearrange <- data.frame(
    name=sample_name,
    fitmu = fitmu,
    fitmu_candidate = fitmu_candidate,
    fitcell = fitcell,
    fits = fits,
    fitp = fitp,
    intermu = intermu,
    intermu_candidate = intermu_candidate,
    intertrust = intertrust,
    intercell = intercell,
    inters = inters,
    interp = interp,
    backp = backp,
    stringsAsFactors = FALSE
  )

  #print(data.rearrange)

  write.table(data.rearrange, output.file.name, sep='\t', row.names=F, quote=F);

}

#' Main cluster vaf update if we can not find peak directly
#' @return updated Peak
maincluster.candidate.update<-function(){
  main.vaf<-TEATIME$main_cluster_vaf
  depth=TEATIME$depth
  beta=TEATIME$beta


  temp.vaf<-(0.5+min(main.vaf))/2

  iter=0
  repeat {
    old.vaf <- if (iter == 0) temp.vaf else new.vaf
    inita=depth*c(0.5,old.vaf)
    initb=depth-inita
    probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
    df<-data.frame(prob=probs,vaf=main.vaf)
    df<-beta_reassign(df)
    df<-df[df$cluster>1,]
    new.vaf<-mean(df$vaf)
    diff<-abs(new.vaf-old.vaf)

    if (diff <= 1 / depth || iter > 100) {break}
    iter=iter+1

  }



  return(new.vaf)
}


#' Find S values from prediction result
#'
#' @param predict.result Data frame with prediction results
#' @return Data frame with selected S values
Find.s.from.predict<-function(predict.result){
  depth=TEATIME$depth
  beta=TEATIME$beta
  result=TEATIME$vafdata

  p.list=unique(predict.result$p)
  s.results <- lapply(p.list, function(p) {
    df<-Get.S( p, p_thre=1e-6)
    #print(df)
    #df<-df[!is.na(df$s),]
    df <- na.omit(df)
    if(nrow(df)>0){
    df$p_value <- p
    data<-pick_s(df)
    #print(data)
    murange=Simulate_ratio_peak(p,data$s)
    data$minmu=min(murange)
    data$maxmu=max(murange)
    return(data)}else{
     return(NULL)  
    }
  })
  #cat('*****************S**********\n');flush.console();
  #print(s.results)
  fit.data<-data.frame()
  valid_indices <- which(!sapply(s.results, is.null))
  s.results <- s.results[valid_indices]
  p.list <- p.list[valid_indices]

  combined_df <- do.call(rbind, s.results)
  for(try in 1:length(p.list)){
    p.pick=p.list[try]
    fit.check <- predict.result[which(predict.result$p == p.pick), ]
    if(nrow(fit.check)>0){
      minmu <- combined_df[which(combined_df$p_value == p.pick), ]$minmu
      maxmu <- combined_df[which(combined_df$p_value == p.pick), ]$maxmu
      fit.check.sub <- fit.check[which(fit.check$mu >= minmu & fit.check$mu <= maxmu), ]
      fit.check.sub=na.omit(fit.check.sub)

      if(nrow(fit.check.sub)==0){
        minmu=minmu-1.96*1.5
        fit.check.sub<-fit.check[fit.check$mu >=minmu & fit.check$mu<=maxmu, ]
        fit.check.sub=na.omit(fit.check.sub)
      }
      if(nrow(fit.check.sub)>0){
        best_fit <- fit.check.sub[which(fit.check.sub$bic == min(fit.check.sub$bic)), ]
        fit.data <- rbind(fit.data, best_fit)
      }
    }
  }

  # Find rows with minimum bic
  if(nrow(fit.data)>0){
    # Find rows with minimum bic
    min_bic <- min(fit.data$bic)
    min_bic_indices <- which(fit.data$bic == min_bic)
    # Randomly select one index if multiple rows have the same minimum BIC
    selected_index <- sample(min_bic_indices, 1)
    # Assign the corresponding 's' value from combined_df
    selected_row <- fit.data[selected_index, ]
    selected_row$s <- combined_df[which(combined_df$p_value == selected_row$p), ]$s


    return(selected_row)
  }else{
    selected_row<-combined_df
    return(combined_df)
  }
}


###Final process####
#' Determine the final mu pick based on given conditions
#'
#' @param fitdiff Numeric, difference value for fit.
#' @param interdiff Numeric, difference value for inter.
#' @param fitdiff2 Numeric, secondary difference value for fit.
#' @param interdiff2 Numeric, secondary difference value for inter.
#' @param fitmu Numeric, mu value for fit.
#' @param intermu Numeric, mu value for inter.
#' @return Numeric, the chosen mu value based on the given conditions.
determine_mupick <- function(fitdiff, interdiff, fitdiff2, interdiff2, fitmu, intermu) {
  case_when(
    is.na(fitdiff) & !is.na(interdiff) ~ intermu,
    !is.na(fitdiff) & is.na(interdiff) ~ fitmu,
    is.na(fitdiff) & is.na(interdiff) ~ NA_real_,
    fitdiff2 < interdiff2 ~ fitmu,
    TRUE ~ intermu
  )
}

#' Adjust mu values based on given conditions
#'
#' @param mu Numeric, the original mu value.
#' @param mu_candidate Numeric, the candidate mu value for adjustment.
#' @param times Numeric, the threshold value for adjustment.
#' @return Numeric, the adjusted mu value.
adjust_mu <- function(mu, mu_candidate, times) {
  adjusted_mu <- ifelse(mu_candidate != 1 &
                          (mu / mu_candidate > times | mu / mu_candidate < 1 / times),
                        0.5 * mu_candidate + 0.5 * mu,
                        mu)
  return(adjusted_mu)
}


adjust_p<-function(data,magosp,beta,cut=0.3){
  data$goodp<-magosp
  data<-data[!(is.na(data$interp)& is.na(data$fitp)),]

  data$pickp.alt <- ifelse(
          is.na(data$intermu),  # If intermu is NA, directly use fitp
          data$fitp,
          ifelse(
            data$mupick == data$intermu,
            coalesce(data$interp, data$fitp),  # Prioritize interp if available
            coalesce(data$fitp, data$interp)   # Otherwise, use fitp first
          )
        )

        data$pickp.close <- apply(data, 1, function(row) {
          candidates <- as.numeric(c(row["fitp"], row["interp"],row["backp"]))
          candidates <- candidates[!is.na(candidates)]
          if (length(candidates) == 0) return(NA)
          freq1 <- as.numeric(row["goodp"]) # Access freq1 directly from the current row
          closest <- candidates[which.min(abs(candidates - freq1))]
          closest
        })

        data<-data[!is.na(data$pickp.close),]
        data$pickp.close<-as.numeric(data$pickp.close)

        data$mupick.new <- ifelse(
          is.na(data$fitp),  # If fitp is NA, directly assign intermu
          data$intermu,
          ifelse(
            abs(data$pickp.close - data$fitp) < 0.01,
            coalesce(data$fitmu, data$intermu),  # Prioritize fitmu if available
            coalesce(data$intermu, data$fitmu)  # Otherwise, use intermu first
          )
        )
        data$mupick.choose<-ifelse(round(abs(data$pickp.alt-data$goodp),2)>=cut,data$mupick.new,data$mupick)
        data$mupick<-data$mupick.choose
        data$pickp<-ifelse(data$mupick.choose==data$mupick.new,data$pickp.close,data$pickp.alt)

        data$pickp<-ifelse(round(abs(data$goodp-data$pickp),2)>=cut,data$goodp/2+data$pickp/2,data$pickp)
        data$picks <- ifelse(
          is.na(data$intermu),
          data$fits,
          ifelse(
            data$mupick == data$intermu,
            data$inters,data$fits
          )
        )
        data$pickt1 <- ifelse(
          is.na(data$intermu),  # If intermu is NA, directly use fitp
          data$fit_len_diff/data$mupick,
          ifelse(
            data$mupick == data$intermu,
            data$inter_len_diff/data$mupick,data$fit_len_diff/data$mupick
          )
        )
        #data$pickp<-ifelse(data$mupick==data$intermu,data$interp,data$fitp)
        data$picktend=((log(data$pickp)-log(1-data$pickp))/(log(2)*beta)+(1+data$picks)*data$pickt1)/data$picks
        data$picktend=ifelse(data$picktend>0,data$picktend,NA)
        valid_idx <- which(data$picktend >= data$pickt1)


        invalid_idx <- setdiff(seq_len(nrow(data)), valid_idx)

        data[invalid_idx, c("picks", "mupick", "pickp", "pickt1", "picktend")] <- NA
        return(data)
}



#' Process and finalize the data
#'
#' @param data.rearrange Data frame
#' @param Rbest.data Data frame, Rbest file
#' This function modifies and writes the final data to a file.
#' @return None
final.process<-function(data.rearrange,Rbest.data){
  depth=TEATIME$depth
  beta=TEATIME$beta
  output.folder=TEATIME$outputfolder
  output.prefix=TEATIME$outputprefix
  sample_name= TEATIME$id
  magosp=TEATIME$magosp
  write_final<-TEATIME$write_final
  missample <- Rbest.data$samplename[which(Rbest.data$len_adj == 1)]

  working.dir <- paste(output.folder, '/', sep='');
  final.file.name <- paste(working.dir, output.prefix, '.final.txt', sep='');

  final.data=data.frame(
    name=sample_name,
    mu=NA,
    s=NA,
    t1=NA,
    tend=NA,
    p=NA
  )
  if(data.rearrange$label=="inter"){
    data=data.rearrange
    times=2
    data <- data %>%
      mutate(
        intermu = adjust_mu(intermu, intermu_candidate, times),
        fitmu = adjust_mu(fitmu, fitmu_candidate, times)
      )
    data$mupick <- ifelse(
    is.na(data$fitdiff) & !is.na(data$interdiff), data$intermu, # fit_len_ratio is NA
    ifelse(
      !is.na(data$fitdiff) & is.na(data$interdiff), data$fitmu, # inter_len_ratio is NA
      ifelse(
        is.na(data$fitdiff) & is.na(data$interdiff), NA, # Both are NA
        ifelse(data$fitdiff2 < data$interdiff2, data$fitmu, data$intermu) # Neither are NA
      )
    )
  )

    if(Rbest.data$len_adj==1){
      data$mupick<-ifelse(data$name %in% missample,data$intermu,data$mupick)
      }

  
    data <- data %>%
      mutate(mupick_low_depth = determine_mupick(fitdiff, interdiff, fitdiff2, interdiff2, fitmu, intermu))

    data<-adjust_p(data,magosp,beta)

    final.data=data.frame(
      name=sample_name,
      mu=data$mupick,
      s=data$picks,
      t1=data$pickt1,
      tend=data$picktend,
      p=data$pickp
    )

  }

  if(write_final){
  write.table(final.data, final.file.name, sep='\t', row.names=F, quote=F);
  }
  return(final.data)



}

#' Compare three case with Rbest
#'
#' @return None
Post_process<-function(){
  depth=TEATIME$depth
  beta=TEATIME$beta
  output.folder=TEATIME$outputfolder
  output.prefix=TEATIME$outputprefix
  sample_name= TEATIME$id
  second.vaf<-TEATIME$second_cluster_vaf$normal
  main.vaf<-TEATIME$main_cluster_vaf
  vaf.all<-c(main.vaf,second.vaf)

  working.dir <- paste(output.folder, '/', sep='');
  update.guess.file.name <- paste(working.dir, output.prefix, '.all.guess.update.txt', sep='');
  Rbest.file.name <- paste(working.dir, output.prefix, '.rbest.info.txt', sep='');


  data.rearrange<-read.table(update.guess.file.name, sep='\t', header=T, stringsAsFactors=F);
  Rbest.data<-read.table(Rbest.file.name, sep='\t', header=T, stringsAsFactors=F);

  m<-automixfit(vaf.all, type = "beta",Nc =2:10,thresh=0,k = 6,Niter.max=10000)
  a=m["a",]
  b=m["b",]
  mean.a.b.new=a/(a+b)
  close_05_vaf=mean.a.b.new[which.min(abs(mean.a.b.new - 0.5))]
  mean_a_b_filtered <- mean.a.b.new[which(mean.a.b.new > min(main.vaf) & mean.a.b.new < close_05_vaf)]

  if(length(mean_a_b_filtered) == 0){
    insert.vaf<-maincluster.candidate.update()
    mean.a.b<-c(mean.a.b.new[which(mean.a.b.new<close_05_vaf)],insert.vaf)
    data.rearrange$clonallen=1
  }else{
    mean.a.b<-mean.a.b.new[which(mean.a.b.new<close_05_vaf)]
    data.rearrange$clonallen=2
  }

  data.rearrange$minvaf=min(mean.a.b)
  data.rearrange$maxvaf=max(mean.a.b)
  data.rearrange$closevaf=mean.a.b[which.min(abs(mean.a.b - 0.5))]

  # For fit
  fitp=data.rearrange$fitp
  fits=data.rearrange$fits
  if(!is.na(fitp)) {
    fit1 = fitp/2 + (1-fitp)/(2 * exp(log(2) * beta * 1))
    data.rearrange$closefitvaf=mean.a.b[which.min(abs(mean.a.b - fit1))]
    data.rearrange$fitdiff=abs(mean.a.b[which.min(abs(mean.a.b - fit1))] - fit1)

    fit2= fitp/2
    data.rearrange$closefitvaf2=mean.a.b[which.min(abs(mean.a.b - fit2))]
    data.rearrange$fitdiff2=abs(mean.a.b[which.min(abs(mean.a.b - fit2))] - fit2)



    inita=depth*c(0.5,fit1,fit2)
    initb=depth-inita
    probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
    df<-data.frame(prob=probs,vaf=main.vaf)
    df<-beta_reassign(df)

    data.rearrange$fit_len_diff=data.rearrange$fitmu*data.rearrange$fitcell
    data.rearrange$fit_len_ratio=(data.rearrange$fitmu*data.rearrange$fitcell)/nrow(df[df$cluster>1,])

  } else {
    data.rearrange$closefitvaf=NA
    data.rearrange$fitdiff=NA
    data.rearrange$closefitvaf2=NA
    data.rearrange$fitdiff2=NA
    data.rearrange$fit_len_diff=NA
    data.rearrange$fit_len_ratio=NA
  }

  # For inter
  interp=data.rearrange$interp
  inters=data.rearrange$inters
  if(!is.na(interp)) {
    inter1 = interp/2 + (1-interp)/(2 * exp(log(2) * beta * 1))
    data.rearrange$closeintervaf=mean.a.b[which.min(abs(mean.a.b - inter1))]
    data.rearrange$interdiff=abs(mean.a.b[which.min(abs(mean.a.b - inter1))] - inter1)

    inter2=interp/2
    data.rearrange$closeintervaf2=mean.a.b[which.min(abs(mean.a.b - inter2))]
    data.rearrange$interdiff2=abs(mean.a.b[which.min(abs(mean.a.b - inter2))] - inter2)



    inita=depth*c(0.5,inter1,inter2)
    initb=depth-inita
    probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
    df<-data.frame(prob=probs,vaf=main.vaf)
    df<-beta_reassign(df)
    intermu<-data.rearrange$intermu
    #intermu=ifelse(data.rearrange[i_index,"intermu_candidate"]!=1 & data.rearrange[i_index,"intermu"]/data.rearrange[i_index,"intermu_candidate"]>10,0.5*data.rearrange[i_index,"intermu_candidate"]+0.5*data.rearrange[i_index,"intermu"],data.rearrange[i_index,"intermu"])
    #intermu=ifelse(data.rearrange[i_index,"intermu_candidate"]!=1 & data.rearrange[i_index,"intermu"]/data.rearrange[i_index,"intermu_candidate"]>10,data.rearrange[i_index,"intermu_candidate"],data.rearrange[i_index,"intermu"])
    data.rearrange$inter_len_diff=intermu*data.rearrange$intercell

    data.rearrange$inter_len_ratio=(intermu*data.rearrange$intercell)/nrow(df[df$cluster>1,])


  } else {
    #inter1 = NA # or any other fallback value or action
    data.rearrange$closeintervaf=NA
    data.rearrange$interdiff=NA
    data.rearrange$closeintervaf2=NA
    data.rearrange$interdiff2=NA
    data.rearrange$inter_len_diff=NA
    data.rearrange$inter_len_ratio=NA

  }

  # For bac
  backp=data.rearrange$backp
  if(!is.na(backp)) {
    bac1 = backp/2 + (1-backp)/(2 * exp(log(2) * beta * 1))
    data.rearrange$closebacvaf=mean.a.b[which.min(abs(mean.a.b - bac1))]
    data.rearrange$bacdiff=abs(mean.a.b[which.min(abs(mean.a.b - bac1))] - bac1)

    bac2 = (1-backp)/(2 * exp(log(2) * beta * 1))
    data.rearrange$closebacvaf2=mean.a.b[which.min(abs(mean.a.b - bac2))]
    data.rearrange$bacdiff2=abs(mean.a.b[which.min(abs(mean.a.b - bac2))] - bac2)

  } else {
    data.rearrange$closebacvaf=NA # or any other fallback value or action
    data.rearrange$bacdiff=NA
    data.rearrange$closebacvaf2=NA
    data.rearrange$bacdiff2=NA
  }



  #Rbest.data only use len_adj
  inter_sample <- Rbest.data$samplename[which(Rbest.data$len_adj > 1)]
  missample <- Rbest.data$samplename[which(Rbest.data$len_adj == 1)]

  cat('Rbest case infer and miss:...\n');
  cat('   ',inter_sample,missample,'...'); flush.console();

  fit_sum <- sum(data.rearrange[c("fitdiff2", "fitdiff2")], na.rm = TRUE)
  if (any(is.na(data.rearrange[c("fitdiff2", "fitdiff2")]))) {
    fit_sum <- Inf
  }

  inter_sum <- sum(data.rearrange[c("interdiff2", "interdiff2")], na.rm = TRUE)
  #inter_sum <- min(data.rearrange[i_index, "interdiff"],data.rearrange[i_index, "interdiff2"])
  if (any(is.na(data.rearrange[c("interdiff2", "interdiff2")]))) {
    inter_sum <- Inf
  }

  bac_sum <- sum(data.rearrange[c("bacdiff2", "bacdiff2")], na.rm = TRUE)
  if (any(is.na(data.rearrange[c("bacdiff2", "bacdiff2")]))) {
    bac_sum <- Inf
  }




  cat('Our case...\n');
  cat('   ',data.rearrange$name,'...'); flush.console();

  data.rearrange$label <- NA
  if(data.rearrange$name %in% inter_sample){
    data.rearrange$label<- "inter"
  }

  if(data.rearrange$name %in% missample){
    if(inter_sum >bac_sum) {
      data.rearrange$label <- "bac"
    } else{
      data.rearrange$label<- "inter"
    }
  }

  final.data=final.process(data.rearrange,Rbest.data)
  return(final.data)

}



#' Perform Evolutionary parameter estimation from sequencing data.
#' @param input.file Three options: 1. VCF like input for magos if steps has 0; 2. magos result;
#' 3.manually, it should be dataframe with columns: vaf.1, depth.1,colors
#' colors refers to the cluster corresponding to each mutation
#' @param beta Numeric, survival rate. Default is 0.9
#' @param depth Numeric,Only enter if you know exactly what the value is.
#' @param p_thre p value cutoff used in most function. Default is 0.01
#' @param magos_object boolean value Check if input is from magos result
#' @param output.folder the folder to write the temporary and final prediction files.
#' @param output.prefix prefix used to name the temporary and final prediction files.
#' @param id Assign unique id for each analysis
#' @param steps A vector of integers indicating which functions to execute. 0: run MAGOS with input; 1: prepare VAF files; 2. Rbest estimate. 3. mutation rate estimate; 4. Fitness estimate. 5 Final prediction.
#' @param write_final if want to write final result in .txt file
#' @param debug_mode If want to show error
#' @param seed An optional seed for reproducibility.
#' @return NULL
#' @export
TEATIME.run <- function(input.file,beta=0.9,depth=NA,p_thre=0.01,magos_object=T,output.folder="./", output.prefix="TEATIME",id='T01',steps=1:5,write_final=T,debug_mode=F,purity_set=0,seed=NA) {
  output.folder <- paste(output.folder, '/', sep='');
  if(!is.na(seed)){
    set.seed(seed)
  }
  if(0 %in% steps) {
    ## filter CNV
    # Ensure the last column is numeric
    input.file[, ncol(input.file)] <- as.numeric(input.file[, ncol(input.file)])

    # Filter rows where the last column (CN) is 2 (normal)
    input.file <- input.file[input.file[, ncol(input.file)] == 2, ]

    # Remove the last column (CN) after filtering
    input.file <- input.file[, -ncol(input.file)]

    input.file = mag.single.run(input.file,fold= T)
    magos_object=T

    #prepare.vaf.data(vafdata=input.file, beta=beta,depth=depth, magos_object=magos_object,output.folder=output.folder, output.prefix=output.prefix,id=id);
  }
  if(!is.na(seed)){
    seed=seed+1
    set.seed(seed)
  }
  if(1 %in% steps) {

    #prepare.vaf.data = function(vafdata, beta,depth,magos_object,output.folder, output.prefix,id,write_final=write_final,debug_mode=debug_mode)
    prepare.vaf.data(vafdata=input.file, beta=beta,depth=depth, magos_object=magos_object,output.folder=output.folder, output.prefix=output.prefix,id=id,write_final=write_final,debug_mode=debug_mode,purity_set=purity_set);

    }
  if(2 %in% steps) {
    ##Rbest
    Rbest_classify();
  }
  if(3 %in% steps) {
    Run.para.estimate.maincluster(p_thre=p_thre)
  }
  if(4 %in% steps) {
    All.guess.update();
  }
  if(5 %in% steps) {
    final.data=Post_process();
    return(final.data)
  }
}
