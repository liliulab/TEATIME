# Production estimation pipeline. Functions are isolated in `.prod_defs` (a local
# environment) and copied into a per-run dispatch environment by `.load_production`.
.prod_defs <- local({
  ## ---- optimizeS.R ----
#source("Script/optimizeM.R")
#source("Script/loss.R")
round2 = function(x, digits) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}
find.best<-function(pick.s.range,M.data,sub.vaf,vaf.t1.main,p.f,mu.simu,beta,depth){
  r.data<-data.frame(vaf=numeric(),
                     r=numeric(),
                     mu=numeric(),
                     s=numeric()
                     
  )
  
  current.sub.vaf<-c(M.data[M.data<=vaf.t1.main],sub.vaf[sub.vaf<=vaf.t1.main])
  
  G1.sub<-data.frame(current.sub.vaf)
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
  if(pick.s.range$max==pick.s.range$min){
    vaf.list=c(pick.s.range$max)
  }else{
    #test.range<-G1.sub[G1.sub$vaf<=pick.s.range$max & G1.sub$vaf>=pick.s.range$min,]
    vaf.list<-seq(pick.s.range$min,pick.s.range$max,0.01)
    n_steps <- ceiling((pick.s.range$max - pick.s.range$min) / 0.01)
    if (vaf.list[n_steps] != pick.s.range$max) {
      vaf.list <- c(vaf.list, pick.s.range$max)
    }
  }
  datas=data.frame(
    vaf=numeric(),
    p=numeric()
  )
  for(i in seq(1,length(vaf.list),1)){
    
    cut=vaf.list[i]
    
    s=log(p.f/(cut*2))/(log(2)*beta)-1
    if(s>0){
    div=round2(log((1-p.f)/(2*cut))/(log(2)*beta),0)
    if(div>0){
      for (j in 1:div) { 
        cell.num<-exp(log(2)*beta*j)
        
        vec<-(1-p.f)/(2*cell.num)
        
        if(cut>vec){
          div=j-1
          break
        }
        
      }
    
    }
    div<-ifelse(div>0,div,0)
    if(div==0){
      # Create an empty vector to store the results
      
      cell.num<-exp(log(2)*beta*(1+s))
      if(cell.num <=exp(log(2)*beta*2)){
        vec<-c(vaf.t1.main,vaf.t1.main/(exp(log(2)*beta*1)),cut)
        vec<-data.frame(vec)
        vec$a<-depth*vec$vec
        vec$b<-depth-vec$a
        
        sub.include<-current.sub.vaf[current.sub.vaf>=min(vec$vec)]
        
        probs <- sapply(1:nrow(vec), function(i) dbeta(sub.include, vec$a[i], vec$b[i]))
        max_rows <- apply(probs, 1, which.max)
    
        idx <- which(max_rows == 2)
        count1 <- length(idx)
        mu.2=count1/exp(log(2)*beta*1)
        
        
        #p<-ifelse(count_max_rows<expect.mu,pnorm(count_max_rows,mean =expect.mu, sd=sqrt(expect.mu)),1-pnorm(count_max_rows,mean =expect.mu, sd=sqrt(expect.mu)))
        datas[i,"vaf"]=cut
        datas[i,"p"]=lossfun(mu.simu,mu.2,0,0,1)
        #datas[i,"p"]=(((mu.simu+1)-(mu.2+1)))^2
      }else{
        p=p.f
        vec <- c(vaf.t1.main)
        cell.count<-floor(s+1)
        for (j in 1:cell.count) { 
          
          cell.num<-exp(log(2)*beta*j)
          
          vec[j+1]<-vaf.t1.main/cell.num
          
          
          
        }
        vec[cell.count+2]<-cut
        
        vec<-data.frame(vec)
        vec$a<-depth*vec$vec
        vec$b<-depth-vec$a
        
        sub.include<-current.sub.vaf
        
        probs <- sapply(1:nrow(vec), function(i) dbeta(sub.include, vec$a[i], vec$b[i]))
        max_rows <- apply(probs, 1, which.max)
        
        clustercount <- data.frame(table(max_rows))
        clustercount$max_rows<-as.numeric(as.character(clustercount$max_rows))
        missing_value <- setdiff(1:cell.count+2, clustercount$max_rows)
        
        # Add the missing value to the dataframe with a frequency of 0
        if(length(missing_value) > 0) {
          clustercount <- rbind(clustercount, data.frame(max_rows = missing_value, Freq = 0))
        }
        
        clustercount<-clustercount[order(clustercount$max_rows),]
        
        
        clustercount<-cbind(clustercount,vec)
        clustercount<-clustercount[-c(1,nrow(clustercount)),]
        clustercount<-clustercount[order(clustercount$vec,decreasing = T),]
        #clustercount<-clustercount[-cell.div,]
        clustercount<-clustercount[clustercount$Freq>0,]
        mu.1=mu.simu
        clustercount$cumulative <- cumsum(clustercount$Freq)
        clustercount$x<-1/clustercount$vec
        fm2<-lm(cumulative~x,data=clustercount)
        slope<-summary(fm2)$coefficients[2,1]
        mu.2=slope/p
        datas[i,"vaf"]=cut
        datas[i,"p"]=lossfun(mu.simu,mu.2,0,0,1)
        
      }
      
      
      
      
 
    }
    if(div>=1){
      
      
      # Create an empty vector to store the results
      #expect.mu<-mu.simu*((p.f/(cut*2))+1)
      
      fit.div<-floor(1+s)
      close.vaf<-min((1-p.f)/(2*exp(log(2)*beta*div)),p.f/(2*exp(log(2)*beta*fit.div)))
      vec<-c(cut,close.vaf)
      vec<-data.frame(vec)
      vec$a<-depth*vec$vec
      vec$b<-depth-vec$a
      
      ##count1
      sub.include<-current.sub.vaf[current.sub.vaf>=min(vec$vec)]
     
      probs <- sapply(1:nrow(vec), function(i) dbeta(sub.include, vec$a[i], vec$b[i]))
      # Assign each number to the largest beta distribution
      max_rows <- apply(probs, 1, which.max)
      #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
      
      # Get the indices of rows where max_rows is equal to 2
      idx <- which(max_rows == 1)
      
      # Count the number of rows where max_rows is equal to 2
      count1 <- length(idx)
      
      
      
      
      ##count2
      
      idivg1=1
        pbeta_value <- pbeta(current.sub.vaf[current.sub.vaf < min(vec$vec)], vec$a[idivg1], vec$b[idivg1])
        vecdata<-data.frame(vaf=current.sub.vaf[current.sub.vaf < min(vec$vec)],p=pbeta_value)
        vecdata$p<-ifelse(vecdata$vaf>vec$a[idivg1]/(vec$a[idivg1]+vec$b[idivg1]),1-vecdata$p,vecdata$p)
        count2<- nrow(vecdata[vecdata$p>0.45,])
  
  
        mu.2=(count1+count2)/exp(log(2)*beta*(1+s))
      
      
      #p<-ifelse(count_max_rows<expect.mu,pnorm(count_max_rows,mean =expect.mu, sd=sqrt(expect.mu)),1-pnorm(count_max_rows,mean =expect.mu, sd=sqrt(expect.mu)))
      datas[i,"vaf"]=cut
      datas[i,"p"]=lossfun(mu.simu,mu.2,0,0,1)
      #datas[i,"p"]=(((mu.simu+1)-(mu.2+1)))^2
    }
     
      

      
    }else{
      datas[i,"vaf"]=cut
      datas[i,"p"]=Inf
    }
    datas[i,"s"]=s
    }
    
    
  
  
  datas1<-na.omit(datas)
  
  
  if(nrow(datas1)>0){
    index<-which(datas1$p ==min(datas1$p))[1]
    
    pick.s.vaf<-datas1[index,]
  }else{
    pick.s.vaf<-datas[1,]
  }
  return(pick.s.vaf)
  
}
beta_reassign<-function(df){
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
  
  # Create final VAF list based on frequency for each cluster
  final_df <- proportion_df %>%
    pivot_longer(cols = starts_with("prob"), names_to = "cluster", values_to = "freq") %>%
    mutate(cluster = as.integer(gsub("prob\\.", "", cluster))) %>% # Convert cluster to integer
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
generate_sequence <- function(n, s,set_max=F,set_max_num=0) {
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
    if (i %% n == 0) {
      # Update the current value according to the rule (1 + s) * (n - 1)
      current_val <- (1 + s) * (i/n)
    } else {
      # Otherwise, simply increment the current value by 1
      current_val <- current_val + 1
    }
    cell.list <- c(cell.list, current_val)
  }
  
  return(cell.list)
}
s_update_process<-function(give_n,df,cluster.result,vaf_set,min.s.detect,vaf.t1,depth,beta,p,p_thre=0.05){
  #print(cluster.result)
  s=cluster.result[cluster.result$cluster==give_n,]$new_s
  print(s)
  i=1
  #print(give_n)
  while(i<=100 ){
  round.cell.max<-exp(log(2)*beta*(1+s))
  max_n=give_n+2
  if(give_n>1){
  cell.list <-generate_sequence(give_n,s)
  }else{
    cell.list<-(1:max_n)*(1+s)
  }
  
  
  
  result_vector <- sapply(cell.list, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  mix_check<-sapply(result_vector, function(i) {
    # Replace this with your actual calculation
    calculated_value <- log((1-p)/2/i)/(log(2) * beta)
    calculated_value<-ifelse(calculated_value>0,1,0)
    return(calculated_value)
  })
  cluster_idx <- 1:length(mix_check)
  
  # Create a named vector for easier matching
  named_mix_vec <- setNames(mix_check, cluster_idx)
  
  inita=depth*c(vaf.t1,result_vector)
  initb=depth-inita
  probs <- sapply(1:length(inita), function(i) dbeta(vaf_set, inita[i], initb[i]))
  df<-data.frame(
    prob=probs,
    vaf=vaf_set
  )
  
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  df<-beta_reassign(df)
  df<-df[df$cluster>1,]
  # Sampling function
  #print(df)
  sample_or_not <- function(data, cluster_id,fraction, mix_vec) {
    #print(mix_vec)
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
  
  ##if s close to cell 2, 4, 8,(higher less possible)
  cluster.result<- sampled_df %>%
    group_by(cluster) %>%
    summarise(
      count = n(),
      mean_vaf = mean(vaf),
      min_vaf=min(vaf)
    )
  cluster.result<-cluster.result[cluster.result$cluster>1,]
  
  cluster.result$new_s<-log(p/(cluster.result$mean_vaf*2))/(log(2)*beta)-1
  #set_mins<-ifelse(min.s.detect<0.1,min.s.detect,0.1)
  cluster.result$alt_s<-log(p/(cluster.result$min_vaf*2))/(log(2)*beta)-1
  #set_mins<-ifelse(min.s.detect<0.1,min.s.detect,0.1)
  
  cluster.result$new_s<-ifelse(cluster.result$new_s<0,cluster.result$alt_s,cluster.result$new_s)
  cluster.result<-cluster.result[cluster.result$new_s<=min.s.detect,]
  cluster.result$cluster<-cluster.result$cluster-1
  update_s=cluster.result[cluster.result$cluster==give_n,]$new_s
  s=update_s
  if(abs(update_s-s)<p_thre){
    break
  }
  i=i+1
  }
  
  return(update_s)
  
}

estimate_s<-function(svalue.list,result,vaf.t1,p,depth,beta){
  n <- length(svalue.list)+1
  end.vaf<-p/(2*exp(log(2)*beta*(1+max(svalue.list))))
  vaf_set<-result[result$vaf.1<vaf.t1 & result$vaf.1 > end.vaf,]$vaf.1
  bic.list<-c()
  vaf.list<-c()
  start.list<-c()
  past_diff<-data.frame()
  
  for(i in 1:length(svalue.list)){
      s=svalue.list[i]
      s=ifelse(s>0,s,0)
      give_n=i
      if(give_n>1){
        cell.list <-generate_sequence(give_n,s,set_max=T,set_max_num=length(svalue.list))
      }else{
        cell.list<-(1:length(svalue.list))*(1+s)
      }
      #print(cell.list)
      result_vector <- sapply(cell.list, function(i) {
        # Replace this with your actual calculation
        calculated_value <- p/(2 * exp(log(2) * beta * i))
        return(calculated_value)
      })
      mix_check<-sapply(result_vector, function(i) {
        # Replace this with your actual calculation
        calculated_value <- log((1-p)/2*i)/(log(2) * beta)
        calculated_value<-ifelse(calculated_value>0,1,0)
        return(calculated_value)
      })
      cluster_idx <- 1:length(mix_check)
      
      # Create a named vector for easier matching
      named_mix_vec <- setNames(mix_check, cluster_idx)
      
      inita=depth*c(vaf.t1,result_vector)
      p_vec=c(vaf.t1,result_vector)
      initb=depth-inita
      probs <- sapply(1:length(inita), function(i) dbeta(vaf_set, inita[i], initb[i]))
      df<-data.frame(
        prob=probs,
        vaf=vaf_set
      )
      
      # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
      df<-beta_reassign(df)
      df <- df %>%
        mutate(p = p_vec[cluster])
      
      
      
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
      
      #sampled_df<-sampled_df[df$cluster<=max(i,3),]
      
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
      cluster.result$try<-i
      #print(cluster.result)
      # If you want to get the total log likelihood, you would sum these up
      total_log_likelihood <- sum(sampled_df$log_likelihood)
      bic<- compute_BIC(total_log_likelihood, 1, length(vaf_set))
      aic <- compute_AIC(total_log_likelihood, 1)
      bic.list<-c(bic.list,bic)
      vaf.list<-c(vaf.list,sum(abs(cluster.result$realvaf-cluster.result$mean_vaf))/nrow(cluster.result))
      if(i>1){
        end_row=i-1
        down=0
        for(iter in 1:end_row){
         past_diff.pick<-past_diff[past_diff$try==iter,]
         #current_row=iter
         diff_past <- past_diff.pick[1:iter,]$realvaf - past_diff.pick[1:iter,]$mean_vaf
         diff_cluster <- cluster.result[1:iter,]$realvaf - cluster.result[1:iter,]$mean_vaf
         
         # Compare the differences and get a binary vector (1 if true, 0 otherwise)
         binary_result <- ifelse(diff_past - diff_cluster >= 0, 1, 0)
         
         
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

Get_second_peak_ratio<-function(simu_vaf,vaf_list,cell_list,depth,pick.ratio=1){
  inita=depth*vaf_list
  initb=depth-inita
  probs <- sapply(1:length(inita), function(i) dbeta(simu_vaf, inita[i], initb[i]))
  df<-data.frame(
    prob=probs,
    vaf=simu_vaf
  )
  
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  df<-beta_reassign(df)
  ##if s close to cell 2, 4, 8,(higher less possible)
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
}

simulate_peak<-function(vaf_list,cell_list,vaf_set,depth){
  # Sample data
  num_decimal <- nchar(as.character(depth))
  # Assume mu and depth and num_decimal are defined in your environment
  #sizes <- generate_sizes(length(vaf_list))
  ratio=0
  lowerratio=0
  lower.pick=0.1
  for(i in 1:10){
   vaf.div<-c(1:length(vaf_list))
  all_simulated_vafs_list <- sapply(vaf.div, function(idx) {
    vaf <- vaf_list[idx]
    size <- cell_list[idx]*100
    round(rbinom(size, depth, vaf) / depth, num_decimal)
  })
  all_simulated_vafs <- unlist(all_simulated_vafs_list)
  all_simulated_vafs<-all_simulated_vafs[all_simulated_vafs>min(vaf_set)]
  ratio=ratio+Get_second_peak_ratio(all_simulated_vafs,vaf_list,cell_list,depth,pick.ratio=1)
  
  sampled_vafs_list <- lapply(vaf.div, function(idx) {
    vaf_vec=all_simulated_vafs_list[idx]
    sample_size <- ceiling(length(vaf_vec)/cell_list[idx])
    sample(vaf_vec, size = sample_size)
  })
  
  # If you want the result to be a vector (flattening the list), you can use:
  sampled_vafs_vector <- unlist(sampled_vafs_list)
  sampled_vafs_vector<-sampled_vafs_vector[sampled_vafs_vector>min(vaf_set)]
  lowerratio=lowerratio+Get_second_peak_ratio(sampled_vafs_vector,vaf_list,cell_list,depth,pick.ratio=1/cell_list[2])
  
  
  
  
  
  
  
  }
  ratio=ratio/10
  lowerratio=lowerratio/10
  return(c(ratio,lowerratio))
}

Simulate_ratio_peak<-function(result,p,beta,depth,s){
  ##note we have p, mu
  vaf.t1=p/2
  n <- floor(1+s)# Calculate the length of the sequence
  # Generate the sequence
  
  give_n=n
  if(give_n>1){
    cell.list <-generate_sequence(give_n,s,set_max=T,set_max_num=3)
  }else{
    cell.list<-(1:3)*(1+s)
  }
  
  vaf_set<-result[result$vaf.1<vaf.t1,]$vaf.1
  
  result_vector <- sapply(cell.list, function(i) {
      # Replace this with your actual calculation
      calculated_value <- p/(2 * exp(log(2) * beta * i))
      return(calculated_value)
    })
    
  
  vaf_list<-c(vaf.t1,result_vector)
  cell_list<-c(1,exp(log(2) * beta * cell.list))
  simu_ratio<-simulate_peak(vaf_list,cell_list,vaf_set,depth)
  
    mix_check<-sapply(result_vector, function(i) {
      # Replace this with your actual calculation
      calculated_value <- log((1-p)/2/i)/(log(2) * beta)
      calculated_value<-ifelse(calculated_value>0,1,0)
      return(calculated_value)
    })
    cluster_idx <- 1:length(mix_check)
    
    # Create a named vector for easier matching
    named_mix_vec <- setNames(mix_check, cluster_idx)
    
    
    #sUPPOSE it is from 2 4 6
    inita=depth*c(vaf.t1,result_vector)
    initb=depth-inita
    probs <- sapply(1:length(inita), function(i) dbeta(vaf_set, inita[i], initb[i]))
    df<-data.frame(
      prob=probs,
      vaf=vaf_set
    )
    
    # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
    df<-beta_reassign(df)
    ##if s close to cell 2, 4, 8,(higher less possible)
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
    
    cluster.result<- sampled_df %>%
      group_by(cluster) %>%
      summarise(
        count = n(),
        mean_vaf = mean(vaf),
        min_vaf=min(vaf)
      )
    mu.range<-cluster.result[cluster.result$cluster==2,]$count/simu_ratio/exp(log(2) * beta * cell.list[1])
  return(mu.range)
}

Get.S<-function(result,p,beta,depth,p_thre=0.05){
  
  ##note we have p, mu
  vaf.t1=p/2
  
  
  ##s would not be 2 unrealistic
  
  min.s.detect<-log(p/(quantile(result$vaf.1,0.01)*2))/(log(2)*beta)-1
  round.cell.max<-exp(log(2)*beta*(1+min.s.detect))
  n <- max(2,floor(1+min.s.detect)) # Calculate the length of the sequence
  # Generate the sequence
  
  cell.list <-1:n
  
  vaf_set<-result[result$vaf.1<vaf.t1 & result$vaf.1 > quantile(result$vaf.1,0.05),]$vaf.1
  #print(vaf_set)
  if(length(vaf_set)==0){
    vaf_set<-result[result$vaf.1<vaf.t1,]$vaf.1
    
  }
  if(length(vaf_set)>0){
  #p.now/(vafmin*2)>exp(log(2)*beta)
  
  # Check if the maximum value is not a power of 2, then append it
 
  result_vector <- sapply(cell.list, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  
  mix_check<-sapply(result_vector, function(i) {
    # Replace this with your actual calculation
    calculated_value <- log((1-p)/2/i)/(log(2) * beta)
    calculated_value<-ifelse(calculated_value>0,1,0)
    return(calculated_value)
  })
  cluster_idx <- 1:length(mix_check)
  
  # Create a named vector for easier matching
  named_mix_vec <- setNames(mix_check, cluster_idx)
  
  
  #sUPPOSE it is from 2 4 6
  inita=depth*c(vaf.t1,result_vector)
  initb=depth-inita
  probs <- sapply(1:length(inita), function(i) dbeta(vaf_set, inita[i], initb[i]))
  df<-data.frame(
    prob=probs,
    vaf=vaf_set
  )
  #print(df)
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  if(length(vaf_set)>1){
  df<-beta_reassign(df)}else{
#if(length(vaf_set)==0){
   df$cluster<-1  #??????
   }
  ##if s close to cell 2, 4, 8,(higher less possible)
  if(nrow(df[df$cluster>1,])>0){
  df<-df[df$cluster>1,]
  #}else{

  #print(df)
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
  
  cluster.result<- sampled_df %>%
    group_by(cluster) %>%
    summarise(
      count = n(),
      mean_vaf = mean(vaf),
      min_vaf=min(vaf)
    )
  cluster.result<-cluster.result[cluster.result$cluster>1,]
  cluster.result$new_s<-log(p/(cluster.result$mean_vaf*2))/(log(2)*beta)-1
  cluster.result$alt_s<-log(p/(cluster.result$min_vaf*2))/(log(2)*beta)-1
  #set_mins<-ifelse(min.s.detect<0.1,min.s.detect,0.1)
  
  cluster.result$new_s<-ifelse(cluster.result$new_s<0,cluster.result$alt_s,cluster.result$new_s)
  cluster.result<-cluster.result[cluster.result$new_s<=min.s.detect,]
  cluster.result$cluster<-cluster.result$cluster-1
  check_div<-1:nrow(cluster.result)
  check_div<-check_div[check_div %in% cluster.result$cluster]
  #if(nrow(cluster.result)==1){
  # check_div<-min(cluster.result$cluster)
  #}
  #print(check_div)

  svalue.list<-sapply(check_div, function(give_n) {
    # Replace this with your actual calculation
    svalue<-s_update_process(give_n,df,cluster.result,vaf_set,min.s.detect,vaf.t1,depth,beta,p,p_thre)
    return(svalue)
  })
  #print(svalue.list)
  if(length(svalue.list[svalue.list>0])>0){
  #svalue<-s_update_process(give_n,df,result_vector,vaf_set,min.s.detect,vaf.t1)
  estimate.s.data<-estimate_s(svalue.list,result,vaf.t1,p,depth,beta)
  estimate.s.data<-estimate.s.data[estimate.s.data$s>0,]
  }else{
    print('unable to process due to precision')
    estimate.s.data<-data.frame(
      s=NA,
      bic=NA,
      vaf=NA,
      start=NA
    )
  }
  }else{
    print('unable to process due to precision')
    estimate.s.data<-data.frame(
      s=NA,
      bic=NA,
      vaf=NA,
      start=NA
    )

    }
  }else{
    print('unable to process due to precision')                                 
    estimate.s.data<-data.frame(                                                
      s=NA,                                                                     
      bic=NA,                                                                   
      vaf=NA,                                                                   
      start=NA                                                                  
    )         
}
  return(estimate.s.data)
  
  
}

  ## ---- optimizeM.R ----
# Bundled reference: production estimation helpers. Loaded into a private
# environment by the dispatcher in zzz.R. All library/source lines stripped
# — `optimizeS.R` is sourced into the same env by the dispatcher BEFORE this
# file. Package dependencies (RBesT, dplyr, tidyr, likelihoodExplore) are
# declared in Imports.  `cliffDelta` (was rcompanion) is replaced by the
# package's `.cliffs_delta_abs` helper via an env-level alias the dispatcher
# installs.
round2 = function(x, digits) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
First_run_test<-function(first_vaf_list,second.vaf,give.vaf,dynamic,start.div,end.div,p,beta,depth,num_decimal,Min_Sample_size=6){
  
  #start.div=20
  #end.div=35
  
  collect.data=data.frame()
  idx<-1
  cell.div<-start.div:end.div
  vaf_set=c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  mu_list<-round(length(vaf_set)/cell.div)
  # Create a data frame to keep track of mu_list and corresponding cell.div
  mu_data <- data.frame(mu_list, cell.div)
  
  # Remove duplicates in mu_list and randomly select one corresponding cell.div
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
  div.list<-unique_df$cell.div
  uq=0.75
  lq=0.25
  if(length(unique_df)>=4){
  div.list<-div.list[div.list<=quantile(div.list, uq) & div.list >= quantile(div.list, lq)]
  }
  #print(unique_df$cell.div)
  for (cell.div in .cap_candidates(div.list)){
    #print(cell.div)
    possibleError <- tryCatch({
      if(dynamic){
        
        left_before_t1_vaf<-p/2+(1-p)/(2*exp(log(2)*beta*cell.div))
        a<-c(depth*give.vaf,depth*left_before_t1_vaf)
        b<-depth-a
        second.vaf.require.assign<-second.vaf[second.vaf<=left_before_t1_vaf]
        probs <- sapply(1:length(a), function(i) dbeta(second.vaf.require.assign, a[i], b[i]))
        
        #a<-c(give.vaf,left_before_t1_vaf)
        #b<-c(depth,depth)
        #probs <- sapply(1:length(a), function(i) rbinom(second.vaf.require.assign*depth, size = b[i], prob = a[i]))
        
        df<-data.frame(
          prob=probs,
          vaf=second.vaf.require.assign
        )
        
        df<-beta_reassign(df)
        last.div.set<-df[df$cluster==2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
        
        vaf_set=c(first_vaf_list,second.vaf[second.vaf>left_before_t1_vaf],last.div.set)
      }else{
        vaf_set=c(first_vaf_list,second.vaf[second.vaf>give.vaf])
      }
      mu_est=length(vaf_set)/cell.div
      if(mu_est>3){
        i_values <- 1:cell.div
        
        result_vector <- sapply(i_values, function(i) {
          # Replace this with your actual calculation
          calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
          return(calculated_value)
        })
        
        a<-depth*result_vector
        b<-depth-a
        probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
        df<-data.frame(
          prob=probs,
          vaf=vaf_set
        )
        
        
        df<-beta_reassign(df)
        df<-df[!duplicated(df),]
        
        mulist<-c()
        plist<-c()
        for(try in 1:3){
          mu_from_real<-normal_slope_get(df,p,beta,result_vector,num_decimal)
          #print(mu_from_real)
          #print(mu_est)
          mu_from_simu<-normal_slope_simu(cell.div,mu_est,p,depth,beta,num_decimal)
          #print(mu_from_simu)
          mean_list <- mean(mu_from_simu)
          sd_list <- max(sd(mu_from_simu),1.5)
          
          z_score <- (mu_from_real - mean_list) / sd_list
          mulist<-c(mulist,mu_from_real)
          plist<-c(plist,abs(z_score))
        }
        closest_index <- which.min(plist)
        mu_from_real <- mulist[closest_index]
        #closest_mu1 <- mulist1[which.min(abs(mulist1 - median(mulist1)))]
        z_score<-plist[closest_index]
      
      
    }},error=function(e){
      e
    })
    
    
    
    if(!inherits(possibleError, "error")){
      if(mu_est>3){
      collect.data[idx,"mu"]=mu_est
      collect.data[idx,"mu_real"]=mu_from_real
      collect.data[idx,"mu_simu"]=mean_list
      collect.data[idx,"cell.div"]=cell.div
      collect.data[idx,"z_score"]=z_score
      idx<-idx+1
      }else{
    collect.data[idx,"mu"]=NA                                                 
      collect.data[idx,"mu_real"]=NA                                            
      collect.data[idx,"mu_simu"]=NA                                            
      collect.data[idx,"cell.div"]=NA                                           
      collect.data[idx,"z_score"]=NA                                            
      idx<-idx+1  
        

    }
    }else{
      collect.data[idx,"mu"]=NA
      collect.data[idx,"mu_real"]=NA
      collect.data[idx,"mu_simu"]=NA
      collect.data[idx,"cell.div"]=NA
      collect.data[idx,"z_score"]=NA
      idx<-idx+1
    }
    
  }
  
  
  
  collect.data<-na.omit(collect.data)
  collect.data$p<-p
 
 
  return(collect.data)
}
Second_run_test<-function(first_vaf_list,second.vaf,give.vaf,sec.div.vaf,right_df,dynamic,celldivlist,p,beta,depth,num_decimal,Min_Sample_size=6,p_thre=0.05){

  collect.data2<-data.frame()
  idx<-1
  vaf_set=c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  for (cell.div in .cap_candidates(celldivlist)){
    
    collect.data2[idx,"cell.div"]=cell.div
    
    n_simulations <- getOption("teatime.peak_nsim", 100000)
    if(dynamic){
      
      left_before_t1_vaf<-p/2+(1-p)/(2*exp(log(2)*beta*cell.div))
      a<-c(depth*give.vaf,depth*left_before_t1_vaf)
      b<-depth-a
      second.vaf.require.assign<-second.vaf[second.vaf<=left_before_t1_vaf]
      probs <- sapply(1:length(a), function(i) dbeta(second.vaf.require.assign, a[i], b[i]))
      
      #a<-c(give.vaf,left_before_t1_vaf)
      #b<-c(depth,depth)
      #probs <- sapply(1:length(a), function(i) rbinom(second.vaf.require.assign*depth, size = b[i], prob = a[i]))
      
      df<-data.frame(
        prob=probs,
        vaf=second.vaf.require.assign
      )
      
      df<-beta_reassign(df)
      last.div.set<-df[df$cluster==2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
      
      vaf_set=c(first_vaf_list,second.vaf[second.vaf>left_before_t1_vaf],last.div.set)
    }else{
      vaf_set=c(first_vaf_list,second.vaf[second.vaf>give.vaf])
    }
    
    mu_suppose=length(vaf_set)/cell.div
    
    k<-round(mu_suppose)
    # Monte Carlo simulation to generate VAFs using the binomial distribution
    
    simulated_vafs.right <- round(rbinom(n_simulations, depth, sec.div.vaf) / depth,num_decimal)
    
    # For left_df
    if(k<Min_Sample_size){
      
      largest_right_values <- generate_bootstrap_samples(right_df[1:k,]$vaf,Min_Sample_size,num_decimal)
    }else{
      # For right_df
      
      largest_right_values <- right_df[1:k,]$vaf
    }
    largest_right_values <- na.omit(largest_right_values)
    largest_right_values <- as.vector(largest_right_values)
    # Perform a two-sample Kolmogorov-Smirnov test to get a p-value
    
    w.right.result <- wilcox.test(simulated_vafs.right, largest_right_values) 
    
    data1_long <- data.frame(value = simulated_vafs.right, group = "data1")
    data2_long <- data.frame(value = largest_right_values, group = "data2")
    combined_data <- rbind(data1_long, data2_long)
    
    cd.right=cliffDelta(value ~ group, data = combined_data)
    # Print p-value
    collect.data2[idx,"simu.right.mean"]=mean(simulated_vafs.right)
    collect.data2[idx,"right.mean"]=mean(largest_right_values)
    collect.data2[idx,"right.cd"]=abs(cd.right)
    collect.data2[idx,"right.p"]=w.right.result$p.value
    collect.data2[idx,"mu_est"]=mu_suppose
    idx<-idx+1
    
  }
  #print(collect.data2)
  collect.data2<-na.omit(collect.data2)
  Temp_keep<-collect.data2
  collect.data2<-Temp_keep[Temp_keep$right.p>p_thre,]
  #collect.data<-collect.data[collect.data$mu_est>3,]
  if(nrow(collect.data2)==0){
    collect.data2=Temp_keep  }
  
  # Convert specific columns to ranks and sum those ranks
  
  collect.data2$right_rank<-rank(collect.data2$right.cd)
  if(nrow(collect.data2)>2){
  top_10_percent=floor(0.5 * nrow(collect.data2))
    
    
  sorted_df <- collect.data2[order(collect.data2$right_rank), ]
  top_right <- head(sorted_df, top_10_percent)
  top_pick<-top_right
  }else{
    top_pick<-collect.data2
  }
  #print(top_pick)
  #compare_set,cell.div,mu,p,depth,beta=0.8,num_decimal
  top_pick_p <- apply(top_pick, 1, function(row) {
    cell_div_value <- row['cell.div']
    mu_value <- row['mu_est']
    i_values <- 1:cell_div_value
    result_vector <- sapply(i_values, function(i) {
      # Replace this with your actual calculation
      calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
      return(calculated_value)
    })
    result_vector<-c(0.5,p/2,result_vector)
    
    a<-depth*result_vector
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=vaf_set
    )
    
    df<-beta_reassign(df)
    data<-df[df$cluster>=2,]$vaf
    simulate_from_estimation(data, cell_div_value, mu_value, p, depth, beta,num_decimal)
  })
  
  top_pick_p <- as.data.frame(t(top_pick_p))
  top_pick_p$V1<-ifelse(top_pick_p$V1<0.05,0,top_pick_p$V1)
  top_pick_p$V2<-ifelse(top_pick_p$V2<0.05,0,top_pick_p$V2)
  top_pick$r1<-top_pick_p$V1
  top_pick$r2<-top_pick_p$V2
  
  top_pick$loglike<-top_pick_p$V3
  top_pick$aic<-top_pick_p$V4
  top_pick$bic<-top_pick_p$V5
  
  
  top_pick$top_pick_score<-top_pick$r2+1/top_pick$right_rank
  #top_pick<- top_pick[order(top_pick$top_pick_score), ]
  if(nrow(top_pick)>1){
  top_pick<- top_pick[top_pick$top_pick_score>median(top_pick$top_pick_score), ]
  }
  top_pick <- top_pick[order(top_pick$right_rank), ]
  pick_mu_cell.div<-top_pick[1:min(10,nrow(top_pick)),]
  
  #print(top_pick[abs(top_pick$mu_est-16)==min(abs(top_pick$mu_est-16)),])
  #print(pick_mu_cell.div)
  return(pick_mu_cell.div)

}
Iterate_P_optimize<-function(clear,give.vaf,upper_clonal_vaf,clonal.vaf.left,second.vaf,depth,beta,num_decimal,dynamic=F,p_thre=0.05){
  p=give.vaf*2
  reliable=1
  suppose_right_vaf=p/2+(1-p)/(2*exp(log(2)*beta))
  if(clear){
    first_vaf_list<-c()
    #first_vaf_list<-second.vaf[second.vaf>give.vaf]
  }else{
  if(suppose_right_vaf > upper_clonal_vaf){
    #a<-c(depth*suppose_right_vaf)
    #b<-c(depth-a)
    #probs <-pbeta(clonal.vaf.left, a, b)
    #combined_df <- data.frame(vaf = clonal.vaf.left, probs = probs)
    #freq_df <- combined_df %>% 
    #  group_by(vaf, probs) %>% 
    #  summarise(count = n()) %>% 
    #  mutate(freq_to_assign = round(count * probs)) %>% 
    #  ungroup()
    
    # Expand the dataset based on calculated frequencies
   # first_vaf_list <- freq_df %>% 
    #  uncount(freq_to_assign, .remove = FALSE) %>% 
    #  pull(vaf)
    
    a<-c(depth*0.5,depth*suppose_right_vaf)
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(clonal.vaf.left, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=clonal.vaf.left
    )
    
    df<-beta_reassign(df)
    first_vaf_list<-df[df$cluster==2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
    
  }else{
    a<-c(depth*upper_clonal_vaf,depth*suppose_right_vaf)
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(clonal.vaf.left, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=clonal.vaf.left
    )
    
    df<-beta_reassign(df)
    first_vaf_list<-df[df$cluster==2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
    
  }
  }
  vaf_set<-c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  
  Min_Sample_size=6
  i_values <- 1:20
  result_vector <- sapply(i_values, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  
  a<-depth*result_vector
  b<-depth-a
  probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
  df<-data.frame(
    prob=probs,
    vaf=vaf_set
  )
 
  df<-beta_reassign(df)
  df<-df[!duplicated(df),]
  df_count_freq <- df %>%
    group_by(cluster) %>%
    summarise(count = sum(freq), .groups = 'drop')
  if(min(df_count_freq$cluster)>1){
    new_row <- data.frame(cluster = 1, count = 0)
    df_count_freq <- rbind(new_row, df_count_freq)
    
    # Sort the dataframe by cluster
    df_count_freq <- df_count_freq[order(df_count_freq$cluster), ]
  }
  #print(df_count_freq)
  if(nrow(df_count_freq)>2){
  max_mu=df_count_freq[1,]$count+df_count_freq[2,]$count
  #print(max_mu)
  
  if(nrow(df_count_freq)>6 & df_count_freq[2,"count"]/df_count_freq[3,"count"]<=4){
    df_count_freq<-df_count_freq[2:3,]
  }else{
    df_count_freq<-df_count_freq[2,]
  }
  
  dynamic=F
  total_count_temp=length(first_vaf_list)+length(second.vaf[second.vaf>give.vaf])
  min_mu=max(3,min(df_count_freq$count/4))
  #min_mu=max(3,min(df_count_freq$count/2))
  #print(min_mu)
  }else{
   ###??????
    total_count_temp=length(first_vaf_list)+length(second.vaf[second.vaf>give.vaf])
   max_mu=df_count_freq[1,]$count
   min_mu=max(3,min(df_count_freq$count/4)) 
  }
  if(max_mu < min_mu){
    temp=min_mu
    min_mu=max_mu
    max_mu=temp
  }
  start.div=round(total_count_temp/max_mu)
  end.div<-round(total_count_temp/min_mu)
  #print(c(start.div,end.div))
  if(length(first_vaf_list)==0){
  first_vaf_list<-second.vaf[second.vaf>give.vaf]
   }
  collect.data<-First_run_test(first_vaf_list,second.vaf,give.vaf,dynamic,start.div,end.div,p,beta,depth,num_decimal,Min_Sample_size=6)
  if(nrow(collect.data[collect.data$mu_real>3,])>0){
  collect.data<-collect.data[collect.data$mu_real>3,]
  }else{
   collect.data<-collect.data[collect.data$mu_real>1,]
  }
  #print('first done')
  #collect.data<-collect.data[abs(collect.data$z_score)<=1.96,]
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
  possibleError <- tryCatch({
  collect.data2<-Second_run_test(first_vaf_list,second.vaf,give.vaf,sec.div.vaf,right_df,dynamic,celldivlist,p,beta,depth,num_decimal,Min_Sample_size=6,p_thre)
   },error=function(e){
      e
    })
      
if(!inherits(possibleError, "error")){
  #print('Second done')
  if(nrow(collect.data2)==1){
    range1=collect.data2$mu_est-1
    range2=collect.data2$mu_est+1
  }else{
    range1=min(collect.data2$mu_est)
    range2=max(collect.data2$mu_est)
  }
  }else{
   range1=0
   range2=max(collect.data$mu_real)+1
  }
  over12<-collect.data[collect.data$mu_real>=range1 & collect.data$mu_real<=range2,]
  overlap.pick<-NULL
  if(nrow(over12)>0){
    overlap.pick<-over12
  }
  
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

 possibleError <- tryCatch({
  collect.data3<-Second_run_test(first_vaf_list,second.vaf,give.vaf,sec.div.vaf,right_df,dynamic,celldivlist,p,beta,depth,num_decimal,Min_Sample_size=6,p_thre)
  },error=function(e){
      e
    })

  if(!inherits(possibleError, "error")){
  if(nrow(collect.data2)==1){
    range1=collect.data3$mu_est-1
    range2=collect.data3$mu_est+1
  }else{
    range1=min(collect.data3$mu_est)
    range2=max(collect.data3$mu_est)
  }
  }else{
range1=0
range2=max(collect.data$mu_real)+1
}
  #print(range1)
  #print(range2)
  over23<-collect.data[collect.data$mu_real>=range1 & collect.data$mu_real<=range2,]
  overlap_values <- intersect(over12$mu_real, over23$mu_real)
  #print(over23)
  #print(over12)
  #print(overlap_values)
  if(nrow(over23)>0 | nrow(over12)>0){
  if(length(overlap_values) > 0){
    overlap.pick<-over23[over23$mu_real %in% overlap_values, ]
  }else{
    overlap.pick<-rbind(over12,over23)
  }
  }
  }

  #print(overlap.pick)
  ##check three overlap
  if(is.null(overlap.pick)){
    overlap.pick<-collect.data
  }else if(nrow(na.omit(overlap.pick))==0){
   overlap.pick<-collect.data
  }
  overlap.pick<-na.omit(overlap.pick)
  }else{
    overlap.pick<-collect.data
    reliable=0
  }
  #print(overlap.pick)
  #print('Overcheck done')
  #when mu <=3
  vaf_set=c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  collect.data4<-mu_find_small(vaf_set,p,depth, beta,num_decimal,test_mu=min_mu)
  
  mu_small_selection=collect.data4[collect.data4$p_value>0.05 & collect.data4$wx_p>0.05 ,]
  
  vaf_set<-c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  
  pick.cell.div=overlap.pick$cell.div
  pick.mu=overlap.pick$mu_real
  z_score=overlap.pick$z_score
  
  #print(overlap.pick)
  top_pick_p <- apply(overlap.pick, 1, function(row) {
    cell_div_value <- row['cell.div']
    mu_value <- row['mu_real']
    i_values <- 1:cell_div_value
    result_vector <- sapply(i_values, function(i) {
      # Replace this with your actual calculation
      calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
      return(calculated_value)
    })
    result_vector<-c(0.5,p/2,result_vector)
    
    a<-depth*result_vector
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=vaf_set
    )
    
    df<-beta_reassign(df)
    data<-df[df$cluster>=2,]$vaf
    simulate_from_estimation(data, cell_div_value, mu_value, p, depth, beta,num_decimal)
  })
  top_pick <- as.data.frame(t(top_pick_p))
  #overlap.pick$loglike<-
  overlap.pick$loglike<-top_pick$V3
  overlap.pick$aic<-top_pick$V4
  overlap.pick$bic<-top_pick$V5
  pick.log<-overlap.pick$loglike
  pick.bic<-overlap.pick$bic
  pick.aic<-overlap.pick$aic

  if(nrow(mu_small_selection)>0){
    small_mu_pick <- apply(mu_small_selection, 1, function(row) {
      mu_value <- row['mu_est']
      cell_div_value <- round(length(vaf_set)/mu_value)
      simulate_from_estimation(vaf_set, cell_div_value, mu_value, p, depth, beta,num_decimal)
    })
    small_mu_pick <- as.data.frame(t(small_mu_pick))
    small_mu_pick$V1<-ifelse(small_mu_pick$V1<0.05,0,top_pick_p$V1)
    small_mu_pick$V2<-ifelse(small_mu_pick$V2<0.05,0,top_pick_p$V2)
    
    mu_small_selection$r1<-small_mu_pick$V1
    mu_small_selection$r2<-small_mu_pick$V2
    
    mu_small_selection$loglike<-small_mu_pick$V3
    mu_small_selection$aic<-small_mu_pick$V4
    mu_small_selection$bic<-small_mu_pick$V5
    mu_small_selection<-mu_small_selection[mu_small_selection$bic==min(mu_small_selection$bic),]
    if(mu_small_selection$bic>min(pick_mu_cell.div$bic)*2){
      pick.mu=mu_small_selection$mu_est
      pick.bic=mu_small_selection$bic
      pick.log=mu_small_selection$loglike
      pick.aic<-mu_small_selection$aic
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
calculate_mu_high <- function(G1, beta) {
 
  possibleError <- tryCatch({
    mu.turn <-breakpoints(cumsum~x,data=G1,h=3/nrow(G1))
    
  },error=function(e){
    e
  })
  if(!inherits(possibleError, "error")){
    #REAL WORK
    mu.turn <-breakpoints(cumsum~x,data=G1,h=3/nrow(G1))
    #fm1 <- lm(cumsum ~ breakfactor(mu.turn), data = G1)
    # Extract the coefficients from the model
    bf <- breakfactor(mu.turn)
    
    # Initialize a vector to store AIC values
    AIC_values <- c()
    BIC_values<- c()
    slopes<-c()
    # Loop over each segment and calculate the AIC
    for(i in unique(bf)) {
      segment_data <- subset(G1, bf == i)
      fm_segment <- lm(cumsum ~ x, data = segment_data)
      coeffs <- coef(fm_segment)
      slope <- coeffs[-1]*(-1)*beta*log(2)
      slopes<-c(slopes,slope)
      AIC_values <-c(AIC_values, AIC(fm_segment))
      BIC_values<-c(BIC_values, BIC(fm_segment))
    }
    
    
    #coeffs <- coef(fm1)
    #slopes <- coeffs[-1]*(-1)*beta*log(2)
    
    best_segment <- which.min(AIC_values)
    best_segment1 <- which.min(BIC_values)
    mu_1=(slopes[best_segment]+slopes[best_segment1])/2
  }else{
    fm1<-lm(cumsum~x,data=G1)
    mu_1=summary(fm1)$coefficients[2,1]*(-1)*beta*log(2)
  }
  return(mu_1)  # Return mu_high
}

distribute_remainder <- function(row) {
  total <- sum(row[starts_with("prob")])
  freq_rounded <- floor(row[starts_with("prob")] / total * row$n)
  
  if(sum(freq_rounded) == 0) {
    freq_rounded[order(-row[starts_with("prob")], decreasing = TRUE)[seq_len(row$n)]] <- 1
  } else {
    diff <- row$n - sum(freq_rounded)
    freq_to_add <- replace(rep(0, length(freq_rounded)), order(-row[starts_with("prob")], decreasing = TRUE)[seq_len(diff)], 1)
    freq_rounded <- freq_rounded + freq_to_add
  }
  
  return(freq_rounded)
}

normal_slope_get<-function(df,p,beta,result_vector,num_decimal){
  
  G1<-df%>% 
    group_by(cluster) %>%
    summarise(count = sum(freq))
  
  G1<-G1[order(G1$cluster),]
  
  G1 <- G1 %>% 
    arrange(cluster) %>% # Sort by 'cluster' if it's not already sorted
    mutate(vaf = result_vector[cluster]) # Assuming 'result_vector' and 'cluster' align in terms of index.
  
  
  G1.sub <- G1 %>%
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
  total_rows <- nrow(G1)
  G1<-G1[order(G1$cluster),]
  if(nrow(G1)>=6){
  
  #G1<-G1[floor(total_rows * 0.1) + 1:floor(total_rows * 0.9),]
  G1<-G1[2:floor(total_rows * 0.5),]
  }
  
  #else{
  #  G1<-G1[floor(total_rows * 0.05) + 1:floor(total_rows * 0.95),]
  #}
  if(nrow(G1)>2){
  mu_1=calculate_mu_high(G1, beta)
  }else{
  #print('too few point, not reliable')
  ####???????
  mu_1=3
  }
  return(mu_1)
  
}
normal_slope_simu<-function(cell.div,mu,p,depth,beta,num_decimal){
  i_values <- 1:cell.div
  
  result_vector <- sapply(i_values, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  a<-depth*result_vector
  b<-depth-a
   # Pre-allocate storage
  n_try <- 3
  mu_est_list <- numeric(n_try)
  
  # This is constant within the loop, so compute once
  result_matrix <- matrix(result_vector, nrow=length(result_vector), ncol=round(mu), byrow=TRUE)
  depth_matrix <- matrix(depth, nrow=length(result_vector), ncol=round(mu))
  
  # Loop to generate simulated values
  for (try in 1:n_try) {
    # Use matrix operations for speed
    all_simulated_vafs_list <- sapply(result_vector, function(vaf) {
           round(rbinom(round(mu), depth, vaf) / depth, num_decimal)
          }) 
    #all_simulated_vafs <- as.vector(t(all_simulated_vafs_matrix))
    all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
    probs <- sapply(1:length(a), function(i) dbeta(all_simulated_vafs, a[i], b[i]))
    
    df <- data.frame(prob = probs, vaf = all_simulated_vafs)
    #print(head(df))
    # Assuming beta_reassign works with data.frame
    df <- beta_reassign(df)
    df <- unique(df)  # removes duplicates
    
    mu_est_list[try] <- normal_slope_get(df, p, beta, result_vector, num_decimal)
    
  }
  
#  for (try in 1:10){
    # Simulate VAFs for each element in vaf_list and round them
#    all_simulated_vafs_list <- sapply(result_vector, function(vaf) {
#      round(rbinom(round(mu), depth, vaf) / depth, num_decimal)
#    })
#    all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
#    probs <- sapply(1:length(a), function(i) dbeta(all_simulated_vafs, a[i], b[i]))
#    df<-data.frame(
#      prob=probs,
#      vaf=all_simulated_vafs
#    )
##    
#    df<-beta_reassign(df)
 #   df<-df[!duplicated(df),]
    # Concatenate all vectors into a single vector
 #   mu_est_list<-c(mu_est_list,normal_slope_get(df,p,beta,result_vector,num_decimal))
 # }
  #print(mu_est_list)
  return(mu_est_list)
}
Get_slope<-function(vaf_set,p,beta=0.8){
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
  
  #G1.sub$cumsum<-G1.sub$cumsum-init.mu.simu
  #G1.sub<-G1.sub[G1.sub$cumsum>0,]
  #print(vaf.list)
  G1.sub$vaf<-as.numeric(G1.sub$vaf)
  
  G1.sub$lnf<-2*G1.sub$vaf-p
  G1<-G1.sub[G1.sub$lnf>0,]
  G1$x<-log(G1$lnf)
  
 
  num_rows_to_remove <- round(nrow(G1) * 0.05)
  
  # Calculate the index of the rows to remove
  rows_to_remove <- c(1:num_rows_to_remove, (nrow(G1) - num_rows_to_remove + 1):nrow(G1))
  
  
  # Subset the data frame to remove the rows
  if(nrow(G1[-rows_to_remove,]) > 5){
    G1 <- G1[-rows_to_remove,]
  }else{
    G1=G1
  }
  mu_high=calculate_mu_high(G1, beta)
  

  return(mu_high)
}
beta_reassign<-function(df){
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
Get.mean<-function(numbers,cutlen){
  # Create an empty list to hold the subsets
  subsets <- list()
  
  # Create a temporary list to hold each subset of 4 numbers
  temp_list <- c()
  
  # Loop through the sorted numbers
  for (i in 1:length(numbers)) {
    # Add the current number to the temporary list
    temp_list <- c(temp_list, numbers[i])
    
    # Check if the length of the temporary list is 4
    if (length(temp_list) == cutlen) {
      # Calculate the mean of the temporary list and append it to the list of subsets
      subset_mean <- mean(temp_list)
      subsets <- c(subsets, subset_mean)
      
      # Clear the temporary list
      temp_list <- c()
    }
  }
  
  # If there are any remaining numbers in the temporary list, calculate the mean and append it to the list of subsets
  if (length(temp_list) > 0) {
    subset_mean <- mean(temp_list)
    subsets <- c(subsets, subset_mean)
  }
  
  # Print the list of subset means
  return(subsets)
}
Get.result<-function(G1.sub,p,M.data,sub.vaf,mago.result,depth,checkp,beta){
  print("*******getting mu*********")
  if(p>0.4 & p<0.5){
    M.data<-c(M.data[M.data>=p/2],sub.vaf[sub.vaf>=p/2])
  }
  G1.sub<-data.frame(M.data)
  colnames(G1.sub)<-"vaf"
  
  G1.sub$count<-1
  
  G1.sub <- G1.sub %>%
    mutate(vaf = format(vaf, nsmall = 3)) %>%
    arrange(desc(vaf), desc(count)) %>%
    mutate(cumsum = cumsum(count)) %>%
    group_by(vaf) %>%
    slice_max(cumsum) %>%
    ungroup()
  
  #G1.sub$cumsum<-G1.sub$cumsum-init.mu.simu
  #G1.sub<-G1.sub[G1.sub$cumsum>0,]
  #print(vaf.list)
  G1.sub$vaf<-as.numeric(G1.sub$vaf)
    
    G1.sub$lnf<-2*G1.sub$vaf-p
    G1<-G1.sub[G1.sub$lnf>0,]
    G1$x<-log(G1$lnf)
    #plot(G1$x,G1$cumsum)
    vaf.first<-p/2 + (1-p)/(2*2^beta)
    compare.vec<-c(0.5,vaf.first)
    compare.vec<-sort(compare.vec)
    #print(vec)
    # Find the indices of the closest values in "vec" for each element in "assgin.vaf"
    compare.vec<-data.frame(compare.vec)
    compare.vec$a<-depth*compare.vec$compare.vec
    compare.vec$b<-depth-compare.vec$a
    probs <- sapply(1:nrow(compare.vec), function(i) dbeta(M.data, compare.vec$a[i], compare.vec$b[i]))
    
    # Assign each number to the largest beta distribution
    max_rows <- apply(probs, 1, which.max)
    #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
    
    cluster.count <- table(max_rows)[nrow(compare.vec)]
    M.cul<-length(M.data)-cluster.count
    
    initial.vaf.cut<-G1[G1$cumsum==cluster.count,]
    assgin.vaf<-M.data[M.data<initial.vaf.cut$vaf]
    #print(length(assgin.vaf))
    assgin.vaf<-sort(assgin.vaf,decreasing = T)
    G1<-G1[G1$cumsum>cluster.count,]
    # Calculate the number of rows to remove from each end
    num_rows_to_remove <- round(nrow(G1) * 0.05)
    
    # Calculate the index of the rows to remove
    rows_to_remove <- c(1:num_rows_to_remove, (nrow(G1) - num_rows_to_remove + 1):nrow(G1))
    
    
    # Subset the data frame to remove the rows
    if(nrow(G1[-rows_to_remove,]) > 2){
    G1 <- G1[-rows_to_remove,]
    }else{
      G1=G1
    }
    # Subset the data frame to remove the rows
    fm2<-lm(cumsum~x,data=G1)
    #plot(G1$x,G1$cumsum)
    r.n<-summary(fm2)$r.squared
    count=1
    #plot(G1$x,G1$cumsum)
    tryCatch({
      mu.turn <-breakpoints(cumsum~x,data=G1,h=3/nrow(G1))
      
      #fm1 <- lm(cumsum ~ breakfactor(mu.turn), data = G1)
      # Extract the coefficients from the model
      
      coeffs <- coef(fm1)
      slopes <- coeffs[-1]*(-1)*beta*log(2)
      bf <- breakfactor(mu.turn)
      
      # Initialize a vector to store AIC values
      AIC_values <- c()
      BIC_values<- c()
      slopes<-c()
      # Loop over each segment and calculate the AIC
      for(i in unique(bf)) {
        segment_data <- subset(G1, bf == i)
        fm_segment <- lm(cumsum ~ x, data = segment_data)
        coeffs <- coef(fm_segment)
        slope <- coeffs[-1]*(-1)*beta*log(2)
        slopes<-c(slopes,slope)
        AIC_values <-c(AIC_values, AIC(fm_segment))
        BIC_values<-c(BIC_values, BIC(fm_segment))
      }
      
      
      #coeffs <- coef(fm1)
      #slopes <- coeffs[-1]*(-1)*beta*log(2)
      
      best_segment <- which.min(AIC_values)
      best_segment1 <- which.min(BIC_values)
      mu_high=(slopes[best_segment]+slopes[best_segment1])/2
      count=length(slopes)
    },error=function(e){
      count=1
      fm1<-lm(cumsum~x,data=G1)
      mu_high=summary(fm1)$coefficients[2,1]*(-1)*beta*log(2)
      #mu_high=
    })
    
    max_count=max(count+3,10)
    min_count=1
    
    # Store the slopes for each breakpoint
    
    min.error<-Inf
    magocheck=nrow(mago.result[mago.result$vaf >0.25,])
    if(magocheck<=2 | p>0.8 ){
      print("*****Special optimize for Overwhlem Fitness")
      #stun=1
      for (cell.div.c in min_count:max_count){
        #print(cell.div.c)
        cell.div=cell.div.c
        mu.1=M.cul/cell.div
        if(mu.1>cluster.count){
          ptest<-pnorm(mu.1,mean = cluster.count,sd=sqrt(cluster.count))}else{
            ptest<-0.01
          }
        vec <- numeric(cell.div)
        if(cell.div<=1){
          cell.num<-exp(log(2)*beta*1)
          vafmut<-p/2 + (1-p)/(2*cell.num)
          
          vec<-data.frame(vafmut)
          
          
          vec$real<-mean(assgin.vaf)
          
          error=abs(vec$vafmut-vec$real)/vec$vafmut+ptest
        }else{
          for (j in 1:cell.div) { 
            cell.num<-exp(log(2)*beta*j)
            if(p/2 + (1-p)/(2*cell.num) <min(assgin.vaf)){
              vec[j]<-p/2 + (1-p)/(2*cell.num)
            }else{
              vec[j] <- p/2 + (1-p)/(2*cell.num)  # Calculate the value for the current division number and store it in the vector
            }     
          }
          if((p/2) >=min(assgin.vaf)){
            vec[j]<-p/2
            #error=Inf
            #mu.1=round2(M.cul/cell.div,0)
            
          }
          vec<-vec[which(vec != 0)]
          if((1-p)/(2*2) >=min(assgin.vaf)){
            #vec[j+1]<-(1-p)/(2*2^beta)
            error=Inf
          }else{
            #vec<-sort(vec)
            if(length(vec)<cell.div){
              error=Inf
            }else{
              vec<-data.frame(vec)
              
              
              realmean<-Get.mean(assgin.vaf,round2(mu.1,0))
              realmean<-unlist(realmean)
              diff.real.exp<-abs(length(realmean)-nrow(vec))
              if(diff.real.exp>0 & diff.real.exp<=2){
                if(length(realmean)>nrow(vec)){
                  realmean[nrow(vec)]<-mean(realmean[nrow(vec):length(realmean)])
                  realmean<-realmean[1:nrow(vec)]
                }else{
                  realmean[length(realmean):nrow(vec)]<-0
                }
                
              }
              
              if(diff.real.exp>2){
                realmean[1:nrow(vec)]<-0
              }
              vec$real<-realmean
              vec$diff<-abs(vec$vec-vec$real)/vec$vec
              #print(vec)
              
              error=sum(vec$diff)/cell.div+ptest
            }
            
            
          }
          
        }
        #print(error)
        
        if(error<min.error & is.finite(error)){
          min.error=error
          mu=mu.1
        }
        
        
        
      }
      # Store the slopes for each breakpoint
      
     
    }else{ 
      r.list<-c()
      Fail.cell<-c()
      contp=F
      if(p>=0.4 & p< checkp & magocheck >2 & depth>=50){
        #print("*******Using magos results*******")
        
        magomu=mago.result[2,"sum"]
        
        Temp.G1<-mago.result[2:magocheck,c("vaf","sum")]
        Temp.G1$vaf<-as.numeric(Temp.G1$vaf)
        Temp.G1$cumsum<-cumsum(Temp.G1$sum)
        Temp.G1$lnf<-2*Temp.G1$vaf-p
        Temp.G1<-Temp.G1[Temp.G1$lnf>0,]
        Temp.G1$x<-log(Temp.G1$lnf)
        if(nrow(Temp.G1)>1){
        
        fm2<-lm(cumsum~x,data=Temp.G1)
        mu.2=summary(fm2)$coefficients[2,1]*(-1)*beta*log(2)
        mu.2=min(mu.2,cluster.count)
        
        min.error=lossfun(magomu,mu.2,0,1,magocheck-1)
          mu=magomu
          
        }else{
         
          min.error=0.5*M.cul+0.5
            mu=magomu
        }
        contp=F
        
        
      }
    for (cell.div.c in min_count:max_count){
      #print(cell.div.c)
      cell.div=cell.div.c
      mu.1=M.cul/cell.div
      vec <- numeric(cell.div)
      if(cell.div<=1){
        cell.num<-exp(log(2)*beta*1)
        vafmut<-p/2 + (1-p)/(2*cell.num)
        a<-depth*vafmut
        b<-depth-a
        
        ##count1
        pbeta_value <- pbeta(assgin.vaf, a,b)
        vecdata<-data.frame(vaf=assgin.vaf,p=pbeta_value)
        vecdata$p<-ifelse(vecdata$vaf>a/(a+b),1-vecdata$p,vecdata$p)
        vecdata$p<-p.adjust(vecdata$p,method='BH')
        mu.3=0
        r.n=0.1
        mu.2=min(nrow(vecdata[vecdata$p>0.1,]),cluster.count)
        if(contp){
          ptest=ifelse(mu.1>magomu,1-pnorm(mu.1,mean=magomu,sd=sqrt(magomu)),pnorm(mu.1,mean=magomu,sd=sqrt(magomu)))
        }else{
          ptest<-0
        }
        error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)+1-ptest
      }else{
        for (j in 1:cell.div) { 
          cell.num<-exp(log(2)*beta*j)
          if(p/2 + (1-p)/(2*cell.num) <min(assgin.vaf)){
            vec[j]<-p/2 + (1-p)/(2*cell.num)
          }else{
            vec[j] <- p/2 + (1-p)/(2*cell.num)  # Calculate the value for the current division number and store it in the vector
          }     
        }
        if((p/2) >=min(assgin.vaf)){
          vec[j]<-p/2
          #error=Inf
          #mu.1=round2(M.cul/cell.div,0)
          
        }
        vec<-vec[which(vec != 0)]
        if((1-p)/(2*2) >=min(assgin.vaf)){
          #vec[j+1]<-(1-p)/(2*2^beta)
          error=Inf
        }else{
          #vec<-sort(vec)
          if(length(vec)<cell.div){
            error=Inf
          }else{
            vec<-data.frame(vec)
            vec$a<-depth*vec$vec
            vec$b<-depth-vec$a
            #result<-test_distribution(assgin.vaf,vec)
            probs <- sapply(1:nrow(vec), function(i) dbeta(assgin.vaf, vec$a[i], vec$b[i]))
            
            # Assign each number to the largest beta distribution
            max_rows <- apply(probs, 1, which.max)
            #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
            
            clustercount <- data.frame(table(max_rows))
            clustercount$max_rows<-as.numeric(as.character(clustercount$max_rows))
            missing_value <- setdiff(1:cell.div, clustercount$max_rows)
            
            # Add the missing value to the dataframe with a frequency of 0
            if(length(missing_value) > 0) {
              clustercount <- rbind(clustercount, data.frame(max_rows = missing_value, Freq = 0))
            }
            clustercount<-clustercount[order(clustercount$max_rows),]
            clustercount<-cbind(clustercount,vec)
            clustercount<-clustercount[order(clustercount$vec,decreasing = T),]
            clustercount<-clustercount[clustercount$Freq>0,]
            mu.1=M.cul/nrow(clustercount)
            clustercount$cumulative <- cumsum(clustercount$Freq)
            clustercount$x<-2*clustercount$vec-p
            clustercount<-clustercount[clustercount$x>0,]
            clustercount$x<-log(clustercount$x)
            
            if(nrow(clustercount)<=1){
              
              error=Inf
            }
            if(nrow(clustercount)==2){
              Fail.cell<-c(Fail.cell,cell.div.c)
              r.n=0
              fm2<-lm(cumulative~x,data=clustercount)
              mu.2=summary(fm2)$coefficients[2,1]*(-1)*beta*log(2)
              mu.2=min(mu.2,cluster.count)
              mu.3=summary(fm2)$coefficients[1,1]*log(2)/log(1-p)
              if(contp){
                ptest=ifelse(mu.1>magomu,1-pnorm(mu.1,mean=magomu,sd=sqrt(magomu)),pnorm(mu.1,mean=magomu,sd=sqrt(magomu)))
              }else{
                ptest<-0
              }
              error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)+1-ptest
            }
            if(nrow(clustercount)>2){
              fm2<-lm(cumulative~x,data=clustercount)
              #plot(G1$x,G1$cumsum)
              mu.2=summary(fm2)$coefficients[2,1]*(-1)*beta*log(2)
              mu.2=min(mu.2,cluster.count)
              mu.3=summary(fm2)$coefficients[1,1]*log(2)/log(1-p)
              
              r.n<-summary(fm2)$r.squared
              r.list<-c(r.list,r.n)
              if(contp){
                ptest=ifelse(mu.1>magomu,1-pnorm(mu.1,mean=magomu,sd=sqrt(magomu)),pnorm(mu.1,mean=magomu,sd=sqrt(magomu)))
              }else{
                ptest<-0
              }
              error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)+1-ptest
            }
            #print(mu.2)
            #print(r.n)
            #error=((log(mu.1+1)-log(mu.2+1)))^2+(1-r.n)^2
            #print(cell.div.c)
            
            #error=(((mu.1+1)-(mu.2+1)))^2
          }
          
        }
        
      }
      
      #print(mu.1,mu.2,r.n)
     #print(mu.2)
      #print(r.n)
      #print(error)
      #
      if(error<min.error & is.finite(error)){
        min.error=error
        mu=mu.1
      }
      
      
      
    }
    # Store the slopes for each breakpoint
    if(length(r.list)<1){
      r.list<-c(r.list,0.8)
    }
    if(length(Fail.cell)>0){
      for (cell.div.c in Fail.cell){
        cell.div=cell.div.c
        mu.1=M.cul/cell.div
        vec <- numeric(cell.div)
        for (j in 1:cell.div) { 
          cell.num<-exp(log(2)*beta*j)
          if(p/2 + (1-p)/(2*cell.num) <min(assgin.vaf)){
            vec[j]<-p/2 + (1-p)/(2*cell.num)
          }else{
            vec[j] <- p/2 + (1-p)/(2*cell.num)  # Calculate the value for the current division number and store it in the vector
          }     
        }
        if((p/2) >=min(assgin.vaf)){
          vec[j]<-p/2
          #error=Inf
          mu.1=M.cul/cell.div
          
        }
        vec<-vec[which(vec != 0)]
        if((1-p)/(2*2) >=min(assgin.vaf)){
          #vec[j+1]<-(1-p)/(2*2^beta)
          error=Inf
        }else{
          #vec<-sort(vec)
          if(length(vec)<cell.div){
            error=Inf
          }else{
            vec<-data.frame(vec)
            vec$a<-depth*vec$vec
            vec$b<-depth-vec$a
            #result<-test_distribution(assgin.vaf,vec)
            probs <- sapply(1:nrow(vec), function(i) dbeta(assgin.vaf, vec$a[i], vec$b[i]))
            
            # Assign each number to the largest beta distribution
            max_rows <- apply(probs, 1, which.max)
            #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
            
            clustercount <- data.frame(table(max_rows))
            clustercount$max_rows<-as.numeric(as.character(clustercount$max_rows))
            missing_value <- setdiff(1:cell.div, clustercount$max_rows)
            
            # Add the missing value to the dataframe with a frequency of 0
            if(length(missing_value) > 0) {
              clustercount <- rbind(clustercount, data.frame(max_rows = missing_value, Freq = 0))
            }
            clustercount<-clustercount[order(clustercount$max_rows),]
            clustercount<-cbind(clustercount,vec)
            clustercount<-clustercount[order(clustercount$vec,decreasing = T),]
            clustercount<-clustercount[clustercount$Freq>0,]
            mu.1=M.cul/nrow(clustercount)
            clustercount$cumulative <- cumsum(clustercount$Freq)
            clustercount$x<-2*clustercount$vec-p
            clustercount<-clustercount[clustercount$x>0,]
            clustercount$x<-log(clustercount$x)
            fm2<-lm(cumulative~x,data=clustercount)
            #plot(clustercount$x,clustercount$cumulative)
            mu.2=summary(fm2)$coefficients[2,1]*(-1)*beta*log(2)
            mu.2=min(mu.2,cluster.count)
            mu.3=summary(fm2)$coefficients[1,1]*log(2)/log(1-p)
            r.n<-max(r.list)
            if(contp){
              ptest=ifelse(mu.1>magomu,1-pnorm(mu.1,mean=magomu,sd=sqrt(magomu)),pnorm(mu.1,mean=magomu,sd=sqrt(magomu)))
            }else{
              ptest<-0
            }
            error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)+1-ptest
            #error=((log(mu.1+1)-log(mu.2+1)))^2+(1-r.n)^2
            #error=(((mu.1+1)-(mu.2+1)))^2
          }
        }
        if(error<min.error){
          min.error=error
          mu=mu.1
        }
      }
      
    }
    }
    
   
    
    if(p>=0.4 & magocheck >1 & depth>=100){
    print("*******Using magos results*******")
     
     mu=mu_high
     cluster.count1=mago.result[1,"sum"]
     
   }

    result.data<-list(pf = p, mu = mu,initalcount=cluster.count,initalvaf=initial.vaf.cut$vaf,rf=r.n)
    
  return(result.data)
}
Opt.fun<-function(x,G1.sub,M.data) {
  p=x
  
 
  G1.sub<-data.frame(M.data)
  colnames(G1.sub)<-"vaf"
  
  G1.sub$count<-1
  
  G1.sub <- G1.sub %>%
    mutate(vaf = format(vaf, nsmall = 3)) %>%
    arrange(desc(vaf), desc(count)) %>%
    mutate(cumsum = cumsum(count)) %>%
    group_by(vaf) %>%
    slice_max(cumsum) %>%
    ungroup()
  
  #G1.sub$cumsum<-G1.sub$cumsum-init.mu.simu
  #G1.sub<-G1.sub[G1.sub$cumsum>0,]
  #print(vaf.list)
  G1.sub$vaf<-as.numeric(G1.sub$vaf)
  #print(p)
  G1.sub$lnf<-2*G1.sub$vaf-p
  G1<-G1.sub[G1.sub$lnf>0,]
  G1$x<-log(G1$lnf)
  vaf.first<-p/2 + (1-p)/(2*2^beta)
  compare.vec<-c(0.5,vaf.first)
  
  compare.vec<-sort(compare.vec)
  #print(vec)
  # Find the indices of the closest values in "vec" for each element in "assgin.vaf"
  compare.vec<-data.frame(compare.vec)
  compare.vec$a<-depth*compare.vec$compare.vec
  compare.vec$b<-depth-compare.vec$a
  probs <- sapply(1:nrow(compare.vec), function(i) dbeta(M.data, compare.vec$a[i], compare.vec$b[i]))
  
  # Assign each number to the largest beta distribution
  max_rows <- apply(probs, 1, which.max)
  #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
  
  cluster.count <- table(max_rows)[nrow(compare.vec)]
  
  M.cul<-length(M.data)-cluster.count
  initial.vaf.cut<-G1[G1$cumsum==cluster.count,]
  assgin.vaf<-M.data[M.data<initial.vaf.cut$vaf]
  assgin.vaf<-sort(assgin.vaf,decreasing = T)
  G1<-G1[G1$cumsum>cluster.count,]
  num_rows_to_remove <- round(nrow(G1) * 0.05)
  
  # Calculate the index of the rows to remove
  rows_to_remove <- c(1:num_rows_to_remove, (nrow(G1) - num_rows_to_remove + 1):nrow(G1))
  
  # Subset the data frame to remove the rows
  G1 <- G1[-rows_to_remove,]
  count=1
  tryCatch({
    mu.turn <-breakpoints(cumsum~x,data=G1,h=3/nrow(G1))
    fm1 <- lm(cumsum ~ breakfactor(mu.turn), data = G1)
    # Extract the coefficients from the model
    
    coeffs <- coef(fm1)
    slopes <- coeffs[-1]*(-1)*beta
    count=length(slopes)
  },error=function(e){
    count=1
  })
  
  
  max_count=max(count+3,10)
  min_count=1
  
  # Store the slopes for each breakpoint
  
  error.list<-c()
  Fail.cell=c()
  r.list<-c()
  for (cell.div.c in min_count:max_count){
    cell.div=cell.div.c
    #print(cell.div)
    mu.1=M.cul/cell.div
    vec <- numeric(cell.div)
    if(cell.div<=1){
      cell.num<-exp(log(2)*beta*1)
      vafmut<-p/2 + (1-p)/(2*cell.num)
      a<-depth*vafmut
      b<-depth-a
      
      ##count1
      pbeta_value <- pbeta(assgin.vaf, a,b)
      vecdata<-data.frame(vaf=assgin.vaf,p=pbeta_value)
      vecdata$p<-ifelse(vecdata$vaf>a/(a+b),1-vecdata$p,vecdata$p)
      mu.3=0
      r.n=0.1
      mu.2=min(nrow(vecdata[vecdata$p>0.1,]),cluster.count)
      error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)
      error.list<-c(error.list,error)
      #error=(((mu.1+1)-(mu.2+1)))^2
    }else{
      for (j in 1:cell.div) { 
        cell.num<-exp(log(2)*beta*j)
        if(p/2 + (1-p)/(2*cell.num) <min(assgin.vaf)){
          vec[j]<-p/2 + (1-p)/(2*cell.num)
        }else{
          vec[j] <- p/2 + (1-p)/(2*cell.num)  # Calculate the value for the current division number and store it in the vector
        }     
      }
      if((p/2) >=min(assgin.vaf)){
        vec[j]<-p/2
        #error=Inf
        #cell.div<-cell.div+1
        mu.1=M.cul/cell.div
        
      }
      vec<-vec[which(vec != 0)]
      if((1-p)/(2*2) >=min(assgin.vaf)){
        #vec[j+1]<-(1-p)/(2*2^beta)
        error=Inf
      }else{
        vec<-sort(vec)
        if(length(vec)<cell.div){
          error=Inf
        }else{
          vec<-data.frame(vec)
          vec$a<-depth*vec$vec
          vec$b<-depth-vec$a
          #result<-test_distribution(assgin.vaf,vec)
          probs <- sapply(1:nrow(vec), function(i) dbeta(assgin.vaf, vec$a[i], vec$b[i]))
          
          
          # Assign each number to the largest beta distribution
          max_rows <- apply(probs, 1, which.max)
          #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
          
          
          clustercount <- data.frame(table(max_rows))
          clustercount$max_rows<-as.numeric(as.character(clustercount$max_rows))
          missing_value <- setdiff(1:cell.div, clustercount$max_rows)
          
          # Add the missing value to the dataframe with a frequency of 0
          if(length(missing_value) > 0) {
            clustercount <- rbind(clustercount, data.frame(max_rows = missing_value, Freq = 0))
          }
          
          clustercount<-clustercount[order(clustercount$max_rows),]
          
          clustercount<-cbind(clustercount,vec)
          clustercount<-clustercount[order(clustercount$vec,decreasing = T),]
          #clustercount<-clustercount[-cell.div,]
          clustercount<-clustercount[clustercount$Freq>0,]
          mu.1=M.cul/nrow(clustercount)
          clustercount$cumulative <- cumsum(clustercount$Freq)
          clustercount$x<-2*clustercount$vec-p
          clustercount<-clustercount[clustercount$x>0,]
          clustercount$x<-log(clustercount$x)
          if(nrow(clustercount)==1){
        
            error=Inf
          }
          if(nrow(clustercount)==2){
            Fail.cell<-c(Fail.cell,cell.div.c)
            r.n=0
            fm2<-lm(cumulative~x,data=clustercount)
            mu.2=summary(fm2)$coefficients[2,1]*(-1)*beta*log(2)
            mu.3=summary(fm2)$coefficients[1,1]*log(2)/log(1-p)
            error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)
          }
          if(nrow(clustercount)>2){
            fm2<-lm(cumulative~x,data=clustercount)
            #plot(G1$x,G1$cumsum)
            mu.2=summary(fm2)$coefficients[2,1]*(-1)*beta*log(2)
            mu.3=summary(fm2)$coefficients[1,1]*log(2)/log(1-p)
            r.n<-summary(fm2)$r.squared
            r.list<-c(r.list,r.n)
            error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)
          }
         

        }
        
      }
      
    
    
    error.list<-c(error.list,error)
    
    
  }
  }
  
  if(length(r.list)<1){
    r.list<-c(r.list,0.8)
  }
  if(length(Fail.cell)>0){
    for (cell.div.c in Fail.cell){
      cell.div=cell.div.c
      mu.1=M.cul/cell.div
      vec <- numeric(cell.div)
      for (j in 1:cell.div) { 
        cell.num<-exp(log(2)*beta*j)
        if(p/2 + (1-p)/(2*cell.num) <min(assgin.vaf)){
          vec[j]<-p/2 + (1-p)/(2*cell.num)
        }else{
          vec[j] <- p/2 + (1-p)/(2*cell.num)  # Calculate the value for the current division number and store it in the vector
        }     
      }
      if((p/2) >=min(assgin.vaf)){
        vec[j]<-p/2
        #error=Inf
        mu.1=M.cul/cell.div
        
      }
      vec<-vec[which(vec != 0)]
      if((1-p)/(2*2) >=min(assgin.vaf)){
        #vec[j+1]<-(1-p)/(2*2^beta)
        error=Inf
      }else{
        #vec<-sort(vec)
        if(length(vec)<cell.div){
          error=Inf
        }else{
          vec<-data.frame(vec)
          vec$a<-depth*vec$vec
          vec$b<-depth-vec$a
          #result<-test_distribution(assgin.vaf,vec)
          probs <- sapply(1:nrow(vec), function(i) dbeta(assgin.vaf, vec$a[i], vec$b[i]))
          
          
          # Assign each number to the largest beta distribution
          max_rows <- apply(probs, 1, which.max)
          #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
          
          clustercount <- data.frame(table(max_rows))
          clustercount$max_rows<-as.numeric(as.character(clustercount$max_rows))
          missing_value <- setdiff(1:cell.div, clustercount$max_rows)
          
          # Add the missing value to the dataframe with a frequency of 0
          if(length(missing_value) > 0) {
            clustercount <- rbind(clustercount, data.frame(max_rows = missing_value, Freq = 0))
          }
          clustercount<-clustercount[order(clustercount$max_rows),]
          clustercount<-cbind(clustercount,vec)
          clustercount<-clustercount[order(clustercount$vec,decreasing = T),]
          clustercount<-clustercount[clustercount$Freq>0,]
          mu.1=M.cul/nrow(clustercount)
          clustercount$cumulative <- cumsum(clustercount$Freq)
          clustercount$x<-2*clustercount$vec-p
          clustercount<-clustercount[clustercount$x>0,]
          clustercount$x<-log(clustercount$x)
          fm2<-lm(cumulative~x,data=clustercount)
          #plot(clustercount$x,clustercount$cumulative)
          mu.2=summary(fm2)$coefficients[2,1]*(-1)*beta*log(2)
          mu.3=summary(fm2)$coefficients[1,1]*log(2)/log(1-p)
          r.n<-max(r.list)
          error=lossfun(mu.1,mu.2,mu.3,r.n,cell.div.c)
          #error=(((mu.1+1)-(mu.2+1)))^2
        }
      }
      error.list<-c(error.list,error)
    }
    
  }
  
  error=min(error.list)
  
  return(error)
}
# Function to generate bootstrap samples
generate_bootstrap_samples <- function(original_data, n_samples,num_decimal) {
  # Initialize an empty vector to store bootstrap samples
  bootstrap_samples <- numeric(n_samples)
  
  # Generate n_samples bootstrap samples
  for (i in 1:n_samples) {
    resample <- sample(original_data, size = length(original_data), replace = TRUE)
    bootstrap_samples[i] <- mean(resample)  # Or any other statistic you're interested in
  }
  
  return(bootstrap_samples)
}
# Function to compute the log likelihood for a mixture of p's
beta_likelihood<-function(vaf,shape1,shape2){
    ll <- sum(((shape1-1)*log(vaf)+(shape2-1)*log(1-vaf)-log(beta(shape1,shape2))))
    return(ll)
  
}
compute_AIC <- function(log_likelihood, num_params) {
  return(-2 * log_likelihood + 2 * num_params)
}

compute_BIC <- function(log_likelihood, num_params, sample_size) {
  return(-2 * log_likelihood + num_params * log(sample_size))
}

log_likelihood_bac_fit_compare <- function(data, p_vec,depth) {
  
  
  a<-depth*p_vec
  b<-depth-a
  probs <- sapply(1:length(a), function(i) dbeta(data, a[i], b[i]))
  df<-data.frame(
    prob=probs,
    vaf=data
  )
  
  if(length(p_vec)>1){
    # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
    df<-beta_reassign(df)
    
    df <- df %>%
      mutate(p = p_vec[cluster])
  }else{
    df$p<-p_vec
  }
  
 
  df <- df %>%
   mutate(log_likelihood = likbinom(x = round(vaf*depth),
                                  size = depth, prob = p))
  
  # If you want to get the total log likelihood, you would sum these up
  total_log_likelihood <- sum(df$log_likelihood)
  
  return(total_log_likelihood)
}

log_likelihood_mixture <- function(data, p_vec,depth) {
  
  
  a<-depth*p_vec
  b<-depth-a
  probs <- sapply(1:length(a), function(i) dbeta(data, a[i], b[i]))
  df<-data.frame(
    prob=probs,
    vaf=data
  )
 
  if(length(p_vec)>1){
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  df<-beta_reassign(df)
  
  df <- df %>%
    mutate(p = p_vec[cluster])
  }else{
    df$p<-p_vec
  }

  df <- df %>%
  mutate(log_likelihood = likbeta(x = vaf,
                                    shape1 = p*depth, shape2 = depth-p))
  #df <- df %>%
  # mutate(log_likelihood = likbinom(x = round(vaf*depth),
  #                                size = depth, prob = p))
  
  # If you want to get the total log likelihood, you would sum these up
  total_log_likelihood <- sum(df$log_likelihood)
  
  return(total_log_likelihood)
}




peak_test<-function(df,Min_Sample_size,vaf,k,depth,num_decimal){

  if(nrow(df)<k){
    return(c(NA,NA,NA,NA))
  }else{
    #Define start.div and end.div
    n_simulations <- getOption("teatime.peak_nsim", 100000)
    #n_simulations <- 1000
    simulated_vafs.right <- round(rbinom(n_simulations, depth, vaf) / depth,num_decimal)
    
    # For left_df
    if(k<Min_Sample_size){
      largest_right_values <- generate_bootstrap_samples(df[1:k,]$vaf,Min_Sample_size,num_decimal)
    }else{

      largest_right_values <- df[1:k,]$vaf
    }
    largest_right_values <- na.omit(largest_right_values)
    largest_right_values <- as.vector(largest_right_values)
    # Perform a two-sample Kolmogorov-Smirnov test to get a p-value
    #print(length(simulated_vafs.right))
   # print(length(largest_right_values))
    w.right.result <- wilcox.test(simulated_vafs.right, largest_right_values) 
    
    data1_long <- data.frame(value = simulated_vafs.right, group = "data1")
    data2_long <- data.frame(value = largest_right_values, group = "data2")
    combined_data <- rbind(data1_long, data2_long)
    
    cd.right=cliffDelta(value ~ group, data = combined_data)
    return(c(mean(simulated_vafs.right),mean(largest_right_values),abs(cd.right),w.right.result$p.value))
  }
  
  
  
}
simulate_from_estimation<-function(compare_set,cell.div,mu,p,depth,beta=0.8,num_decimal){
  vaf_list<-c(p/2)
  #print(cell.div)
  for(i in 1:cell.div){
  cell.num<-exp(log(2)*beta*i)

  vaf_list<-c(vaf_list,p/2 + (1-p)/(2*cell.num))

  }
  
  ##perform 10 times
  p_v_1_values = c()
  p_v_2_values = c()
  
  for (try in 1:getOption("teatime.sim_ntry", 50L)){
  # Simulate VAFs for each element in vaf_list and round them
  all_simulated_vafs_list <- sapply(vaf_list, function(vaf) {
    round(rbinom(round(mu), depth, vaf) / depth, num_decimal)
  })
  
  # Concatenate all vectors into a single vector
  all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
  p_v_1=ks.test(compare_set, all_simulated_vafs)$p.value
  p_v_2=wilcox.test(compare_set, all_simulated_vafs)$p.value
  p_v_1_values = c(p_v_1_values, p_v_1)
  p_v_2_values = c(p_v_2_values, p_v_2)
  }
  mean_p_v_1 = mean(p_v_1_values)
  mean_p_v_2 = mean(p_v_2_values)
  
  log_likelihood = log_likelihood_mixture(compare_set, vaf_list,depth) # example
  num_params = length(vaf_list)
  
  # Sample size
  sample_size = length(compare_set)
  AIC_value = compute_AIC(log_likelihood, num_params)
  BIC_value = compute_BIC(log_likelihood, num_params, sample_size)
  return(c(mean_p_v_1,mean_p_v_2,log_likelihood,AIC_value,BIC_value))
}

#If mu is reasonable large, try to find good-fit parameter
fit_domi_slope_check<-function(p,vaf_set,start.div,depth,beta,num_decimal,dynamic=F){

  collect.data=data.frame()
  idx<-1
  min_mu=3
  end.div<-round(length(vaf_set)/min_mu)
  cell.div<-start.div:end.div
  
  #vaf_set=c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  mu_list<-round(length(vaf_set)/cell.div)
  # Create a data frame to keep track of mu_list and corresponding cell.div
  mu_data <- data.frame(mu_list, cell.div)
  
  # Remove duplicates in mu_list and randomly select one corresponding cell.div
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
  div.list<-unique_df$cell.div
  uq=0.75
  lq=0.25
  if(length(div.list)>=4){
  div.list<-div.list[div.list<=quantile(div.list, uq) & div.list >= quantile(div.list, lq)]
  }
  #dynamic=T
  #print(unique_df$cell.div)
  main.vaf<-vaf_set
  for (cell.div in .cap_candidates(div.list)){
    #print(cell.div)
    if(dynamic){
      
      i_values <- 1:cell.div
      
      result_vector <- sapply(i_values, function(i) {
        # Replace this with your actual calculation
        calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
        return(calculated_value)
      })
      a<-depth*c(0.5,p/2,result_vector)
      b<-depth-a
      vaf.require.assign<-main.vaf
      probs <- sapply(1:length(a), function(i) dbeta(vaf.require.assign, a[i], b[i]))
      
    
      df<-data.frame(
        prob=probs,
        vaf=vaf.require.assign
      )
      
      df<-beta_reassign(df)
      
      vaf_set<-df[df$cluster>2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
      
      #vaf_set=c(first_vaf_list,second.vaf[second.vaf>left_before_t1_vaf],last.div.set)
    }else{
      vaf_set=vaf_set[vaf_set>p/2]
    }
    mu_est=length(vaf_set)/cell.div
    #print(cell.div)
    #print(mu_est)
    if(mu_est>3){
      i_values <- 1:cell.div
      
      result_vector <- sapply(i_values, function(i) {
        # Replace this with your actual calculation
        calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
        return(calculated_value)
      })
      
      a<-depth*result_vector
      b<-depth-a
      probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
      df<-data.frame(
        prob=probs,
        vaf=vaf_set
      )
      
      
      df<-beta_reassign(df)
      df<-df[!duplicated(df),]
      
      mulist<-c()
      plist<-c()
      for(try in 1:3){
        mu_from_real<-normal_slope_get(df,p,beta,result_vector,num_decimal)
        #print(mu_from_real)
        #print(mu_est)
        mu_from_simu<-normal_slope_simu(cell.div,mu_est,p,depth,beta,num_decimal)
        #print(mu_from_simu)
        mean_list <- mean(mu_from_simu)
        sd_list <- max(sd(mu_from_simu),1.5)
        
        z_score <- (mu_from_real - mean_list) / sd_list
        mulist<-c(mulist,mu_from_real)
        plist<-c(plist,abs(z_score))
      }
      closest_index <- which.min(plist)
      mu_from_real <- mulist[closest_index]
      #closest_mu1 <- mulist1[which.min(abs(mulist1 - median(mulist1)))]
      z_score<-plist[closest_index]
      
      collect.data[idx,"mu"]=mu_est
      collect.data[idx,"mu_real"]=mu_from_real
      collect.data[idx,"mu_simu"]=mean_list
      collect.data[idx,"cell.div"]=cell.div
      collect.data[idx,"z_score"]=z_score
      idx<-idx+1
    }
  }
  
  
  
  collect.data<-na.omit(collect.data)
  collect.data$p<-p
  
  
  return(collect.data)
  
}


mu_find_large<-function(vaf_set,celldivlist,p,depth, beta=0.8,num_decimal,p_thre=0.05){
  Min_Sample_size=6
  min_mu=3
  

  collect.data<-data.frame()
  idx<-1
  for (cell.div in .cap_candidates(celldivlist)){
    
    collect.data[idx,"cell.div"]=cell.div
    i_values <- 1:cell.div
    result_vector <- sapply(i_values, function(i) {
      # Replace this with your actual calculation
      calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
      return(calculated_value)
    })
    result_vector<-c(0.5,p/2,result_vector)
    
    a<-depth*result_vector
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=vaf_set
    )
    
    df<-beta_reassign(df)
    
    mu1<-max(nrow(df[df$cluster>=3,])/cell.div,1)
    left_df <- df[df$cluster==2+cell.div,]
    
    right_df <-df[df$cluster==3,]
    if(nrow(right_df)<nrow(left_df)){
      right_df<-df[df$cluster==2,]
    }
    if(nrow(right_df)>mu1 & mu1>3){
      left_most_vaf=p/2 + (1-p)/(2 * exp(log(2) * beta * cell.div))
      right_most.vaf=p/2 + (1-p)/(2 * exp(log(2) * beta * 1))
      
      
      #df<-df[df$cluster==3,]
      
      left_df$abs_diff <- abs(left_df$vaf - left_most_vaf)
      
      # Sort the data frame by absolute differences
      #left_df <- left_df[order(left_df$abs_diff), ]

  
      left_df$prob<-pbeta(left_df$vaf,depth*left_most_vaf,depth-depth*left_most_vaf)
      left_df$score<-left_df$prob/max(left_df$prob)-left_df$abs_diff/max(left_df$abs_diff)
      left_df <- left_df[order(-left_df$score), ]
      
      right_df$abs_diff <- abs(right_df$vaf - right_most.vaf)
      
      #right_df <- right_df[order(right_df$abs_diff), ]
      right_df$prob<-pbeta(right_df$vaf,depth*right_most.vaf,depth-depth*right_most.vaf)
      
      right_df$score<-right_df$prob/max(right_df$prob)-right_df$abs_diff/max(right_df$abs_diff)
      right_df <- right_df[order(-right_df$score), ]
      
    
    k<-round(mu1)
    # Monte Carlo simulation to generate VAFs using the binomial distribution
    left_result<-peak_test(left_df,Min_Sample_size,left_most_vaf,k,depth,num_decimal)
    right_result<-peak_test(right_df,Min_Sample_size,right_most.vaf,k,depth,num_decimal)
    # Print p-value
    collect.data[idx,"simu.left.mean"]=left_result[1]
    collect.data[idx,"simu.right.mean"]=right_result[1]
    collect.data[idx,"left.mean"]=left_result[2]
    collect.data[idx,"right.mean"]=right_result[2]
    collect.data[idx,"left.cd"]=left_result[3]
    collect.data[idx,"right.cd"]=right_result[3]
    collect.data[idx,"left.p"]=left_result[4]
    collect.data[idx,"right.p"]=right_result[4]
    
    collect.data[idx,"mu_est"]=mu1
    idx<-idx+1
    }
  }
  
  
  Temp_keep<-collect.data
  collect.data<-Temp_keep[Temp_keep$left.p>p_thre | Temp_keep$right.p>p_thre,]
  collect.data<-collect.data[complete.cases(collect.data[ , 3]),]
  collect.data$left.cd<-ifelse(is.na(collect.data$left.cd),10,collect.data$left.cd)
  #collect.data<-collect.data[collect.data$mu_est>3,]
  if(nrow(collect.data)==0){
    collect.data=Temp_keep
    collect.data<-collect.data[complete.cases(collect.data[ , 3]),]
    collect.data$left.cd<-ifelse(is.na(collect.data$left.cd),10,collect.data$left.cd)
  }
    
  # Convert specific columns to ranks and sum those ranks
  if (all(collect.data$left.cd == 10)){
    collect.data$right_rank<-rank(collect.data$right.cd)
    
    if(nrow(collect.data)>2){
    top_10_percent=floor(0.5 * nrow(collect.data))
    
    sorted_df <- collect.data[order(collect.data$right_rank), ]
    top_right <- head(sorted_df, top_10_percent)
    
    
    top_pick<-top_right}else{
      top_pick<-collect.data
    }
   
  }else{
    collect.data$left_rank<-rank(collect.data$left.cd)
    collect.data$right_rank<-rank(collect.data$right.cd)
    if(nrow(collect.data)>2){
    top_10_percent=floor(0.5 * nrow(collect.data))
    
    sorted_df <- collect.data[order(collect.data$left_rank), ]
    top_left <- head(sorted_df, top_10_percent)
    
    sorted_df <- collect.data[order(collect.data$right_rank), ]
    top_right <- head(sorted_df, top_10_percent)
    
    
    top_pick<-rbind(top_left,top_right)
    top_pick<-top_pick[!duplicated(top_pick),]}else{
      top_pick<-collect.data
    }
  }
  
  
  
  
  #compare_set,cell.div,mu,p,depth,beta=0.8,num_decimal
  top_pick_p <- apply(top_pick, 1, function(row) {
    cell_div_value <- row['cell.div']
    mu_value <- row['mu_est']
    i_values <- 1:cell_div_value
    result_vector <- sapply(i_values, function(i) {
      # Replace this with your actual calculation
      calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
      return(calculated_value)
    })
    result_vector<-c(0.5,p/2,result_vector)
    
    a<-depth*result_vector
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=vaf_set
    )
    
    df<-beta_reassign(df)
    data<-df[df$cluster>=2,]$vaf
    simulate_from_estimation(data, cell_div_value, mu_value, p, depth, beta,num_decimal)
  })
  
  top_pick_p <- as.data.frame(t(top_pick_p))
  top_pick_p$V1<-ifelse(top_pick_p$V1<0.05,0,top_pick_p$V1)
  top_pick_p$V2<-ifelse(top_pick_p$V2<0.05,0,top_pick_p$V2)
  top_pick$r1<-top_pick_p$V1
  top_pick$r2<-top_pick_p$V2
  
  top_pick$loglike<-top_pick_p$V3
  top_pick$aic<-top_pick_p$V4
  top_pick$bic<-top_pick_p$V5

  
  top_pick$top_pick_score<-top_pick$r2+1/top_pick$right_rank
  #top_pick<- top_pick[order(top_pick$top_pick_score), ]
  if(nrow(top_pick)>1){
    top_pick<- top_pick[top_pick$top_pick_score>median(top_pick$top_pick_score), ]
  }
  top_pick <- top_pick[order(top_pick$right_rank), ]
  pick_mu_cell.div<-top_pick[1:min(10,nrow(top_pick)),]
  
  #print(top_pick[abs(top_pick$mu_est-16)==min(abs(top_pick$mu_est-16)),])
  #print(pick_mu_cell.div)
  return(pick_mu_cell.div)

}

mu_find_large_bac<-function(vaf_set,start.div,p,depth, beta=0.8,num_decimal,p_thre=0.05){
  Min_Sample_size=6
  min_mu=3
  collect.data<-data.frame()
  right_most.vaf <- p/2 + (1-p)/(2 * exp(log(2) * beta ))
  
  abs_diff <- abs(vaf_set - right_most.vaf)
  
  
  right_df<-data.frame(
    vaf=vaf_set,
    abs_diff=abs_diff
  )
  
  #right_df <- right_df[order(right_df$abs_diff), ]
  right_df$prob<-pbeta(vaf_set,depth*right_most.vaf,depth-depth*right_most.vaf)
  
  right_df$score<-right_df$prob/max(right_df$prob)-right_df$abs_diff/max(right_df$abs_diff)
  right_df <- right_df[order(-right_df$score), ]
  idx<-1
  for (mu in min_mu:length(vaf_set)){
    
    collect.data[idx,"mu_est"]=mu
    

    k<-mu
    # Monte Carlo simulation to generate VAFs using the binomial distribution
    #print(k)
    right_result<-peak_test(right_df,Min_Sample_size,right_most.vaf,k,depth,num_decimal)
      # Print p-value
     
      collect.data[idx,"simu.right.mean"]=right_result[1]
      collect.data[idx,"right.mean"]=right_result[2]
      collect.data[idx,"right.cd"]=right_result[3]
      collect.data[idx,"right.p"]=right_result[4]
      #collect.data[idx,"mu_est"]=mu1
      idx<-idx+1
    
  }
  
  
  Temp_keep<-collect.data
  collect.data<-Temp_keep[Temp_keep$right.p>p_thre,]
  if(nrow(collect.data)==0){
    collect.data=Temp_keep
   
  }
  
  # Convert specific columns to ranks and sum those ranks

  collect.data$right_rank<-rank(collect.data$right.cd)
  #print(collect.data) 
  if(nrow(collect.data)>4){
  top_10_percent=floor(0.25 * nrow(collect.data))
    
  sorted_df <- collect.data[order(collect.data$right_rank), ]
  top_right <- head(sorted_df, top_10_percent)
  }else if(nrow(collect.data)>2){
    top_10_percent=floor(0.5 * nrow(collect.data))
    
    sorted_df <- collect.data[order(collect.data$right_rank), ]
    top_right <- head(sorted_df, top_10_percent)
  }else{
    top_right<-collect.data
  } 
    
  top_pick<-top_right
    
  
  
  #compare_set,cell.div,mu,p,depth,beta=0.8,num_decimal
  top_pick_p <- apply(top_pick, 1, function(row) {
    mu_value <- row['mu_est']
    largest_right_values <- right_df[1:mu_value,]$vaf
    ##perform 10 times
    p_v_1_values = c()
    p_v_2_values = c()
    
    for (try in 1:getOption("teatime.sim_ntry", 50L)){
      # Simulate VAFs for each element in vaf_list and round them
      all_simulated_vafs_list <- round(rbinom(round(mu_value), depth, right_most.vaf) / depth, num_decimal)
      
      # Concatenate all vectors into a single vector
      all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
      p_v_1=ks.test(largest_right_values, all_simulated_vafs)$p.value
      p_v_2=wilcox.test(largest_right_values, all_simulated_vafs)$p.value
      p_v_1_values = c(p_v_1_values, p_v_1)
      p_v_2_values = c(p_v_2_values, p_v_2)
    }
    mean_p_v_1 = mean(p_v_1_values)
    mean_p_v_2 = mean(p_v_2_values)
    
    log_likelihood = log_likelihood_mixture(largest_right_values, right_most.vaf,depth) # example
    #log_likelihood=log_likelihood/length(largest_right_values)
    num_params = 1
    
    # Sample size
    sample_size = length(largest_right_values)
    AIC_value = compute_AIC(log_likelihood, num_params)
    BIC_value = compute_BIC(log_likelihood, num_params, sample_size)
    return(c(mean_p_v_1,mean_p_v_2,log_likelihood,AIC_value,BIC_value))
  })
  
  top_pick_p <- as.data.frame(t(top_pick_p))
  top_pick_p$V1<-ifelse(top_pick_p$V1<0.05,0,top_pick_p$V1)
  top_pick_p$V2<-ifelse(top_pick_p$V2<0.05,0,top_pick_p$V2)
  top_pick$r1<-top_pick_p$V1
  top_pick$r2<-top_pick_p$V2
  
  top_pick$loglike<-top_pick_p$V3
  top_pick$aic<-top_pick_p$V4
  top_pick$bic<-top_pick_p$V5
  
  top_pick$top_pick_score<-top_pick$r2+1/top_pick$right_rank
  top_pick<- top_pick[order(-top_pick$top_pick_score), ]
  
  pick_mu_cell.div<-top_pick[1:min(10,nrow(top_pick)),]
  
  return(pick_mu_cell.div)
  
}
mu_find_small<-function(vaf_set,p,depth, beta,num_decimal,test_mu=3){
  ##Now if mu is very small, then it is a line
  mu_collect=c()
  
  # Calculate slopes
  slope1 <- Get_slope(vaf_set, p)
  slope2 <- Get_slope(vaf_set, p + 0.01)
  slope3 <- Get_slope(vaf_set, p - 0.01)
  
  # Add small random noise in case the slopes are the same
  if (slope1 == slope2 || slope1 == slope3 || slope2 == slope3) {
    noise <- runif(3, min = -1e-5, max = 1e-5) # Change the range of noise as needed
    slope1 <- slope1 + noise[1]
    slope2 <- slope2 + noise[2]
    slope3 <- slope3 + noise[3]
  }
  
  # Add to mu_collect
  mu_collect <- c(mu_collect, slope1, slope2, slope3)
  
  # Means to test against
  means_to_test <- 1:test_mu
  
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
    wilcox.test(diff2, diff1)$p.value
  }
  )
  # Combine results into a data frame
  result_df <- data.frame(mu_est = means_to_test, p_value = test_results, wx_p = test_results2)
  
  ###mu small selection
  return(result_df)
  
}
fit.over.check<-function(main.vaf,depth,beta,p_thre=0.05){
  num_decimal <- nchar(as.character(depth))
  right_most.vaf<-0
  a<-depth*0.5
  b<-depth-a
  
  probs <- pbeta(main.vaf, a, b)
  combined_df <- data.frame(main.vaf = main.vaf, probs = probs)
  #print(combined_df)
 # border_vaf<-max(combined_df[combined_df$probs<=0.05,"main.vaf"]) ##first range vaf
    if(nrow(combined_df[combined_df$probs<=0.05,]) >0){
  border_vaf<-max(combined_df[combined_df$probs<=0.05,"main.vaf"]) ##first range vaf
  }else{
    border_vaf<-min(combined_df$main.vaf) ##first range vaf
  }
  ##First find left most vaf
  tryCatch({
      # First attempt
      m <- automixfit(main.vaf, type = "beta", Nc = 2:10, thresh = 0, Ninit = min(50, round(length(main.vaf) / 5)), k = 6, Niter.max = 10000)
    }, error = function(e) {
      tryCatch({
        # Second attempt if the first fails
        m <- automixfit(main.vaf, type = "beta", Nc = 2:10, thresh = 0, Ninit = 3, k = 6, Niter.max = 10000)
      }, error = function(e) {
        # Third attempt if the second also fails
        m <- automixfit(main.vaf, type = "beta", Nc = 2:10, thresh = 0, Niter.max = 10000)
      })
    })
  #m<-automixfit(main.vaf, type = "beta",Nc =2:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf)/5)), Niter.max=10000)
  #print(m)
  #print(main.vaf)
  a=m["a",]
  b=m["b",]
  mean.a.b=a/(a+b)
  approx_vaf=min(mean.a.b)
  approx_vaf_index=which.min(mean.a.b)
  
  left_most_vaf<-round(approx_vaf,num_decimal)

  
  ##Now find right vaf
  inita=c(0.5*depth,a[approx_vaf_index])
  initb=c(0.5*depth,b[approx_vaf_index])
  probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
  df<-data.frame(
    prob=probs,
    vaf=main.vaf
  )
  
        # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  df<-beta_reassign(df)
  #a[1]=
  main.vaf_update<-df[df$cluster<2,]$vaf
  
  right_most.vaf<-approx_vaf
  
  right_save<-right_most.vaf
  ra_save<-a[approx_vaf_index]
  rb_save<-b[approx_vaf_index]
  iter=1
  while (TRUE & iter<=100){
      tryCatch({
      # First attempt
      m <- automixfit(main.vaf_update, type = "beta", Nc = 1:10, thresh = 0, Ninit = min(50, round(length(main.vaf_update) / 5)), k = 6, Niter.max = 10000)
    }, error = function(e) {
      tryCatch({
        # Second attempt if the first fails
        m <- automixfit(main.vaf_update, type = "beta", Nc = 1:10, thresh = 0, Ninit = 3, k = 6, Niter.max = 10000)
      }, error = function(e) {
        # Third attempt if the second also fails
        m <- automixfit(main.vaf_update, type = "beta", Nc = 1:10, thresh = 0, Niter.max = 10000)
      })
    })
  #m<-automixfit(main.vaf_update, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf_update)/5)), Niter.max=10000)
  updatea=m["a",]
  updateb=m["b",]
  mean.a.b=updatea/(updatea+updateb)
  approx_vaf=min(mean.a.b)
  #print(approx_vaf)
  if(approx_vaf>=0.5 | length(mean.a.b)==1){
   break
  }
  approx_vaf_index=which.min(mean.a.b)
  
  currenta=c(0.5*depth,updatea[approx_vaf_index])
  currentb=c(0.5*depth,updateb[approx_vaf_index])
  ra<-currenta[2]
  rb<- currentb[2]
  probs <- sapply(1:length(currenta), function(i) dbeta(main.vaf_update, currenta[i], currentb[i]))

  df<-data.frame(
    prob=probs,
    vaf=main.vaf_update
  )
  
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability

  df<-beta_reassign(df)
  main.vaf_update<-df[df$cluster<2,]$vaf
  #print(diff)
  if(approx_vaf<0.5){
    right_most.vaf<-approx_vaf
    if(approx_vaf>right_save){
      right_save<-approx_vaf
    }
  }
  iter<-iter+1
  }

  ##consider theratical first peak
  right.th.vaf<-left_most_vaf+(0.5-left_most_vaf)/exp(log(2)*beta)
  right_most.vaf=max(right_most.vaf,border_vaf,right.th.vaf,right_save)
  right_most.vaf=round(right_most.vaf,num_decimal)
  
  p<-(right_most.vaf*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)
  
  
 
  
  
  
  #Clonal vaf to analysis 
  vaf_set<-main.vaf
  
  
  
  
  #If mu is reasonable large, try to find good-fit parameter
  if(right_most.vaf == left_most_vaf){
    start.div<-1
  }else{
    start.div<-2
  }
  #print(vaf_set)
  first_check<-fit_domi_slope_check(p,vaf_set,start.div,depth,beta,num_decimal,dynamic=F)
  if(nrow(first_check[first_check$mu_real>3,])>0){   # keep cells above the minimum mu
  first_check<-first_check[first_check$mu_real>3,]
  }else{
  first_check<-first_check[first_check$mu_real>1,]
  }
  
  
  if(nrow(first_check[abs(first_check$z_score)<=1.96,])>0){
    ##if over 20 ( very hard to distinguish)
    collect.data<-first_check[abs(first_check$z_score)<=1.96,]
  }else{
    collect.data<-first_check
  }
    celldivlist<-collect.data$cell.div
    # Sort the data frame by absolute differences
    collect.data2<- mu_find_large(vaf_set,celldivlist,p,depth, beta=beta,num_decimal,p_thre)
    #print('Second done')
    if(nrow(collect.data2)==1){
      range1=collect.data2$mu_est-1
      range2=collect.data2$mu_est+1
    }else{
      range1=min(collect.data2$mu_est)
      range2=max(collect.data2$mu_est)
    }
    if(range2>median(collect.data$mu)){
      
    over12<-collect.data[collect.data$mu_real>=range1 & collect.data$mu_real<=range2,]
    }else{
      over12<-collect.data[collect.data$mu_real>=range1,]
    }
    overlap.pick<-NULL
    if(nrow(over12)>0){
      overlap.pick<-over12
    }
    ##check three overlap
    if(is.null(overlap.pick)){
      overlap.pick<-collect.data
    }
    
    overlap.pick<-na.omit(overlap.pick)
  
  #print(overlap.pick)
  #print('Overcheck done')
  #when mu <=3
  

  
  
  ##Now if mu is very small, then it is a line
  i_values <- 1:20
  result_vector <- sapply(i_values, function(i) {
    # Replace this with your actual calculation
    calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
    return(calculated_value)
  })
  result_vector<-c(0.5,p/2,result_vector)
  
  a<-depth*result_vector
  b<-depth-a
  probs <- sapply(1:length(a), function(i) dbeta(main.vaf, a[i], b[i]))
  df<-data.frame(
    prob=probs,
    vaf=main.vaf
  )
  
  df<-beta_reassign(df)
  
  result_df<-mu_find_small(df[df$cluster>=2,]$vaf,p,depth, beta,num_decimal)
  ###mu small selection
  
  mu_small_selection=result_df[result_df$p_value>0.05 & result_df$wx_p>0.05 ,]


  vaf_set<-main.vaf
  

  #print(overlap.pick)
  top_pick_p <- apply(overlap.pick, 1, function(row) {
    cell_div_value <- row['cell.div']
    mu_value <- row['mu_real']
    i_values <- 1:cell_div_value
    result_vector <- sapply(i_values, function(i) {
      # Replace this with your actual calculation
      calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
      return(calculated_value)
    })
    result_vector<-c(0.5,p/2,result_vector)
    
    a<-depth*result_vector
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=vaf_set
    )
    
    df<-beta_reassign(df)
    data<-df[df$cluster>=2,]$vaf
    simulate_from_estimation(data, cell_div_value, mu_value, p, depth, beta,num_decimal)
  })
  top_pick <- as.data.frame(t(top_pick_p))
  #overlap.pick$loglike<-
  overlap.pick$loglike<-top_pick$V3
  overlap.pick$aic<-top_pick$V4
  overlap.pick$bic<-top_pick$V5
  pick.cell.div=overlap.pick$cell.div
  pick.mu=overlap.pick$mu_real
  z_score=overlap.pick$z_score
  
  pick.log<-overlap.pick$loglike
  pick.bic<-overlap.pick$bic
  pick.aic<-overlap.pick$aic
  
  if(nrow(mu_small_selection)>0){
    
    small_mu_pick <- apply(mu_small_selection, 1, function(row) {
      mu_value <- row['mu_est']
      cell_div_value <- round(length(df[df$cluster>=2,]$vaf)/mu_value)
      
      simulate_from_estimation(df[df$cluster>=2,]$vaf, cell_div_value, mu_value, p, depth, beta,num_decimal)
    })
    small_mu_pick <- as.data.frame(t(small_mu_pick))
    small_mu_pick$V1<-ifelse(small_mu_pick$V1<0.05,0,top_pick_p$V1)
    small_mu_pick$V2<-ifelse(small_mu_pick$V2<0.05,0,top_pick_p$V2)
   
    mu_small_selection$r1<-small_mu_pick$V1
    mu_small_selection$r2<-small_mu_pick$V2
    
    mu_small_selection$loglike<-small_mu_pick$V3
    mu_small_selection$aic<-small_mu_pick$V4
    mu_small_selection$bic<-small_mu_pick$V5

    mu_small_selection<-mu_small_selection[mu_small_selection$bic==min(mu_small_selection$bic),]
    
    if(mu_small_selection$bic>min(pick_mu_cell.div$bic)*2){
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
  #all.p.data$clear<-clear
  all.p.data<-all.p.data[order(-all.p.data$score),]
  lower_bound<-max(collect.data2$mu_est)
  lower_bound2<-min(collect.data2$mu_est)
  all.p.data$lowerbound1<-lower_bound
  all.p.data$lowerbound2<-lower_bound2
  #print(fit_pick)
  return(all.p.data)
  
}

bac.over.check<-function(main.vaf,second.vaf,depth,beta,p_thre=0.05){
  num_decimal <- nchar(as.character(depth))
  ##First find left most vaf
      tryCatch({
      # First attempt
      m <- automixfit(second.vaf, type = "beta", Nc = 1:10, thresh = 0, Ninit = min(50, round(length(second.vaf) / 5)), k = 6, Niter.max = 10000)
    }, error = function(e) {
      tryCatch({
        # Second attempt if the first fails
        m <- automixfit(second.vaf, type = "beta", Nc = 1:10, thresh = 0, Ninit = 3, k = 6, Niter.max = 10000)
      }, error = function(e) {
        # Third attempt if the second also fails
        m <- automixfit(second.vaf, type = "beta", Nc = 1:10, thresh = 0, Niter.max = 10000)
      })
    })
  #m<-automixfit(second.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(second.vaf)/3)), Niter.max=10000)
  a=m["a",]
  b=m["b",]
  mean.a.b=a/(a+b)
  left_most_vaf=min(mean.a.b)
  right_most.vaf=max(mean.a.b)
  right_vaf_index=which.max(mean.a.b)
  keep.right.vaf<-right_most.vaf
  probs <- sapply(1:length(a), function(i) dbeta(second.vaf, a[i], b[i]))
  df<-data.frame(
    prob=probs,
    vaf=second.vaf
  )
  
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  if(length(mean.a.b)>1){
  df<-beta_reassign(df)
  second.update.vaf<-df[df$cluster==right_vaf_index,]$vaf
  }else{
    second.update.vaf<-df$vaf
  }
  #a[1]=
  iter<-1
  while (iter<=100 & length(mean.a.b)>1){
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
    if(length(mean.a.b)<2){
      break
    }
    right_vaf_index=which.max(mean.a.b)
    probs <- sapply(1:length(a), function(i) dbeta(second.update.vaf, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=second.update.vaf
    )
    
    # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
    df<-beta_reassign(df)
    #a[1]=
    second.update.vaf<-df[df$cluster==right_vaf_index,]$vaf
    if (length(second.update.vaf) == 0 || (length(unique(second.update.vaf)) == 1)) {
      break
    }
    }else{
      break
    }
  iter=iter+1
  }
  right_most.vaf<-max(right_most.vaf,keep.right.vaf)
  #print(right_most.vaf)
  p=(right_most.vaf*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)
  
  if(p<0){
    p=1-2*exp(log(2)*beta)*right_most.vaf
    if(p<0){
      p<-0.01
    }
    
  }
  righta=c(0.5*depth,right_most.vaf*depth)
  rightb=depth-righta
  #mu_est=Get_slope(second.vaf,p)
  probs <- sapply(1:length(righta), function(i) dbeta(main.vaf, righta[i], rightb[i]))
  df<-data.frame(
    prob=probs,
    vaf=main.vaf
  )
  
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  df<-beta_reassign(df)
  
  vaf_set<-c(df[df$cluster==2,]$vaf,second.update.vaf)
  
  #If mu is reasonable large, try to find good-fit parameter
  if(left_most_vaf == right_most.vaf){
    start.div<-1
  }else{
    start.div<-2
  }
  #print(start.div)
  pick_mu_cell.div=mu_find_large_bac(vaf_set,start.div,p,depth, beta=beta,num_decimal,p_thre)
  
  #print('curr')
  vaf_set=c(df[df$cluster==2,]$vaf,second.vaf)
  ##Now if mu is very small, then it is a line
  result_df<-mu_find_small(vaf_set,p,depth, beta,num_decimal)
  ###mu small selection
  
  mu_small_selection=result_df[result_df$p_value>0.05 & result_df$wx_p>0.05 ,]
  
  

  pick.mu=pick_mu_cell.div$mu_est
  pick.bic=pick_mu_cell.div$bic
  pick.aic=pick_mu_cell.div$aic
  pick.log=pick_mu_cell.div$loglike
  
  
  if(nrow(mu_small_selection)>0){
    
    small_mu_pick <- apply(mu_small_selection, 1, function(row) {
      mu_value <- row['mu_est']
      cell_div_value <- round(length(vaf_set)/mu_value)
      simulate_from_estimation(vaf_set, cell_div_value, mu_value, p, depth, beta,num_decimal)
    })
    small_mu_pick <- as.data.frame(t(small_mu_pick))
    small_mu_pick$V1<-ifelse(small_mu_pick$V1<0.05,0,top_pick_p$V1)
    small_mu_pick$V2<-ifelse(small_mu_pick$V2<0.05,0,top_pick_p$V2)
    
    mu_small_selection$r1<-small_mu_pick$V1
    mu_small_selection$r2<-small_mu_pick$V2
    
    mu_small_selection$loglike<-small_mu_pick$V3
    mu_small_selection$aic<-small_mu_pick$V4
    mu_small_selection$bic<-small_mu_pick$V5
    
    
    mu_small_selection<-mu_small_selection[mu_small_selection$bic==min(mu_small_selection$bic),]
    
    if(mu_small_selection$bic>min(pick_mu_cell.div$bic)*2){
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

Find_p_process<-function(second.vaf,mean.a.b,upper_clonal_vaf,clonal.vaf,clonal.vaf.left,beta){
 
    tryCatch({
      # First attempt
      m <- automixfit(second.vaf, type = "beta", Nc = 1:10, thresh = 0, Ninit = min(50, round(length(second.vaf) / 3)), k = 6, Niter.max = 10000)
    }, error = function(e) {
      tryCatch({
        # Second attempt if the first fails
        m <- automixfit(second.vaf, type = "beta", Nc = 1:10, thresh = 0, Ninit = 3, k = 6, Niter.max = 10000)
      }, error = function(e) {
        # Third attempt if the second also fails
        m <- automixfit(second.vaf, type = "beta", Nc = 1:10, thresh = 0, Niter.max = 10000)
      })
    })
  #m<-automixfit(second.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(second.vaf)/3)), Niter.max=10000)
  a1=m["a",]
  b1=m["b",]
  mean.a.b1=a1/(a1+b1)
  second_clone_left_vaf=min(mean.a.b1)
  clear=F
  if(length(mean.a.b)>1){
    left_limit=(min(second.vaf)+second_clone_left_vaf)/2
    right_limit=max(((upper_clonal_vaf+max(clonal.vaf.left))/2*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)/2,(max(second.vaf)+second_clone_left_vaf)/2)
    
  }else{
    #Try from min to a left_limit to right_limit
    left_limit=(min(second.vaf)+second_clone_left_vaf)/2
    right_limit=((max(second.vaf)+upper_clonal_vaf)/2*2*exp(log(2)*beta)-1)/(exp(log(2)*beta)-1)/2
    #print(left_limit)
    #print(right_limit)
    if(right_limit>0.45 | right_limit<left_limit){
      right_limit=(max(second.vaf)+second_clone_left_vaf)/2
    }
    clear=T
  }
  
  if(length(second.vaf[second.vaf>right_limit])<=6){
    sorted_vaf <- sort(second.vaf, decreasing = TRUE)
    right_limit <- sorted_vaf[max(round(0.25*length(sorted_vaf)),6)]
  }
  normal.return<-NULL
  z_score_min<-19.6
  left <- left_limit
  middle <- (left_limit + right_limit) / 2
  right <- right_limit
  left1 <- left_limit + 0.25 * (right_limit - left_limit)
  right1<-middle+0.25 * (right_limit - left_limit)
  # Quantile points sequence
  #my_seq <- c(left, left1, middle, right1,right)
  my_seq <- c(left,middle,right)
  my_seq<-my_seq[my_seq>0.5/(1+exp(log(2)*beta))]
  return(list(my_seq = my_seq, clear = clear))
}

Update_P_process<-function(second.vaf,give.vaf,upper_clonal_vaf,clonal.vaf,clonal.vaf.left){
  
  iter<-1
  while(iter<=100){
  
  suppose_right_vaf=p/2+(1-p)/(2*exp(log(2)*beta))
  if(suppose_right_vaf > upper_clonal_vaf){
    #a<-c(depth*suppose_right_vaf)
    #b<-c(depth-a)
    #probs <-pbeta(clonal.vaf.left, a, b)
    #combined_df <- data.frame(vaf = clonal.vaf.left, probs = probs)
    #freq_df <- combined_df %>% 
    #  group_by(vaf, probs) %>% 
    #  summarise(count = n()) %>% 
    #  mutate(freq_to_assign = round(count * probs)) %>% 
    #  ungroup()
    
    # Expand the dataset based on calculated frequencies
    # first_vaf_list <- freq_df %>% 
    #  uncount(freq_to_assign, .remove = FALSE) %>% 
    #  pull(vaf)
    
    a<-c(depth*0.5,depth*suppose_right_vaf)
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(clonal.vaf.left, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=clonal.vaf.left
    )
    
    df<-beta_reassign(df)
    first_vaf_list<-df[df$cluster==2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
    
  }else{
    a<-c(depth*upper_clonal_vaf,depth*suppose_right_vaf)
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(clonal.vaf.left, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=clonal.vaf.left
    )
    
    df<-beta_reassign(df)
    first_vaf_list<-df[df$cluster==2,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
    
  }
  vaf_set<-c(first_vaf_list,second.vaf[second.vaf>give.vaf])
  
    # Replace this with your actual calculation
    i_values <- 1:20
    result_vector <- sapply(i_values, function(i) {
      # Replace this with your actual calculation
      calculated_value <- p/2 + (1-p)/(2 * exp(log(2) * beta * i))
      return(calculated_value)
    })
    a<-depth*result_vector
    b<-depth-a
    probs <- sapply(1:length(a), function(i) dbeta(vaf_set, a[i], b[i]))
    df<-data.frame(
      prob=probs,
      vaf=vaf_set
    )
    
    df<-beta_reassign(df)
    #df<-df[!duplicated(df),]
    if(min(df$cluster)>1){
      df$cluster<-df$cluster-(min(df$cluster))+1
    }
    cluster.result<- df %>%
      group_by(cluster) %>%
      summarise(
        count = n(),
        mean_vaf = mean(vaf)
      )
    
    cluster.result$realp<-((2 * exp(log(2) * beta * cluster.result$cluster)) * cluster.result$mean_vaf-1)/(exp(log(2) * beta * cluster.result$cluster)-1)

    updatep<-min(cluster.result[1:3,]$realp)
   print(updatep)
   print(p)
   if(abs(give.vaf*2-updatep)>0.05){
     give.vaf<-updatep/2
     p<-updatep
   }else{
     break
   }
   iter<-iter+1 
  }
  print(updatep)
  
}

normal.check<-function(mago.result.filter,result,second.vaf,maxvaf,minvaf,depth,beta,p_thre=0.05){
  num_decimal <- nchar(as.character(depth))
  maxvaf<-maxvaf[order(maxvaf$x,decreasing = T),]
  minvaf<-minvaf[order(minvaf$x,decreasing = T),]
  value.max<-mago.result.filter[1,"max"]
  value.min<-mago.result.filter[1,"min"]
  color<-c(maxvaf[maxvaf$x == value.max, 1],minvaf[minvaf$x == value.min, 1])
  clonal.vaf<-result[result$colors %in% color,]$vaf.1
  tryCatch({
      # First attempt
      m <- automixfit(clonal.vaf, type = "beta", Nc = 1:10, thresh = 0, Ninit = min(50, round(length(clonal.vaf) / 3)), k = 6, Niter.max = 10000)
    }, error = function(e) {
      tryCatch({
        # Second attempt if the first fails
        m <- automixfit(clonal.vaf, type = "beta", Nc = 1:10, thresh = 0, Ninit = 3, k = 6, Niter.max = 10000)
      }, error = function(e) {
        # Third attempt if the second also fails
        m <- automixfit(clonal.vaf, type = "beta", Nc = 1:10, thresh = 0, Niter.max = 10000)
      })
    })
  #m<-automixfit(clonal.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(clonal.vaf)/3)), Niter.max=10000)
  a=m["a",]
  b=m["b",]
  mean.a.b=a/(a+b)
  
  upper_clonal_vaf=min(mean.a.b)
  
  approx_vaf_index=which.min(mean.a.b)
  probs <- sapply(1:length(a), function(i) dbeta(clonal.vaf, a[i], b[i]))
  df<-data.frame(
    prob=probs,
    vaf=clonal.vaf
  )
  
  df<-beta_reassign(df)
  clonal.vaf.left<-df[df$cluster==approx_vaf_index,]$vaf ##Need clonal.vaf.left  upper_clonal_vaf
  
  
  if(nrow(mago.result.filter)>1){
    #value.max<-mago.result.filter[2,"max"]
    #value.min<-mago.result.filter[2,"min"]
    #color<-c(maxmaf[maxmaf$x == value.max, 1],minvaf[minvaf$x == value.min, 1])
    #second.vaf<-result[result$colors %in% color,]$vaf.1
    second.vaf=second.vaf
    second.try=T
  }else{
    temp_vaf<-(upper_clonal_vaf+0.5)/2
  second.vaf<-clonal.vaf[clonal.vaf<temp_vaf]
  second.try=F
  }
  seq_data<-Find_p_process(second.vaf,mean.a.b,upper_clonal_vaf,clonal.vaf,clonal.vaf.left,beta)
  my_seq<-seq_data$my_seq
  clear<-seq_data$clear
  if(length(my_seq)<=1 & second.try==T ){
    print('Main clone too big')
    temp_vaf<-(upper_clonal_vaf+0.5)/2
    second.vaf<-clonal.vaf[clonal.vaf<temp_vaf]
    seq_data<-Find_p_process(second.vaf,mean.a.b,upper_clonal_vaf,clonal.vaf,clonal.vaf.left,beta)
    my_seq<-seq_data$my_seq
    clear<-seq_data$clear
  }
  #print(my_seq)
  #print(clear)
  
  
  if(length(my_seq)<=0){
    print('Can not find point for estimate P')
    all.p.data <- data.frame(
      cell.div = NA,
      mu = NA,
      loglike = NA,
      bic = NA,
      aic = NA,
      p = NA,
      z_score = NA,
      reliable = NA,
      z_score1 = NA,
      bicrank = NA,
      score = 0,
      clear = NA
    )
     }else{
  
  all.p.data<-data.frame()
  for(give.vaf in my_seq){
    #print(give.vaf)
    opt.data<-NULL
    start_time <- Sys.time()
   tryCatch({
      opt.data<-Iterate_P_optimize(clear,give.vaf,upper_clonal_vaf,clonal.vaf.left,second.vaf,depth,beta,num_decimal,p_thre)
      if(nrow(opt.data[abs(opt.data$z_score)<=1.96,])<=0){
        opt.data<-opt.data[opt.data$mu>=median(opt.data$mu),]
      }
    },error=function(e){
      print(e)
    })
    
    #opt.data<-Iterate_P_optimize(clear,give.vaf,upper_clonal_vaf,clonal.vaf.left,second.vaf,depth,beta,num_decimal)
    #start_time <- Sys.time()
     end_time <- Sys.time()
    time_durations<-end_time-start_time
    print(time_durations)
    if(is.null(opt.data)){
      print('similaution failed')
    }else{
      all.p.data<-rbind(all.p.data,opt.data)
    }
  }
  #print(all.p.data)
  if(nrow(all.p.data)>0){
  all.p.data$mu<-round(all.p.data$mu,num_decimal)
  #print(all.p.data)
  if(nrow(all.p.data[abs(all.p.data$z_score)<=1.96,])>=2){
    all.p.data=all.p.data[abs(all.p.data$z_score)<=1.96,]
  }
  #print(all.p.data)
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
  
  #print(all.p.data)
  }else{
    all.p.data <- data.frame(
      cell.div = NA,
      mu = NA,
      loglike = NA,
      bic = NA,
      aic = NA,
      p = NA,
      z_score = NA,
      reliable = NA,
      z_score1 = NA,
      bicrank = NA,
      score = 0,
      clear = NA
    )
    
  }}
  return(all.p.data)
  
}
Optimize.M<-function(depth,beta,result,force=F,cut.start=1,auto=T){
  penalty=0.1
  r.data<-data.frame(vaf=numeric(),
                     r=numeric(),
                     mu=numeric()
                     
  )
  # Calculate the minimum accepted VAF based on depth
  result$vaf.1 <- result[result$vaf.1>=1/depth,]
  
  # Determine the number of decimal places to keep based on the min_vaf
  num_decimal <- nchar(as.character(depth))
 
  # Round the vaf to the determined number of decimal places
  result$vaf.1  <- round(result$vaf.1, num_decimal)
  
  if(depth <100){
    
  }
  
  if(depth>=100){
    result$count<-1
    meanvaf<-aggregate(result$vaf.1, list(result$colors), mean)
    maxvaf<-aggregate(result$vaf.1, list(result$colors), max)
    minvaf<-aggregate(result$vaf.1, list(result$colors), min)
    
    sumvaf<-aggregate(result$count, list(result$colors), sum)
    mago.result<-data.frame(max=maxvaf$x,min=minvaf$x,vaf=meanvaf$x,sum=sumvaf$x)
    mago.result$count<-1
    mago.result<-mago.result[order(mago.result$max,decreasing = T),]

    # print the updated data
    mago.result<-mago.result[mago.result$count>0,]
    
    
    print(mago.result)
    mago.result.filter<-mago.result[mago.result$max>0.25 | mago.result$min >0.25,]
    mago.result.filter$count<-1
    for (i in 1:nrow(mago.result.filter)) {
      
      
      # check if the previous range should be merged
      if (mago.result.filter$sum[i] <= 6) {
        
        # merge the previous range
        mago.result.filter$min[i-1] <- mago.result.filter$min[i]
        
        mago.result.filter$count[i] <- 0
      }
      
      
    }
    mago.result.filter<-mago.result.filter[mago.result.filter$count>0,]
    mago.result.filter<-mago.result.filter[order(mago.result.filter$max,decreasing = T),]
    
    ##First check if overwhlem bac or fitness
    main.vaf=result[result$colors %in% maxvaf[maxvaf$x ==max(maxvaf$x),1],]$vaf.1
    #Suppose fitness >90%
    
    fit.check<-fit.over.check(main.vaf,depth,beta)
    
   
    #Suppose bac over
    # Sort maxmaf$x in ascending order
    
    # Get the second highest value
    second_highest_value <-mago.result.filter[2,"max"]
    
    # Update main.vaf to consider rows corresponding to the second highest value
    second.vaf <- result[result$colors %in% maxvaf[maxvaf$x == second_highest_value, 1],]$vaf.1
    
    bac.check<-bac.over.check(main.vaf,second.vaf,depth,beta)
    
    #check magos by far max bound
    if(nrow(mago.result.filter)==1){
      value<-mago.result.filter[1,"max"]
      mix.vaf<-result[result$colors %in% maxvaf[maxvaf$x == value, 1],]$vaf.1
      pk=1
      
    }else{
      maxvaf<-maxvaf[order(maxvaf$x,decreasing = T),]
      mix.vaf<-result[result$colors %in% maxvaf[1:2, 1],]$vaf.1
      pk=2
    }
   
    
    
    
    vaf.main<-result[result$colors %in% meanvaf[meanvaf$x>0.25,1],]
    vaf.supp<-result[result$colors %in% maxvaf[maxvaf$x>0.25,1],]
    vaf.second<-result[result$colors %in% meanvaf[which(meanvaf$x<=0.25),1],]
    #minvaf[which(minvaf$x>min(minvaf$x) & minvaf$x<max(minvaf$x)),1]
    main.vaf=vaf.main$vaf.1
    main.supp=vaf.supp$vaf.1
    sub.vaf=vaf.second$vaf.1
    sub.vaf=sub.vaf[sub.vaf>1/depth]
    M.data<-main.vaf
    
    max.value=sort(mago.result$max, decreasing = TRUE)[1]
    min.value=sort(mago.result$min, decreasing = TRUE)[2]
    #mean.value=sort(mago.result$vaf, decreasing = TRUE)[2]
    
    p1=min((2^(1+beta)*max.value-1)/(2^beta-1),0.95)
    p2=max((2^(1+beta)*min.value-1)/(2^beta-1),0.01)
  
    
    
    }
  
  
  G1.sub<-data.frame(M.data)
  colnames(G1.sub)<-"vaf"
  
  G1.sub$count<-1
  
  G1.sub <- G1.sub %>%
    mutate(vaf = format(vaf, nsmall = 3)) %>%
    arrange(desc(vaf), desc(count)) %>%
    mutate(cumsum = cumsum(count)) %>%
    group_by(vaf) %>%
    slice_max(cumsum) %>%
    ungroup()
  
  #G1.sub$cumsum<-G1.sub$cumsum-init.mu.simu
  #G1.sub<-G1.sub[G1.sub$cumsum>0,]
  #print(vaf.list)
  G1.sub$vaf<-as.numeric(G1.sub$vaf)
  
  #vaf.list<-vaf.list$vaf
  
  
  
  
  
  
  
  if(p3>0.4){
    x <- seq(p3-0.05,min(0.95,p1),0.01)
    n_steps <- ceiling((min(0.95,p1) - p3) / 0.01)
    if (x[n_steps] != min(0.95,p1)) {
      x <- c(x, min(0.95,p1))
    }
  }else{
    
    
    x <- seq(max(0,p2),min(0.95,p3),0.01)
    n_steps <- ceiling((min(0.95,p3) - max(0,p2)) / 0.01)
    if (x[n_steps] != min(0.95,p3)) {
      x <- c(x, min(0.95,p3))
    }
  }
 
  
  if(force & p2<=0 & p1>=0.4){
    x <- seq(max(p3-0.05,0),min(0.95,p1),0.01)
    n_steps <- ceiling((min(0.95,p1) - max(p3-0.05,0)) / 0.01)
    if (x[n_steps] != min(0.95,p1)) {
      x <- c(x, min(0.95,p1))
    }
    
  }
  temp.test <- sapply(x, function(x) {
    tryCatch(Opt.fun(x,G1.sub,M.data), error = function(e) {
      # If an error occurs, store the value as Inf
      return(Inf)
    })
  })
  # Apply the square function to each element of the vector using lapply()
  temp.test <- unlist(temp.test)
  minvalue<-min(temp.test)
  if(p3>0.4){
    min_indices <- which(temp.test == min(temp.test))
    # find the maximum index among the minimum indices
    max_index <- max(min_indices)
    # find the corresponding value in x
    p <- x[max_index]
   
  }else{
    p<-x[which.min(temp.test)]
  }
  
  
  if(p>0.4){
    print('**** deep optimize*****')
    p.temp<-p3-0.01
    if(p3>0.4){
      if(p1>0.8){
        x <- seq(p.temp,min(0.95,p1+0.1),0.002)
        n_steps <- ceiling((min(0.95,p1+0.1) - p.temp) / 0.002)
        if (x[n_steps] != min(0.95,p1+0.1)) {
          x <- c(x, min(0.95,p1+0.1))
        }
      }else{
        x <- seq(p.temp,min(0.8,p1+0.05),0.002)
        n_steps <- ceiling((min(0.8,p1+0.05) - p.temp) / 0.002)
        if (x[n_steps] != min(0.8,p1+0.05)) {
          x <- c(x, min(0.8,p1+0.05))
        }
      }
      
     
    }else{
      
      
      x <- seq(max(0,p2),min(0.95,p*1.2),0.01)
      n_steps <- ceiling((min(0.95,p*1.2) - max(0,p2)) / 0.01)
      if (x[n_steps] != min(0.95,p*1.2)) {
        x <- c(x, min(0.95,p*1.2))
      }
    }

    if(force & p2<=0 & p1>=0.4){
      x <- seq(max(p3-0.01,0),p1+0.05,0.002)
      n_steps <- ceiling((p1+0.05 - max(p3-0.01,0)) / 0.002)
      if (x[n_steps] != p1+0.05) {
        x <- c(x, p1+0.05)
      }
      
    }
    temp.test <- sapply(x, function(x) {
      tryCatch(overall.opt(M.data,x,depth,beta,mago.result,sub.vaf,checkp), error = function(e) {
        # If an error occurs, store the value as Inf
        return(Inf)
      })
    })
    # Apply the square function to each element of the vector using lapply()
    temp.test <- unlist(temp.test)
    minvalue<-min(temp.test)
    if(p3>0.4){
      min_indices <- which(temp.test == min(temp.test))
      # find the maximum index among the minimum indices
      max_index <- max(min_indices)
      # find the corresponding value in x
      p <- x[max_index]
      
    }else{
      p<-x[which.min(temp.test)]
    }
    
  }
  
  print("****optimize finish****")
  if(p>0.4){
    result.data<-Get.result(G1.sub,p,M.data,sub.vaf,mago.result,depth,checkp,beta)
  }else{
    G1.sub$lnf<-1/G1.sub$vaf
    G1<-G1.sub[G1.sub$lnf>0,]
    G1$x<-G1$lnf
    #plot(G1$x,G1$cumsum)
    vaf.first<-p/2 + (1-p)/(2*2^beta)
    compare.vec<-c(0.5,vaf.first)
    compare.vec<-sort(compare.vec)
    #print(vec)
    # Find the indices of the closest values in "vec" for each element in "assgin.vaf"
    compare.vec<-data.frame(compare.vec)
    compare.vec$a<-depth*compare.vec$compare.vec
    compare.vec$b<-depth-compare.vec$a
    probs <- sapply(1:nrow(compare.vec), function(i) dbeta(M.data, compare.vec$a[i], compare.vec$b[i]))
    
    # Assign each number to the largest beta distribution
    max_rows <- apply(probs, 1, which.max)
    #cluster <- sapply(M.data, function(x) which.min(abs(x - compare.vec)))
    
    cluster.count <- table(max_rows)[nrow(compare.vec)]
    M.cul<-length(M.data)-cluster.count
    
    initial.vaf.cut<-G1[G1$cumsum==cluster.count,]
    assgin.vaf<-M.data[M.data<initial.vaf.cut$vaf]
    assgin.vaf<-sort(assgin.vaf,decreasing = T)
    G1<-G1[G1$cumsum>cluster.count,]
    plot(G1$x,G1$cumsum)
    # Calculate the number of rows to remove from each end
    num_rows_to_remove <- round(nrow(G1) * 0.05)
    
    # Calculate the index of the rows to remove
    rows_to_remove <- c(1:num_rows_to_remove, (nrow(G1) - num_rows_to_remove + 1):nrow(G1))
    
    # Subset the data frame to remove the rows
    G1 <- G1[-rows_to_remove,]
    fm2<-lm(cumsum~x,data=G1)
    mu<-summary(fm2)$coefficients[2,1]
    r.n<-summary(fm2)$r.squared
    if (depth>=100){
    tryCatch({
      mu.turn <-breakpoints(cumsum~x,data=G1,h=3/nrow(G1))
      
      #fm1 <- lm(cumsum ~ breakfactor(mu.turn), data = G1)
      # Extract the coefficients from the model
      # Create a break factor
      bf <- breakfactor(mu.turn)
      
      # Initialize a vector to store AIC values
      AIC_values <- c()
      BIC_values<- c()
      slopes<-c()
      # Loop over each segment and calculate the AIC
      for(i in unique(bf)) {
        segment_data <- subset(G1, bf == i)
        fm_segment <- lm(cumsum ~ x, data = segment_data)
        #coeffs <- coef(fm_segment)
        slope <- summary(fm_segment)$coefficients[2,1]
        slopes<-c(slopes,slope)
        AIC_values <-c(AIC_values, AIC(fm_segment))
        BIC_values<-c(BIC_values, BIC(fm_segment))
      }
      

      min_aic <- which.min(AIC_values)
      min_Bic <- which.min(BIC_values)
      print(slopes)
      print(is_zigzag(slopes))
      print(AIC_values)
      print(BIC_values)
      mu=(slopes[min_aic]+slopes[min_Bic])/2
      #print(mu)
      count=length(slopes)
    },error=function(e){
      count=1
      fm1<-lm(cumsum~x,data=G1)
      mu=summary(fm1)$coefficients[2,1]*(-1)*beta*log(2)
      #mu_high=
    })
    }
    # Subset the data frame to remove the rows
    print(mu)
    result.data<-list(pf = p, mu = mu,initalcount=cluster.count,initalvaf=initial.vaf.cut$vaf,rf=r.n)
  }
 
  return(result.data)
}




  ## ---- run_est.R ----
# Bundled reference: production estimation pipeline. Loaded into a private
# environment by the dispatcher in zzz.R; `library()` and `source()` lines
# stripped — `optimizeM.R` and `optimizeS.R` are sourced into the same env
# by the dispatcher before this file. `parallel` is in the package Imports.
beta=0.9
backrun<-function(main.vaf,second.vaf,depth,beta,p_thre=0.05,filter=T){

  mulist<-c()
  mulist1<-c()
  mulist2<-c()
  biclist<-c()
  plist<-c()
  plist1<-c()
  data.all<-data.frame()
  for(try in 1:3){
    fit.check<-bac.over.check(main.vaf,second.vaf,depth,beta,p_thre)
    fit.check <- fit.check[order(fit.check$bic),]
    fit.check$try<-try
    data.all<-rbind(data.all,fit.check)
    ##not to remove out
    if(nrow(fit.check[fit.check$mu>3,])>1){
      fit.check<-fit.check[fit.check$mu>3,]
    }
    ##remove out
    if(nrow(fit.check)>3){
      
      # Calculate the mean and standard deviation without the maximum value
      mean_filtered <- median(fit.check$cell.div)
      sd_filtered <-mad(fit.check$cell.div)
      fit.check$z_score <- abs((fit.check$cell.div - mean_filtered) / sd_filtered)
      fit.check <- fit.check[fit.check$z_score < 2,]
      
    }
    mulist1<-c(mulist1,median(fit.check$mu))
    mulist2<-c(mulist1,max(fit.check$mu))
    biclist<-c(biclist,fit.check[fit.check$mu==max(fit.check$mu),]$bic)
    fit.check <- fit.check[1,]
    p=fit.check[1,]$p
    mu.select1=fit.check$mu
    mulist<-c(mulist,mu.select1)
    
    plist<-c(plist,fit.check[1,]$p)
    plist1<-c(plist1,fit.check[which.min(fit.check$mu-median(fit.check$mu)),]$p)
    #print(fit.check)
  }
  if(length(mulist)>0){
    
    mu.select2<-mulist[which.max(mulist)]
    mu.select3<- mulist1[which.max(mulist1)]
    
    mu.select4<-mulist2[which.min(biclist)]
    
    p1<-plist[which.min(mulist-mu.select2)]
    p2<-plist1[which.min(mulist1-mu.select3)]
  }else{
    mu.select2<-NA
    mu.select3<-NA
    mu.select4<-NA
    p1<-NA
    p2<-NA
  }
  
  upper_bound<-length(second.vaf)/3
  #print(fit.check)
  select.data<-data.frame(
    mu1=mu.select2,
    mu2=mu.select3,
    mu3=mu.select4,
    upper_bound=upper_bound,
    p1=p1,
    p2=p2
  )
  if(filter){
    result_list <- list(select.data = select.data, data.all = data.all)
    return(result_list)}else{
    return(data.all)
  }
}
fitrun<-function(main.vaf,depth,beta,p_thre=0.05,filter=T){
  mulist<-c()
  mulist1<-c()
  mulist2<-c()
  biclist<-c()
  plist<-c()
  plist1<-c()
  upperlist1<-c()
  upperlist2<-c()
  data.all<-data.frame()
  for(try in 1:3){
    fit.check<-fit.over.check(main.vaf,depth,beta,p_thre)
    fit.check <- fit.check[order(fit.check$bic),]
    fit.check$try<-try
    data.all<-rbind(data.all,fit.check)
    ##not to remove out
    if(nrow(fit.check[fit.check$mu>3,])>1){
      fit.check<-fit.check[fit.check$mu>3,]
    }
    ##remove out
    if(nrow(fit.check)>3){
      
      # Calculate the mean and standard deviation without the maximum value
      mean_filtered <- median(fit.check$cell.div)
      sd_filtered <-mad(fit.check$cell.div)
      fit.check$z_score <- abs((fit.check$cell.div - mean_filtered) / sd_filtered)
      fit.check <- fit.check[fit.check$z_score < 2,]
     
    }
    mulist1<-c(mulist1,median(fit.check$mu))
    mulist2<-c(mulist1,max(fit.check$mu))
    biclist<-c(biclist,fit.check[fit.check$mu==max(fit.check$mu),]$bic)
    fit.check <- fit.check[1,]
    p=fit.check[1,]$p
    mu.select1=fit.check$mu
    mulist<-c(mulist,mu.select1)
    upperlist1<-c(upperlist1,fit.check$lowerbound1)
    upperlist2<-c(upperlist2,fit.check$lowerbound2)
    plist<-c(plist,fit.check[1,]$p)
    plist1<-c(plist1,fit.check[which.min(fit.check$mu-median(fit.check$mu)),]$p)
    #print(fit.check)
  }
  if(length(mulist)>0){
    
    mu.select2<-mulist[which.max(mulist)]
    mu.select3<- mulist1[which.max(mulist1)]
    
    mu.select4<-mulist2[which.min(biclist)]
    up1<-upperlist1[which.max(upperlist1)]
    up2<-upperlist1[which.max(upperlist1)]
    p1<-plist[which.min(mulist-mu.select2)]
    p2<-plist1[which.min(mulist1-mu.select3)]
  }else{
    mu.select2<-NA
    mu.select3<-NA
    mu.select4<-NA
    up1<-NA
    up2<-NA
    p1<-NA
    p2<-NA
  }

  #print(fit.check)
  select.data<-data.frame(
    mu1=mu.select2,
    mu2=mu.select3,
    mu3=mu.select4,
    up1=up1,
    up2=up2,
    p1=p1,
    p2=p2
  )
  if(filter){
    result_list <- list(select.data = select.data, data.all = data.all)
    return(result_list)}else{
      return(data.all)
    }
  
  
  
}

normalrun<-function(mago.result.filter,result,second.vaf,maxvaf,minvaf,depth,beta,p_thre=0.05,filter=T){
  mulist<-c()
  mulist1<-c()
  mulist2<-c()
  biclist<-c()
  plist<-c()
  plist1<-c()
  data.all<-data.frame()
  for(try in 1:3){
    fit.check<-normal.check(mago.result.filter,result,second.vaf,maxvaf,minvaf,depth,beta,p_thre)
    fit.check <- fit.check[order(-fit.check$score),]
    fit.check$try<-try
    data.all<-rbind(data.all,fit.check)
    if(!is.na(fit.check$mu[1])){   # R>=4.2: if() needs length-1; [1] = pre-4.2 first-element behavior (top-scored row)
      # Pick the row with the maximum score
      mulist1<-c(mulist1,median(fit.check$mu))
      mulist2<-c(mulist1,max(fit.check$mu))
      biclist<-c(biclist,fit.check[fit.check$mu==max(fit.check$mu),]$bic)
      
      fit.check <- fit.check[1,]
      p=fit.check[1,]$p
      mu.select1=fit.check$mu
      mulist<-c(mulist,mu.select1)
      
      plist<-c(plist,fit.check[1,]$p)
      plist1<-c(plist1,fit.check[which.min(fit.check$mu-median(fit.check$mu)),]$p)
    }
  }
  if(length(mulist)>0){
    
    mu.select2<-mulist[which.max(mulist)]
    mu.select3<- mulist1[which.max(mulist1)]
    mu.select4<-mulist2[which.min(biclist)]
    p1<-plist[which.min(mulist-mu.select2)]
    p2<-plist1[which.min(mulist1-mu.select3)]
  }else{
    mu.select2<-NA
    mu.select3<-NA
    mu.select4<-NA
    p1<-NA
    p2<-NA
  }
  #cell.div=fit.check[fit.check$score==mu.select1,]$cell.div
  select.data<-data.frame(
    mu1=mu.select2,
    mu2=mu.select3,
    mu3=mu.select4,
    p1=p1,
    p2=p2
  )
  if(filter){
    result_list <- list(select.data = select.data, data.all = data.all)
    return(result_list)
    }else{
      return(data.all)
    }
  
}
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

pick_s<-function(df){
  df$bicrank<-rank(df$bic)
  df$score<-df$vaf
  if(max(df$start)>1){
    df<-keep_longest_consecutive_rows(df)
  }
  
  #   #drop
  #   
  
  
  if(max(df$start)<1 ){
    #if(max(data.pick$start)<1 | (any(data.pick$bic > 0) && any(data.pick$bic < 0)) ){
    
    #data.pick<-data.pick[data.pick$bicrank <=ceiling(nrow(data.pick) * 0.5),]
    pickdata<-df[df$bic==min(df$bic),]
  }else{
    data.pick.sub<-df
    data.pick.sub$truescore <- 0:(nrow(data.pick.sub)-1)
    data.pick.sub<-data.pick.sub[data.pick.sub$truescore<=data.pick.sub$start & data.pick.sub$start>0,]
    if(nrow(data.pick.sub)>0){
      #data.pick.sub<-data.pick.sub[data.pick$bicrank <=ceiling(nrow(data.pick) * 0.5),]
      pickdata<-df[df$bic==min(data.pick.sub$bic),]
    }else{
      #data.pick<-data.pick[data.pick$bicrank <= ceiling(nrow(data.pick) * 0.5),]
      pickdata<-df[df$bic==min(df$bic),]
    }
  }
  
  
  
  return(pickdata)
  
}

fit_update<-function(main.vaf,depth,beta){
  temp.vaf<-(0.5+min(main.vaf))/2
  inita=depth*c(0.5,temp.vaf)
  initb=depth-inita
  probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
  df<-data.frame(
    prob=probs,
    vaf=main.vaf
  )
  
  # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
  df<-beta_reassign(df)
  ##if s close to cell 2, 4, 8,(higher less possible)
  df<-df[df$cluster>1,]
  new.vaf<-mean(df$vaf)
  old.vaf<-temp.vaf
  diff<-abs(new.vaf-old.vaf)
  iter=0
  while(diff>1/depth){
    old.vaf<-new.vaf
    inita=depth*c(0.5,old.vaf)
    initb=depth-inita
    probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
    df<-data.frame(
      prob=probs,
      vaf=main.vaf
    )
    
    # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
    df<-beta_reassign(df)
    ##if s close to cell 2, 4, 8,(higher less possible)
    df<-df[df$cluster>1,]
    new.vaf<-mean(df$vaf)
    diff<-abs(new.vaf-old.vaf)
    iter=iter+1
    if(iter>100){
      break
    }
    
  }
  
  return(new.vaf)
}
fit_step<-function(predict.result,result,depth=120,beta=0.9){
  p.list=unique(predict.result$p)
  #fit.check<-fit.over.check(main.vaf,depth,beta)
  s.results <- lapply(p.list, function(p) {
    df<-Get.S(result, p, beta, depth, p_thre=1e-6)
    df<-df[!is.na(df$s), ]
    if(nrow(df)>0){
    df$p_value <- p
    data<-pick_s(df)
    murange=Simulate_ratio_peak(result,p,beta,depth,data$s)
    data$minmu=min(murange)
    data$maxmu=max(murange)
    return(data)
    }
  })
  fit.data<-data.frame()
  combined_df <- do.call(rbind, s.results)
  for(try in 1:length(p.list)){
    p.pick=p.list[try]
    fit.check=predict.result[predict.result$p==p.pick,]
    
    if(nrow(fit.check)>0){
      minmu=combined_df[combined_df$p_value==p.pick,]$minmu
      maxmu=combined_df[combined_df$p_value==p.pick,]$maxmu
      
      fit.check.sub<-fit.check[fit.check$mu >=minmu & fit.check$mu<=maxmu, ]
      fit.check.sub=na.omit(fit.check.sub)
      #print(fit.check)
      if(nrow(fit.check.sub)==0){
        minmu=minmu-1.96*1.5
        fit.check.sub<-fit.check[fit.check$mu >=minmu & fit.check$mu<=maxmu, ]
        fit.check.sub=na.omit(fit.check.sub)
      }
      if(nrow(fit.check.sub)>0){
        fit.data<-rbind(fit.data,fit.check.sub[fit.check.sub$bic==min(fit.check.sub$bic),])
      }
    }
  }
  
  # Find rows with minimum bic
  if(nrow(fit.data)>0){
    # Find rows with minimum bic
    min_bic_rows <- fit.data[fit.data$bic == min(fit.data$bic), ]
    
    # If multiple rows have the same minimum bic, randomly select one
    selected_row <- min_bic_rows[sample(nrow(min_bic_rows), 1), ]
    
    selected_row$s<-combined_df[combined_df$p_value==selected_row$p,]$s
    return(selected_row)
  }else{
    selected_row<-combined_df
    return(selected_row)
  }
}

inter_step<-function(predict.result,result,depth=120,beta=0.9){
  p.list=unique(predict.result$p)
  #fit.check<-fit.over.check(main.vaf,depth,beta)
  s.results <- lapply(p.list, function(p) {
    #print(p)
    df<-Get.S(result, p, beta, depth, p_thre=1e-6)
   df<-df[!is.na(df$s), ]                                                      
    if(nrow(df)>0){ 
    df$p_value <- p
    data<-pick_s(df)
    murange=Simulate_ratio_peak(result,p,beta,depth,data$s)
    data$minmu=min(murange)
    data$maxmu=max(murange)
    return(data)
    }
  })
  fit.data<-data.frame()
  combined_df <- do.call(rbind, s.results)
  for(try in 1:length(p.list)){
    p.pick=p.list[try]
    fit.check=predict.result[predict.result$p==p.pick,]
    
    if(nrow(fit.check)>0){
      minmu=combined_df[combined_df$p_value==p.pick,]$minmu
      maxmu=combined_df[combined_df$p_value==p.pick,]$maxmu
      
      fit.check.sub<-fit.check[fit.check$mu >=minmu & fit.check$mu<=maxmu, ]
      fit.check.sub=na.omit(fit.check.sub)
      #print(fit.check)
      if(nrow(fit.check.sub)==0){
        minmu=minmu-1.96*1.5
        fit.check.sub<-fit.check[fit.check$mu >=minmu & fit.check$mu<=maxmu, ]
        fit.check.sub=na.omit(fit.check.sub)
      }
      if(nrow(fit.check.sub)>0){
        fit.data<-rbind(fit.data,fit.check.sub[fit.check.sub$bic==min(fit.check.sub$bic),])
      }
    }
  }
  
  
  if(nrow(fit.data)>0){
    # Find rows with minimum bic
    min_bic_rows <- fit.data[fit.data$bic == min(fit.data$bic), ]
    
    # If multiple rows have the same minimum bic, randomly select one
    selected_row <- min_bic_rows[sample(nrow(min_bic_rows), 1), ]
    
    selected_row$s<-combined_df[combined_df$p_value==selected_row$p,]$s
    return(selected_row)
  }else{
    selected_row<-combined_df
    return(selected_row)
  }
  #return(selected_row)
}

Get_simu_para<-function(pathname,filename){
  simulations = readRDS(paste0(pathname,filename))
  simulation_name=names(simulations)
  para.frame=data.frame()
  for(i in 1:length(simulation_name)){
    para.frame[i,'id']=simulation_name[i]
    content_simu=simulations[[simulation_name[i]]]
    clone_parameters=content_simu$clone_parameters[1:2,]
    #mu=mean(clone_parameters$mutation_rate)
    mu=mean(clone_parameters$mutation_rates)
    #start_time=max(clone_parameters$start_time)
    start_time=max(clone_parameters$clone_start_times)
    #s=(max(clone_parameters$birth_rate)-min(clone_parameters$birth_rate))/min(clone_parameters$birth_rate)
    s=(max(clone_parameters$birthrates)-min(clone_parameters$birthrates))/min(clone_parameters$birthrates)
    #end_time=as.integer(content_simu$simulation_parameters[4])
    end_time=as.integer(content_simu$simulation_parameters[1])
    para.frame[i,'murate']=mu
    para.frame[i,'s']=s
    para.frame[i,'t1']=start_time
    para.frame[i,'tend']=end_time
    para.frame[i,'fitp']=content_simu$cell_numbers[2]/sum(content_simu$cell_numbers)
  }
  return(para.frame)
}

fit.all.run<-function(sid,depth,beta,savepath,dpcode){

  
    sample_name=sid
    wgs.ex<-readRDS(paste0(savepath,sample_name,"/All_",sample_name,".rds"))
    data<-data.frame()
    data2<-data.frame()
    result<-wgs.ex$result
    result$count<-1
    meanvaf<-aggregate(result$vaf.1, list(result$colors), mean)
    maxvaf<-aggregate(result$vaf.1, list(result$colors), max)
    minvaf<-aggregate(result$vaf.1, list(result$colors), min)
    
    sumvaf<-aggregate(result$count, list(result$colors), sum)
    mago.result<-data.frame(max=maxvaf$x,min=minvaf$x,vaf=meanvaf$x,sum=sumvaf$x)
    mago.result$count<-1
    mago.result<-mago.result[order(mago.result$max,decreasing = T),]
    
    # print the updated data
    mago.result<-mago.result[mago.result$count>0,]
    
    mago.result.filter<-mago.result[mago.result$max>0.25 | mago.result$min >0.25,]
    mago.result.filter$count<-1
    for (i in 1:nrow(mago.result.filter)) {
      
      
      # check if the previous range should be merged
      if (mago.result.filter$sum[i] <= 6) {
        
        # merge the previous range
        mago.result.filter$min[i-1] <- mago.result.filter$min[i]
        
        mago.result.filter$count[i] <- 0
      }
      
      
    }
    mago.result.filter<-mago.result.filter[mago.result.filter$count>0,]
    mago.result.filter<-mago.result.filter[order(mago.result.filter$max,decreasing = T),]
    
    
    
    ##First check if overwhlem bac or fitness
    main.vaf=result[result$colors %in% maxvaf[maxvaf$x ==max(maxvaf$x),1],]$vaf.1
    #Suppose fitness >90%

    #fit.check<-fit.over.check(main.vaf,depth,beta)
    second_highest_value <-mago.result.filter[2,"max"]
    
    # Update main.vaf to consider rows corresponding to the second highest value
    second.vaf <- result[result$colors %in% maxvaf[maxvaf$x == second_highest_value, 1],]$vaf.1
    
    
    ##fitness run
    possibleError <- tryCatch({
      fit.check<-fitrun(main.vaf,depth,beta,p_thre=0.01)
      
      
    },error=function(e){
      e
    })
    if(!inherits(possibleError, "error")){
      #REAL WORK
      select_data_result <- fit.check$select.data  # or result[[1]]
      
      # To access data.all
      data_all_result <- fit.check$data.all  # or result[[2]]
      select_data_result$pt<-0.01
      data_all_result$pt<-0.01
      select_data_result$name<-sample_name
      data_all_result$name<-sample_name
      print('binding1')
      data<-rbind(data,select_data_result)
      print('binding2')
      data2<-rbind(data2,data_all_result)
      #print()
    }else{
      # Create data frame with the same structure as select_data_result
      new_select_data <- data.frame(
        mu1 = NA, mu2 = NA, mu3 = NA, up1 = NA, up2 = NA, 
        p1 = NA, p2 = NA, pt = NA, name = sample_name, stringsAsFactors = FALSE
      )
      
      # Create data frame with the same structure as data_all_result
      new_data_all <- data.frame(
        cell.div = NA, mu = NA, loglike = NA, bic = NA, 
        aic = NA, z_score = NA, p = NA, z_score1 = NA, 
        bicrank = NA, score = NA, lowerbound1 = NA, 
        lowerbound2 = NA, try = NA, pt = NA, name = sample_name, stringsAsFactors = FALSE
      )
      
      # Append NA row to 'data'
      data <- rbind(data, new_select_data)
      
      # Create NA row for 'data2' with just the 'name' column filled
      
      # Append NA row to 'data2'
      data2 <- rbind(data2, new_data_all)
    }
    

  

  saveRDS(data,paste0(savepath,sample_name,'/normal.result.fitdomi.depth',dpcode,'.rds'))
  saveRDS(data2,paste0(savepath,sample_name,'/normal.result.fitdomi.allpoint.depth',dpcode,'.rds'))
}
  
  


bac.all.run<-function(sid,depth,beta,savepath,dpcode){
    
    sample_name=sid
    wgs.ex<-readRDS(paste0(savepath,sample_name,"/All_",sample_name,".rds"))
    result<-wgs.ex$result
    data<-data.frame()
    data2<-data.frame()
    #plot(result$vaf.1, result$depth.1, col= result$colors, xlab= 'VAF', ylab='Depth') 
    
    result$count<-1
    meanvaf<-aggregate(result$vaf.1, list(result$colors), mean)
    maxvaf<-aggregate(result$vaf.1, list(result$colors), max)
    minvaf<-aggregate(result$vaf.1, list(result$colors), min)
    
    sumvaf<-aggregate(result$count, list(result$colors), sum)
    mago.result<-data.frame(max=maxvaf$x,min=minvaf$x,vaf=meanvaf$x,sum=sumvaf$x)
    mago.result$count<-1
    mago.result<-mago.result[order(mago.result$max,decreasing = T),]
    
    # print the updated data
    mago.result<-mago.result[mago.result$count>0,]
    
    mago.result.filter<-mago.result[mago.result$max>0.25 | mago.result$min >0.25,]
    
    mago.result.filter$count<-1
    for (i in 1:nrow(mago.result.filter)) {
      
      
      # check if the previous range should be merged
      if (mago.result.filter$sum[i] <= 6) {
        
        # merge the previous range
        mago.result.filter$min[i-1] <- mago.result.filter$min[i]
        
        mago.result.filter$count[i] <- 0
      }
      
      
    }
    mago.result.filter<-mago.result.filter[mago.result.filter$count>0,]
    mago.result.filter<-mago.result.filter[order(mago.result.filter$max,decreasing = T),]
    
    if(nrow(mago.result.filter)==1){
      print('warning sample!')
      mago.result.filter<-rbind(mago.result.filter,mago.result[2,])
    }
    
    
    ##First check if overwhlem bac or fitness
    main.vaf=result[result$colors %in% maxvaf[maxvaf$x ==max(maxvaf$x),1],]$vaf.1
    #Suppose fitness >90%
    #fit.check<-fit.over.check(main.vaf,depth,beta)
    second_highest_value <-mago.result.filter[2,"max"]
    
    # Update main.vaf to consider rows corresponding to the second highest value
    second.vaf <- result[result$colors %in% maxvaf[maxvaf$x %in% second_highest_value, 1],]$vaf.1
    
    
    
    ##fitness run
    
   
    possibleError <- tryCatch({
      fit.check<-backrun(main.vaf,second.vaf,depth,beta,p_thre=0.01)
      
      
    },error=function(e){
      e
    })
    if(!inherits(possibleError, "error")){
    select_data_result <- fit.check$select.data  # or result[[1]]
    
    # To access data.all
    data_all_result <- fit.check$data.all  # or result[[2]]
    select_data_result$pt<-0.01
    data_all_result$pt<-0.01
    select_data_result$name<-sample_name
    data_all_result$name<-sample_name
    data<-rbind(data,select_data_result)
    data2<-rbind(data2,data_all_result)
    }else{
      new_select_data <- data.frame(
        mu1 = NA, mu2 = NA, mu3 = NA, upper_bound = NA, 
        p1 = NA, p2 = NA, pt = NA, name = sample_name, stringsAsFactors = FALSE
      )
      
      # Create data frame with the same structure as data_all_result
      new_data_all <- data.frame(
        cell.div = NA, mu = NA, loglike = NA, bic = NA, 
        aic = NA, p = NA, try = NA, pt = NA, name = sample_name, stringsAsFactors = FALSE
      )
      
      data <- rbind(data, new_select_data)
      
      # Append NA row to 'data2'
      data2 <- rbind(data2, new_data_all)
    }
    

  saveRDS(data,paste0(savepath,sample_name,'/normal.result.bacdomi.depth',dpcode,'.rds'))
  saveRDS(data2,paste0(savepath,sample_name,'/normal.result.bacdomi.allpoint.depth',dpcode,'.rds'))
  
  
  
}

normal.all.run<-function(sid,depth,beta,savepath,dpcode){
 
  data<-data.frame()
  data2<-data.frame()
 
   sample_name=sid
  wgs.ex<-readRDS(paste0(savepath,sample_name,"/All_",sample_name,".rds"))
  result<-wgs.ex$result
    
    #plot(result$vaf.1, result$depth.1, col= result$colors, xlab= 'VAF', ylab='Depth') 
    
    result$count<-1
    meanvaf<-aggregate(result$vaf.1, list(result$colors), mean)
    maxvaf<-aggregate(result$vaf.1, list(result$colors), max)
    minvaf<-aggregate(result$vaf.1, list(result$colors), min)
    
    sumvaf<-aggregate(result$count, list(result$colors), sum)
    mago.result<-data.frame(max=maxvaf$x,min=minvaf$x,vaf=meanvaf$x,sum=sumvaf$x)
    mago.result$count<-1
    mago.result<-mago.result[order(mago.result$max,decreasing = T),]
    
    # print the updated data
    mago.result<-mago.result[mago.result$count>0,]
    
    mago.result.filter<-mago.result[mago.result$max>0.25 | mago.result$min >0.25,]
    mago.result.filter$count<-1
    for (i in 1:nrow(mago.result.filter)) {
      
      
      # check if the previous range should be merged
      if (mago.result.filter$sum[i] <= 6) {
        
        # merge the previous range
        mago.result.filter$min[i-1] <- mago.result.filter$min[i]
        
        mago.result.filter$count[i] <- 0
      }
      
      
    }
    mago.result.filter<-mago.result.filter[mago.result.filter$count>0,]
    mago.result.filter<-mago.result.filter[order(mago.result.filter$max,decreasing = T),]
    
    
    
    
    ##First check if overwhlem bac or fitness
    main.vaf=result[result$colors %in% maxvaf[maxvaf$x ==max(maxvaf$x),1],]$vaf.1
    #Suppose fitness >90%
    #fit.check<-fit.over.check(main.vaf,depth,beta)
    #second_highest_value <-mago.result.filter[2:nrow(mago.result.filter),"max"]
    if(nrow(mago.result.filter)==1){
      print('warning sample!')
      second_highest_value <-mago.result[2,"max"]
      # Update main.vaf to consider rows corresponding to the second highest value
      second.vaf <- result[result$colors %in% maxvaf[maxvaf$x %in% second_highest_value, 1],]$vaf.1
      
    }else{
      second_highest_value <-mago.result.filter[2,"max"]
      # Update main.vaf to consider rows corresponding to the second highest value
      second.vaf <- result[result$colors %in% maxvaf[maxvaf$x %in% second_highest_value, 1],]$vaf.1
      
    }
    
    ##fitness run
    
    possibleError <- tryCatch({
      fit.check<-normalrun(mago.result.filter,result,second.vaf,maxvaf,minvaf,depth,beta,p_thre=0.01)
      
      
    },error=function(e){
      e
    })
    if(!inherits(possibleError, "error")){
      #REAL WORK
      select_data_result <- fit.check$select.data  # or result[[1]]
      
      # To access data.all
      data_all_result <- fit.check$data.all  # or result[[2]]
      select_data_result$pt<-0.01
      data_all_result$pt<-0.01
      select_data_result$name<-sample_name
      data_all_result$name<-sample_name
      data<-rbind(data,select_data_result)
      data2<-rbind(data2,data_all_result)
    }else{
      new_select_data <- data.frame(
        mu1 = NA, mu2 = NA, mu3 = NA, p1 = NA, p2 = NA, 
        pt = NA, name = sample_name, stringsAsFactors = FALSE
      )
        
      # Create data frame with the same structure as data_all_result
      new_data_all <- data.frame(
        cell.div = NA, mu = NA, loglike = NA, bic = NA, 
        aic = NA, p = NA, z_score = NA, reliable = NA, 
        z_score1 = NA, bicrank = NA, score = NA, clear = NA, 
        try = NA, pt = NA, name = sample_name, stringsAsFactors = FALSE
      )
      
      # Append NA row to 'data'
      data <- rbind(data, new_select_data)
      
      
      # Append NA row to 'data2'
      data2 <- rbind(data2, new_data_all)
    }
    
    
    
  
  
  saveRDS(data,paste0(savepath,sample_name,'/normal.result.fitnormalnew.depth',dpcode,'.rds'))
  saveRDS(data2,paste0(savepath,sample_name,'/normal.result.fitnormalnew.allpoint.depth',dpcode,'.rds'))
  
}


All.guess.update<-function(sid,savepath,depth,beta,dpcode){
  
    sample_name<-sid
   fit.predict.all<-readRDS(paste0(savepath,sample_name,'/normal.result.fitdomi.depth',dpcode,'.rds'))
   fit.all<-readRDS(paste0(savepath,sample_name,'/normal.result.fitdomi.allpoint.depth',dpcode,'.rds'))

   bac.all<-readRDS(paste0(savepath,sample_name,'/normal.result.bacdomi.allpoint.depth',dpcode,'.rds'))

  #data.inter<-readRDS('/Users/hchen/R_project/Model_Caner/normal.result.fitnormalnew.allpoint.nomode.rds')
  #data.inter.other<-readRDS('/Users/hchen/R_project/Model_Caner/Other.sample.normal.result.fitnormalnew.allpoint.nomode.rds')
  inter.all<-readRDS(paste0(savepath,sample_name,'/normal.result.fitnormalnew.allpoint.depth',dpcode,'.rds'))
 
  bac.predict.all<-readRDS(paste0(savepath,sample_name,'/normal.result.bacdomi.depth',dpcode,'.rds'))
  
  inter.predict.all<-readRDS(paste0(savepath,sample_name,'/normal.result.fitnormalnew.depth',dpcode,'.rds'))
 
 
 # lower<-Suppose_para
  data.rearrange<-data.frame()
  


   
    sample_name=sid
    i_index=1
    data.rearrange[i_index,"name"]=sample_name
   
    #wgs.ex<-readRDS(paste0("/Users/hchen/R_project/Model_Caner/mob_simulation_magos/result/",sample_name,"/All_",sample_name,".rds"))
    wgs.ex<-readRDS(paste0(savepath,sample_name,"/All_",sample_name,".rds"))
    result<-wgs.ex$result
    
    fit.result<-fit.all[fit.all$name==sample_name &fit.all$pt==0.01,]
    inter.result<-inter.all[inter.all$name==sample_name &inter.all$pt==0.01,]
    bac.result<-bac.all[bac.all$name==sample_name &bac.all$pt==0.01,]
    bac.upper<-bac.predict.all[bac.predict.all$name==sample_name &bac.predict.all$pt==0.01,]
    
    if(nrow(na.omit(fit.result))>0){
      fit.result=na.omit(fit.result)
      fit.data<-fit_step(fit.result,result,depth=depth)
      
      if(!("mu" %in% colnames(fit.data))) {
        mu.fit.pick=fit.predict.all[fit.predict.all$name==sample_name &fit.predict.all$pt==0.01, ]$up1
        
        fit.select.row<-fit.result[fit.result$lowerbound1==mu.fit.pick,]
        fit.select.row<-fit.select.row[sample(nrow(fit.select.row), 1),]
        p<-fit.select.row$p
        mu.s.update<-fit.data[fit.data$p_value==p,]
        mu.select2<-mu.s.update[sample(nrow(mu.s.update), 1),]$minmu
        
        df<-Get.S(result, p, beta=0.9, depth=120, p_thre=1e-6)
        df$p_value <- p
        data<-pick_s(df)
        data.rearrange[i_index,"fitmu"]=mu.fit.pick
        data.rearrange[i_index,"fitmu_candidate"]=mu.select2
        data.rearrange[i_index,"fitcell"]=fit.select.row$cell.div
        data.rearrange[i_index,"fits"]=data$s
        data.rearrange[i_index,"fitp"]=p
        
      }else{
        data.rearrange[i_index,"fitmu"]=fit.data$mu
        data.rearrange[i_index,"fitmu_candidate"]=1
        data.rearrange[i_index,"fitcell"]=fit.data$cell.div
        data.rearrange[i_index,"fits"]=fit.data$s
        data.rearrange[i_index,"fitp"]=fit.data$p
      }
    }else{
      data.rearrange[i_index,"fitmu"]=NA
      data.rearrange[i_index,"fitmu_candidate"]=1
      data.rearrange[i_index,"fitcell"]=NA
      data.rearrange[i_index,"fits"]=NA
      data.rearrange[i_index,"fitp"]=NA
    }
    
    if(nrow(na.omit(inter.result))>0){
      inter.result=na.omit(inter.result)
      inter.data<-inter_step(inter.result,result,depth=depth)
      if(!("mu" %in% colnames(inter.data))) {
        mu.inter.pick=inter.predict.all[inter.predict.all$name==sample_name &inter.predict.all$pt==0.01,]$mu1
        
        inter.select.row<-inter.result[inter.result$mu==mu.inter.pick,]
        inter.select.row<-inter.select.row[sample(nrow(inter.select.row), 1),]
        p<-inter.select.row$p
        mu.s.update<-inter.data[inter.data$p_value==p,]
        mu.select2<-mu.s.update[sample(nrow(mu.s.update), 1),]$minmu
        
        df<-Get.S(result, p, beta=0.9, depth=120, p_thre=1e-6)
        df$p_value <- p
        
        data<-pick_s(df)
        data.rearrange[i_index,"intermu"]=mu.select2
        data.rearrange[i_index,"intermu_candidate"]=mu.inter.pick
        data.rearrange[i_index,"intertrust"]=(2*abs(mu.select2-mu.inter.pick))/(mu.inter.pick+mu.select2)
        #data.rearrange[i_index,"intermu"]<-ifelse( data.rearrange[i_index,"intertrust"]<=1,mu.select2,(mu.inter.pick+mu.select2)/2)
        data.rearrange[i_index,"intercell"]=inter.select.row$cell.div
        data.rearrange[i_index,"inters"]=data$s
        data.rearrange[i_index,"interp"]=p
        
        
      }else{
        #inter.data$cell.div
        data.rearrange[i_index,"intermu"]=inter.data$mu
        data.rearrange[i_index,"intermu_candidate"]=1
        data.rearrange[i_index,"intertrust"]=0
        
        data.rearrange[i_index,"intercell"]=inter.data$cell.div
        data.rearrange[i_index,"inters"]=inter.data$s
        data.rearrange[i_index,"interp"]=inter.data$p
      }
    }else{
      data.rearrange[i_index,"intermu"]=NA
      data.rearrange[i_index,"intermu_candidate"]=NA
      data.rearrange[i_index,"intertrust"]=NA
      data.rearrange[i_index,"intercell"]=NA
      data.rearrange[i_index,"inters"]=NA
      data.rearrange[i_index,"interp"]=NA
    }
    
    data.rearrange[i_index,"backp"]=bac.upper$p1
    
    
  

  saveRDS(data.rearrange,paste0(savepath,sample_name,'/all.sample.all.guess.depth',dpcode,'.rds'))
  
}


Rbest_classify<-function(sample_all,savepath){
  data<-data.frame()
  for(i_index in 1:length(sample_all)){
    
    
   
    sample_name=sample_all[i_index]
    print(sample_name)
   
    wgs.ex<-readRDS(paste0(savepath,sample_name,"/All_",sample_name,".rds"))
    
    result<-wgs.ex$result
    result$count<-1
    meanvaf<-aggregate(result$vaf.1, list(result$colors), mean)
    maxvaf<-aggregate(result$vaf.1, list(result$colors), max)
    minvaf<-aggregate(result$vaf.1, list(result$colors), min)
    
    sumvaf<-aggregate(result$count, list(result$colors), sum)
    mago.result<-data.frame(max=maxvaf$x,min=minvaf$x,vaf=meanvaf$x,sum=sumvaf$x)
    mago.result$count<-1
    mago.result<-mago.result[order(mago.result$max,decreasing = T),]
    
    # print the updated data
    mago.result<-mago.result[mago.result$count>0,]
    main.vaf=result[result$colors %in% maxvaf[maxvaf$x ==max(maxvaf$x),1],]$vaf.1
    possibleError <- tryCatch({
      m<-automixfit(main.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf)/5)), Niter.max=10000)
      
    },error=function(e){
      e
    })
    if(!inherits(possibleError, "error")){
    m<-automixfit(main.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf)/5)), Niter.max=10000)
    a=m["a",]
    b=m["b",]
    mean.a.b=a/(a+b)
    data[i_index,"samplename"]=sample_name
    data[i_index,"minvaf"]=min(mean.a.b)
    data[i_index,"maxvaf"]=max(mean.a.b)
    #must_have_one
    
    mean.a.b.sub<-mean.a.b[abs(mean.a.b-0.5)>0.005]
    #
    #
    data[i_index,"len_2"]=length(mean.a.b)
    data[i_index,"len"]=length(mean.a.b.sub)+1
    }else{
      data[i_index,"samplename"]=sample_name
      data[i_index,"minvaf"]=mean(main.vaf)
      data[i_index,"maxvaf"]=mean(main.vaf)
       #
      data[i_index,"len_2"]=1
      data[i_index,"len"]=1
    }
    
    ##adjust to 0.5
    current.mean<-mean(main.vaf)
    adjustment<-0.5-current.mean
    main.vaf<-main.vaf+adjustment
    possibleError <- tryCatch({
      m<-automixfit(main.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf)/5)), Niter.max=10000)
      
      
    },error=function(e){
      e
    })
    if(!inherits(possibleError, "error")){
    #m<-automixfit(main.vaf, type = "beta",Nc =1:10,thresh=0,k = 6,Ninit=min(50,round(length(main.vaf)/5)), Niter.max=10000)
    a=m["a",]
    b=m["b",]
    mean.a.b=a/(a+b)
    data[i_index,"minvaf_adj"]=min(mean.a.b)
    data[i_index,"maxvaf_adj"]=max(mean.a.b)
    data[i_index,"len_adj"]=length(mean.a.b)
    mean.a.b.sub<-mean.a.b[abs(mean.a.b-0.5)>0.005]
    data[i_index,"len_adj_2"]=length(mean.a.b)
    }else{
      data[i_index,"minvaf_adj"]=mean(main.vaf)
      data[i_index,"maxvaf_adj"]=mean(main.vaf)
      data[i_index,"len_adj"]=1
      data[i_index,"len_adj_2"]=1
    }
    
    
    
    
  }
  saveRDS(data,paste0(savepath,'/Rbest.classify.rds'))
}


Post_process<-function(sid,savepath,depth,beta,dpcode){
  
  sample_name<-sid
 
  data.rearrange<-readRDS(paste0(savepath,sample_name,'/all.sample.all.guess.depth',dpcode,'.rds'))
  #row <- data.rearrange[1, ]

# Check if all values are NA except for 'backp'
  #) {
  #data.rearrange$label <- "bac"
# All columns except 'backp' are NA, and 'backp' is not NA
  #}else{
  for(i_index in 1:nrow(data.rearrange)){
    
    
    #print(lower[i_index,"id"])
    sample_name=data.rearrange[i_index,"name"]
    #wgs.ex<-readRDS(paste0("/Users/hchen/R_project/Model_Caner/mob_simulation_magos/result/",sample_name,"/All_",sample_name,".rds"))
    wgs.ex<-readRDS(paste0(savepath,sample_name,"/All_",sample_name,".rds"))
    result<-wgs.ex$result
    
    result$count<-1
    meanvaf<-aggregate(result$vaf.1, list(result$colors), mean)
    maxvaf<-aggregate(result$vaf.1, list(result$colors), max)
    minvaf<-aggregate(result$vaf.1, list(result$colors), min)
    
    sumvaf<-aggregate(result$count, list(result$colors), sum)
    mago.result<-data.frame(max=maxvaf$x,min=minvaf$x,vaf=meanvaf$x,sum=sumvaf$x)
    mago.result$count<-1
    mago.result<-mago.result[order(mago.result$max,decreasing = T),]
    
    # print the updated data
    mago.result<-mago.result[mago.result$count>0,]
    #data.rearrange[i_index,"clonalmean"]=mago.result[1,"vaf"]
    #data.rearrange[i_index,"secondmean"]=mago.result[2,"vaf"]
    main.vaf=result[result$colors %in% maxvaf[maxvaf$x ==max(maxvaf$x),1],]$vaf.1
    second_highest_value <-mago.result[2,"max"]
    
    # Update main.vaf to consider rows corresponding to the second highest value
    second.vaf <- result[result$colors %in% maxvaf[maxvaf$x == second_highest_value, 1],]$vaf.1
    vaf.all<-c(main.vaf,second.vaf)
    m<-automixfit(vaf.all, type = "beta",Nc =2:10,thresh=0,k = 6,Niter.max=10000)
    #m<-automixfit(vaf.all, type = "beta",Nc =1:10,thresh=0,k = 6, Ninit=min(50,round(length(main.vaf)/5)),Niter.max=10000)
    a=m["a",]
    b=m["b",]
    #m<-automixfit(second.vaf, type = "beta",Nc =1:10,thresh=0,k = 6, Ninit=min(50,round(length(second.vaf)/5)),Niter.max=10000)
    #a1=m["a",]
    #b1=m["b",]
    #mean.a.b.new=c(a/(a+b),a1/(a1+b1))
    mean.a.b.new=a/(a+b)
    
    close_05_vaf=mean.a.b.new[which.min(abs(mean.a.b.new - 0.5))]
    if(length(mean.a.b.new[mean.a.b.new>min(main.vaf) & mean.a.b.new<close_05_vaf])==0){
      insert.vaf<-fit_update(main.vaf,depth,beta)
      mean.a.b<-c(mean.a.b.new[mean.a.b.new<close_05_vaf],insert.vaf)
      data.rearrange[i_index,"clonallen"]=1
    }else{
      mean.a.b<-mean.a.b.new[mean.a.b.new<close_05_vaf]
      data.rearrange[i_index,"clonallen"]=2
    }
    
    
    
    #data.rearrange[i_index,"name"]
    data.rearrange[i_index,"minvaf"]=min(mean.a.b)
    data.rearrange[i_index,"maxvaf"]=max(mean.a.b)
    data.rearrange[i_index,"closevaf"]=mean.a.b[which.min(abs(mean.a.b - 0.5))]
    # For fit1
    if(!is.na(data.rearrange[i_index,"fitp"])) {
      fit1 = data.rearrange[i_index,"fitp"]/2 + (1-data.rearrange[i_index,"fitp"])/(2 * exp(log(2) * beta * 1))
      data.rearrange[i_index,"closefitvaf"]=mean.a.b[which.min(abs(mean.a.b - fit1))]
      data.rearrange[i_index,"fitdiff"]=abs(mean.a.b[which.min(abs(mean.a.b - fit1))] - fit1)
      fit2= data.rearrange[i_index,"fitp"]/2
      data.rearrange[i_index,"closefitvaf2"]=mean.a.b[which.min(abs(mean.a.b - fit2))]
      data.rearrange[i_index,"fitdiff2"]=abs(mean.a.b[which.min(abs(mean.a.b - fit2))] - fit2)
      fit3= data.rearrange[i_index,"fitp"]/(2 * exp(log(2) * beta * 1))
      fit4= data.rearrange[i_index,"fitp"]/(2 * exp(log(2) * beta * (1+data.rearrange[i_index,"fits"])))
      if(fit4>min(second.vaf)){
        data.rearrange[i_index,"closefitvaf3"]=mean.a.b[which.min(abs(mean.a.b - fit4))]
        data.rearrange[i_index,"fitdiff3"]=abs(mean.a.b[which.min(abs(mean.a.b - fit4))] - fit4)
        
      }else{
        data.rearrange[i_index,"closefitvaf3"]=mean.a.b[which.min(abs(mean.a.b - fit3))]
        data.rearrange[i_index,"fitdiff3"]=abs(mean.a.b[which.min(abs(mean.a.b - fit3))] - fit3)
        
      }
      
      inita=depth*c(0.5,fit1,fit2)
      initb=depth-inita
      probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
      df<-data.frame(
        prob=probs,
        vaf=main.vaf
      )
      
      # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
      df<-beta_reassign(df)   
      data.rearrange[i_index,"fit_len_diff"]=data.rearrange[i_index,"fitmu"]*data.rearrange[i_index,"fitcell"]
      data.rearrange[i_index,"fit_len_ratio"]=(data.rearrange[i_index,"fitmu"]*data.rearrange[i_index,"fitcell"])/nrow(df[df$cluster>1,])
      
    } else {
      data.rearrange[i_index,"closefitvaf"]=NA
      data.rearrange[i_index,"fitdiff"]=NA
      data.rearrange[i_index,"closefitvaf2"]=NA
      data.rearrange[i_index,"fitdiff2"]=NA
      data.rearrange[i_index,"closefitvaf3"]=NA
      data.rearrange[i_index,"fitdiff3"]=NA
      data.rearrange[i_index,"fit_len_diff"]=NA
      data.rearrange[i_index,"fit_len_ratio"]=NA
    }
    
    # For inter1
    if(!is.na(data.rearrange[i_index,"interp"])) {
      inter1 = data.rearrange[i_index,"interp"]/2 + (1-data.rearrange[i_index,"interp"])/(2 * exp(log(2) * beta * 1))
      data.rearrange[i_index,"closeintervaf"]=mean.a.b[which.min(abs(mean.a.b - inter1))]
      data.rearrange[i_index,"interdiff"]=abs(mean.a.b[which.min(abs(mean.a.b - inter1))] - inter1)
      inter2=data.rearrange[i_index,"interp"]/2
      data.rearrange[i_index,"closeintervaf2"]=mean.a.b[which.min(abs(mean.a.b - inter2))]
      data.rearrange[i_index,"interdiff2"]=abs(mean.a.b[which.min(abs(mean.a.b - inter2))] - inter2)
      fit3= data.rearrange[i_index,"interp"]/(2 * exp(log(2) * beta * (1+data.rearrange[i_index,"inters"])))
      data.rearrange[i_index,"closeintervaf3"]=mean.a.b[which.min(abs(mean.a.b - fit3))]
      data.rearrange[i_index,"interdiff3"]=abs(mean.a.b[which.min(abs(mean.a.b - fit3))] - fit3)
      
      inita=depth*c(0.5,inter1,inter2)
      initb=depth-inita
      probs <- sapply(1:length(inita), function(i) dbeta(main.vaf, inita[i], initb[i]))
      df<-data.frame(
        prob=probs,
        vaf=main.vaf
      )
      
      # Apply the get_cluster function to each row of probs to get the cluster number for each maximum probability
      df<-beta_reassign(df)   
      
      #intermu<-ifelse(data.rearrange[i_index,"intertrust"]<=1,data.rearrange[i_index,"intermu"],(data.rearrange[i_index,"intermu"]+data.rearrange[i_index,"intermu_candidate"])/2)
      #intermu<-ifelse(data.rearrange[i_index,"intertrust"]<=1,data.rearrange[i_index,"intermu"],data.rearrange[i_index,"intermu_candidate"])
      intermu<-data.rearrange[i_index,"intermu"]
      #intermu=ifelse(data.rearrange[i_index,"intermu_candidate"]!=1 & data.rearrange[i_index,"intermu"]/data.rearrange[i_index,"intermu_candidate"]>10,0.5*data.rearrange[i_index,"intermu_candidate"]+0.5*data.rearrange[i_index,"intermu"],data.rearrange[i_index,"intermu"])
      #intermu=ifelse(data.rearrange[i_index,"intermu_candidate"]!=1 & data.rearrange[i_index,"intermu"]/data.rearrange[i_index,"intermu_candidate"]>10,data.rearrange[i_index,"intermu_candidate"],data.rearrange[i_index,"intermu"])
      data.rearrange[i_index,"inter_len_diff"]=intermu*data.rearrange[i_index,"intercell"]
      
      #data.rearrange[i_index,"inter_len_diff"]=abs(nrow(df[df$cluster>1,])-intermu*data.rearrange[i_index,"intercell"])
      data.rearrange[i_index,"inter_len_ratio"]=(intermu*data.rearrange[i_index,"intercell"])/nrow(df[df$cluster>1,])
      
      
    } else {
      #inter1 = NA # or any other fallback value or action
      data.rearrange[i_index,"closeintervaf"]=NA
      data.rearrange[i_index,"interdiff"]=NA
      data.rearrange[i_index,"closeintervaf2"]=NA
      data.rearrange[i_index,"interdiff2"]=NA
      data.rearrange[i_index,"closeintervaf3"]=NA
      data.rearrange[i_index,"interdiff3"]=NA
      data.rearrange[i_index,"inter_len_diff"]=NA
      data.rearrange[i_index,"inter_len_ratio"]=NA
      
    }
    
    # For bac1
    if(!is.na(data.rearrange[i_index,"backp"])) {
      bac1 = data.rearrange[i_index,"backp"]/2 + (1-data.rearrange[i_index,"backp"])/(2 * exp(log(2) * beta * 1))
      data.rearrange[i_index,"closebacvaf"]=mean.a.b[which.min(abs(mean.a.b - bac1))]
      data.rearrange[i_index,"bacdiff"]=abs(mean.a.b[which.min(abs(mean.a.b - bac1))] - bac1)
      bac2 = (1-data.rearrange[i_index,"backp"])/(2 * exp(log(2) * beta * 1))
      data.rearrange[i_index,"closebacvaf2"]=mean.a.b[which.min(abs(mean.a.b - bac2))]
      data.rearrange[i_index,"bacdiff2"]=abs(mean.a.b[which.min(abs(mean.a.b - bac2))] - bac2)
      
    } else {
      data.rearrange[i_index,"closebacvaf"]=NA # or any other fallback value or action
      data.rearrange[i_index,"bacdiff"]=NA
      data.rearrange[i_index,"closebacvaf2"]=NA
      data.rearrange[i_index,"bacdiff2"]=NA
    }
    
    
    
    
    
    
    
  }
  
  saveRDS(data.rearrange,paste0(savepath,sample_name,'/all.sample.all.guess.start.choose.depth',dpcode,'.rds'))
  
  
  
  data<-readRDS(paste0(savepath,'/Rbest.classify.rds'))
  inter_sample=data[data$len_adj>1,]$samplename
  missample=data[data$len_adj==1,]$samplename
  for(i_index in 1:nrow(data.rearrange)){
    
    
    
    
    
    fit_sum <- sum(data.rearrange[i_index, c("fitdiff2", "fitdiff2")], na.rm = TRUE)
    if(any(is.na(data.rearrange[i_index, c("fitdiff2","fitdiff2")]))) fit_sum <- Inf
    
    inter_sum <- sum(data.rearrange[i_index, c("interdiff2", "interdiff2")], na.rm = TRUE)
    #inter_sum <- min(data.rearrange[i_index, "interdiff"],data.rearrange[i_index, "interdiff2"])
    if(any(is.na(data.rearrange[i_index, c("interdiff", "interdiff2")]))) inter_sum <- Inf
    
    bac_sum <- sum(data.rearrange[i_index, c("bacdiff2", "bacdiff2")], na.rm = TRUE)
    #bac_sum <- min(data.rearrange[i_index, "bacdiff"], data.rearrange[i_index, "bacdiff2"])
    if(any(is.na(data.rearrange[i_index, c("bacdiff", "bacdiff2")]))) bac_sum <- Inf
    
    # Determine the label for the row based on the smallest sum
    #min_sum <- min(inter_sum, bac_sum)
    #min_sum2<-min(inter_sum, bac_sum)
    min_2<-min(inter_sum,bac_sum)
    data.rearrange$label[i_index] <- NA
    if(data.rearrange[i_index,"name"] %in% inter_sample){
      
      data.rearrange$label[i_index] <- "inter"
      
    }
    
    if(data.rearrange[i_index,"name"] %in% missample){
      #data.rearrange$label[i_index]<-NA
      #bac_sum1 == min_2 & data.rearrange[i_index,"backp"]<0.4
      if(inter_sum >bac_sum) {
        data.rearrange$label[i_index] <- "bac"
      } else{
        data.rearrange$label[i_index] <- "inter"
      }
    }
    
    
    
    
    
  }
  
  #data.rearrange[i_index,"clonallen"]
  return(data.rearrange)
}


Inter_post_process<-function(sid,data.rearrange,savepath,depth,beta,dpcode){
  data=data.rearrange[data.rearrange$label=="inter",]
  if(nrow(data)>0){
  times=2
  
  data$intermu <- with(data, ifelse(intermu_candidate != 1 & 
                                      (intermu/intermu_candidate > times | intermu/intermu_candidate < 1/times), 
                                    0.5*intermu_candidate + 0.5*intermu, 
                                    intermu))
  data$fitmu <- with(data, ifelse(fitmu_candidate != 1 & 
                                    (fitmu/fitmu_candidate > times | fitmu/fitmu_candidate < 1/times), 
                                  0.5*fitmu_candidate + 0.5*fitmu, 
                                  fitmu))
  
  
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
  Rbest_data<-readRDS(paste0(savepath,'/Rbest.classify.rds'))
  missample=Rbest_data[Rbest_data$len_adj==1,]$samplename
  data$mupick<-ifelse(data$name %in% missample,data$intermu,data$mupick)

  
  
  data$picks<-ifelse(data$mupick==data$intermu,data$inters,data$fits)  
  #data<-merge(data,Suppose_para,by.x=c("name"),by.y=c("id"),all.x=T)   
  
  data$pickt1<-ifelse(data$mupick==data$intermu,data$inter_len_diff/data$mupick,data$fit_len_diff/data$mupick)  
  #data$realt1<-log(data$t1)/(log(2)*beta)
  
  saveRDS(data,paste0(savepath,sid,'/final.result.depth',dpcode,'.rds'))
  }else{
    
    saveRDS(data,paste0(savepath,sid,'/final.result.depth',dpcode,'.rds'))
  } 
}
  environment()
})
