## S4 classes for a machine learning method to identify gold standards 
## of plant complexes by integrating the information of CORUM and 
## CF-MSdata from the size exclusion chromatography (SEC) experiment.
##
## object oriented programming in R:
##

#####################################################################
##     1. Define function interfaces in S4
#####################################################################


find_super_sub_nonduplicate_set <- function(x, evaluation_data, cutoff, decision){
  # This function is used to 
  #   1. find supper and sub sets from the given sets.
  #   2. find non-duplicate upper/sub sets that satisfying some criteria.
  # The input, x, should be a list containing vectors (sets) with distinct names.
  if (typeof(x) != "list"){
    stop("argument type error: x should be a list")
  }else if (length(levels(factor(names(x)))) < length(x)){
    stop(cat("duplicate set names error: there are dulplicate names in the list"))
  }
  
  # super set is a set that cannot be proper subset of any other sets
  super_names <- c()
  for (c1 in names(x)){
    set_c1 <- x[[c1]]
    is.subset_set_c1 <- F
    for (c2 in names(x)) {
      set_c2 <- x[[c2]]
      if(length(set_c2)>length(set_c1) & all(set_c1 %in% set_c2)){
        is.subset_set_c1 <- T
        break
      }
    }
    if(is.subset_set_c1 == F){
      super_names <- c(super_names, c1)
    }
  }
  
  # find the equivalent sets (duplicates)
  count <- 0
  duplicate_of_set <- list()
  names_super <- super_names
  while(length(names_super)>0){
    count <- count + 1
    name_1 <- names_super[1]
    temp_names <- c(name_1)
    for (s in names_super){
      if (! s==name_1){
        if (all(c(sum(x[[name_1]] %in% x[[s]])==length(x[[name_1]]),
                  length(x[[name_1]]==length(x[[s]]))))){
          temp_names <- c(temp_names, s)
          
        }
      }
    }
    names_super <- names_super[! names_super %in% temp_names]
    duplicate_of_set[[count]] <- temp_names
  }
  
  good_set_selected_nonduplicate <- c()
  for (i in 1:length(duplicate_of_set)){
    c_t <- c()
    v_t <- c()
    for (c in duplicate_of_set[[i]]){
      v <- evaluation_data[[c]]
      if (decision=="max" & v >= cutoff){
        c_t <- c(c_t, c)
        v_t <- c(v_t, v)
      }else if (decision=="min" & v <= cutoff){
        c_t <- c(c_t, c)
        v_t <- c(v_t, v)
      }else{
        c_t <- c(c_t, c)
        v_t <- c(v_t, v)
      }
    }
    if (decision=="max"){
      good_set_selected_nonduplicate <- c(good_set_selected_nonduplicate, c_t[which.max(v_t)])
    }else if(decision=="min"){
      good_set_selected_nonduplicate <- c(good_set_selected_nonduplicate, c_t[which.min(v_t)])
    }else{
      good_set_selected_nonduplicate <- c(good_set_selected_nonduplicate, c_t[length(c_t)])
    } 
    
  }
  sub_names <- names(x)[!(names(x) %in% super_names)]
  
  return(list(superset_names=super_names, subset_names=sub_names, selected_nonduplicate_set_names=good_set_selected_nonduplicate))
}

Mcalc_calculation <- function(complex_name, Mcomplex_to_Plant_mapping, ortholog_set, Mmnono_data, stoichiometry_data){
  # This function is used to calculate the sum of mono mass of proteins in a given set
  #   Mcalc = sum(average(mono masses of orthologs for M_protein_i) * stoichiometry of M_protein_i)
  c <- complex_name
  Mcalc_set_u <- list()
  for (u in names(Mcomplex_to_Plant_mapping)){
    Mcalc_subset_u <- c()
    for (orth in ortholog_set){
      if (orth %in% Mcomplex_to_Plant_mapping[[u]]){
        Mcalc_subset_u <- c(Mcalc_subset_u, orth)
      }
    }
    if (length(Mcalc_subset_u)>1){
      Mcalc_set_u[[u]] <- Mcalc_subset_u
    }
  }
  ortholog_set_left <-ortholog_set[!(ortholog_set %in% unlist(Mcalc_set_u))]
  Mcalc_set_p <- list()
  for (orth in ortholog_set_left){
    Mcalc_subset_p <- c()
    for (u in names(Mcomplex_to_Plant_mapping)){
      if (orth %in% Mcomplex_to_Plant_mapping[[u]]){
        Mcalc_subset_p <- c(Mcalc_subset_p, u)
      }
    }
    Mcalc_set_p[[orth]] <- Mcalc_subset_p
  }
  
  Mcalc <-c()
  for (u in names(Mcalc_set_u)){
    g <- c()
    for (i in Mcalc_set_u[[u]]){
      g <- c(g, i)
    }
    if (c %in% names(stoichiometry_data)){
      if(u %in% names(stoichiometry_data[[c]])){
        sto_num <- stoichiometry_data[[c]][[u]][['stoichiometry']]
      }else{
        sto_num <- 1
      }
      
    }else{
      sto_num <- 1
    }
    mc <- c()
    for (j in g){
      mono_j <- Mmnono_data[[j]]
      mc <- c(mc, mono_j*sto_num)
    }
    Mcalc <- c(Mcalc, mean(mc))
  }
  for (orth in names(Mcalc_set_p)){
    mc <- c()
    for (u in Mcalc_set_p[[orth]]){
      if (c %in% names(stoichiometry_data)){
        if(u %in% names(stoichiometry_data[[c]])){
          sto_num <- stoichiometry_data[[c]][[u]][['stoichiometry']]
        }else{
          sto_num <- 1
        }
        
      }else{
        sto_num <- 1
      }
      mc <- c(mc, Mmnono_data[[orth]]*sto_num)
    }
    Mcalc <- c(Mcalc, mean(mc))
  }
  return(list(Mcalc_set=list(Mcalc_set_u, Mcalc_set_p), Mcalc=sum(Mcalc)))
}


small_world_analysis <- function(names_CORUM_complexes, ortholog_map, prot_peak_map, prots_for_small_world, 
                                 Gaussian_peak_fraction_1, Gaussian_peak_fraction_2, 
                                 Matrix_SEC_B1, Matrix_SEC_B2, parameters_SOM, parameters_APC,
                                 cross_correlated_r){
  Similarity_Matrix_Wcc_Euclidean <- function(prots, M_peakprofile_1, M_peakprofile_2, GPF1, GPF2, weight_v, l){
    x <- cbind(M_peakprofile_1[prots,], M_peakprofile_2[prots,])
    d <- length(x[1,])/2
    w <- diag(d)
    for (i in 1:d){
      for (j in 1:d){
        if (abs(i-j)<l & abs(i-j)>0){
          w[i,j] <- 1 - abs(i-j)/l
        }
      }
    }
    
    w_cross_corr_1 <- x[,1:d] %*% w %*% t(x[,1:d])
    s1 <- w_cross_corr_1
    w_cross_corr_2 <- x[,(d+1):length(x[1,])] %*% w %*% t(x[,(d+1):length(x[1,])])
    s2 <- w_cross_corr_2
    for (i in 1:nrow(w_cross_corr_1)){
      for (j in 1:nrow(w_cross_corr_1)){
        s1[i,j] <- w_cross_corr_1[i,j]/(sqrt(w_cross_corr_1[i,i])*sqrt(w_cross_corr_1[j,j]))
        s2[i,j] <- w_cross_corr_2[i,j]/(sqrt(w_cross_corr_2[i,i])*sqrt(w_cross_corr_2[j,j]))
      }
    }
    s_sim <- weight_v[1]*s1 + weight_v[2]*s2
    
    m <- matrix(rep(GPF1[1,1],length(GPF1)), nrow=length(GPF1))
    for (i in 1:length(GPF1)){
      mt <- matrix(rep(GPF1[i,1],length(GPF1)), nrow=length(GPF1))
      m <- cbind(m, abs(GPF1-mt))
    }
    m <- m[,-1]
    colnames(m) <- rownames(m)
    
    s_GPF1 <- 1-m
    
    m <- matrix(rep(GPF2[1,1],length(GPF2)), nrow=length(GPF2))
    #abs(GF1-m)
    for (i in 1:length(GPF2)){
      mt <- matrix(rep(GPF2[i,1],length(GPF2)), nrow=length(GPF2))
      m <- cbind(m, abs(GPF2-mt))
    }
    m <- m[,-1]
    colnames(m) <- rownames(m)
    
    s_GPF2 <- 1-m
    
    s_GPF <- weight_v[3]*s_GPF1 + weight_v[4]*s_GPF2
    
    s_both <- s_sim + s_GPF
    
    return(list(sim_Wcc=s_sim, sim_GPF=s_GPF, s_both=s_both))
    
  }
  weighted_sim <- function(x,l){
    d <- length(x[1,])/2
    w <- diag(d)
    for (i in 1:d){
      for (j in 1:d){
        if (abs(i-j)<l & abs(i-j)>0){
          w[i,j] <- 1 - abs(i-j)/l
        }
      }
    }
    
    w_cross_corr_1 <- x[,1:d] %*% w %*% t(x[,1:d])
    s1 <- w_cross_corr_1
    w_cross_corr_2 <- x[,(d+1):length(x[1,])] %*% w %*% t(x[,(d+1):length(x[1,])])
    s2 <- w_cross_corr_2
    for (i in 1:nrow(w_cross_corr_1)){
      for (j in 1:nrow(w_cross_corr_1)){
        s1[i,j] <- w_cross_corr_1[i,j]/(sqrt(w_cross_corr_1[i,i])*sqrt(w_cross_corr_1[j,j]))
        s2[i,j] <- w_cross_corr_2[i,j]/(sqrt(w_cross_corr_2[i,i])*sqrt(w_cross_corr_2[j,j]))
      }
    }
    s <- 0.5*s1 + 0.5*s2
    return(s)
  }
  
  data_m <- cbind(Matrix_SEC_B1[prots_for_small_world,],Matrix_SEC_B2[prots_for_small_world,])
  sim_matrix <- weighted_sim(x=data_m, l=cross_correlated_r) 
  rownames(sim_matrix) <- rownames(Matrix_SEC_B1)
  colnames(sim_matrix) <- rownames(Matrix_SEC_B1)
  sim_space <- sim_matrix[upper.tri(sim_matrix, diag = FALSE)]
  GPF_std_1 <- Gaussian_peak_fraction_1/max(Gaussian_peak_fraction_1)
  GPF_std_2 <- Gaussian_peak_fraction_2/max(Gaussian_peak_fraction_2)
  GPF_mean <- (Gaussian_peak_fraction_1+Gaussian_peak_fraction_2)/2
  GPF_mean_std <- GPF_mean/max(GPF_mean)
  G_1 <- GPF_mean_std[,1][1]
  G_matrix <- matrix(abs(rep(G_1, nrow(GPF_mean_std))-GPF_mean_std), nrow=nrow(GPF_mean_std), dimnames=list(rownames(GPF_mean_std), names(G_1)))
  for (i in 2:nrow(GPF_mean_std)){
    G_i <- GPF_mean_std[,1][i]
    M_temp <- matrix(abs(rep(G_i, nrow(GPF_mean_std))-GPF_mean_std), nrow=nrow(GPF_mean_std), dimnames=list(rownames(GPF_mean_std), names(G_i)))
    G_matrix <- cbind(G_matrix, M_temp)
  }
  G_space <- G_matrix[upper.tri(G_matrix, diag = FALSE)]
  
  subcomplex_prediction <- list()
  pdf(paste(paste("Small_World_Analysis_SOMandAPC_Output", Sys.Date(), sep='_'), "pdf", sep='.'))
  count_t <- 0
  for (c in names_CORUM_complexes){
    count_t <- count_t + 1
    print(paste(count_t, paste('Processing complex', c, sep='-'), sep=':'))
    find_name <- function(x, orth_map){
      names_need <- c()
      for (p in levels(factor(x))){
        for (u in names(orth_map)){
          if (p %in% orth_map[[u]]){
            names_need <- c(names_need, u)
          }
        }
      }
      return(levels(factor(names_need)))
    }
    
    ortholog_set <- unname(unlist(ortholog_map[[c]]))
    prot_id_c <- levels(factor(ortholog_set[!is.na(ortholog_set)]))
    u_names_in_c <- find_name(x=prot_id_c, orth_map=ortholog_map[[c]])
    prot_id_temp <- c()
    for (p in prot_id_c){
      if (! p%in% names(prot_peak_map)){
        prot_id_temp <- c(prot_id_temp, p)
      }else{
        prot_id_temp <- c(prot_id_temp, prot_peak_map[[p]])
      }
      
    }
    prot_id_c <- prot_id_temp
    if (length(levels(factor(u_names_in_c)))<2){
      next # only one mammalian id has ortholog
    }
    if (sum(prot_id_c %in% prots_for_small_world)>1){
      prot_id_rep_c <- prot_id_c[prot_id_c %in% prots_for_small_world]
    }else{
      next
    }
    u_names_in_rep_c <- find_name(x=prot_id_rep_c, orth_map=ortholog_map[[c]])
    if (length(levels(factor(u_names_in_rep_c)))<2){
      next # only one mammalian id has detected ortholog
    }
    
    if (length(prot_id_rep_c)==2){
      subcomplex_prediction[[c]][['only_2_prots_before_clustering']] <- prot_id_rep_c
      next
    }
    
    Filtered_Data_raw_SEC <- list(Filtered_Raw_SEC1=Matrix_SEC_B1[prot_id_rep_c,],
                                  Filtered_Raw_SEC2=Matrix_SEC_B2[prot_id_rep_c,],
                                  Filtered_Gaussian_peak_fraction_1=GPF_std_1[prot_id_rep_c,],
                                  Filtered_Gaussian_peak_fraction_2=GPF_std_2[prot_id_rep_c,])
    
    sample_G_matrix <- G_matrix[prot_id_rep_c, prot_id_rep_c]
    G_set <- sample_G_matrix[upper.tri(sample_G_matrix, diag = FALSE)]
    if (sum(G_set<=quantile(G_space, 0.05))==length(G_set)){
      subcomplex_prediction[[c]][['all_prots_before_clustering_GPFdist']] <- prot_id_rep_c
      next
    }
    
    sample_sim_matrix <- sim_matrix[prot_id_rep_c, prot_id_rep_c]
    sim_set <- sample_sim_matrix[upper.tri(sample_sim_matrix, diag = FALSE)]
    if (sum(sim_set>=quantile(sim_space, 0.95))==length(sim_set)){
      subcomplex_prediction[[c]][['all_prots_before_clustering_sim']] <- prot_id_rep_c
      next
    }
    
    num_p <- length(prot_id_rep_c)
    
    if (num_p<6){
      x_limit <- 1
      y_limit <- num_p
    }
    
    if (num_p>=6 & num_p<10){
      x_limit <- 2
      y_limit <- 3
    }
    
    if (num_p>=10 & num_p<15){
      x_limit <- 3
      y_limit <- 3
    }
    
    if (num_p>=15 & num_p<20){
      x_limit <- 3
      y_limit <- 4
    }
    if (num_p>=20 & num_p<35){
      x_limit <- 4
      y_limit <- 4
    }
    
    if (num_p>=35 & num_p<60){
      x_limit <- 4
      y_limit <- 5
    }
    if (num_p>=60 & num_p<80){
      x_limit <- 5
      y_limit <- 5
    }
    
    if (num_p>=80 & num_p<100){
      x_limit <- 5
      y_limit <- 6
    }
    if (num_p>=100){
      x_limit <- 5
      y_limit <- 10
    }
    
    set.seed(123)
    
    system.time(SOM_wcc_eucd <- supersom(Filtered_Data_raw_SEC, 
                                         grid=somgrid(xdim = x_limit, ydim = y_limit, topo = parameters_SOM$topo_shape, 
                                                      neighbourhood.fct = "gaussian", toroidal = parameters_SOM$toroidal_choice),  
                                         dist.fcts = parameters_SOM$dist_measure, user.weights = parameters_SOM$user_weight,
                                         rlen=parameters_SOM$num_of_multiple_of_sample_size,
                                         normalizeDataLayers = FALSE))
    
    label_of_sample <- SOM_wcc_eucd$unit.classif
    cluster_ID <- levels(factor(label_of_sample))
    
    sample_in_each_unit <- list()
    for (i in c(1:(x_limit*y_limit))){
      postion_index <- c(label_of_sample==i)
      name <- paste0('V',as.character(i), sep="")
      if (length(label_of_sample[postion_index])>0){
        sample_in_each_unit[[name]] <- prot_id_rep_c[postion_index]
      }else{
        sample_in_each_unit[[name]] <- NULL
      }
    }
    code_som <- SOM_wcc_eucd$codes
    s <- Similarity_Matrix_Wcc_Euclidean(prots=rownames(code_som$Filtered_Raw_SEC1), 
                                         M_peakprofile_1=code_som$Filtered_Raw_SEC1,
                                         M_peakprofile_2=code_som$Filtered_Raw_SEC2, 
                                         GPF1=code_som$Filtered_Gaussian_peak_fraction_1,
                                         GPF2=code_som$Filtered_Gaussian_peak_fraction_2, 
                                         weight_v=parameters_SOM$user_weight, 
                                         l=cross_correlated_r)$s_both
    if (num_p<4){
      apc <- apcluster(s=s, q=parameters_APC$q_1, maxits=parameters_APC$maxits, convits=parameters_APC$convits, lam=parameters_APC$lam, detail=parameters_APC$detail)
    }else if (num_p>3 & num_p<6){
      apc <- apcluster(s=s, q=parameters_APC$q_2, maxits=parameters_APC$maxits, convits=parameters_APC$convits, lam=parameters_APC$lam, detail=parameters_APC$detail)
    }else{
      apc <- apcluster(s=s, q=parameters_APC$q_3, maxits=parameters_APC$maxits, convits=parameters_APC$convits, lam=parameters_APC$lam, detail=parameters_APC$detail)
    }
    
    
    if (length(apc)==0){
      print(paste('warning: no apc in', c, sep=' '))
      next
    }
    
    unit_labels <- rownames(s)
    group_units <-c()
    for (u in unit_labels){
      for (i in 1:length(apc)){
        if (u %in% names(apc[[i]])){
          group_units[u] <- i
        }
      }
    }
    
    count_sub <- 0
    count_orthoparalog <- 0
    count_single <- 0
    for (g in levels(factor(unname(group_units)))){
      #count_t <- count_t + 1
      pos <-unname(group_units) %in% g
      V <- names(group_units)[pos]
      s_g <- c()
      for (v in V){
        s_g <- c(s_g, sample_in_each_unit[[v]])
      }
      if (length(s_g)>1){
        names_in_s_g <- find_name(x=substr(s_g,1,nchar(rownames((Matrix_SEC_B1))[1])), orth_map=ortholog_map[[c]])
        if (length(levels(factor(names_in_s_g)))>1){
          count_sub <- count_sub + 1
          subcomplex_prediction[[c]][[paste('sub', count_sub, sep='_')]] <- s_g
        }else{
          count_orthoparalog <- count_orthoparalog + 1
          subcomplex_prediction[[c]][[paste('sub_orthoparalog', count_orthoparalog, sep='_')]] <- s_g
        }
      }
      if (length(s_g)==1){
        count_single <- count_single + 1
        subcomplex_prediction[[c]][[paste('singleton', count_single, sep='_')]] <- s_g
      }
    }
    
    plot(SOM_wcc_eucd, type="mapping",labels=prot_id_rep_c, bgcol = pretty_palette(length(apc))[group_units],cex=0.25,heatkeywidth=0.65, 
         main = paste("subcomplexes of", c, sep=' ')) 
    if (x_limit>1 & y_limit>1){
      add.cluster.boundaries(SOM_wcc_eucd, group_units)
    }
  }
  dev.off()
  
  sim_matrix <- Similarity_Matrix_Wcc_Euclidean(prots=prots_for_small_world, 
                                                M_peakprofile_1=Matrix_SEC_B1, 
                                                M_peakprofile_2=Matrix_SEC_B2, 
                                                GPF1=GPF_std_1, GPF2=GPF_std_2, 
                                                weight_v=parameters_SOM$user_weight, 
                                                l=cross_correlated_r)$s_both 
  subprediction <- subcomplex_prediction
  
  print('It needs some time for random sampling and P-value calculation, please wait for ...')
  num_prot_loc <- c()
  for (c in names(subprediction)){
    if(length(subprediction[[c]])>0){
      for (i in names(subprediction[[c]])){
        if (length(subprediction[[c]][[i]])>1){
          num_prot_loc <- c(num_prot_loc, paste(length(unlist(subprediction[[c]])), length(subprediction[[c]][[i]]), collapse = '-'))
        }
      }
    }
  }
  
  num_prot_loc <- levels(factor(num_prot_loc))
  
  sampling_random_com <- list()
  for (ind in num_prot_loc){
    pair_n <- as.numeric(strsplit(ind, split = ' ')[[1]])
    mean_rc <- c()
    for (i in 1:(length(prots_for_small_world)*50)){
      rc_large <- sample(prots_for_small_world, size=pair_n[1], replace=FALSE)
      rc_small <- sample(rc_large, size=pair_n[2], replace=FALSE)
      m <- sim_matrix[rc_small, rc_small]
      s_wccd <- (sum(m)- pair_n[2])/2
      
      avg <- s_wccd/(pair_n[2]*(pair_n[2]-1)/2)
      mean_rc <- c(mean_rc, avg)
    }
    sampling_random_com[[ind]] <- mean_rc
  }
  
  subcomplex_cutoff <- list()
  for (c in names(subprediction)){
    if(length(subprediction[[c]])>0){
      subcomplex_cutoff[[c]] <- list()
      for (i in names(subprediction[[c]])){
        if (length(subprediction[[c]][[i]])>1){
          num_loc <- paste(length(unlist(subprediction[[c]])), length(subprediction[[c]][[i]]), collapse = '-')
          src <- sampling_random_com[[num_loc]]
          sim_com <- sim_matrix[subprediction[[c]][[i]], subprediction[[c]][[i]]]
          s_com <- (sum(sim_com)- length(subprediction[[c]][[i]]))/2
          sim_mean <- s_com/(length(subprediction[[c]][[i]])*(length(subprediction[[c]][[i]])-1)/2)
          cut_off <- quantile(src, 0.95)
          Pvalue <- (sum(src>sim_mean)+1)/(length(src)+1)
          subcomplex_cutoff[[c]][[i]] <- c(sim_mean, cut_off, Pvalue)
        }
      }
    }
  }
  
  subcomplex_overlap <- list()
  for (c in names(subprediction)){
    if(length(subprediction[[c]])>0){
      subcomplex_overlap[[c]] <- list()
      for (i in names(subprediction[[c]])){
        if (length(subprediction[[c]][[i]])>1){
          overlap <- c()
          for (sc in names_CORUM_complexes){
            if (all(subprediction[[c]][[i]] %in% unlist(ortholog_map[[sc]])) & sc != c){
              overlap <- c(overlap, sc)
            }
          }
          if (length(overlap)>0){
            subcomplex_overlap[[c]][[i]] <- overlap
          }else{
            subcomplex_overlap[[c]][[i]] <- "Not available"
          }
          
        }
      }
    }
  }
  print('small world analysis done!')
  return(list(subcomplex_prediction=subcomplex_prediction, subcomplex_cutoff=subcomplex_cutoff, subcomplex_overlap=subcomplex_overlap))
  
}

#############################################################################
##     2. S4 classes for a machine learning method to identify gold standards
#############################################################################
setClass("Protein.Complexes.Gold.Standard.Prediction", 
         slots=c(mammalPlantSpecies="character",
                 mammalCORUM.Data="data.frame",
                 mammalPlantOrtholog.Data="data.frame",
                 mammalCORUM.StoichiometryData="list",
                 plantProtnames.Data="list",
                 plantMmono.Data="list",
                 peakProfileMatrixB1="matrix",
                 peakProfileMatrixB2="matrix",
                 reproduciblePeakFeatures.Data="list", 
                 singlePeakProts="character",
                 multiplePeakProts="character",
                 mostCommonProtIDExample="character",
                 genomeCoverageCutoff="numeric",
                 RappCutoffNonmono = "numeric",
                 RappCutoffChoosingThreshold="numeric",
                 parametersSOM="list",
                 parametersAPC="list",
                 crossCorrelationRadius="numeric",
                 exogenousClusteringResult="data.frame",
                 targetCluster="character",
                 is.withData="logical",
                 axisTextsize="numeric",
                 axisTitleize="numeric",
                 sctterPontSize="numeric",
                 Set.Manipulation.func.find.supper_sub_noduplicate="function",
                 Mcalc.computing.Func="function",
                 Small.World.Analysis.Func="function"),
         prototype=list(Set.Manipulation.func.find.supper_sub_noduplicate=find_super_sub_nonduplicate_set,
                        Small.World.Analysis.Func=small_world_analysis,
                        Mcalc.computing.Func=Mcalc_calculation,
                        is.withData=T,
                        axisTextsize=10,
                        axisTitleize=14,
                        sctterPontSize=1)
)


setGeneric("ProteinComplexGoldStandardPredictionWithinCORUMOrthologousComplexes",
           function(Protein.Complexes.Gold.Standard.Prediction){
             standardGeneric("ProteinComplexGoldStandardPredictionWithinCORUMOrthologousComplexes")
           })
setMethod("ProteinComplexGoldStandardPredictionWithinCORUMOrthologousComplexes", signature("Protein.Complexes.Gold.Standard.Prediction"),
          function(Protein.Complexes.Gold.Standard.Prediction){
            cat( "LOG OF COMPUTATION\n\n", file="LOG OF COMPUTATION.txt" )
            cat("Started\n\n")
            cat("Step 0: To load data, parameters and functions.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 0: To load data, parameters and functions.\n")
            mammalPlantSpecies <- Protein.Complexes.Gold.Standard.Prediction@mammalPlantSpecies
            mammalCORUM.Data <- Protein.Complexes.Gold.Standard.Prediction@mammalCORUM.Data
            mammalPlantOrtholog.Data <- Protein.Complexes.Gold.Standard.Prediction@mammalPlantOrtholog.Data
            mammalCORUM.StoichiometryData <- Protein.Complexes.Gold.Standard.Prediction@mammalCORUM.StoichiometryData
            plantProtnames.Data <- Protein.Complexes.Gold.Standard.Prediction@plantProtnames.Data
            plantMmono.Data <- Protein.Complexes.Gold.Standard.Prediction@plantMmono.Data
            peakProfileMatrixB1 <- Protein.Complexes.Gold.Standard.Prediction@peakProfileMatrixB1
            peakProfileMatrixB2 <- Protein.Complexes.Gold.Standard.Prediction@peakProfileMatrixB2
            reproduciblePeakFeatures.Data <- Protein.Complexes.Gold.Standard.Prediction@reproduciblePeakFeatures.Data
            singlePeakProts <- Protein.Complexes.Gold.Standard.Prediction@singlePeakProts
            multiplePeakProts <- Protein.Complexes.Gold.Standard.Prediction@multiplePeakProts
            genomeCoverageCutoff <- Protein.Complexes.Gold.Standard.Prediction@genomeCoverageCutoff
            RappCutoffNonmono <- Protein.Complexes.Gold.Standard.Prediction@RappCutoffNonmono
            RappCutoffChoosingThreshold <- Protein.Complexes.Gold.Standard.Prediction@RappCutoffChoosingThreshold
            is.withData <- Protein.Complexes.Gold.Standard.Prediction@is.withData
            axisTextsize <- Protein.Complexes.Gold.Standard.Prediction@axisTextsize
            axisTitleize <- Protein.Complexes.Gold.Standard.Prediction@axisTextsize
            sctterPontSize <- Protein.Complexes.Gold.Standard.Prediction@sctterPontSize
            M_name <- mammalPlantSpecies[1]
            P_name <- mammalPlantSpecies[2]
            mostCommonProtIDExample=Protein.Complexes.Gold.Standard.Prediction@mostCommonProtIDExample
            ncharProtID <- nchar(mostCommonProtIDExample)
            parameters_SOM <- Protein.Complexes.Gold.Standard.Prediction@parametersSOM
            parameters_APC <- Protein.Complexes.Gold.Standard.Prediction@parametersAPC
            cross_correlated_r <- Protein.Complexes.Gold.Standard.Prediction@crossCorrelationRadius
            exogenousClusteringResult <- Protein.Complexes.Gold.Standard.Prediction@exogenousClusteringResult
            targetCluster <- Protein.Complexes.Gold.Standard.Prediction@targetCluster
            
            
            cat("Step 0: Done.\n\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 0: Done.\n\n")
            
            cat( "Step 1: To construct reproducible, resolvable and Rapp mean >",  RappCutoffChoosingThreshold, "peaks from",  P_name, "protein SEC data for clustering proteins in the small world, within a CORUM orthologous complex.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat( "Step 1: To construct reproducible, resolvable and Rapp mean >",  RappCutoffChoosingThreshold, "peaks from",  P_name, "protein SEC data for clustering proteins in the small world, within a CORUM orthologous complex.\n")
            
            # find resolvable single peaks and the single peaks with Rapp.avg > Rapp threshold
            resolvablePeaks.FromSinglePeakProt <- c()
            resolvablePeaks.FromSinglePeakProt.Chosen <- c()
            RappMappAvg.resolvablePeaks.FromSinglePeakProt <- list()
            RappMappAvg.resolvablePeaks.FromSinglePeakProt.Chosen <- list()
            for (p in singlePeakProts){
              indexMaxB1 <- which.max(peakProfileMatrixB1[p,])
              indexMaxB2 <- which.max(peakProfileMatrixB2[p,])
              isResolvable <- !all(c(indexMaxB1==1, indexMaxB2==1))
              if (isResolvable){
                resolvablePeaks.FromSinglePeakProt <- c(resolvablePeaks.FromSinglePeakProt, p)
                RappMean <- (reproduciblePeakFeatures.Data[[p]][['B1']][['Rapp1']] + reproduciblePeakFeatures.Data[[p]][['B2']][['Rapp1']])/2
                MappMean <- (reproduciblePeakFeatures.Data[[p]][['B1']][['Mapp1']] + reproduciblePeakFeatures.Data[[p]][['B2']][['Mapp1']])/2
                RappMappAvg.resolvablePeaks.FromSinglePeakProt[[p]][['Rapp_mean']] <- RappMean
                RappMappAvg.resolvablePeaks.FromSinglePeakProt[[p]][['Mapp_mean']] <- MappMean
                isChosen <- RappMean > RappCutoffChoosingThreshold
                if (isChosen){
                  resolvablePeaks.FromSinglePeakProt.Chosen <- c(resolvablePeaks.FromSinglePeakProt.Chosen, p)
                  RappMappAvg.resolvablePeaks.FromSinglePeakProt.Chosen[[p]][['Rapp_mean']] <- RappMean
                  RappMappAvg.resolvablePeaks.FromSinglePeakProt.Chosen[[p]][['Mapp_mean']] <- MappMean
                }
              }
            }
            cat("\tstep 1.1 has filtered out", length(resolvablePeaks.FromSinglePeakProt.Chosen), 
                "reproducible, resolvable and Rapp mean >", RappCutoffChoosingThreshold, "peaks from", 
                length(singlePeakProts), P_name, "single-peak proteins.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tstep 1.1 has filtered out", length(resolvablePeaks.FromSinglePeakProt.Chosen), 
                "reproducible, resolvable and Rapp mean >", RappCutoffChoosingThreshold, "peaks from", 
                length(singlePeakProts), P_name)
            
            
            # peak profile matrix for the chosen peaks from single peak prots
            peakProfileMatrixB1.FromSinglePeakProt.Chosen <- peakProfileMatrixB1[resolvablePeaks.FromSinglePeakProt.Chosen,]
            peakProfileMatrixB2.FromSinglePeakProt.Chosen <- peakProfileMatrixB2[resolvablePeaks.FromSinglePeakProt.Chosen,]
            cat( "\tstep 1.2 has produced peak profile matrices for the peaks chosen by step 1.1\n")
            cat( "\tstep 1.2 has produced peak profile matrices for the peaks chosen by step 1.1\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            
            lengthPeakprofileVec <- ncol(peakProfileMatrixB1)
            fakeMultipeakProts <-c()
            peakProfileMatrixB1.FromMultiplePeakProt <- matrix(numeric(lengthPeakprofileVec), nrow = 1, byrow = TRUE)
            peakProfileMatrixB2.FromMultiplePeakProt <- matrix(numeric(lengthPeakprofileVec), nrow = 1, byrow = TRUE)
            cout_unresolvable_first_peak <- 0
            for(p in multiplePeakProts){
              index_max_1 <- which.max(peakProfileMatrixB1[p,])
              index_max_2 <- which.max(peakProfileMatrixB2[p,])
              Mapp_1_mean <- (reproduciblePeakFeatures.Data[[p]][['B1']][['Mapp1']] + reproduciblePeakFeatures.Data[[p]][['B2']][['Mapp1']])/2
              Mapp_2_mean <- (reproduciblePeakFeatures.Data[[p]][['B1']][['Mapp2']] + reproduciblePeakFeatures.Data[[p]][['B2']][['Mapp2']])/2
              Mapp_3_mean <- (reproduciblePeakFeatures.Data[[p]][['B1']][['Mapp3']] + reproduciblePeakFeatures.Data[[p]][['B2']][['Mapp3']])/2
              num_peak <- sum(c(!is.na(Mapp_1_mean), !is.na(Mapp_2_mean), !is.na(Mapp_3_mean)))
              
              peak_1_B1 <- reproduciblePeakFeatures.Data[[p]][['B1']][['Gaussian_start1']]
              peak_1_B2 <- reproduciblePeakFeatures.Data[[p]][['B2']][['Gaussian_start1']]
              peak_2_B1 <- reproduciblePeakFeatures.Data[[p]][['B1']][['Gaussian_start2']]
              peak_2_B2 <- reproduciblePeakFeatures.Data[[p]][['B2']][['Gaussian_start2']]
              peak_3_B1 <- reproduciblePeakFeatures.Data[[p]][['B1']][['Gaussian_start3']]
              peak_3_B2 <- reproduciblePeakFeatures.Data[[p]][['B2']][['Gaussian_start3']]
              
              Gaussian_peak_1_B1 <- reproduciblePeakFeatures.Data[[p]][['B1']][['Gaussian_peak1']]
              Gaussian_peak_1_B2 <- reproduciblePeakFeatures.Data[[p]][['B2']][['Gaussian_peak1']]
              Gaussian_peak_2_B1 <- reproduciblePeakFeatures.Data[[p]][['B1']][['Gaussian_peak2']]
              Gaussian_peak_2_B2 <- reproduciblePeakFeatures.Data[[p]][['B2']][['Gaussian_peak2']]
              Gaussian_peak_3_B1 <- reproduciblePeakFeatures.Data[[p]][['B1']][['Gaussian_peak3']]
              Gaussian_peak_3_B2 <- reproduciblePeakFeatures.Data[[p]][['B2']][['Gaussian_peak3']]
              
              if(sum(!is.na(c(peak_1_B1,peak_2_B1,peak_3_B1)))>num_peak | sum(!is.na(c(peak_1_B2,peak_2_B2,peak_3_B2)))>num_peak){
                if((sum(!is.na(c(peak_1_B1,peak_2_B1,peak_3_B1)))>num_peak)==F | (sum(!is.na(c(peak_1_B2,peak_2_B2,peak_3_B2)))>num_peak)==F){
                  fakeMultipeakProts <-c(fakeMultipeakProts, p)
                  next
                }
              }
              
              if(sum(!is.na(c(peak_1_B1,peak_2_B1,peak_3_B1)))>num_peak & sum(!is.na(c(peak_1_B2,peak_2_B2,peak_3_B2)))>num_peak){
                peak_1_B1 <- peak_2_B1
                peak_1_B2 <- peak_2_B2
                peak_2_B1 <- peak_3_B1
                peak_2_B2 <- peak_3_B2
                peak_3_B1 <- NA
                peak_3_B2 <- NA
              }
              
              if(all(c(!is.na(peak_1_B1), !is.na(peak_2_B1)))){
                trough_1_B1 <- peak_1_B1 - 1 + which.min(peakProfileMatrixB1[p, peak_1_B1:peak_2_B1])
              }
              if(all(c(!is.na(peak_1_B2), !is.na(peak_2_B2)))){
                trough_1_B2 <- peak_1_B2 - 1 + which.min(peakProfileMatrixB2[p, peak_1_B2:peak_2_B2])
              }
              if(all(c(!is.na(peak_2_B1), !is.na(peak_3_B1)))){
                trough_2_B1 <- peak_2_B1 - 1 + which.min(peakProfileMatrixB1[p, peak_2_B1:peak_3_B1])
              }
              if(all(c(!is.na(peak_2_B2), !is.na(peak_3_B2)))){
                trough_2_B2 <- peak_2_B2 - 1 + which.min(peakProfileMatrixB2[p, peak_2_B2:peak_3_B2])
              }
              
              if (any(c(index_max_1==1, index_max_2==1, Mapp_1_mean>850))){
                cout_unresolvable_first_peak <- cout_unresolvable_first_peak + 1
                if(num_peak==2){
                  if (abs(reproduciblePeakFeatures.Data[[p]][['B1']][[paste0('Gaussian_peak', 2)]]-
                          reproduciblePeakFeatures.Data[[p]][['B2']][[paste0('Gaussian_peak', 2)]])<=2){
                    M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:lengthPeakprofileVec])), nrow = 1)
                    rownames(M_p2_B1) <- paste(p, 22, sep='_')
                    peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                    M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:lengthPeakprofileVec])), nrow = 1)
                    rownames(M_p2_B2) <- paste(p, 22, sep='_')
                    peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                  }
                }else if(num_peak==3){
                  for (i in 2:3){
                    for (j in 2:3){
                      if (abs(reproduciblePeakFeatures.Data[[p]][['B1']][[paste0('Gaussian_peak', i)]]-
                              reproduciblePeakFeatures.Data[[p]][['B2']][[paste0('Gaussian_peak', j)]])<=2){
                        if (all(c(i==2,j==2))){
                          M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:trough_2_B1]), numeric(lengthPeakprofileVec-trough_2_B1)), nrow = 1)
                          rownames(M_p2_B1) <- paste(p, 22, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                          M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:trough_2_B2]), numeric(lengthPeakprofileVec-trough_2_B2)), nrow = 1)
                          rownames(M_p2_B2) <- paste(p, 22, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                        }
                        if (all(c(i==2,j==3))){
                          M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:trough_2_B1]), numeric(lengthPeakprofileVec-trough_2_B1)), nrow = 1)
                          rownames(M_p2_B1) <- paste(p, 23, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                          M_p3_B2 <- matrix(c(numeric(trough_2_B2-1), as.numeric(peakProfileMatrixB2[p, trough_2_B2:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B2) <- paste(p, 23, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p3_B2)
                        }
                        if (all(c(i==3, j==2))){
                          M_p3_B1 <- matrix(c(numeric(trough_2_B1-1), as.numeric(peakProfileMatrixB1[p, trough_2_B1:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B1) <- paste(p, 32, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p3_B1)
                          M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:trough_2_B2]), numeric(lengthPeakprofileVec-trough_2_B2)), nrow = 1)
                          rownames(M_p2_B2) <- paste(p, 32, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                        }
                        if (all(c(i==3, j==3))){
                          M_p3_B1 <- matrix(c(numeric(trough_2_B1-1), as.numeric(peakProfileMatrixB1[p, trough_2_B1:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B1) <- paste(p, 33, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p3_B1)
                          M_p3_B2 <- matrix(c(numeric(trough_2_B2-1), as.numeric(peakProfileMatrixB2[p, trough_2_B2:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B2) <- paste(p, 33, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p3_B2)
                        }
                        
                      }
                    }
                  }
                }
              }else{
                if(num_peak==2){
                  for (i in 1:2){
                    for (j in 1:2){
                      if (abs(reproduciblePeakFeatures.Data[[p]][['B1']][[paste0('Gaussian_peak', i)]]-
                              reproduciblePeakFeatures.Data[[p]][['B2']][[paste0('Gaussian_peak', j)]])<=2){
                        if (all(c(i==1, j==1))){
                          M_p1_B1 <- matrix(c(as.numeric(peakProfileMatrixB1[p, 1:trough_1_B1]), numeric(lengthPeakprofileVec-trough_1_B1)), nrow = 1)
                          rownames(M_p1_B1) <- paste(p, 11, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p1_B1)
                          M_p1_B2 <- matrix(c(as.numeric(peakProfileMatrixB2[p, 1:trough_1_B2]), numeric(lengthPeakprofileVec-trough_1_B2)), nrow = 1)
                          rownames(M_p1_B2) <- paste(p, 11, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p1_B2)
                        }
                        if (all(c(i==1, j==2))){
                          M_p1_B1 <- matrix(c(as.numeric(peakProfileMatrixB1[p, 1:trough_1_B1]), numeric(lengthPeakprofileVec-trough_1_B1)), nrow = 1)
                          rownames(M_p1_B1) <- paste(p, 12, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p1_B1)
                          M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p2_B2) <- paste(p, 12, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                        }
                        if (all(c(i==2, j==1))){
                          M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p2_B1) <- paste(p, 21, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                          
                          M_p1_B2 <- matrix(c(as.numeric(peakProfileMatrixB2[p, 1:trough_1_B2]), numeric(lengthPeakprofileVec-trough_1_B2)), nrow = 1)
                          rownames(M_p1_B2) <- paste(p, 21, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p1_B2)
                        }
                        if (all(c(i==2, j==2))){
                          M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p2_B1) <- paste(p, 22, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                          M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p2_B2) <- paste(p, 22, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                        }
                      }
                    }
                  }
                  
                }else if(num_peak==3){
                  for (i in 1:3){
                    for (j in 1:3){
                      if (abs(reproduciblePeakFeatures.Data[[p]][['B1']][[paste0('Gaussian_peak', i)]]-
                              reproduciblePeakFeatures.Data[[p]][['B2']][[paste0('Gaussian_peak', j)]])<=2){
                        if (all(c(i==1, j==1))){
                          M_p1_B1 <- matrix(c(as.numeric(peakProfileMatrixB1[p, 1:trough_1_B1]), numeric(lengthPeakprofileVec-trough_1_B1)), nrow = 1)
                          rownames(M_p1_B1) <- paste(p, 11, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p1_B1)
                          M_p1_B2 <- matrix(c(as.numeric(peakProfileMatrixB2[p, 1:trough_1_B2]), numeric(lengthPeakprofileVec-trough_1_B2)), nrow = 1)
                          rownames(M_p1_B2) <- paste(p, 11, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p1_B2)
                        }
                        if (all(c(i==1, j==2))){
                          M_p1_B1 <- matrix(c(as.numeric(peakProfileMatrixB1[p, 1:trough_1_B1]), numeric(lengthPeakprofileVec-trough_1_B1)), nrow = 1)
                          rownames(M_p1_B1) <- paste(p, 12, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p1_B1)
                          M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:trough_2_B2]), numeric(lengthPeakprofileVec-trough_2_B2)), nrow = 1)
                          rownames(M_p2_B2) <- paste(p, 12, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                        }
                        if (all(c(i==1, j==3))){
                          M_p1_B1 <- matrix(c(as.numeric(peakProfileMatrixB1[p, 1:trough_1_B1]), numeric(lengthPeakprofileVec-trough_1_B1)), nrow = 1)
                          rownames(M_p1_B1) <- paste(p, 13, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p1_B1)
                          M_p3_B2 <- matrix(c(numeric(trough_2_B2-1), as.numeric(peakProfileMatrixB2[p, trough_2_B2:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B2) <- paste(p, 13, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p3_B2)
                        }
                        if (all(c(i==2, j==1))){
                          M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:trough_2_B1]), numeric(lengthPeakprofileVec-trough_2_B1)), nrow = 1)
                          rownames(M_p2_B1) <- paste(p, 21, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                          M_p1_B2 <- matrix(c(as.numeric(peakProfileMatrixB2[p, 1:trough_1_B2]), numeric(lengthPeakprofileVec-trough_1_B2)), nrow = 1)
                          rownames(M_p1_B2) <- paste(p, 21, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p1_B2)
                          
                        }
                        if (all(c(i==2, j==2))){
                          M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:trough_2_B1]), numeric(lengthPeakprofileVec-trough_2_B1)), nrow = 1)
                          rownames(M_p2_B1) <- paste(p, 22, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                          M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:trough_2_B2]), numeric(lengthPeakprofileVec-trough_2_B2)), nrow = 1)
                          rownames(M_p2_B2) <- paste(p, 22, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                        }
                        if (all(c(i==2, j==3))){
                          M_p2_B1 <- matrix(c(numeric(trough_1_B1-1), as.numeric(peakProfileMatrixB1[p, trough_1_B1:trough_2_B1]), numeric(lengthPeakprofileVec-trough_2_B1)), nrow = 1)
                          rownames(M_p2_B1) <- paste(p, 23, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p2_B1)
                          M_p3_B2 <- matrix(c(numeric(trough_2_B2-1), as.numeric(peakProfileMatrixB2[p, trough_2_B2:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B2) <- paste(p, 23, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p3_B2)
                        }
                        if (all(c(i==3, j==1))){
                          M_p3_B1 <- matrix(c(numeric(trough_2_B1-1), as.numeric(peakProfileMatrixB1[p, trough_2_B1:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B1) <- paste(p, 31, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p3_B1)
                          M_p1_B2 <- matrix(c(as.numeric(peakProfileMatrixB2[p, 1:trough_1_B2]), numeric(lengthPeakprofileVec-trough_1_B2)), nrow = 1)
                          rownames(M_p1_B2) <- paste(p, 31, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p1_B2)
                        }
                        if (all(c(i==3, j==2))){
                          M_p3_B1 <- matrix(c(numeric(trough_2_B1-1), as.numeric(peakProfileMatrixB1[p, trough_2_B1:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B1) <- paste(p, 32, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p3_B1)
                          M_p2_B2 <- matrix(c(numeric(trough_1_B2-1), as.numeric(peakProfileMatrixB2[p, trough_1_B2:trough_2_B2]), numeric(lengthPeakprofileVec-trough_2_B2)), nrow = 1)
                          rownames(M_p2_B2) <- paste(p, 32, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p2_B2)
                        }
                        if (all(c(i==3, j==3))){
                          M_p3_B1 <- matrix(c(numeric(trough_2_B1-1), as.numeric(peakProfileMatrixB1[p, trough_2_B1:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B1) <- paste(p, 33, sep='_')
                          peakProfileMatrixB1.FromMultiplePeakProt <- rbind(peakProfileMatrixB1.FromMultiplePeakProt,M_p3_B1)
                          M_p3_B2 <- matrix(c(numeric(trough_2_B2-1), as.numeric(peakProfileMatrixB2[p, trough_2_B2:lengthPeakprofileVec])), nrow = 1)
                          rownames(M_p3_B2) <- paste(p, 33, sep='_')
                          peakProfileMatrixB2.FromMultiplePeakProt <- rbind(peakProfileMatrixB2.FromMultiplePeakProt,M_p3_B2)
                        }
                      }
                    }
                  }
                }
              }
            }
            
            peakProfileMatrixB1.FromMultiplePeakProt <- peakProfileMatrixB1.FromMultiplePeakProt[-1,]
            colnames(peakProfileMatrixB1.FromMultiplePeakProt) <- colnames(peakProfileMatrixB1.FromSinglePeakProt.Chosen)
            peakProfileMatrixB2.FromMultiplePeakProt <- peakProfileMatrixB2.FromMultiplePeakProt[-1,]
            colnames(peakProfileMatrixB2.FromMultiplePeakProt) <- colnames(peakProfileMatrixB2.FromSinglePeakProt.Chosen)
            
            
            resolvablePeaks.FromMultiplePeakProt <- rownames(peakProfileMatrixB1.FromMultiplePeakProt)
            RappMappAvg.resolvablePeaks.FromMultiplePeakProt <- list()
            RappMappAvg.resolvablePeaks.FromMultiplePeakProt.Chosen <- list()
            
            for (p in resolvablePeaks.FromMultiplePeakProt){
              real_p <- substr(p, 1, ncharProtID)
              peak_id_1 <- substr(p, ncharProtID+2, ncharProtID+2)
              peak_id_2 <- substr(p, ncharProtID+3, ncharProtID+3)
              MappMean <- (reproduciblePeakFeatures.Data[[real_p]][['B1']][[paste0("Mapp",peak_id_1)]]+
                             reproduciblePeakFeatures.Data[[real_p]][['B2']][[paste0("Mapp",peak_id_2)]])/2
              RappMean <- (reproduciblePeakFeatures.Data[[real_p]][['B1']][[paste0("Rapp",peak_id_1)]]+
                             reproduciblePeakFeatures.Data[[real_p]][['B2']][[paste0("Rapp",peak_id_2)]])/2
              RappMappAvg.resolvablePeaks.FromMultiplePeakProt[[p]][['Rapp_mean']] <- RappMean
              RappMappAvg.resolvablePeaks.FromMultiplePeakProt[[p]][['Mapp_mean']] <- MappMean
              isChosen <- RappMean > RappCutoffChoosingThreshold
              if (isChosen){
                RappMappAvg.resolvablePeaks.FromMultiplePeakProt.Chosen[[p]][['Rapp_mean']] <- RappMean
                RappMappAvg.resolvablePeaks.FromMultiplePeakProt.Chosen[[p]][['Mapp_mean']] <- MappMean
              }
            }
            
            resolvablePeaks.FromMultiplePeakProt.Chosen <- names(RappMappAvg.resolvablePeaks.FromMultiplePeakProt.Chosen)
            cat("\tstep 1.3 has filtered out", length(resolvablePeaks.FromMultiplePeakProt.Chosen), 
                "reproducible, resolvable and Rapp mean >", RappCutoffChoosingThreshold, "peaks from", 
                length(multiplePeakProts), P_name, "multiple-peak proteins.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tstep 1.3 has filtered out", length(resolvablePeaks.FromMultiplePeakProt.Chosen), 
                "reproducible, resolvable and Rapp mean >", RappCutoffChoosingThreshold, "peaks from", 
                length(multiplePeakProts), P_name, "multiple-peak proteins.\n")
            
            
            peakProfileMatrixB1.FromMultiplePeakProt.Chosen <- peakProfileMatrixB1.FromMultiplePeakProt[resolvablePeaks.FromMultiplePeakProt.Chosen, ]
            peakProfileMatrixB2.FromMultiplePeakProt.Chosen <- peakProfileMatrixB2.FromMultiplePeakProt[resolvablePeaks.FromMultiplePeakProt.Chosen, ]
            cat("\tstep 1.4 has produced peak profile matrices for the peaks chosen by step 1.3\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tstep 1.4 has produced peak profile matrices for the peaks chosen by step 1.3\n")
            
            peakProfileMatrixB1.Chosen <- rbind(peakProfileMatrixB1.FromSinglePeakProt.Chosen, peakProfileMatrixB1.FromMultiplePeakProt.Chosen)
            peakProfileMatrixB2.Chosen <- rbind(peakProfileMatrixB2.FromSinglePeakProt.Chosen, peakProfileMatrixB2.FromMultiplePeakProt.Chosen)
            cat("\tstep 1.5 has produced peak profile matrices for the peaks chosen by steps 1.1 and 1.2, both single-and-multiple-peak proteins\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tstep 1.5 has produced peak profile matrices for the peaks chosen by steps 1.1 and 1.2, both single-and-multiple-peak proteins\n")
            
            # mapping: protein ID -> peak
            multiplePeakProts.Chosen <- substr(resolvablePeaks.FromMultiplePeakProt.Chosen, 1, ncharProtID)
            multiplePeakProts.Chosen <- levels(factor(multiplePeakProts.Chosen, levels = unique(multiplePeakProts.Chosen)))
            multiplePeakProtsToPeaksMap <- list()
            for (p in multiplePeakProts.Chosen){
              multiplePeakProtsToPeaksMap[[p]] <- resolvablePeaks.FromMultiplePeakProt.Chosen[substr(resolvablePeaks.FromMultiplePeakProt.Chosen, 1, ncharProtID) %in% p]
            }
            cat("\tstep 1.6 has produced a mapping from multiple-peak protein ID to peaks generated from multiple-peak proteins.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tstep 1.6 has produced a mapping from multiple-peak protein ID to peaks generated from multiple-peak proteins.\n")

            GaussianPeakFractionB1FromSinglePeakProt <- 
              matrix(reproduciblePeakFeatures.Data[[rownames(peakProfileMatrixB1.FromSinglePeakProt.Chosen)[1]]][['B1']][['Gaussian_peak1']], nrow=1, dimnames=list(c(rownames(peakProfileMatrixB1.FromSinglePeakProt.Chosen)[1]), c("Gaussian_peak_fraction")))
            GaussianPeakFractionB2FromSinglePeakProt <- 
              matrix(reproduciblePeakFeatures.Data[[rownames(peakProfileMatrixB2.FromSinglePeakProt.Chosen)[1]]][['B2']][['Gaussian_peak1']], nrow=1, dimnames=list(c(rownames(peakProfileMatrixB2.FromSinglePeakProt.Chosen)[1]), c("Gaussian_peak_fraction")))
            for (p in rownames(peakProfileMatrixB1.FromSinglePeakProt.Chosen)[-1]){
              m1 <- matrix(reproduciblePeakFeatures.Data[[p]][['B1']][['Gaussian_peak1']], nrow=1, dimnames=list(c(p), c("Gaussian_peak_fraction")))
              GaussianPeakFractionB1FromSinglePeakProt <- rbind(GaussianPeakFractionB1FromSinglePeakProt, m1)
              m2 <- matrix(reproduciblePeakFeatures.Data[[p]][['B2']][['Gaussian_peak1']], nrow=1, dimnames=list(c(p), c("Gaussian_peak_fraction")))
              GaussianPeakFractionB2FromSinglePeakProt <- rbind(GaussianPeakFractionB2FromSinglePeakProt, m2)
            }
            
            for (p in rownames(peakProfileMatrixB1.FromMultiplePeakProt.Chosen)){
              real_p <- substr(p, 1, ncharProtID)
              peak_index_1 <- substr(p, ncharProtID+2, ncharProtID+2)
              peak_index_2 <- substr(p, ncharProtID+3, ncharProtID+3)
              m1 <- matrix(reproduciblePeakFeatures.Data[[real_p]][['B1']][[paste0('Gaussian_peak', peak_index_1)]], nrow=1, dimnames=list(c(p), c("Gaussian_peak_fraction")))
              GaussianPeakFractionB1FromSinglePeakProt <- rbind(GaussianPeakFractionB1FromSinglePeakProt, m1)
              m2 <- matrix(reproduciblePeakFeatures.Data[[real_p]][['B2']][[paste0('Gaussian_peak', peak_index_2)]], nrow=1, dimnames=list(c(p), c("Gaussian_peak_fraction")))
              GaussianPeakFractionB2FromSinglePeakProt <- rbind(GaussianPeakFractionB2FromSinglePeakProt, m2)
            }
            
            GaussianPeakFractionB1.Chosen <- GaussianPeakFractionB1FromSinglePeakProt
            GaussianPeakFractionB2.Chosen <- GaussianPeakFractionB2FromSinglePeakProt
            cat( "\tstep 1.7 has produced Gaussian peak franction matrices for the peaks chosen by steps 1.1 and 1.2, both single-and-multiple-peak proteins.\n")
            cat( "\tstep 1.7 has produced Gaussian peak franction matrices for the peaks chosen by steps 1.1 and 1.2, both single-and-multiple-peak proteins.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            
            write.csv(as.data.frame(peakProfileMatrixB1.Chosen, row.names=rownames(peakProfileMatrixB1.Chosen)), 
                      paste(paste("Table 1-1 Peak Profile B1", Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
            write.csv(as.data.frame(GaussianPeakFractionB1.Chosen, row.names=rownames(GaussianPeakFractionB1.Chosen)), 
                      paste(paste("Table 1-2 Gaussian Peak Fraction B1", Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
            write.csv(as.data.frame(peakProfileMatrixB2.Chosen, row.names=rownames(peakProfileMatrixB2.Chosen)), 
                      paste(paste("Table 1-3 Peak Profile B12", Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
            write.csv(as.data.frame(GaussianPeakFractionB2.Chosen, row.names=rownames(GaussianPeakFractionB1.Chosen)), 
                      paste(paste("Table 1-4 Gaussian Peak Fraction B2", Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
            cat( "\tstep 1.8 has produced Table 1: Peak Profile and Gaussian Peak Fraction for the peaks chosen by steps 1.1 and 1.2, both single-and-multiple-peak proteins.\n")
            cat( "\tstep 1.8 has produced Table 1: Peak Profile and Gaussian Peak Fraction for the peaks chosen by steps 1.1 and 1.2, both single-and-multiple-peak proteins.\n", file="LOG OF COMPUTATION.txt", append=TRUE)

            protsSelectedForSmallWorld <- rownames(peakProfileMatrixB1.Chosen)
            cat("Step 1: Done. The peaks for small world analysis are the peaks chosen by steps 1.1 and 1.2, from both single-and-multiple-peak proteins.\n\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 1: Done. The peaks for small world analysis are the peaks chosen by steps 1.1 and 1.2, from both single-and-multiple-peak proteins.\n\n")
            
            peaksConstructedForSmallWorldClustering <- list(peakProfileMatrixB1.FromSinglePeakProt.Chosen=peakProfileMatrixB1.FromSinglePeakProt.Chosen,
                                                            peakProfileMatrixB1.FromSinglePeakProt.Chosen=peakProfileMatrixB1.FromSinglePeakProt.Chosen,
                                                            peakProfileMatrixB1.FromMultiplePeakProt.Chosen=peakProfileMatrixB1.FromMultiplePeakProt.Chosen,
                                                            peakProfileMatrixB2.FromMultiplePeakProt.Chosen=peakProfileMatrixB2.FromMultiplePeakProt.Chosen,
                                                            peakProfileMatrixB1.Chosen=peakProfileMatrixB1.Chosen,
                                                            peakProfileMatrixB2.Chosen=peakProfileMatrixB2.Chosen,
                                                            GaussianPeakFractionB1.Chosen=GaussianPeakFractionB1.Chosen,
                                                            GaussianPeakFractionB2.Chosen=GaussianPeakFractionB2.Chosen,
                                                            multiplePeakProtsToPeaksMap=multiplePeakProtsToPeaksMap,
                                                            protsSelectedForSmallWorld=protsSelectedForSmallWorld)
            
            
            cat("Step 2: To construct CORUM-Inparanoid Orthologous Complexes.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 2: To construct CORUM-Inparanoid Orthologous Complexes.\n")

            # set up M-P ortholog mapping
            M_species <- mammalPlantOrtholog.Data$Species[1]
            P_species <- mammalPlantOrtholog.Data$Species[2]
            M_uniprot_space_with_ortholog <- levels(factor(filter(mammalPlantOrtholog.Data, Species==M_species)$ProtID))
            P_ortholog_space <- levels(factor(filter(mammalPlantOrtholog.Data, Species==P_species)$ProtID))
            if(length(M_uniprot_space_with_ortholog) == length(filter(mammalPlantOrtholog.Data, Species==M_species)$ProtID)){
              print("No duplicate in Mammallian protein ID's in the ortholog data.")
            }else{
              warning("Duplicates exist mammalian protein ID's in the ortholog data.")
            }
            
            MprotToPprotMap <- list() # map: each mammalian uniprot -> plant orthologs
            for (i in min(mammalPlantOrtholog.Data$GroupID):max(mammalPlantOrtholog.Data$GroupID)){
              g <- filter(mammalPlantOrtholog.Data, GroupID==i) #g is index of each group
              u_set <- filter(g, Species==M_species)$ProtID #mammalian Uniprots in the group
              p_set <- filter(g, Species==P_species)$ProtID # plantID in the group
              for (u in u_set){
                MprotToPprotMap[[u]] <- p_set # correspondence between mammalian and plant proteins
              }
            }
            
            if(length(P_ortholog_space) == length(filter(mammalPlantOrtholog.Data, Species==P_species)$ProtID)){
              print("No duplicate in plant protein ID's in the ortholog data")
            }else{
              warning("Duplicates exist plant protein ID's in the ortholog data.")
            }
            
            PprotToMprotMap <- list() # map: each  plant ortholog -> mammalian uniprots
            for (i in min(mammalPlantOrtholog.Data$GroupID):max(mammalPlantOrtholog.Data$GroupID)){
              g <- filter(mammalPlantOrtholog.Data, GroupID==i) #g is index of each group
              p_set <- filter(g, Species==P_species)$ProtID # plant IDs in the group
              u_set <- filter(g, Species==M_species)$ProtID # mammalian Uniprots in the group
              for (p in p_set){
                PprotToMprotMap[[p]] <- u_set # correspondence mammalian mammalian and plant proteins
              }
            }
            cat("\tstep 2.1 has set up ortholog correspondence between", length(M_uniprot_space_with_ortholog), M_name, 
                "proteins and", length(P_ortholog_space), P_name, "proteins, based on Inparanoid ortho-paralog information.\n",
                file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tstep 2.1 has set up ortholog correspondence between", length(M_uniprot_space_with_ortholog), M_name, 
                "proteins and", length(P_ortholog_space), P_name, "proteins, based on Inparanoid ortho-paralog information.\n")
            
            # set up orthologous complexes
            CORUMComplexNames <- levels(factor(mammalCORUM.Data$ComplexName))
            genomeCoverage <- list()
            MprotToPprotMap.inEachCORUMComplex <- list()
            CORUMComplexSize <- list()
            numCoveredByOrtholog <- list()
            orthologComplexSize <- list()
            
            for (c in CORUMComplexNames){
              MprotsInComplex <- levels(factor(filter(mammalCORUM.Data, ComplexName==c)$Mammalian_Prot_ID))
              c_size <- length(MprotsInComplex)
              howManyCoveredByOrtholog <- MprotsInComplex %in% names(MprotToPprotMap) #covered by ortholog assignment
              if (any(howManyCoveredByOrtholog)){
                genomeCoverage[[c]] <- sum(howManyCoveredByOrtholog)/c_size
                CORUMComplexSize[[c]] <- c_size
                `numCoveredByOrtholog`[[c]] <- sum(howManyCoveredByOrtholog)
                MprotToPprotMap.inEachCORUMComplex[[c]] <- list()
                for (u in MprotsInComplex){
                  if (length(MprotToPprotMap[[u]])>0){
                    MprotToPprotMap.inEachCORUMComplex[[c]][[u]] <- MprotToPprotMap[[u]]
                  }else{
                    MprotToPprotMap.inEachCORUMComplex[[c]][[u]]<- NA
                  }
                }
                orthologComplexSize[[c]] <- sum(!is.na(unlist(MprotToPprotMap.inEachCORUMComplex[[c]])))
              }
            }
            # table of orthologous complexes
            tableOfOrthoCORUMComplexes <- NULL
            for (c in names(genomeCoverage)){
              table_c <- filter(mammalCORUM.Data, ComplexName==c)
              rownames(table_c) <- table_c$Mammalian_Prot_ID
              for (u in rownames(table_c)){
                table_c[u, 'Ortholog_ID'] <- paste(MprotToPprotMap[[u]], collapse=',')
              }
              table_c[, 'GenomeCoverage'] <- numeric(length(rownames(table_c))) + genomeCoverage[[c]]
              table_c[, 'Number_of_ortholog_subunits_in_a_given_complex'] <- numeric(length(rownames(table_c))) +orthologComplexSize[[c]]
              table_c[, 'Number_subunits'] <- numeric(length(rownames(table_c))) + CORUMComplexSize[[c]]
              tableOfOrthoCORUMComplexes <- rbind(tableOfOrthoCORUMComplexes, table_c)
            }
            write.csv(tableOfOrthoCORUMComplexes, paste(paste(paste0(paste0(paste0(paste0("Table 2. ", P_name), " orthologous CORUM "), M_name), " complexes"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = F)
            cat("\tstep 2.2 has set up orthologous CORUM complexes and generated Table 2:",  P_name, " orthologous CORUM ", M_name, " complexes.\n")
            cat("\tstep 2.2 has set up orthologous CORUM complexes and generated Table 2:",  P_name, " orthologous CORUM ", M_name, " complexes.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            
            orthoComplexesAboveGcoverageCutoff <- c()
            orthoComplexesAboveGcoverageCutoff.withMultiorthologs <- c()
            orthoComplexesBelowGcoverageCutoff <- c()
            orthoComplexesBelowGcoverageCutoff.withMultiorthologs <- c()
            orthoComplexesFullGcoverage <- c()
            
            for (c in names(genomeCoverage)){
              if (genomeCoverage[[c]]>=genomeCoverageCutoff){
                orthoComplexesAboveGcoverageCutoff <- c(orthoComplexesAboveGcoverageCutoff, c)
                ortholog_set <- unname(unlist(MprotToPprotMap.inEachCORUMComplex[[c]]))
                if (length(ortholog_set[!is.na(ortholog_set)])>1){
                  orthoComplexesAboveGcoverageCutoff.withMultiorthologs <- c(orthoComplexesAboveGcoverageCutoff.withMultiorthologs, c)
                }
              }else{
                orthoComplexesBelowGcoverageCutoff <- c(orthoComplexesBelowGcoverageCutoff, c)
                ortholog_set <- unname(unlist(MprotToPprotMap.inEachCORUMComplex[[c]]))
                if (length(ortholog_set[!is.na(ortholog_set)])>1){
                  orthoComplexesBelowGcoverageCutoff.withMultiorthologs <- c(orthoComplexesBelowGcoverageCutoff.withMultiorthologs, c)
                }
              }
              if (genomeCoverage[[c]]>0.9999999){
                orthoComplexesFullGcoverage <- c(orthoComplexesFullGcoverage, c)
              }
            }
            orthoCOMRUMCOmplex.sepratedByGcoverageCutoff <- list(orthoComplexesAboveGcoverageCutoff=orthoComplexesAboveGcoverageCutoff, 
                                                                 orthoComplexesAboveGcoverageCutoff.withMultiorthologs=orthoComplexesAboveGcoverageCutoff.withMultiorthologs,
                                                                 orthoComplexesBelowGcoverageCutoff=orthoComplexesBelowGcoverageCutoff,
                                                                 orthoComplexesBelowGcoverageCutoff.withMultiorthologs=orthoComplexesBelowGcoverageCutoff.withMultiorthologs)
            cat( "\tstep 2.3 has separated orthologous CORUM complexes according to if above or below the genome coverage value =", genomeCoverageCutoff, 
                 "and if assigned with multiple orthologs.\n\t\t 1) CORUM has", length(CORUMComplexNames), "distinct", M_name, "complexes.", length(genomeCoverage), "of them have positive genome coverage.\n", 
                 "\t\t 2) There are",
                 length(orthoComplexesAboveGcoverageCutoff.withMultiorthologs), "out of", length(genomeCoverage), "complexes that respectively and correspondingly have their own multiple", 
                 P_name, "orthologs and high genome coverages.\n", 
                 "\t\t 3) There are",
                 length(orthoComplexesFullGcoverage), "out of", length(CORUMComplexNames), "complexes that have 100% genome coverages.\n",
                 "\t\t 4) Note that the CORUM complexes are identified by their complex name and CORUM ID which were used as the ComplexName in this project.\n",
                 file="LOG OF COMPUTATION.txt", append=TRUE)

            # plot number of subunits of Mammal complex
            num_subunits <- c()
            num_type <- c()
            for (c in levels(factor(tableOfOrthoCORUMComplexes$ComplexName))){
              c_t <- filter(tableOfOrthoCORUMComplexes, ComplexName==c)
              num_subunits <- c(num_subunits, c_t$Number_subunits[1])
              num_type <- c(num_type, as.character(c_t$Number_subunits[1]))
            }
            num_subunit_data <- data.frame(num_s=num_subunits,n_type=num_type)
            peak_subunit_all <- ggplot(num_subunit_data, aes(x=factor(num_s)))+
              geom_bar(stat="count", width=0.8, fill="steelblue")+
              labs(x=paste0(paste0("Number of subunits in ", M_name), " complexes"), y = "Count")+
              theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                    axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))
            
            peaks_numsubunits <- c(1/3, 1/2, 2/3, 1)
            peak_subunit_plot <- list()
            j <- 0
            for (i in peaks_numsubunits){
              j <- j + 1
              num_subunits_i <- c()
              for (c in levels(factor(tableOfOrthoCORUMComplexes$ComplexName))){
                c_t <- filter(tableOfOrthoCORUMComplexes, ComplexName==c)
                if (c_t$GenomeCoverage[1]==as.numeric(i)){
                  num_subunits_i <- c(num_subunits_i, c_t$Number_subunits[1])
                }
              }
              num_subunit_data_partial <- data.frame(num_subunits=num_subunits_i)
              peak_subunit_plot[[j]] <- ggplot(num_subunit_data_partial, aes(x=factor(num_subunits)))+
                geom_bar(stat="count", width=0.5, fill="steelblue")+
                labs(title=paste(paste0(paste0("Subunit coverage of CORUM complexes by ", P_name), " genome ="), fractions(i), sep=" "),x=paste0(paste0("Number of subunits in ", P_name), " orthocomplexes"), y = "Count")+
                theme(plot.title=element_text(hjust=1, vjust=0.5, margin=margin(t=40,b=-30)), title =element_text(size=10),
                      axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                      axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))
            }
            
            pdf(paste(paste(paste0(paste(paste0("figure file 1. Bar plot - Number of subunits in", M_name), P_name, sep="-"), " orthologous complexes" ), Sys.Date(), sep='_'), "pdf", sep='.'), width=12, height=8)
            print(peak_subunit_all)
            for (i in 1:length( peak_subunit_plot)){
              print( peak_subunit_plot[[i]])
            }
            dev.off()  
            cat("\tstep 2.4 has generated figure file 1: Bar plots - Number of subunits in", M_name, "complexes, based on Table 2.\n")
            cat("\tstep 2.4 has generated figure file 1: Bar plots - Number of subunits in", M_name, "complexes, based on Table 2.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            
            # plot geome-coverage dist for M-P complexes
            corum_cover_P <- unlist(genomeCoverage)
            data_cover_P <- data.frame(cover = corum_cover_P)
            data_cover_round_P <- mutate(data_cover_P, cover_r=round(cover,2))
            pdf(paste(paste(paste0(paste0(paste0(paste0("figure file 2. Genome-coverage disttribution plot for ", M_name), "-"), P_name), " complexes"), Sys.Date(), sep='_'), "pdf", sep='.'), width=8, height=8)
            print(ggplot(data_cover_round_P, aes(x=factor(cover_r)))+
                    geom_bar(stat="count", width=0.8, fill="steelblue")+
                    labs(x=paste0(paste0("Subunit coverage of CORUM complexes by ",  P_name), " genome"), y = "Count")+
                    theme(plot.title=element_text(hjust=1, vjust=0.5, margin=margin(t=40,b=-30)), title =element_text(size=8),
                          axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                          axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))+
                    scale_x_discrete(breaks=c(0, 0.33, 0.5, 0.67, 1), labels=c(0, 0.33, 0.5, 0.67, 1)))
            dev.off() 
            cat("\tstep 2.5 has generated figure file 2: Bar plot - Genome-coverage distribution plot for", M_name,"-", P_name, "orthologous complexes, based on Table 2.\n")
            cat("\tstep 2.5 has generated figure file 2: Bar plot - Genome-coverage distribution plot for", M_name,"-", P_name, "orthologous complexes, based on Table 2.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            
            Mcalc_c_fullcoverage <- list()
            for (c in orthoComplexesFullGcoverage){
              orthologs_c <- unname(unlist(MprotToPprotMap.inEachCORUMComplex[[c]]))
              if (length(orthologs_c)<1){
                next
              }
              Mcalc_c_fullcoverage[[c]] <- Mcalc_calculation(complex_name=c, Mcomplex_to_Plant_mapping=MprotToPprotMap.inEachCORUMComplex[[c]], 
                                                             ortholog_set=orthologs_c,
                                                             Mmnono_data=plantMmono.Data,
                                                             stoichiometry_data=mammalCORUM.StoichiometryData)$Mcalc
            }
            pdf(paste(paste("figure file 3. Scatter plot - Mcalc of Mammal complex v.s. Mcalc of Plant complex (full coverage)", Sys.Date(), sep='_'), "pdf", sep='.'), width=8, height=8)
            Mcalc_c_conserve_P <- c()
            Mcalc_c_conserve_M <- c()
            complexes_for_conserve <- c()
            for (c in orthoComplexesFullGcoverage){
              if(c %in% names(Mcalc_c_fullcoverage)){
                complexes_for_conserve <- c(complexes_for_conserve, c)
                Mcalc_c_conserve_P <- c(Mcalc_c_conserve_P, Mcalc_c_fullcoverage[[c]])
                temp_c_M <- filter(mammalCORUM.Data, ComplexName==c)
                mass_M <- c()
                for (u in temp_c_M$Mammalian_Prot_ID){
                  temp_Mu <- filter(temp_c_M, Mammalian_Prot_ID==u)
                  if (c %in% names(mammalCORUM.StoichiometryData) & u %in% names(mammalCORUM.StoichiometryData[[c]])){
                    sto_n <- mammalCORUM.StoichiometryData[[c]][[u]]$stoichiometry
                  }else{
                    sto_n <- 1
                  }
                  mass_M <- c(mass_M, temp_Mu$Mass_kDa[1]*sto_n)
                }
                Mcalc_c_conserve_M <- c(Mcalc_c_conserve_M, sum(mass_M))
              }
            }
            Mcalc_plant <- Mcalc_c_conserve_P
            Mcalc_mammal <- Mcalc_c_conserve_M
            mcalc_covserve_MP <- data.frame(Mcalc_plant=Mcalc_plant, Mcalc_mammal=Mcalc_mammal)
            limit_xy <- max(c(max(Mcalc_c_conserve_P),max(Mcalc_c_conserve_M)))
            l <- limit_xy %/% 10
            if (l>100){
              length_here <- 100
            }else if (l<=100 && l>50){
              length_here <- 50
            }else{
              length_here <- 25
            }
            print(ggplot(mcalc_covserve_MP, aes(Mcalc_plant, Mcalc_mammal) ) +   geom_point(size=1) + geom_smooth(method=lm, size=0.5, alpha=0.6)+
                    labs(x=paste0(paste0("Mcalc of ", P_name), " orthocomplexes (kDa)"), y = paste0(paste0("Mcalc of ", M_name), " complexes (kDa)"))+
                    theme(plot.title=element_text(hjust=1, vjust=0.5, margin=margin(t=40,b=-30)), title =element_text(size=12),
                          axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                          axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))+
                    scale_x_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy))+
                    scale_y_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy)))
            dev.off()
            LM_covserve_MP <- lm(Mcalc_mammal ~ Mcalc_plant)
            cat("\tstep 2.6 has generated figure file 3: Scatter plot - Mcalc of", M_name,"complexes V.S. Mcalc of", P_name, "orthologous complexes, based on full genome-coverage complexes in Table 2.\n")
            cat("\tstep 2.6 has generated figure file 3: Scatter plot - Mcalc of", M_name,"complexes V.S. Mcalc of", P_name, "orthologous complexes, based on full genome-coverage complexes in Table 2.\n",
                "\t\t linear regression summary of Mcalc of", M_name, "of complexes on Mcalc of", P_name, "orthologous complexes:\n",
                file="LOG OF COMPUTATION.txt", append=TRUE)
            capture.output(summary(LM_covserve_MP), file="LOG OF COMPUTATION.txt", append=TRUE)
            
            pdf(paste(paste(paste(paste0("figure file 4. Scatter plot-mono mass of ", M_name), paste0(" proteins v.s. mmono mass of ", P_name)), Sys.Date(), sep='_'), "pdf", sep='.'), width=8, height=8)
            fullcover_data <- filter(tableOfOrthoCORUMComplexes, GenomeCoverage==1)
            prots_M <- levels(factor(fullcover_data$Mammalian_Prot_ID))
            plant_mass_each_Mprot <- list()
            for (p in prots_M){
              plantid <- MprotToPprotMap[[p]]
              plant_mass <- c()
              for (i in plantid){
                plant_mass <- c(plant_mass, plantMmono.Data[[i]])
              }
              plant_mass_each_Mprot[[p]] <- mean(plant_mass)
            }
            plant_mass <- c()
            mammal_mass <- c()
            for (p in prots_M){
              mammal_mass <- c(mammal_mass, filter(fullcover_data, Mammalian_Prot_ID==p)$`Mass_kDa`[1])
              plant_mass <- c(plant_mass, plant_mass_each_Mprot[[p]])
            }
            limit_xy <- max(c(max(mammal_mass),max(plant_mass)))
            l <- limit_xy %/% 10
            if (l>100){
              length_here <- 100
            }else if (l<=100 && l>50){
              length_here <- 50
            }else{
              length_here <- 25
            }
            p_m_mass <- data.frame(plant_mass=plant_mass,mammal_mass=mammal_mass)
            print(ggplot(p_m_mass, aes(x=plant_mass, y=mammal_mass) ) + geom_point(size=1) + geom_smooth(method=lm, size=0.5, alpha=0.6)+
                    labs(x=paste0(paste0("Mmono of subunits in ", P_name), " orthocomplexes (kDa)"), y = paste0(paste0("Mmono of subunits in ", M_name), " (kDa)"))+
                    theme(plot.title=element_text(hjust=1, vjust=0.5, margin=margin(t=40,b=-30)), title =element_text(size=8),
                          axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                          axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))+
                    scale_x_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy))+
                    scale_y_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy)))
            dev.off()
            mammal_monomass <- mammal_mass
            plant_monomass <- plant_mass
            LM_mono <- lm(mammal_monomass ~ plant_monomass)
            cat("\tstep 2.7 has generated figure file 4: Scatter plot - Mmono masses of", M_name,"proteins V.S. mmono masses of", P_name, "proteins, based on full genome-coverage complexes in Table 2.\n")
            cat("\tstep 2.7 has generated figure file 4: Scatter plot - Mmono masses of", M_name,"proteins V.S. mmono masses of", P_name, "proteins, based on full genome-coverage complexes in Table 2.\n",
                "\t\t linear regression summary of mmono masses of", M_name, "proteins on", P_name, "proteins:\n",
                file="LOG OF COMPUTATION.txt", append=TRUE)
            capture.output(summary(LM_mono), file="LOG OF COMPUTATION.txt", append=TRUE)

            # plot Mapp human v.s. Mcalc rice
            Mclac_c_compare_P <- c()
            Mapp_c_compare_P <- c()
            data_comparision_Mapp_Mcalc <- NULL
            
            for (c in orthoComplexesFullGcoverage){
              ortholog_set <- unname(unlist(MprotToPprotMap.inEachCORUMComplex[[c]]))
              Mappavg_c <- c()
              OrthologID_c <- c()
              Mapp_detected <- c()
              for (p in ortholog_set){
                if (p %in% resolvablePeaks.FromSinglePeakProt.Chosen){
                  OrthologID_c <- c(OrthologID_c, p)
                  Mappavg_p <- RappMappAvg.resolvablePeaks.FromSinglePeakProt.Chosen[[p]][['Mapp_mean']]
                  if(is.null(Mappavg_p)){
                    Mappavg_c <- c(Mappavg_c, "not detected by the chosen peaks")
                  }else{
                    Mappavg_c <- c(Mappavg_c, Mappavg_p)
                    Mapp_detected <- c(Mapp_detected, Mappavg_p)
                  }
                }else if (p %in% multiplePeakProts.Chosen){
                  peaks_p <- multiplePeakProtsToPeaksMap[[p]]
                  for (peak in peaks_p){
                    OrthologID_c <- c(OrthologID_c, peak)
                    Mappavg_p <- RappMappAvg.resolvablePeaks.FromMultiplePeakProt.Chosen[[peak]][['Mapp_mean']]
                    Mappavg_c <- c(Mappavg_c, Mappavg_p)
                    Mapp_detected <- c(Mapp_detected, Mappavg_p)
                  }
                }else{
                  OrthologID_c <- c(OrthologID_c, p)
                  Mappavg_c <- c(Mappavg_c, "not detected by the chosen peaks")
                }
                if (!is.null(Mapp_detected)){
                  Mclac_c <- Mcalc_c_fullcoverage[[c]]
                  Mapp_c_compare_P <- c(Mapp_c_compare_P, Mapp_detected)
                  Mclac_c_compare_P <- c(Mclac_c_compare_P, rep(Mclac_c, length(Mapp_detected)))
                }
              }
              ComplexName_c <-rep(c, length(OrthologID_c))
              Mcalc_c <- rep(Mcalc_c_fullcoverage[[c]], length(OrthologID_c))
              df_t <- data.frame(ComplexName=ComplexName_c, OrthologID=OrthologID_c, `Mapp.avg of subuints`=Mappavg_c, `Mclc of complex`=Mcalc_c)
              data_comparision_Mapp_Mcalc <- rbind(data_comparision_Mapp_Mcalc, df_t)
            }
            
            mcalc_mapp_compare_P <- data.frame(Mclac_c_compare_P=Mclac_c_compare_P, Mapp_c_compare_P=Mapp_c_compare_P)
            limit_x <- max(Mclac_c_compare_P)
            limit_y <- max(Mapp_c_compare_P)
            lx <- limit_x %/% 10
            if (lx>100){
              length_x <- 100
            }else if (lx<=100 && lx>50){
              length_x <- 50
            }else{
              length_x <- 25
            }
            ly <- limit_y %/% 10
            if (ly>100){
              length_y <- 100
            }else if (lx<=100 && ly>50){
              length_y <- 50
            }else{
              length_y <- 25
            }
            
            length_here <- max(c(length_x, length_y))
            p_compare <- ggplot(mcalc_mapp_compare_P, aes(Mclac_c_compare_P, Mapp_c_compare_P)) + geom_point(size=1)+
              labs(x=paste0(paste0("Mcalc of ", P_name), " orthocomplexes (kDa)"), y = paste0(paste0("Mapp.avg of subunits in ", P_name), " orthocomplexes (kDa)"))+
              theme(plot.title=element_text(hjust=1, vjust=0.5, margin=margin(t=40,b=-30)), title =element_text(size=8),
                    axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                    axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))+
              scale_x_continuous(breaks=seq(0,limit_x,length_here),limits = c(0,limit_x))+
              scale_y_continuous(breaks=seq(0,limit_y,length_here),limits = c(0,limit_y))
            pdf(paste(paste(paste0(paste0(paste0("figure file 5. Scatter plot - Mcalc of ", P_name), " orthocomplexes"),paste0(paste0("Mapp.avg of subunits in ", P_name), " orthocomplexes")), Sys.Date(), sep='_'), "pdf", sep='.'), width=8, height=8*limit_y/limit_x)
            print(p_compare)
            dev.off()
            
            cat("\tstep 2.8 has generated figure file 5: Scatter plot - Mcalc of", P_name,"orthocomplexes V.S.Mapp.avg of subunits in", P_name, "orthocomplexes, based on full genome-coverage complexes in Table 2.\n")
            cat("\tstep 2.8 has generated figure file 5: Scatter plot - Mcalc of", P_name,"orthocomplexes V.S.Mapp.avg of subunits in", P_name, "orthocomplexes, based on full genome-coverage complexes in Table 2.\n",
                "\t\t The Mapp.avg values are from the chosen peaks\n",
                file="LOG OF COMPUTATION.txt", append=TRUE)
            
            # find super nonduplicate complex for small world analysis
            Mprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog <- list()
            for (c in orthoComplexesAboveGcoverageCutoff.withMultiorthologs){
              u_c <- levels(factor(names(MprotToPprotMap.inEachCORUMComplex[[c]])))
              Mprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog[[c]] <- u_c
            }
            super_goodcoverage_complex_basedon_Mprots <- find_super_sub_nonduplicate_set(x=Mprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog, evaluation_data=genomeCoverage, cutoff=0, decision="NA")[[1]]
            sub_goodcoverage_complex_basedon_Mprots <- find_super_sub_nonduplicate_set(x=Mprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog, evaluation_data=genomeCoverage, cutoff=0, decision="NA")[[2]]
            super_goodcoverage_complex_basedon_Mprots_nonduplicate <- find_super_sub_nonduplicate_set(x=Mprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog, evaluation_data=genomeCoverage, cutoff=0, decision="NA")[[3]]
            super_sub_goodcoverage_complex_basedon_Mprots <- list(super_goodcoverage_complex_basedon_Mprots=super_goodcoverage_complex_basedon_Mprots,
                                                                  sub_goodcoverage_complex_basedon_Mprots=sub_goodcoverage_complex_basedon_Mprots,
                                                                  super_goodcoverage_complex_basedon_Mprots_nonduplicate=super_goodcoverage_complex_basedon_Mprots_nonduplicate)
            
            # find super-sets and remove duplicates according to P-prots
            #From output produced by step 1.1, searching for super and sub orthologous complexes in terms of experiment-detected plant proteins 
            detected_plantprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog <- list()
            for (c in orthoComplexesAboveGcoverageCutoff.withMultiorthologs){
              plantprot_c <- levels(factor(unname(unlist(MprotToPprotMap.inEachCORUMComplex[[c]]))))
              
              detected_prots <- plantprot_c[plantprot_c %in% substr(protsSelectedForSmallWorld,1,ncharProtID)]
              if(length(detected_prots)>1){
                detected_plantprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog[[c]] <- detected_prots
              }
            }
            
            find.supper_sub_noduplicate_basedon_detected_Pprots <- 
              find_super_sub_nonduplicate_set(x=detected_plantprot_in_Mcomplex_goodGenomcover_with_multiOrhtolog, evaluation_data=genomeCoverage, cutoff=0, decision="NA")
            super_goodcoverage_complex_basedon_detected_Pprots <- find.supper_sub_noduplicate_basedon_detected_Pprots[[1]]
            sub_goodcoverage_complex_basedon_detected_Pprots <- find.supper_sub_noduplicate_basedon_detected_Pprots[[2]]
            super_goodcoverage_complex_basedon_detected_Pprots_nonduplicate <- find.supper_sub_noduplicate_basedon_detected_Pprots[[3]]
            super_sub_goodcoverage_complex_basedon_detected_Pprots <- list(super_goodcoverage_complex_basedon_detected_Pprots=super_goodcoverage_complex_basedon_detected_Pprots,
                                                                           sub_goodcoverage_complex_basedon_detected_Pprots=sub_goodcoverage_complex_basedon_detected_Pprots,
                                                                           super_goodcoverage_complex_basedon_detected_Pprots_nonduplicate=super_goodcoverage_complex_basedon_detected_Pprots_nonduplicate)
            cat("\tstep 2.9 has generated super-sets of orthologous complexes and removed the duplicates according to the", P_name, "protein peaks chosen for small world analysis.\n",
                "\t\t 1) There are", length(super_goodcoverage_complex_basedon_Mprots_nonduplicate), 
                "out of", length(orthoComplexesAboveGcoverageCutoff.withMultiorthologs), "complexes produced in step 2.3 that are super-nonduplicate complexes in terms of", P_name, "orthologs.\n",
                "\t\t 2) There are", length(super_goodcoverage_complex_basedon_detected_Pprots_nonduplicate), 
                "out of", length(super_goodcoverage_complex_basedon_Mprots_nonduplicate), "complexes produced in the step above that are super-nonduplicate complexes in terms of", P_name, "peaks chosen for small world anlysis.\n",
                "\t\t 3) These", length(super_goodcoverage_complex_basedon_detected_Pprots_nonduplicate), 
                "complexes are CORUM complexex to which small world clustering algorithm will be applied in next step.\n",
                file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 2: Done. The orthologous CORUM complexexes are ready for small world analysis.\n\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 2: Done.  The orthologous CORUM complexexes are ready for small world analysis.\n\n")
            cat("Step 3: To implement small world clustering algorithm.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 3: To implement small world clustering algorithm.\n")
            
            small_world_analysis_output <- small_world_analysis(names_CORUM_complexes=super_goodcoverage_complex_basedon_detected_Pprots_nonduplicate, prot_peak_map=multiplePeakProtsToPeaksMap,
                                                                ortholog_map=MprotToPprotMap.inEachCORUMComplex, 
                                                                prots_for_small_world=protsSelectedForSmallWorld,
                                                                Gaussian_peak_fraction_1=GaussianPeakFractionB1.Chosen, Gaussian_peak_fraction_2=GaussianPeakFractionB2.Chosen,
                                                                Matrix_SEC_B1=peakProfileMatrixB1.Chosen, Matrix_SEC_B2=peakProfileMatrixB2.Chosen,
                                                                parameters_SOM=parameters_SOM, parameters_APC=parameters_APC, cross_correlated_r=cross_correlated_r)
            
            subcomplex_prediction <- small_world_analysis_output$subcomplex_prediction
            subcomplex_cutoff <- small_world_analysis_output$subcomplex_cutoff
            subcomplex_overlap <- small_world_analysis_output$subcomplex_overlap
            cat("\tStep 3.1 has finished the small world analysis and the bootstrap P-value calculation.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.1 has finished the small world analysis and the bootstrap P-value calculation.\n")
            
            # summarize the result of small world analysis according to the algorithm
            #searching for the correctly predicted monomers and ortho-paralog-monomers.
            count_sub <- 0
            count_sub_som <- 0
            count_sub_som_Rappsmall <- 0
            count_sub_only2 <- 0
            count_sub_only2_Rappsmall <- 0
            count_homolog <- 0
            count_sub_all <- 0
            count_sub_all_Rappsmall <- 0
            count_single <- 0
            count_homolog_Rappsmall <- 0
            count_single_Rappsmall <- 0
            count_homolog_Rapp_allsmall <- 0
            count_som_Rapp_allsmall <- 0
            monomers <- c()
            comsub_monomers <- c()
            homologs_monomer_found <- c()
            Rappavg_calc <- function(x){
              if (nchar(x)==ncharProtID){
                return(RappMappAvg.resolvablePeaks.FromSinglePeakProt.Chosen[[x]][['Rapp_mean']])
              }else{
                return(RappMappAvg.resolvablePeaks.FromMultiplePeakProt.Chosen[[x]][['Rapp_mean']])
              }
            }
            for (c in names(subcomplex_prediction)){
              sub_names <- names(subcomplex_prediction[[c]])
              sub_names_sep <- strsplit(sub_names, split='_')
              count_sub <- count_sub + length(sub_names)
              for (i in 1:length(sub_names)){
                if ("sub" %in% sub_names_sep[[i]] & !("orthoparalog" %in% sub_names_sep[[i]])){
                  count_sub_som <- count_sub_som + 1
                  som_set <- subcomplex_prediction[[c]][[i]]
                  for (p in som_set){
                    Rappavg <- Rappavg_calc(x=p)
                    if(Rappavg<=RappCutoffNonmono){
                      count_sub_som_Rappsmall <- count_sub_som_Rappsmall + 1
                      break
                    }
                  }
                  all_sub_check <- 0
                  for (p in som_set){
                    Rappavg <- Rappavg_calc(x=p)
                    if(Rappavg<=RappCutoffNonmono){
                      all_sub_check <- all_sub_check + 1
                    }
                  }
                  if (all_sub_check == length(som_set)){
                    count_som_Rapp_allsmall <- count_som_Rapp_allsmall + 1
                  }
                }else if("orthoparalog" %in% sub_names_sep[[i]]){
                  count_homolog <- count_homolog + 1
                  homolog_set <- subcomplex_prediction[[c]][[i]]
                  for (p in homolog_set){
                    Rappavg <- Rappavg_calc(x=p)
                    if(Rappavg<=RappCutoffNonmono){
                      count_homolog_Rappsmall <- count_homolog_Rappsmall + 1
                      break
                    }
                  }
                  all_homo_check <- 0
                  for (p in homolog_set){
                    Rappavg <- Rappavg_calc(x=p)
                    if(Rappavg<=RappCutoffNonmono){
                      all_homo_check <- all_homo_check + 1
                    }
                  }
                  if (all_homo_check == length(homolog_set)){
                    count_homolog_Rapp_allsmall <- count_homolog_Rapp_allsmall + 1
                    homologs_monomer_found <- c(homologs_monomer_found, paste(c,sub_names[i],sep=":"))
                  }
                }else if("only" %in% sub_names_sep[[i]]){
                  count_sub_only2 <- count_sub_only2 + 1
                  only2_set <- subcomplex_prediction[[c]][[i]]
                  for (p in only2_set){
                    Rappavg <- Rappavg_calc(x=p)
                    if(Rappavg<=RappCutoffNonmono){
                      count_sub_only2_Rappsmall <- count_sub_only2_Rappsmall + 1
                      break
                    }
                  }
                }else if ("all" %in% sub_names_sep[[i]]){
                  count_sub_all <- count_sub_all + 1
                  all_set <- subcomplex_prediction[[c]][[i]]
                  for (p in all_set){
                    Rappavg <- Rappavg_calc(x=p)
                    if(Rappavg<=RappCutoffNonmono){
                      count_sub_all_Rappsmall <- count_sub_all_Rappsmall + 1
                    }
                  }
                }else if ("singleton" %in% sub_names_sep[[i]]){
                  count_single <- count_single + 1
                  single_set <- subcomplex_prediction[[c]][[i]]
                  for (p in single_set){
                    Rappavg <- Rappavg_calc(x=p)
                    if(Rappavg<=RappCutoffNonmono){
                      count_single_Rappsmall <- count_single_Rappsmall + 1
                      comsub_monomers <- c(comsub_monomers, paste(paste(c,sub_names[i],sep="-"), p, sep=": "))
                      monomers <- c(monomers, p)
                    }
                  }
                }
              }
            }
            predicted_monomers <- list(monomers=monomers, ortho_paralog_monomers=homologs_monomer_found)
            cat("\tStep 3.2 has found the correctly predicted monomers and ortho-paralog-monomers.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.2 has found the correctly predicted monomers and ortho-paralog-monomers.\n")
        
            # calculate Mcalc and Mapp for the predicted subcomplexes, calculating Mcalc and Mapp for the predicted subcomplexes.
            Mcalc_sub <- list()
            Mapp_mean_sub <- list()
            Mappavg_calc <- function(x){
              if (nchar(x)==ncharProtID){
                return(RappMappAvg.resolvablePeaks.FromSinglePeakProt.Chosen[[x]][['Mapp_mean']])
              }else{
                return(RappMappAvg.resolvablePeaks.FromMultiplePeakProt.Chosen[[x]][['Mapp_mean']])
              }
            }
            for (c in names(subcomplex_prediction)){
              for (sub in names(subcomplex_prediction[[c]])){
                Mcalc_sub[[c]][[sub]] <- Mcalc_calculation(complex_name=c, Mcomplex_to_Plant_mapping=MprotToPprotMap.inEachCORUMComplex[[c]], 
                                                           ortholog_set=substr(subcomplex_prediction[[c]][[sub]],1,ncharProtID),
                                                           Mmnono_data=plantMmono.Data,
                                                           stoichiometry_data=mammalCORUM.StoichiometryData)$Mcalc
                Mapp_sub <- c()
                for (p in subcomplex_prediction[[c]][[sub]]){
                  Mapp_sub <- c(Mapp_sub, Mappavg_calc(x=p))
                }
                Mapp_mean_sub[[c]][[sub]] <- mean(Mapp_sub)
              }
            }
            
            Mcalc_c <- list()
            for (c in names(MprotToPprotMap.inEachCORUMComplex)){
              orthologs_c <- unname(unlist(MprotToPprotMap.inEachCORUMComplex[[c]]))
              orthologs_c <- orthologs_c[!is.na(orthologs_c)]
              Mcalc_c[[c]] <- Mcalc_calculation(complex_name=c, Mcomplex_to_Plant_mapping=MprotToPprotMap.inEachCORUMComplex[[c]], 
                                                ortholog_set=orthologs_c,
                                                Mmnono_data=plantMmono.Data,
                                                stoichiometry_data=mammalCORUM.StoichiometryData)$Mcalc
            }
            
            # generate data table of all information of small world analysis
            subcomplex_prediction_alldata <- NULL
            for (c in names(subcomplex_prediction)){
              subunits <- c()
              plant_id <- c()
              mapp_m <- c()
              mcalc_m <- c()
              mcalc_p <- c()
              mcalc_p_sub <- c()
              prot_name <- c()
              gene_name <- c()
              mmono <- c()
              Mapp1_peak <- c()
              Mapp2_peak <- c()
              Rapp1_peak <- c()
              Rapp2_peak <- c()
              sub_id <- c()
              p_v <- c()
              plantid_identified <- c()
              Mapp_Subcomplex <- c()
              
              for (i in names(subcomplex_prediction[[c]])){
                u_space_sub_i <- c()
                for (p_id in subcomplex_prediction[[c]][[i]]){
                  u_set <- c()
                  for (u in names(MprotToPprotMap.inEachCORUMComplex[[c]])){
                    if (substr(p_id,1,ncharProtID) %in% MprotToPprotMap.inEachCORUMComplex[[c]][[u]]){
                      u_set <- c(u_set, u)
                      mapp_m <- c(mapp_m, filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$`Mass_kDa`[1])
                      mcalc_m <- c(mcalc_m, filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$`Mcalc_kDa`[1])
                      prot_name <- c(prot_name, filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$Protein_names[1])
                      gene_name <- c(gene_name, filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$Gene_names[1])
                    }
                  }
                  u_space_sub_i <- c(u_space_sub_i, u_set)
                  plantid_identified <- c(plantid_identified, p_id)
                  plant_id <- c(plant_id, rep(p_id, length(u_set)))
                  if (nchar(p_id)==ncharProtID){
                    Mapp1_peak <- c(Mapp1_peak, rep(reproduciblePeakFeatures.Data[[p_id]][['B1']][['Mapp1']], length(u_set)))
                    Mapp2_peak <- c(Mapp2_peak, rep(reproduciblePeakFeatures.Data[[p_id]][['B2']][['Mapp1']], length(u_set)))
                    Rapp1_peak <- c(Rapp1_peak, rep(reproduciblePeakFeatures.Data[[p_id]][['B1']][['Rapp1']], length(u_set)))
                    Rapp2_peak <- c(Rapp2_peak, rep(reproduciblePeakFeatures.Data[[p_id]][['B2']][['Rapp1']], length(u_set)))
                  }
                  if (nchar(p_id)>ncharProtID){
                    realp_id <- substr(p_id, 1, ncharProtID)
                    peak_b1 <- substr(p_id, ncharProtID+2, ncharProtID+2)
                    peak_b2 <- substr(p_id, ncharProtID+3, ncharProtID+3)
                    Mapp1_peak <- c(Mapp1_peak, rep(reproduciblePeakFeatures.Data[[realp_id]][['B1']][[paste0('Mapp', peak_b1)]], length(u_set)))
                    Mapp2_peak <- c(Mapp2_peak, rep(reproduciblePeakFeatures.Data[[realp_id]][['B2']][[paste0('Mapp', peak_b2)]], length(u_set)))
                    Rapp1_peak <- c(Rapp1_peak, rep(reproduciblePeakFeatures.Data[[realp_id]][['B1']][[paste0('Rapp', peak_b1)]], length(u_set)))
                    Rapp2_peak <- c(Rapp2_peak, rep(reproduciblePeakFeatures.Data[[realp_id]][['B2']][[paste0('Rapp', peak_b2)]], length(u_set)))
                  }
                }
                sub_id <- c(sub_id, rep(i, length(u_space_sub_i)))
                mcalc_p_sub <- c(mcalc_p_sub, rep(Mcalc_sub[[c]][[i]], length(u_space_sub_i)))
                Mapp_Subcomplex <- c(Mapp_Subcomplex, rep(Mapp_mean_sub[[c]][[i]], length(u_space_sub_i)))
                if (length(subcomplex_cutoff[[c]][[i]])>1){
                  p_v <- c(p_v, rep(subcomplex_cutoff[[c]][[i]][3], length(u_space_sub_i)))
                }else{
                  p_v <- c(p_v, rep(NA, length(u_space_sub_i)))
                }
                subunits <- c(subunits, u_space_sub_i)          
              }
              u_space_c <- names(MprotToPprotMap.inEachCORUMComplex[[c]])
              u_comp <- c()
              for (u in u_space_c){
                pid_u <- MprotToPprotMap.inEachCORUMComplex[[c]][[u]]
                pid_not_identified_in_u <- pid_u[!(pid_u %in% substr(plantid_identified,1,ncharProtID))]
                if (length(pid_not_identified_in_u)>0){
                  plant_id <- c(plant_id, pid_not_identified_in_u)
                  mapp_m <- c(mapp_m, rep(filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$Mass_kDa[1], length(pid_not_identified_in_u)))
                  mcalc_m <- c(mcalc_m, rep(filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$Mcalc_kDa[1], length(pid_not_identified_in_u)))
                  prot_name <- c(prot_name, rep(filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$Protein_names[1], length(pid_not_identified_in_u)))
                  gene_name <- c(gene_name, rep(filter(mammalCORUM.Data, ComplexName==c & Mammalian_Prot_ID==u)$Gene_names[1], length(pid_not_identified_in_u)))
                  u_comp <- c(u_comp, rep(u, length(pid_not_identified_in_u)))
                }
                
              }
              
              Mapp1_peak <- c(Mapp1_peak, rep(NA, length(u_comp)))
              Mapp2_peak <- c(Mapp2_peak, rep(NA, length(u_comp)))
              Rapp1_peak <- c(Rapp1_peak, rep(NA, length(u_comp)))
              Rapp2_peak <- c(Rapp2_peak, rep(NA, length(u_comp)))
              
              sub_id <- c(sub_id, rep(NA, length(u_comp)))
              p_v <- c(p_v, rep(NA, length(u_comp)))
              subunits <- c(subunits, u_comp) 
              
              com_name <- rep(c, length(subunits))
              num_subnuits <- rep(CORUMComplexSize[[c]], length(subunits))
              coverage <- rep(genomeCoverage[[c]], length(subunits))
              mcalc_p <- rep(Mcalc_c[[c]], length(subunits))
              mcalc_p_sub <- c(mcalc_p_sub, rep(NA, length(u_comp)))
              Mapp_Subcomplex <- c(Mapp_Subcomplex, rep(NA, length(u_comp)))
              mapp_avg_p <- rep(NA, length(subunits))
              plant_prot_name <- rep(NA, length(subunits))
              pids <- levels(factor(unname(unlist(MprotToPprotMap.inEachCORUMComplex[[c]]))))
              num_plant_pid <- rep(length(pids), length(subunits))
              mmono <- rep(NA, length(subunits))
              
              df_t <- data.frame(ComplexName=com_name, Number_of_mammalian_subunits_in_a_given_complex=num_subnuits, Coverage=coverage, 
                                 Mammalian_Prot_ID=subunits, Mammalian_Protein_Name=prot_name, 
                                 Mammalian_Gene_Name=gene_name, Mass_kDa=mapp_m, Mcalc_of_CORUM_Complex_kDa=mcalc_m, 
                                 Plant_Ortholog_ID=plant_id, Plant_Protein_Name=plant_prot_name, 
                                 Number_of_plant_subunits_in_a_given_complex=num_plant_pid,Mmono_kDa=mmono, 
                                 Mcalc_of_Plant_Complex_kDa=mcalc_p, Mapp_avg_kDa=mapp_avg_p,
                                 Mapp1_peak=Mapp1_peak, Mapp2_peak=Mapp2_peak, Rapp1_peak=Rapp1_peak, Rapp2_peak=Rapp2_peak,          
                                 SubcomplexID=sub_id,  P_value=p_v, Mcalc_Subcomplex=mcalc_p_sub, Mapp_Subcomplex=Mapp_Subcomplex 
              )
              subcomplex_prediction_alldata <- rbind(subcomplex_prediction_alldata, df_t)
            }
            
            putative_stoichiometry <- c()
            reference_stoichiometry <-c()
            for(i in 1:dim(subcomplex_prediction_alldata)[1]){
              c <- subcomplex_prediction_alldata[i, "ComplexName"]
              u <- subcomplex_prediction_alldata[i, "Mammalian_Prot_ID"]
              if (c %in% names(mammalCORUM.StoichiometryData)){
                if (u %in% names(mammalCORUM.StoichiometryData[[c]])){
                  putative_stoichiometry <- c(putative_stoichiometry, mammalCORUM.StoichiometryData[[c]][[u]][["stoichiometry"]])
                  reference_stoichiometry <-c(reference_stoichiometry, mammalCORUM.StoichiometryData[[c]][[u]][["reference"]])
                }else{
                  putative_stoichiometry <- c(putative_stoichiometry, NA)
                  reference_stoichiometry <-c(reference_stoichiometry, NA)
                }
              }else{
                putative_stoichiometry <- c(putative_stoichiometry, NA)
                reference_stoichiometry <-c(reference_stoichiometry, NA)
              } 
            }

            for(i in 1:dim(subcomplex_prediction_alldata)[1]){
              loc_i <- subcomplex_prediction_alldata$Plant_Ortholog_ID[i]
              if (! is.na(loc_i)){
                if (substr(loc_i,1,1)==substr(mostCommonProtIDExample,1,1)){
                  loc <- substr(loc_i, 1,ncharProtID)
                }else{
                  loc <- loc_i
                }
              }
              
              if (! is.na(loc)){
                if (length(plantMmono.Data[[loc]])<1){
                  warning(paste(paste("No mono mass for", loc, sep = " "), "and use 0 as mono mass", sep = " "))
                  print(paste(paste("No mono mass for", loc, sep = " "), "and use 0 as mono mass", sep = " "))
                  mono_t <- 0
                }else{
                  mono_t <- plantMmono.Data[[loc]]
                }
                subcomplex_prediction_alldata[i, 'Mmono_kDa'] <- mono_t
                if (length(plantProtnames.Data[[loc]])<1){
                  warning(paste(paste("No protname for", loc, sep = " "), "and use NA as protname", sep = " "))
                  print(paste(paste("No protname for", loc, sep = " "), "and use NA as protname", sep = " "))
                  protname_t <- NA
                }else{
                  protname_t <- plantProtnames.Data[[loc]]
                }
                subcomplex_prediction_alldata[i, 'Mmono_kDa'] <- mono_t
                subcomplex_prediction_alldata[i, 'Plant_Protein_Name'] <- protname_t
                subcomplex_prediction_alldata[i, 'Mapp_avg_kDa'] <- (subcomplex_prediction_alldata[i, 'Mapp1_peak']+subcomplex_prediction_alldata[i, 'Mapp2_peak'])/2
              }
            }
            subcomplex_prediction_alldata <- mutate(subcomplex_prediction_alldata,putative_stoichiometry=putative_stoichiometry,reference_stoichiometry=reference_stoichiometry)
            write.csv(subcomplex_prediction_alldata, paste(paste("Table 3 - Subcomplex prediction output of all relevant data", Sys.Date(), sep='_'), "csv", sep='.'), row.names = F)
            cat("\tStep 3.3 has generated Table 3 : Subcomplex prediction output of all relevant data.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.3 has generated Table 3 : Subcomplex prediction output of all relevant data.\n")
            
            # find non-duplicate subcomplexes with p-value<=0.05
            find_certain_p_v_subcomplexes_no_duplicate <- function(subcomplex_prediction, subcomplex_cutoff){
              count <- 0
              good_subcomplex_pred <- list()
              good_subcomplex_comsub <- list()
              for(c in names(subcomplex_prediction)){
                for(sub in names(subcomplex_prediction[[c]])){
                  if (sub %in% names(subcomplex_cutoff[[c]])){
                    P_v <- subcomplex_cutoff[[c]][[sub]][3]
                    if(P_v <= 0.05){
                      good_subcomplex_pred[[c]][[sub]][['LOC']] <- subcomplex_prediction[[c]][[sub]]
                      good_subcomplex_pred[[c]][[sub]][['P_value']] <- P_v
                      good_subcomplex_comsub[[paste(c, sub, sep=':')]] <- subcomplex_prediction[[c]][[sub]]
                      count <- count + 1
                    }
                  }
                }
              }
              
              good_subcomplex_name_set <- names(good_subcomplex_comsub)
              count <- 0
              duplicate_of_subcomplex <- list()
              name_set <- good_subcomplex_name_set
              while(length(name_set)>0){
                count <- count + 1
                name_1 <- name_set[1]
                temp_set <- c(name_1)
                for (n in name_set){
                  if (! n==name_1){
                    if (all(c(sum(good_subcomplex_comsub[[name_1]] %in% good_subcomplex_comsub[[n]])==length(good_subcomplex_comsub[[name_1]]),
                              length(good_subcomplex_comsub[[name_1]]==length(good_subcomplex_comsub[[n]]))))){
                      temp_set <- c(temp_set, n)
                      
                    }
                  }
                }
                name_set <- name_set[! name_set %in% temp_set]
                duplicate_of_subcomplex[[count]] <- temp_set
              }
              
              good_subcomplex_selected_noduplicate <- c()
              for (i in 1:length(duplicate_of_subcomplex)){
                good_subcomplex_selected_noduplicate <- c(good_subcomplex_selected_noduplicate, duplicate_of_subcomplex[[i]][1])
              }
              return(good_subcomplex_selected_noduplicate)
            }
            
            subcomplex_names_noduplicates_goodcover <- find_certain_p_v_subcomplexes_no_duplicate(subcomplex_prediction=subcomplex_prediction,
                                                                                                  subcomplex_cutoff=subcomplex_cutoff)
            cat("\tStep 3.4 has filtered out", length(subcomplex_names_noduplicates_goodcover),  "non-duplicate subcomplexes with p-value<=0.05\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.4 has filtered out", length(subcomplex_names_noduplicates_goodcover),  "non-duplicate subcomplexes with p-value<=0.05\n")
            
            
            subcomplex_goldstd <- NULL
            for (c_sub in subcomplex_names_noduplicates_goodcover){
              c <- strsplit(c_sub, split=':')[[1]][1]
              sub <- strsplit(c_sub, split=':')[[1]][2]
              d_t <- filter(subcomplex_prediction_alldata, `ComplexName`==c & `SubcomplexID`==sub)
              subcomplex_goldstd <- rbind(subcomplex_goldstd, d_t)
            }
            
            subcomplex_goldstd <- mutate(subcomplex_goldstd, 
                                         Ratio_Mapp_Mcalc_subcomplex=as.numeric(Mapp_Subcomplex)/as.numeric(Mcalc_Subcomplex))
            subcomplex_goldstd <- mutate(subcomplex_goldstd,Rapp_avg=(Rapp1_peak+Rapp2_peak)/2)
            
            subcomplex_goldstd_Rapplarge <- filter(subcomplex_goldstd, Rapp_avg>RappCutoffNonmono)
            
            com_sub_Rapplarge <- paste(subcomplex_goldstd_Rapplarge$ComplexName, subcomplex_goldstd_Rapplarge$SubcomplexID, sep=':')
            com_sub_Rapplarge <- levels(factor(com_sub_Rapplarge))
            
            subcomplex_goldstd_Rapplarge_nosingleton <- NULL
            com_sub_Rapplarge_nosingleton <- c()
            singleton_withoutmono <- c()
            for (cs in com_sub_Rapplarge){
              c <- strsplit(cs, split = ':')[[1]][1]
              s <- strsplit(cs, split = ':')[[1]][2]
              data_cs <- filter(subcomplex_goldstd_Rapplarge, ComplexName==c, SubcomplexID==s)
              if(dim(data_cs)[1]>1){
                subcomplex_goldstd_Rapplarge_nosingleton <- rbind(subcomplex_goldstd_Rapplarge_nosingleton, data_cs)
                com_sub_Rapplarge_nosingleton <- c(com_sub_Rapplarge_nosingleton, cs)
              }else{
                singleton_withoutmono <- c(singleton_withoutmono, cs)
              }
            }
            cat("\tStep 3.5 has filtered out", length(com_sub_Rapplarge_nosingleton),  "non-duplicate subcomplexes with p-value<=0.05 and without monomer.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.5 has filtered out", length(com_sub_Rapplarge_nosingleton),  "non-duplicate subcomplexes with p-value<=0.05 and without monomer.\n")
            #plot good sub-complex identification result
            Mcalc_sub_nodup <- c()
            Mapp_sub_nodup <- c()
            Mcalc_std <- list()
            Mapp_std <- list()
            for (comsub in com_sub_Rapplarge_nosingleton){
              com <- strsplit(comsub, split = ':')[[1]][1]
              sub <- strsplit(comsub, split = ':')[[1]][2]
              
              loc_set <- filter(subcomplex_goldstd_Rapplarge, ComplexName==com, SubcomplexID==sub)$Plant_Ortholog_ID
              Mcalc_sub_nodup <- c(Mcalc_sub_nodup, Mcalc_calculation(complex_name=com, Mcomplex_to_Plant_mapping=MprotToPprotMap.inEachCORUMComplex[[com]], 
                                                                      ortholog_set=substr(loc_set,1,ncharProtID),
                                                                      Mmnono_data=plantMmono.Data,
                                                                      stoichiometry_data=mammalCORUM.StoichiometryData)$Mcalc)
              Mcalc_std[[comsub]] <- Mcalc_calculation(complex_name=com, Mcomplex_to_Plant_mapping=MprotToPprotMap.inEachCORUMComplex[[com]], 
                                                       ortholog_set=substr(loc_set,1,ncharProtID),
                                                       Mmnono_data=plantMmono.Data,
                                                       stoichiometry_data=mammalCORUM.StoichiometryData)$Mcalc
              Mapp_detected_com <- c()
              for (p in loc_set){
                Mapp_detected_com <- c(Mapp_detected_com, Mappavg_calc(x=p))
              }
              Mapp_std[[comsub]] <- mean(Mapp_detected_com)
              Mapp_sub_nodup <- c(Mapp_sub_nodup, mean(Mapp_detected_com))
            }
            
            rownames(subcomplex_goldstd_Rapplarge) <- 1:dim(subcomplex_goldstd_Rapplarge)[1]
            subcomplex_goldstd_Rapplarge_nosingleton <- NULL
            for (i in 1:dim(subcomplex_goldstd_Rapplarge)[1]){
              cs <- paste(subcomplex_goldstd_Rapplarge[i,"ComplexName"], subcomplex_goldstd_Rapplarge[i,"SubcomplexID"], sep=":")
              if(cs=="NA:NA"){
              }
              if (cs %in% com_sub_Rapplarge_nosingleton){
                d_ttt <- subcomplex_goldstd_Rapplarge[i,]
                d_ttt[1,"Mcalc_Subcomplex"] <- Mcalc_std[[cs]]
                d_ttt[1,"Mapp_Subcomplex"] <- Mapp_std[[cs]]
                subcomplex_goldstd_Rapplarge_nosingleton <- rbind(subcomplex_goldstd_Rapplarge_nosingleton, d_ttt)
              }
              
            }
            subcomplex_goldstd_Rapplarge_nosingleton <- mutate(subcomplex_goldstd_Rapplarge_nosingleton, 
                                                               Ratio_Mapp_Mcalc_subcomplex=as.numeric(Mapp_Subcomplex)/as.numeric(Mcalc_Subcomplex))
            
            write.csv(subcomplex_goldstd_Rapplarge_nosingleton, paste(paste("Table 4 - Gold standard subcomplexes(complexes)", Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
            data_ratiosmall_Rapplarge_nosingleton <- filter(subcomplex_goldstd_Rapplarge_nosingleton,Ratio_Mapp_Mcalc_subcomplex<=2 & 0.5<=Ratio_Mapp_Mcalc_subcomplex)
            data_ratiosmall_Rapplarge_nosingleton <- filter(subcomplex_goldstd_Rapplarge_nosingleton,Ratio_Mapp_Mcalc_subcomplex<=2 & 0.5<=Ratio_Mapp_Mcalc_subcomplex)
            com_sub_ratiosmall_Rapplarge_nosingleton <- 
              levels(factor(paste(data_ratiosmall_Rapplarge_nosingleton$ComplexName, data_ratiosmall_Rapplarge_nosingleton$SubcomplexID,sep=":")))
            cat("\tStep 3.6 has filtered out", length(com_sub_ratiosmall_Rapplarge_nosingleton),  "non-duplicate subcomplexes with p-value<=0.05 and 0.5<=Mapp_Mcalc_Ratio<=2 without monomer.\n",
                "\t\t has generated Table 4: Gold standard subcomplexes/complexes.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.6 has filtered out", length(com_sub_ratiosmall_Rapplarge_nosingleton),  "non-duplicate subcomplexes with p-value<=0.05 and 0.5<=Mapp_Mcalc_Ratio<=2 without monomer.\n",
                "\t\t has generated Table 4: Gold standard subcomplexes/complexes\n.")
            goldstd_complex_names <- c()
            for (i in com_sub_ratiosmall_Rapplarge_nosingleton){
              cs <- strsplit(i,split=":",fixed=T)[[1]]
              goldstd_complex_names <- c(goldstd_complex_names, cs[1])
            }
            goldstd_complex_names <- levels(factor(goldstd_complex_names))

            mcalc_mapp_sub <- data.frame(x = Mcalc_sub_nodup, y=Mapp_sub_nodup)
            mcalc_sub_x <- mcalc_mapp_sub$x
            mapp_sub_y <- mcalc_mapp_sub$y
            filtered_x <- c()
            filtered_y <- c()
            class_ratio <- c()
            filtered_consub_names <- c()
            for (i in 1:length(mcalc_sub_x)){
              if(mapp_sub_y[i]/mcalc_sub_x[i]<=2 & mapp_sub_y[i]/mcalc_sub_x[i]>=0.5){
                filtered_x <- c(filtered_x, mcalc_sub_x[i])
                filtered_y <- c(filtered_y, mapp_sub_y[i])
                class_ratio <- c(class_ratio, "0.5 <= Ratio <= 2")
                filtered_consub_names <- c(filtered_consub_names, subcomplex_names_noduplicates_goodcover[i])
              }else{
                class_ratio <- c(class_ratio, "Ratio out of [0.5, 2]")
              }
            }
            ratio <- mapp_sub_y/mcalc_sub_x
            mcalc_mapp_sub_class <- mutate(mcalc_mapp_sub, Class=class_ratio)
            data_tt <- filter(mcalc_mapp_sub_class, Class=="0.5 <= Ratio <= 2")
            df_abline <- data.frame(intercept=c(0,0,0),slope=c(max(data_tt$y/data_tt$x), 1, min(data_tt$y/data_tt$x)),linetype=factor(c(1,3,5)))
            limit_xy <- max(c(max(mcalc_mapp_sub_class$x),max(mcalc_mapp_sub_class$y)))
            l <- limit_xy %/% 10
            if (l>100){
              length_here <- 100
            }else if (l<=100 && l>50){
              length_here <- 50
            }else{
              length_here <- 25
            }
            pdf(paste(paste(paste0(paste0(paste0(paste0(paste0("figure file 6. Scatter plot - Mcalc of ", P_name), " subcomplexes"), " V.S. Mapp-avg of subunits in "), P_name), " subcomplexes"), Sys.Date(), sep='_'), "pdf", sep='.'), width=8, height=8)
            print(ggplot(mcalc_mapp_sub_class, aes(x, y)) + geom_point(shape=21,aes(fill= Class),  size = 3,alpha = 0.6)+
              geom_abline(data=df_abline, aes(intercept=intercept,slope=slope, linetype=linetype))+
              labs(x=paste0(paste0("Mcalc of ", P_name), " subcomplexes (kDa)"), y = paste0(paste0("Mapp-avg of subunits in ", P_name), " subcomplexes (kDa)"))+
              theme(plot.title=element_text(hjust=1, vjust=0.5, margin=margin(t=40,b=-30)), title =element_text(size=10),
                    axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                    axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))+
              theme(legend.title = element_text(size = 10,face = 10),legend.text=element_text(size=10),legend.position = "bottom")+
              scale_x_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy))+
              scale_y_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy))+
              scale_linetype_manual(name = "Separation line",
                                    values=c(2,1,3),labels=c(paste0("slpoe=",round(max(data_tt$y/data_tt$x),2)),paste0("slpoe=",1),paste0("slpoe=",round(min(data_tt$y/data_tt$x),2))))
            )
            dev.off()
            cat("\tStep 3.7 has generated figure file 6: Scatter plot - Mcalc of", P_name, "subcomplexes V.S. Mapp-avg of subunits in", P_name, "subcomplexes.\n",
                "\t\t Theis figure is supported by Table 4.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.7 has generated figure file 6: Scatter plot - Mcalc of", P_name, "subcomplexes V.S. Mapp-avg of subunits in", P_name, "subcomplexes.\n",
                "\t\t Theis figure is supported by Table 4.\n")
            
            pdf(paste(paste("figure file 7. The subcomplex Clustering identification for CORUM complexes - Mapp vs Mcalc", Sys.Date(), sep='_'), "pdf", sep='.'))
            for (c in names(subcomplex_prediction)){
              data_plot <- data.frame(prot_id_name=NA, Cluster.ID_P.value=NA, Mcalc=NA, Mapp=NA, Rapp=NA)
              for (sub in names(subcomplex_prediction[[c]])){
                loc_set <- subcomplex_prediction[[c]][[sub]]
                sub_pv_t <- rep(paste(sub,format(subcomplex_cutoff[[c]][[sub]][[3]], scientific = T,digits=3),sep=", "), length(loc_set))
                Mcalc_t <- rep(Mcalc_sub[[c]][[sub]], length(loc_set))
                Mapp_t <- c()
                Rapp_t <- c()
                prot_id_name_t <- c()
                for (p in loc_set){
                  if (nchar(p)>ncharProtID){
                    real_p <- substr(p,1,ncharProtID)
                    peak_id_1 <- substr(p,ncharProtID+2,ncharProtID+2)
                    peak_id_2 <- substr(p,ncharProtID+3,ncharProtID+3)
                    Mapp_p <- (reproduciblePeakFeatures.Data[[real_p]][['B1']][[paste0("Mapp",peak_id_1)]] + 
                                 reproduciblePeakFeatures.Data[[real_p]][['B2']][[paste0("Mapp",peak_id_2)]])/2
                    Rapp_p <- (reproduciblePeakFeatures.Data[[real_p]][['B1']][[paste0("Rapp",peak_id_1)]]+
                                 reproduciblePeakFeatures.Data[[real_p]][['B2']][[paste0("Rapp",peak_id_2)]])/2
                    Mapp_t <- c(Mapp_t, Mapp_p)
                    Rapp_t <- c(Rapp_t, Rapp_p)
                  }else{
                    Mapp_p <- (reproduciblePeakFeatures.Data[[p]][['B1']][["Mapp1"]] + 
                                 reproduciblePeakFeatures.Data[[p]][['B2']][["Mapp1"]])/2
                    Rapp_p <- (reproduciblePeakFeatures.Data[[p]][['B1']][['Rapp1']] + 
                                 reproduciblePeakFeatures.Data[[p]][['B2']][['Rapp1']])/2
                    Mapp_t <- c(Mapp_t, Mapp_p)
                    Rapp_t <- c(Rapp_t, Rapp_p)
                  }
                  if (length(plantProtnames.Data[[substr(p, 1, ncharProtID)]])>0){
                    prot_p_name <- strsplit(plantProtnames.Data[[substr(p, 1, ncharProtID)]], split=' ')[[1]][1]
                  }else{
                    prot_p_name <- NA
                  }
                  prot_id_name_t <- c(prot_id_name_t, paste(p, prot_p_name, sep=", "))
                  
                }
                data_t <- data.frame(prot_id_name=prot_id_name_t, Cluster.ID_P.value=sub_pv_t, Mcalc=Mcalc_t, Mapp=Mapp_t, Rapp=Rapp_t)
                data_plot <- rbind(data_plot, data_t)
              }
              data_plot <- data_plot[-1,]
              Rapp_check <- c()
              Rapp_c_t <- data_plot$Rapp
              for (r in Rapp_c_t){
                if (r <= RappCutoffNonmono){
                  Rapp_check <- c(Rapp_check, paste0("Rapp.avg <= ", RappCutoffNonmono))
                }else{
                  Rapp_check <- c(Rapp_check, paste0("Rapp.avg > ", RappCutoffNonmono))
                }
              }
              data_plot <- mutate(data_plot, Rapp_check=Rapp_check)
              data_c_t <- data_plot
              limit_x <-max(data_c_t$Mcalc)
              limit_y <-max(data_c_t$Mapp)
              limit_all <- max(c(limit_x,limit_y))
              
              print(ggplot(data = data_c_t, aes(Mcalc, Mapp))+
                      geom_point(shape=21,aes(fill= Cluster.ID_P.value), size = 3,alpha = 0.6) +
                      geom_abline(aes(intercept=0,slope=1),linetype=2,alpha=0.7)+
                      labs(x=paste0(paste0("Mcalc of ", P_name),  " subcomplexes/singletons (kDa)"), y = paste0(paste0("Mapp-avg of subunits in ", P_name), " orthocomplexes (kDa)"),title=c)+
                      geom_text_repel(aes(Mcalc, Mapp,label=prot_id_name,color=Cluster.ID_P.value),size=3,
                                      family = 'Times',
                                      fontface = 'bold',
                                      box.padding = unit(0.5, 'lines'),
                                      point.padding = unit(1, 'lines'),
                                      segment.color = '#cccccc',
                                      segment.size = 0.25,
                                      arrow = arrow(length = unit(0.01, 'npc')),
                                      force = 1,
                                      max.iter = 3e3, max.overlaps = Inf
                      )+
                      geom_point(data=data_c_t[data_c_t$Rapp_check == paste0("Rapp.avg <= ", RappCutoffNonmono),],shape=24,size=2,alpha = 0.6) +
                      theme(plot.title=element_text(vjust=-0.5), title =element_text(size=8),
                            axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
                            axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                            legend.key.size=unit(0.5,"cm"),legend.title = element_text(size = 6,face = 6),legend.text=element_text(size=5),legend.position = "bottom")+
                      scale_x_continuous(breaks=seq(0,limit_all+50,50),limits = c(0,limit_all+50))+scale_y_continuous(breaks=seq(0,limit_all+50,50),limits = c(0,limit_all+50))
              )
            }
            dev.off()
            
            pdf(paste(paste("figure file 8. subcomplex proteins peak profile plot", Sys.Date(), sep='_'), "pdf", sep='.'))
            x = 1:dim(peakProfileMatrixB1.Chosen)[2]
            for (c in names(subcomplex_prediction)){
              for (sub in names(subcomplex_prediction[[c]])){
                proteins <- subcomplex_prediction[[c]][[sub]]
                ylimit_1 <- max(peakProfileMatrixB1.Chosen[proteins,])
                ylimit_2 <- max(peakProfileMatrixB2.Chosen[proteins,])
                par(mfrow=c(2,1))
                plot(x, peakProfileMatrixB1.Chosen[proteins[1], ], type="b", col=1, xlab="Fraction", ylab="Peak Profile", main='SEC BIO 1',cex=0.8, ylim=range(0,ylimit_1))
                axis(side=1,at=1:dim(peakProfileMatrixB1.Chosen)[2],labels=1:dim(peakProfileMatrixB1.Chosen)[2], cex.axis=0.5)
                c1 <- 1
                if (length(proteins)>1){
                  for (i in 2:length(proteins)){
                    c1 <- c1 + 1
                    lines(peakProfileMatrixB1.Chosen[proteins[i], ], type="b", col=c1)
                  }
                }
                
                legend(13,ylimit_1-ylimit_1/10, proteins, lwd=numeric(length(proteins))+2, col=c(1:length(proteins)), cex=0.5)
                plot(x, peakProfileMatrixB2.Chosen[proteins[1], ], type="b", col=1, xlab="Fraction", ylab="Peak Profile", main='SEC BIO 2',cex=0.8, ylim=range(0,ylimit_2))
                axis(side=1,at=1:dim(peakProfileMatrixB1.Chosen)[2],labels=1:dim(peakProfileMatrixB1.Chosen)[2], cex.axis=0.5)
                c2 <- 1
                if (length(proteins)>1){
                  for (i in 2:length(proteins)){
                    c2 <- c2 + 1
                    lines(peakProfileMatrixB2.Chosen[proteins[i], ], type="b", col=c2)
                  }
                }
                legend(13,ylimit_2-ylimit_2/10, proteins, lwd=numeric(length(proteins))+2, col=c(1:length(proteins)), cex=0.5)
                mtext(paste(c, sub, sep='--'), side=1,cex = 0.5,line = -1, outer=TRUE)
                
              }
            }
            dev.off()  
            
            cat("\tStep 3.8 has generated figure file 7: Scatter plot - The subcomplex Clustering identification for CORUM complexes - Mapp vs Mcalc.\n",
                "\t\t and figure file 8: Curve plot - subcomplex proteins peak profile plot.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 3.8 has generated figure file 7: Scatter plot - The subcomplex Clustering identification for CORUM complexes - Mapp vs Mcalc.\n",
                "\t\t and figure file 8: Curve plot - subcomplex proteins peak profile plot.\n")
            
            cat("Step 3: Done. The gold standard complexes have been predicted.\n\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 3: Done. The gold standard complexes have been predicted.\n\n")
            cat("Step 4: To calculate intactness and purity for the given exogenous clustering result.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 4: To calculate intactness and purity for the given exogenous clustering result.\n")
            
            subcomplex_goldstd_ratiosmall_Rapplarge_nosingleton <- data_ratiosmall_Rapplarge_nosingleton
            prots_in_clusters <- levels(factor(exogenousClusteringResult$protID))
            prots_in_clusters_real <- levels(factor(substr(prots_in_clusters, 1, ncharProtID)))
            # find gold std that have multiple overlap with clusters
            gold_std_overlap_protsInclusters <- list()
            for (comsub in com_sub_ratiosmall_Rapplarge_nosingleton){
              com <- strsplit(comsub, split = ':')[[1]][1]
              sub <- strsplit(comsub, split = ':')[[1]][2]
              prots_std <- filter(subcomplex_goldstd_ratiosmall_Rapplarge_nosingleton, ComplexName==com, SubcomplexID==sub)$Plant_Ortholog_ID
              prots_real_std <- levels(factor(substr(prots_std, 1, ncharProtID)))
              if (sum(prots_real_std %in% prots_in_clusters_real)>1){
                gold_std_overlap_protsInclusters[[comsub]] <- prots_real_std[prots_real_std %in% prots_in_clusters_real]
              }
            }
            
            find_unit_for_sample <- function(sample_data, samples_in_unit){
              unit_found <- c()
              for (s in sample_data){
                for (id in names(samples_in_unit)){
                  if (s %in% substr(samples_in_unit[[id]],1, ncharProtID)){
                    unit_found <- c(unit_found, id)
                  }
                }
              }
              return (unit_found)
            }
            
            prots_in_each_cluster <- list()
            for (cluster_type in names(exogenousClusteringResult)[2:dim(exogenousClusteringResult)[2]]){
              prots_in_each_cluster[[cluster_type]] <- list()
              sub_data <- exogenousClusteringResult[, c("protID",cluster_type)]
              cluster_id_t <- levels(factor(exogenousClusteringResult[,cluster_type]))
              for(c in cluster_id_t){
                prots_t_c <- filter(sub_data, sub_data[,2]==as.numeric(c))$protID
                prots_in_each_cluster[[cluster_type]][[c]] <- prots_t_c
              }
            }
            prots_smallworld_used<- levels(factor(substr(protsSelectedForSmallWorld, 1, ncharProtID)))
            intact_purity_for_each_clutertype <- list()
            individual_separation_fraction <- c()
            for (type_c in names(prots_in_each_cluster)){
              individual_separation <- 0
              for(comsub in names(gold_std_overlap_protsInclusters)){
                prots_cs <- gold_std_overlap_protsInclusters[[comsub]]
                find_cs <- find_unit_for_sample(sample_data=prots_cs, samples_in_unit=prots_in_each_cluster[[type_c]])
                if (length(find_cs)>0){
                  table_find_cs <- table(factor(find_cs))
                  unit_id <- names(table_find_cs[table_find_cs==max(table_find_cs)])[1]
                  if (max(table_find_cs)==1){
                    intact_purity_for_each_clutertype[[type_c]][[comsub]] <- c(intactness=0, 
                                                                               purity=1/sum(substr(prots_in_each_cluster[[type_c]][[unit_id]], 1,ncharProtID) %in% prots_smallworld_used))
                    individual_separation <- individual_separation + 1
                  }else{
                    intact_purity_for_each_clutertype[[type_c]][[comsub]] <- c(intactness=max(table_find_cs)/length(prots_cs), 
                                                                               purity=max(table_find_cs)/sum(substr(prots_in_each_cluster[[type_c]][[unit_id]], 1,ncharProtID) %in% prots_smallworld_used))
                  }
                }
              }
              individual_separation_fraction <- c(individual_separation_fraction, individual_separation/length(gold_std_overlap_protsInclusters))
            }
            # plot comparison intactness and purity
            complex_namesofallgood <- levels(factor(subcomplex_goldstd_ratiosmall_Rapplarge_nosingleton$ComplexName))
            data_here <- filter(subcomplex_prediction_alldata, !is.na(Mapp_avg_kDa))
            data_here <- mutate(data_here, Rapp_avg =(Rapp1_peak+Rapp2_peak)/2)
            data_here <- filter(data_here, Rapp_avg>RappCutoffChoosingThreshold)
            if (is.na(targetCluster)){
              targetCluster <- 1
            }
            intactness_purity_goodcover_corum <- list()
            intactness_num_corum <- list()
            com_num_corum <- list()
            clusters_in_goodsub_corum <- list()
            clusters_for_each_protein_corum <- list()
            for (c in goldstd_complex_names){
              subcom <- substr(unname(unlist(subcomplex_prediction[[c]])), 1,ncharProtID)
              over_lap <- levels(factor(subcom[subcom %in% prots_in_clusters_real]))
              if (length(over_lap)>1){
                find_c <- find_unit_for_sample(sample_data=over_lap, samples_in_unit=prots_in_each_cluster[[targetCluster]])
                if (length(find_c)>0){
                  for (k in 1:length(over_lap)){
                    clusters_for_each_protein_corum[[c]][[over_lap[k]]] <- find_c[k]
                  }
                  table_find_c <- table(factor(find_c))
                  unit_id <- names(table_find_c[table_find_c==max(table_find_c)])[1]
                  if (max(table_find_c)==1){
                    intactness_purity_goodcover_corum[[c]] <- c(intactness=0, 
                                                                purity=1/sum(substr(prots_in_each_cluster[[targetCluster]][[unit_id]], 1,ncharProtID) %in% prots_smallworld_used))
                  }else{
                    intactness_purity_goodcover_corum[[c]] <- c(intactness=max(table_find_c)/length(over_lap), 
                                                                purity=max(table_find_c)/sum(substr(prots_in_each_cluster[[targetCluster]][[unit_id]], 1,ncharProtID) %in% prots_smallworld_used))
                  }
                  
                  intactness_num_corum[[c]] <- max(table_find_c)
                  #com_num_Y[[c]][[i]] <- length(subcomplex_rice_prediction[[c]][[i]])
                }else{
                  intactness_purity_goodcover_corum[[c]] <- rep("not available", 2)
                }
              }else{
                intactness_purity_goodcover_corum[[c]] <- rep("not available", 2)
              }
            }
            goldstd_ratiosmall_Rapplarge_nosingleton <- mutate(data_ratiosmall_Rapplarge_nosingleton, 
                                                               Inactness_subcomplex_based=rep(NA,dim(data_ratiosmall_Rapplarge_nosingleton)[1]),
                                                               Inactness_corum_based=rep(NA,dim(data_ratiosmall_Rapplarge_nosingleton)[1]),
                                                               Purity_subcomplex_based=rep(NA,dim(data_ratiosmall_Rapplarge_nosingleton)[1]),
                                                               Purity_corum_based=rep(NA,dim(data_ratiosmall_Rapplarge_nosingleton)[1]),
                                                               ClusterID_in_exogenousClusteringResult=rep(NA,dim(data_ratiosmall_Rapplarge_nosingleton)[1]))
            for(i in 1:dim(goldstd_ratiosmall_Rapplarge_nosingleton)[1]){
              c <- goldstd_ratiosmall_Rapplarge_nosingleton[i,"ComplexName"]
              s <- goldstd_ratiosmall_Rapplarge_nosingleton[i,"SubcomplexID"]
              cs <- paste(c,s,sep=":")
              if(cs %in% names(intact_purity_for_each_clutertype[[targetCluster]])){
                if(intact_purity_for_each_clutertype[[targetCluster]][[cs]][1]=="not available"){
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Inactness_subcomplex_based"] <- NA
                }else{
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Inactness_subcomplex_based"] <- intact_purity_for_each_clutertype[[targetCluster]][[cs]][1]
                }
                if(intact_purity_for_each_clutertype[[targetCluster]][[cs]][2]=="not available"){
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Purity_subcomplex_based"] <- NA
                }else{
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Purity_subcomplex_based"] <- intact_purity_for_each_clutertype[[targetCluster]][[cs]][2]
                }
                
              }else{
                goldstd_ratiosmall_Rapplarge_nosingleton[i, "Inactness_subcomplex_based"] <- NA
                goldstd_ratiosmall_Rapplarge_nosingleton[i, "Purity_subcomplex_based"] <- NA
              }
              if (c %in% names(intactness_purity_goodcover_corum)){
                if(intactness_purity_goodcover_corum[[c]][1]=="not available"){
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Inactness_corum_based"] <- NA
                }else{
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Inactness_corum_based"] <- intactness_purity_goodcover_corum[[c]][1]
                }
                if(intactness_purity_goodcover_corum[[c]][2]=="not available"){
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Purity_corum_based"] <- NA
                }else{
                  goldstd_ratiosmall_Rapplarge_nosingleton[i, "Purity_corum_based"] <- intactness_purity_goodcover_corum[[c]][2]
                }
              }else{
                goldstd_ratiosmall_Rapplarge_nosingleton[i, "Inactness_corum_based"] <- NA
                goldstd_ratiosmall_Rapplarge_nosingleton[i, "Purity_corum_based"] <- NA
              }
              
              prot_tt <- substr(goldstd_ratiosmall_Rapplarge_nosingleton[i,"Plant_Ortholog_ID"],1, ncharProtID)
              where_prot <- find_unit_for_sample(sample_data=prot_tt, samples_in_unit=prots_in_each_cluster[[targetCluster]])
              where_prot <- paste(where_prot,collapse = ";")
              if (length(where_prot)>0){
                goldstd_ratiosmall_Rapplarge_nosingleton[i, "ClusterID_in_exogenousClusteringResult"] <- where_prot
              }else{
                goldstd_ratiosmall_Rapplarge_nosingleton[i, "ClusterID_in_exogenousClusteringResult"] <- NA
              }
              
            }
            
            gold_std_result_summary <- goldstd_ratiosmall_Rapplarge_nosingleton
            c_names <- levels(factor(goldstd_ratiosmall_Rapplarge_nosingleton$ComplexName,
                                     levels = unique(goldstd_ratiosmall_Rapplarge_nosingleton$ComplexName)))
            
            sub_l <- c()
            size_of_complex_by_smallworld <- c()
            for (c in c_names){
              
              sub_data <- filter(goldstd_ratiosmall_Rapplarge_nosingleton, ComplexName==c)
              size_of_complex_by_smallworld <- c(size_of_complex_by_smallworld, 
                                                 rep(length(unlist(subcomplex_prediction[[c]])),length(sub_data$Plant_Ortholog_ID)))
              s_names <- levels(factor(sub_data$SubcomplexID,
                                       levels = unique(sub_data$SubcomplexID)))
              for (s in s_names){
                sub_i_data <- filter(sub_data, SubcomplexID==s)
                length_subi <- length(levels(factor(sub_i_data$Plant_Ortholog_ID)))
                sub_l <- c(sub_l, rep(length_subi,length(sub_i_data$Plant_Ortholog_ID)))
              }
              
            }
            
            gold_std_result_summary <- mutate(gold_std_result_summary, Number_of_detected_subunits=sub_l,
                                              Number_of_subunits_detected_by_prots_used_in_complex=size_of_complex_by_smallworld)
            
            fullyAssembledComplexes <- levels(factor(filter(gold_std_result_summary,Inactness_corum_based==1,Inactness_subcomplex_based==1,)$ComplexName))
            
            gold_std_result_summary <- mutate(gold_std_result_summary, fully_assembled_complexes=rep(NA, dim(gold_std_result_summary)[1]))
            
            gold_std_result_summary <- mutate(gold_std_result_summary, class_ratio_Mapp_Mcalc=rep(NA, dim(gold_std_result_summary)[1]))
            for (i in 1:dim(gold_std_result_summary)[1]){
              if (gold_std_result_summary[i,"Ratio_Mapp_Mcalc_subcomplex"]<=2 & 0.5 <= gold_std_result_summary[i,"Ratio_Mapp_Mcalc_subcomplex"]){
                gold_std_result_summary[i,"class_ratio_Mapp_Mcalc"] <- "0.5 <= Ratio <= 2"
              }else{
                gold_std_result_summary[i,"class_ratio_Mapp_Mcalc"] <- "Ratio out of [0.5, 2]"
              }
            }
            breif_summary_gold_std <- NULL
            c_order <- levels(factor(gold_std_result_summary$ComplexName, levels = unique(gold_std_result_summary$ComplexName)))
            for (c in c_order){
              datasub_t <- filter(gold_std_result_summary, ComplexName==c)
              sub_order <- levels(factor(datasub_t$SubcomplexID, levels = unique(datasub_t$SubcomplexID)))
              for (s in sub_order){
                data_tobind <- filter(datasub_t, SubcomplexID==s)[1,]
                fac <- c[c %in% fullyAssembledComplexes]
                if (length(fac)>0){
                  data_tobind[,"fully_assembled_complexes"] <- fac
                }
                breif_summary_gold_std <- rbind(breif_summary_gold_std, data_tobind)
              }
            }
            write.csv(gold_std_result_summary, paste(paste("Table 5 - Summary of gold standard subcomplexes(complexes)", Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
            write.csv(breif_summary_gold_std, paste(paste("Table 5.1 - Brief summary of gold standard subcomplexes(complexes)", Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
            
            data_ratiosmall_Rapplarge_nosingleton <- filter(subcomplex_goldstd_Rapplarge_nosingleton,Ratio_Mapp_Mcalc_subcomplex<=2 & 0.5<=Ratio_Mapp_Mcalc_subcomplex)
            com_sub_ratiosmall_Rapplarge_nosingleton <- 
              levels(factor(paste(data_ratiosmall_Rapplarge_nosingleton$ComplexName, data_ratiosmall_Rapplarge_nosingleton$SubcomplexID,sep=":")))
            cat("\tStep 4.1 has generated Table 5: Summary of gold standard subcomplexes/complexes.\n", 
                "\t\t has generated Table 5.1: Brief summary of gold standard subcomplexes(complexes)\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 4.1 has generated Table 5: Summary of gold standard subcomplexes/complexes\n.",
                "\t\t has generated Table 5.1: Brief summary of gold standard subcomplexes(complexes)\n")
            if (length(fullyAssembledComplexes)>0){
              subcomplex_goldstd_Rapplarge_nosingleton <- mutate(subcomplex_goldstd_Rapplarge_nosingleton, fully_assembled_complexes=rep(NA, dim(subcomplex_goldstd_Rapplarge_nosingleton)[1]))
              
              subcomplex_goldstd_Rapplarge_nosingleton <- mutate(subcomplex_goldstd_Rapplarge_nosingleton, Class=rep(NA, dim(subcomplex_goldstd_Rapplarge_nosingleton)[1]))
              for (i in 1:dim(subcomplex_goldstd_Rapplarge_nosingleton)[1]){
                if (subcomplex_goldstd_Rapplarge_nosingleton[i,"Ratio_Mapp_Mcalc_subcomplex"]<=2 & 0.5<=subcomplex_goldstd_Rapplarge_nosingleton[i,"Ratio_Mapp_Mcalc_subcomplex"]){
                  subcomplex_goldstd_Rapplarge_nosingleton[i,"Class"] <- "0.5 <= Ratio <= 2"
                }else{
                  subcomplex_goldstd_Rapplarge_nosingleton[i,"Class"] <- "Ratio out of [0.5, 2]"
                }
              }
              plot_summary_gold_std <- NULL
              c_order <- levels(factor(subcomplex_goldstd_Rapplarge_nosingleton$ComplexName, levels = unique(subcomplex_goldstd_Rapplarge_nosingleton$ComplexName)))
              for (c in c_order){
                datasub_t <- filter(subcomplex_goldstd_Rapplarge_nosingleton, ComplexName==c)
                sub_order <- levels(factor(datasub_t$SubcomplexID, levels = unique(datasub_t$SubcomplexID)))
                for (s in sub_order){
                  data_tobind <- filter(datasub_t, SubcomplexID==s)[1,]
                  fac <- c[c %in% fullyAssembledComplexes]
                  if (length(fac)>0){
                    data_tobind[,"fully_assembled_complexes"] <- fac
                  }
                  plot_summary_gold_std <- rbind(plot_summary_gold_std, data_tobind)
                }
              }
              limit_xy <- max(c(max(plot_summary_gold_std$Mcalc_Subcomplex),max(plot_summary_gold_std$Mapp_Subcomplex)))
              l <- limit_xy %/% 10
              if (l>100){
                length_here <- 100
              }else if (l<=100 && l>50){
                length_here <- 50
              }else{
                length_here <- 25
              }
              pdf(paste(paste(paste0(paste0(paste0(paste0(paste0("figure file 6_with fully assembled_Scatter plot - Mcalc of ", P_name), " subcomplexes"), " V.S. Mapp-avg of subunits in "), P_name), " subcomplexes"), Sys.Date(), sep='_'), "pdf", sep='.'), width=8, height=8)
              print(ggplot(data=plot_summary_gold_std, aes(x=Mcalc_Subcomplex, y=Mapp_Subcomplex)) + geom_point(shape=21,aes(fill=Class), size = 3,alpha = 0.6)+
                      geom_abline(data=df_abline, aes(intercept=intercept,slope=slope, linetype=linetype))+
                      labs(x=paste0(paste0("Mcalc of ", P_name), " subcomplexes (kDa)"), y = paste0(paste0("Mapp-avg of subunits in ", P_name), " subcomplexes (kDa)"))+
                      geom_text_repel(aes(Mcalc_Subcomplex, Mapp_Subcomplex,label=fully_assembled_complexes),size=2,color="black",
                                      family = 'Times',
                                      fontface = 'bold',
                                      box.padding = unit(0.5, 'lines'),
                                      point.padding = unit(1, 'lines'),
                                      segment.color = '#cccccc',
                                      segment.size = 0.25,
                                      arrow = arrow(length = unit(0.01, 'npc')),
                                      force = 1,
                                      max.iter = 3e3, max.overlaps = Inf
                      )+
                      theme(plot.title=element_text(hjust=1, vjust=0.5, margin=margin(t=40,b=-30)), title =element_text(size=10),
                            axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                            axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize))+
                      theme(legend.title = element_text(size = 10,face = 10),legend.text=element_text(size=10),legend.position = "bottom")+
                      scale_x_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy))+
                      scale_y_continuous(breaks=seq(0,limit_xy,length_here),limits = c(0,limit_xy))+
                      scale_linetype_manual(name = "Separation line",
                                            values=c(2,1,3),labels=c(paste0("slpoe=",round(max(data_tt$y/data_tt$x),2)),paste0("slpoe=",1),paste0("slpoe=",round(min(data_tt$y/data_tt$x),2))))
              )
              dev.off()
              cat("\t\tStep 4.1.1 has generated figure file 6-with fully assembled: Scatter plot - Mcalc of", P_name, "subcomplexes V.S. Mapp-avg of subunits in", P_name, "subcomplexes.\n",
                  "\t\t This figure is supported by Table 5.1.\n", file="LOG OF COMPUTATION.txt", append=TRUE)
              cat("\t\tStep 4.1.1 has generated figure file 6-with fully assembled: Scatter plot - Mcalc of", P_name, "subcomplexes V.S. Mapp-avg of subunits in", P_name, "subcomplexes.\n",
                  "\t\t This figure is supported by Table 5.1.\n")
            }
            
            intactness_corum <- c()
            purity_corum <- c()
            for (c in names(intactness_purity_goodcover_corum)){
              intactness_corum <- c(intactness_corum, intactness_purity_goodcover_corum[[c]][1])
              purity_corum <- c(purity_corum, intactness_purity_goodcover_corum[[c]][2])
            }
            intactness_corum_t <- as.numeric(intactness_corum[intactness_corum != "not available"])
            purity_corum_t <- as.numeric(purity_corum[purity_corum != "not available"])

            purity_sub_t <- c()
            intactness_sub_t <- c()
            for (i in names(intact_purity_for_each_clutertype[[targetCluster]])){
              purity_sub_t <- c(purity_sub_t, intact_purity_for_each_clutertype[[targetCluster]][[i]][2])
              intactness_sub_t <- c(intactness_sub_t, intact_purity_for_each_clutertype[[targetCluster]][[i]][1])
            }

            typeOfIP <- c(rep('Subcomplex-based', length(intactness_sub_t)), rep('CORUM-based', length(intactness_corum_t)))
 
            data_intact_purity <- data.frame(intactness=c(round(intactness_sub_t,3), round(intactness_corum_t,3)),
                                             purity=c(round(purity_sub_t,3), round(purity_corum_t,3)),
                                             Type=typeOfIP)
            pdf(paste(paste("figure file 9. Bar plot - Comparision of intactnesses between subcomplex-baed and CORUM-based", Sys.Date(), sep='_'), "pdf", sep='.'), width=12, height=8)
            print(ggplot(data_intact_purity, aes(x=factor(intactness),fill=Type))+
              geom_bar(stat="count", width=0.8, position=position_dodge(0.9))+geom_text(aes(label=..count..),stat = 'count',position = position_dodge(0.9), vjust=-0.1)+
              labs(x="Intactness", y = "Count")+
              theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                    axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize),
                    legend.title = element_text(size = 10,face = 10),legend.text=element_text(size=10),legend.position = "bottom")
            )
            dev.off()  
            pdf(paste(paste("figure file 10. Bar plot - Comparision of purities between subcomplex-baed and CORUM-based", Sys.Date(), sep='_'), "pdf", sep='.'), width=12, height=8)
            
            print(ggplot(data_intact_purity, aes(x=factor(purity),fill=Type))+
              geom_bar(stat="count", width=0.8, position=position_dodge(0.9))+geom_text(aes(label=..count..),stat = 'count',position = position_dodge(0.9), vjust=-0.1)+
              labs(x="Purity", y = "Count")+
              theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = axisTextsize), axis.text.y = element_text(size = axisTextsize), 
                    axis.title.x = element_text(size = axisTitleize), axis.title.y = element_text(size = axisTitleize),
                    legend.title = element_text(size = 10,face = 10),legend.text=element_text(size=10),legend.position = "bottom")
            )
            dev.off() 
            cat("\tStep 4.2 has generated figure files 9 and 10: Bar plots - Comparision of intactnesses/purities between subcomplex-baed and CORUM-based.\n",
                file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("\tStep 4.2 has generated figure files 9 and 10: Bar plots - Comparision of intactnesses/purities between subcomplex-baed and CORUM-based.\n")
            
            if (dim(exogenousClusteringResult)[2]>2){
              num_cluster_type <- dim(exogenousClusteringResult)[2] - 1
              intactness_variation <- as.data.frame(matrix(1:(num_cluster_type*length(names(gold_std_overlap_protsInclusters))), nrow=num_cluster_type))
              purity_variation <- as.data.frame(matrix(1:(num_cluster_type*length(names(gold_std_overlap_protsInclusters))), nrow=num_cluster_type))
              
              colnames(intactness_variation) <- names(gold_std_overlap_protsInclusters)
              colnames(purity_variation) <- names(gold_std_overlap_protsInclusters)
              
              for(i in 1:dim(intactness_variation)[1]){
                for(colnm in names(intactness_variation)){
                  intactness_variation[i,colnm] <- intact_purity_for_each_clutertype[[i]][[colnm]][1]
                  purity_variation[i,colnm] <- intact_purity_for_each_clutertype[[i]][[colnm]][2]
                }
              }
              
              intactness_variation <- mutate(intactness_variation, num_clusters=as.numeric(substr(names(intact_purity_for_each_clutertype),2,10)))
              purity_variation <- mutate(purity_variation, num_clusters=as.numeric(substr(names(intact_purity_for_each_clutertype),2,10)))
              
              to_drop <-c()
              check_it <- intactness_variation[1,]
              for (i in 1:(dim(intactness_variation)[2]-1)){
                if (check_it[1,i]==0){
                  to_drop <- c(to_drop, i)
                }
              }
              
              purity_variation_plot <- purity_variation[, -to_drop]
              intactness_variation_plot <- intactness_variation[, -to_drop]
              median_intact <- c()
              median_purity <- c()
              for(i in 1:dim(purity_variation_plot)[1]){
                median_intact <- c(median_intact, median(t(intactness_variation_plot[i, 1:(dim(purity_variation_plot)[2]-1)])[,1]))
                median_purity <- c(median_purity, median(t(purity_variation_plot[i, 1:(dim(purity_variation_plot)[2]-1)])[,1]))
              }
              intactness_variation_plot <- mutate(intactness_variation_plot, median_intact=median_intact)
              purity_variation_plot <- mutate(purity_variation_plot, median_purity=median_purity)
              
              dfp <- NULL
              
              for(i in c(1:(dim(intactness_variation_plot)[2]-2))){
                
                temp_df <- data.frame(x=as.numeric(purity_variation_plot[,dim(intactness_variation_plot)[2]-1]), y=purity_variation_plot[,i], 
                                      subcomplex=rep(i, num_cluster_type), known_size=rep(length(gold_std_overlap_protsInclusters[[names(intactness_variation_plot)[i]]]),num_cluster_type))
                dfp <- rbind(dfp,temp_df)
                
              }
              
              
              dfi <- NULL
              for(i in c(1:(dim(intactness_variation_plot)[2]-2))){
                
                temp_df <- data.frame(x=as.numeric(intactness_variation_plot[,dim(intactness_variation_plot)[2]-1]), y=intactness_variation_plot[,i], 
                                      subcomplex=rep(i, num_cluster_type), known_size=rep(length(gold_std_overlap_protsInclusters[[names(intactness_variation_plot)[i]]]),num_cluster_type))
                dfi <- rbind(dfi,temp_df)
              }
              
              p_1 <- ggplot() + geom_line(data=dfp,aes(x=x,y=y,group=subcomplex,colour=known_size),size=1,alpha=0.7)+
                geom_line(data=intactness_variation_plot,aes(x=num_clusters,y=median_purity),size=1,color="black")+
                scale_colour_gradient(low = "green", high = "red")+
                labs(x="Number of clusters", y = "Purity")+
                theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
                      legend.title = element_text(size = 8,face = 6),legend.text=element_text(size=6),legend.position = "bottom")+
                annotate(geom = "text", x = 850, y = 0.56, label = "Median", hjust = "left",size=3)+
                scale_x_continuous(breaks=seq(0,1500,100),limits = c(0,1500))+scale_y_continuous(breaks=seq(0,1,0.1),limits = c(0,1))+
                guides(color=guide_colorbar(title='Number of subunits')) 
              
              
              p_2 <- ggplot() + geom_line(data=dfi,aes(x=x,y=y,group=subcomplex,colour=known_size),size=1,alpha=0.7)+
                geom_line(data=intactness_variation_plot,aes(x=num_clusters,y=median_intact),size=1,color="black")+
                scale_colour_gradient(low = "green", high = "red")+
                labs(x="Number of clusters", y = "Intactness")+
                theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
                      legend.title = element_text(size = 8,face = 6),legend.text=element_text(size=6),legend.position = "bottom")+
                annotate(geom = "text", x = 800, y = 0.645, label = "Median", hjust = "left",size=3)+
                scale_x_continuous(breaks=seq(0,1500,100),limits = c(0,1500))+scale_y_continuous(breaks=seq(0,1,0.1),limits = c(0,1))+
                guides(color=guide_colorbar(title='Number of subunits')) 
              pdf(paste(paste("figure file 11. Curve plots - the intactness and purity variation based on subcomplexes", Sys.Date(), sep='_'), "pdf", sep='.'), width=8, height=6)
              print(p_1)
              print(p_2)
              dev.off()
              cat("\tStep 4.3 has generated figure file 11: Curve plot - the intactness and purity variation based on subcomplexes.\n",
                           file="LOG OF COMPUTATION.txt", append=TRUE)
              cat("\tStep 4.3 has generated figure file 11: Curve plot - the intactness and purity variation based on subcomplexes.\n")

              monomers_real <- levels(factor(substr(monomers, 1, ncharProtID)))
              monomers_real <- monomers_real[monomers_real %in% prots_in_clusters_real]
              #sum(monomers_real %in% prots_in_clusters_real)
              homologs_monomer_found <- c()
              monomer_purity_for_each_clutertype <- list()
              homolog_monomer_purity_for_each_clutertype <- list()
              for (type_c in names(prots_in_each_cluster)){
                for(m in monomers_real){
                  find_cs <- find_unit_for_sample(sample_data=m, samples_in_unit=prots_in_each_cluster[[type_c]])
                  if (length(find_cs)>0){
                    table_find_cs <- table(factor(find_cs))
                    unit_id <- names(table_find_cs[table_find_cs==max(table_find_cs)])[1]
                    monomer_purity_for_each_clutertype[[type_c]][[m]] <- 1/sum(substr(prots_in_each_cluster[[type_c]][[unit_id]], 1,ncharProtID) %in% prots_smallworld_used)
                  }
                }
                for(comsub in homologs_monomer_found){
                  prots_cs <- gold_std_overlap_protsInclusters[[comsub]]
                  find_cs <- find_unit_for_sample(sample_data=prots_cs, samples_in_unit=prots_in_each_cluster[[type_c]])
                  if (length(find_cs)>0){
                    table_find_cs <- table(factor(find_cs))
                    unit_id <- names(table_find_cs[table_find_cs==max(table_find_cs)])[1]
                    if (max(table_find_cs)==1){
                      homolog_monomer_purity_for_each_clutertype[[type_c]][[comsub]] <- 1/sum(substr(prots_in_each_cluster[[type_c]][[unit_id]], 1,ncharProtID) %in% prots_smallworld_used)
                    }else{
                      homolog_monomer_purity_for_each_clutertype[[type_c]][[comsub]] <- max(table_find_cs)/sum(substr(prots_in_each_cluster[[type_c]][[unit_id]], 1,ncharProtID) %in% prots_smallworld_used)
                    }
                  }
                }
              }
              
              monomer_purity_variation <- as.data.frame(matrix(1:(num_cluster_type*length(monomers_real)), nrow=num_cluster_type))
              
              colnames(monomer_purity_variation) <- monomers_real
              for(i in 1:dim(monomer_purity_variation)[1]){
                for(colnm in monomers_real){
                  monomer_purity_variation[i,colnm] <- monomer_purity_for_each_clutertype[[i]][[colnm]]
                }
              }
              
              monomer_purity_variation <- mutate(monomer_purity_variation, num_clusters=as.numeric(substr(names(monomer_purity_for_each_clutertype),2,10)))
              
              dfp_m <- NULL
              
              for(i in c(1:(dim(monomer_purity_variation)[2]-1))){
                
                temp_df <- data.frame(x=as.numeric(monomer_purity_variation[,dim(monomer_purity_variation)[2]]), y=monomer_purity_variation[,i], 
                                      mono_name=rep(names(monomer_purity_variation)[i], dim(monomer_purity_variation)[1]),
                                      subcomplex=rep(i, num_cluster_type))
                dfp_m <- rbind(dfp_m,temp_df)
                
              }
  
              colors_mono <- 1:length(monomers_real)
              P_purity_m <- ggplot() + geom_line(data=dfp_m,aes(x=x,y=y,group=mono_name,colour=factor(mono_name)),size=1,alpha=0.5)+
                labs(x="Number of clusters", y = "Purity")+
                theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
                      axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
                      legend.title = element_text(size = 6,face = 6),legend.text=element_text(size=6),legend.position = "bottom")+
                scale_colour_manual(name="Protein names",values=colors_mono)
              if(length(monomers_real)>0){
                pdf(paste(paste("figure file 12. Curve plot - the monomer variation based on subcomplexes", Sys.Date(), sep='_'), "pdf", sep='.'), width=12, height=8)
                print(P_purity_m)
                dev.off()
                cat("\tStep 4.4 has generated figure file 12. Curve plot - the monomer variation based on subcomplexes.\n",
                    file="LOG OF COMPUTATION.txt", append=TRUE)
                cat("\tStep 4.4 has generated figure file 12. Curve plot - the monomer variation based on subcomplexes.\n")
              }
              
            }
            
            
            if (is.withData){
              return(list(protsSelectedForSmallWorld=protsSelectedForSmallWorld,
                          MprotToPprotMap_inEachCORUMComplex=MprotToPprotMap.inEachCORUMComplex,
                          subcomplex_prediction=subcomplex_prediction,
                          subcomplex_cutoff=subcomplex_cutoff,
                          subcomplex_overlap=subcomplex_overlap,
                          predicted_monomers=predicted_monomers,
                          subcomplex_goldstd_Rapplarge_nosingleton=subcomplex_goldstd_Rapplarge_nosingleton,
                          subcomplex_goldstd_ratiosmall_Rapplarge_nosingleton=data_ratiosmall_Rapplarge_nosingleton,
                          mcalc_mapp_sub_class=mcalc_mapp_sub_class))
            }
            cat("Step 4: Done. All Completed.\n\n", file="LOG OF COMPUTATION.txt", append=TRUE)
            cat("Step 4: Done. All Completed.\n\n")
            })

########################################################################
