#-------------------------------------------------------------------------------------------------------
#' @title ridge_closed_form_result
cv_helper <- function(N, fold){
  ## helper function for generating cross-validation sets
  ## args
  ## N: number of sample size
  ## fold: number of folds
  ## values
  ## perm: a permutation of 1 to N
  ## idx: matrix of fold by 2 with first col being starting index and second col being ending index
  valid_num = floor(N/fold)
  set.seed(123)
  perm = sample(1:N, size = N)
  idx1 = seq(1,N,valid_num)
  idx2 = c(idx1[-1]-1,N)
  list(perm=perm, idx=cbind(idx1,idx2))
}

# #-------------------------------------------------------------------------------------------------------
minmax_lambda <- function(lst){
  ## get the minimum and maximum of lambda searched in cross-validation of an elastic net model
  ## args
  ## lst: an object returned by glmnet
  ## value
  ## min_lam: smallest lambda searched in glmnet cross-validation
  ## max_lam: largest lambda searched in glmnet cross-validation
  max_lam = max(unlist(lapply(lst, function(x){max(x$lambda)})))
  min_lam = min(unlist(lapply(lst, function(x){min(x$lambda)})))
  c(min_lam, max_lam)
}

#-------------------------------------------------------------------------------------------------------
#' @title mean_squared_error
#' @description Calculating MSE
#' @param y_est,y_test estimated y and actual y
#' @return MSE


mean_squared_error <- function(y_est, y_test){
  mse = mean((y_est - y_test)^2)
  return(mse)
}

#----------------------------------------------------------------------------------------------------
TL_estimate_prepare<-function(X,w.src,target_index){

  w.src<-w.src[,-target_index,drop=FALSE]

  if(ncol(w.src)!=1) {
    w.src_unit = apply(w.src, 2, function(x) x/sqrt(sum(x^2)))
    ##ensemble w via eigen values
    B = X%*%w.src_unit
    G = cor(B)
    w_weight = (eigen(G)$vectors[,1])^2
    w.src = w.src_unit%*%w_weight
  }

  return(w.src)
}

#----------------------------------------------------------------------------------------------------
#Note that the X_all is for the target dataset
get_initial_penalty<-function(X_all,n_tis,single_initial_est_all,single_initial_pen_all,n_snps){
  
  initial_pen<-list()
  w_src<-matrix(NA,nrow=n_snps,ncol=n_tis)
  angle_res<-matrix(NA,nrow=(n_tis-1),ncol=n_tis)

  for (t in 1:n_tis){ 
    w.src<-single_initial_est_all[,-t,drop=FALSE]
    beta<-single_initial_est_all[,t,drop=FALSE]
    w.src_unit = apply(w.src, 2, function(x) x/sqrt(sum(x^2)))

    if(ncol(w.src)!=1){
      B = X_all[[t]]%*%w.src_unit
      G = cor(B)
      w_weight = (eigen(G)$vectors[,1,drop=FALSE])^2
      w.src = w.src_unit%*%w_weight
    }
    
    var = c(var(beta),var(w.src))
    rho = as.numeric(cor(beta, w.src))
    
    lam_range<-single_initial_pen_all[,t]
    lam_V = seq(lam_range[1], lam_range[2], length.out = 5) 
    eta_range<-rho * lam_V * sqrt(var[1]/var[2])

    initial_pen[[t]]<-cbind(lam_V,eta_range)
    w_src[,t]<-w.src
  }

  return(list(penalty_initial=initial_pen,w_src=w_src))
}

#Part of the code is modified from https://github.com/yiminghu/CTIMP/blob/master/main.R 
TransTWAS = function(fold = 5, 
                      n_tune = 5, 
                      n_cores=n_cores,
                      gene_name=gene_name,
                      dir_geno=dir_geno,
                      dir_exp=dir_exp,
                      dir_output=dir_output,
                      tar_tis_name=tar_tis_name){
    string<-list.files(paste0(dir_exp))
    string<-str_extract(string, "(?<=Y_).*(?=.txt)")
    n_tis = length(string)

    Y<-list()
    for (t in 1:length(string)){
      Y[[t]]<-read.table(file=paste0(dir_exp,"/Y_",string[t],".txt"),header=FALSE)
    }

    tar_tis_index<-which(string==paste0(tar_tis_name))
    Y<- c(Y[tar_tis_index], Y[-tar_tis_index])  
    string<-c(string[tar_tis_index],string[-tar_tis_index])

    rds_file<-snp_attach(paste0(dir_geno))
    geno<-as.matrix(rds_file$genotypes[,])
    N = nrow(geno)
    n_snps = ncol(geno)

    if(!is.numeric(fold)|fold<0){
      stop("fold should be a non-negative integer.")
    }
    
    if(!is.numeric(n_tune)|n_tune<0){
      stop("n_tune should be a non-negative integer.")
    }
    
    if(!is.numeric(n_cores)|n_cores<0){
      stop("n_cores should be a non-negative integer.")
    }
    
    if(!is.numeric(n_tis)|n_tis<0){
      stop("n_tis should be a non-negative integer.")
    }
    
    if (file.exists(dir_geno)==FALSE){
      stop("Genotype matrix does not exists.") 
    }
    
    if (file.exists(dir_exp)==FALSE){
      stop("Expression file does not exists.")
    }
    
    if (length(string)==0||length(string)==1){
      stop("At least two source tissues are needed for TransTWAS.")
    }

    if(tar_tis_name %in% string==FALSE){
      stop("The target tissue name is not found in current gene's expression files.")
    }
    
    for(i in 1:ncol(geno)){
      geno[is.na(geno[,i]), i] <- mean(geno[,i], na.rm = TRUE)
    }
    
    for(j in 1:ncol(geno)){ 
      geno[,j] = geno[,j] - mean(geno[,j])
    }
    
    rownames(geno)<-rds_file$fam$sample.ID
    sub_id=rownames(geno)
    sub_id_map = list() #get the subject ID in each tissue
    
    for (t in 1:n_tis){
      tmp = rep(0, nrow(Y[[t]]))
      tmp <- Y[[t]][,1] #get the ID of the people
      sub_id_map[[t]] <- tmp
    }
    
    cv_config = cv_helper(N, fold)
    cv_perm = cv_config$perm		
    cv_idx = cv_config$idx
    
    penalty_initial<-list()
    
    test_index = cv_perm[cv_idx[1,1]:cv_idx[1,2]] 
    test_id = sub_id[test_index] # sample id in the testing set
    tuning_index = cv_perm[cv_idx[1%%fold+1,1]:cv_idx[1%%fold+1,2]]
    tuning_id = sub_id[tuning_index]
    
    X_test = list()
    Y_test = list()
    X_tune = list()
    Y_tune = list()
    X_train = list()
    Y_train = list()
    X_train_tune = list()
    Y_train_tune = list()
    X_all = list()
    Y_all = list()
    
    for(t in 1:n_tis){  
      X_all_tmp = sub_id_map[[t]]
      X_train_tmp = sub_id_map[[t]][!(sub_id_map[[t]] %in% c(tuning_id,test_id))]
      Y_train_tmp = !(sub_id_map[[t]] %in% c(tuning_id,test_id))
      
      X_tuning_tmp = sub_id_map[[t]][(sub_id_map[[t]] %in% tuning_id)]
      Y_tuning_tmp = (sub_id_map[[t]] %in% tuning_id)
      
      X_test_tmp = sub_id_map[[t]][(sub_id_map[[t]] %in% test_id)]
      Y_test_tmp = (sub_id_map[[t]] %in% test_id)
      
      X_train[[t]] <- geno[X_train_tmp,]
      Y_train[[t]] = Y[[t]][Y_train_tmp, 2]%>%as.matrix()
      
      X_tune[[t]] <- geno[X_tuning_tmp,]
      Y_tune[[t]] = Y[[t]][Y_tuning_tmp, 2]%>%as.matrix()
      
      X_test[[t]] <- geno[X_test_tmp,]
      Y_test[[t]] = Y[[t]][Y_test_tmp, 2]%>%as.matrix()
      
      X_train_tune[[t]]=rbind(X_train[[t]],X_tune[[t]])
      Y_train_tune[[t]]=c(Y_train[[t]],Y_tune[[t]])
      
      #This should be subsetting w.r.t. the rownames,correct this later.
      X_all[[t]] = geno[X_all_tmp,]
      Y_all[[t]] = Y[[t]][,2]%>% as.matrix()
    }
    
    #This thing only need to do once for now. i.e, just do it once for the target. This is for the ease of use.
    single_initial_est_all = matrix(0, nrow=n_snps, ncol=n_tis) #N of snp by N tissue matrix
    single_summary_all = list()
    
    if(n_cores>=2){
      registerDoParallel(cores=n_cores)
      single_summary_all<-foreach(t=1:n_tis,.inorder=FALSE,.packages="glmnet",.combine="cbind") %dopar% {
        res<-cv.glmnet(X_all[[t]],Y_all[[t]], family = 'gaussian', alpha = 0, nfolds = 5) 
        res_all<-c(as.vector(res$glmnet.fit$beta[,which.min(res$cvm)]),c(min(res$lambda),max(res$lambda)))  #the last two rows are the penalty parameter's range
        return(res_all)
      }
      stopImplicitCluster()
    }else{
      single_summary_all<-foreach(t=1:n_tis,.inorder=FALSE,.packages="glmnet",.combine="cbind")%do% {
        res<-cv.glmnet(X_all[[t]],Y_all[[t]], family = 'gaussian', alpha = 0, nfolds = 5) 
        res_all<-c(as.vector(res$glmnet.fit$beta[,which.min(res$cvm)]),c(min(res$lambda),max(res$lambda)))  #the last two rows are the penalty parameter's range
        return(res_all)
      }
    }
    
    colnames(single_summary_all)<-NULL
    single_initial_est_all<-single_summary_all[-((n_snps+1):(n_snps+2)),]
    single_initial_pen_all<-single_summary_all[((n_snps+1):(n_snps+2)),]
    
    TL_initial<-get_initial_penalty(X_all=X_all,
                                    n_tis=n_tis,
                                    single_initial_est_all=single_initial_est_all,
                                    single_initial_pen_all=single_initial_pen_all,
                                    n_snps=n_snps)
    
    penalty_initial<-TL_initial$penalty_initial
    w_src<-TL_initial$w_src  #the source estimates are all very small
    
    rm(TL_initial)
    
    l_2D <-penalty_initial[[1]][,1]
    e_2D <-penalty_initial[[1]][,2] 
    ridge_res<-ridge_closed_form_cpp_econ(x=X_train[[1]], y=Y_train[[1]], lam=l_2D, w=w_src[,1], eta=e_2D, x_val=X_tune[[1]], y_val=Y_tune[[1]])
    #get the penalty parameter
    temp<-c(ridge_res[1,3],ridge_res[1,4])
    tune_para_2D<-temp
    
    #retrain the model on the whole dataset using the optimal penalty parameter. This practice is observed in other paper
    weight<-ridge_closed_form_best_cpp_econ(x=X_all[[1]], y=Y_all[[1]], lam=as.numeric(tune_para_2D[1]), w=w_src[,1], eta=as.numeric(tune_para_2D[2]))   
    downstream_est = data.frame(rds_file$map, weight)
    remove(weight)
    colnames(downstream_est)[7]<-"weights"
    
    #---------------------------------------------------------------------------
    pred_exp<-as.matrix(X_all[[1]])%*%downstream_est[,7] #get the GreX here
    weight_df<-data.frame(
      gene=gene_name,
      marker_ID=downstream_est$marker.ID,
      chromosome=downstream_est$chromosome,
      allele1=downstream_est$allele1,
      allele2=downstream_est$allele2,  
      weight=downstream_est[,7])
    
    weight_df<-weight_df[weight_df$weight!=0,]
    if(nrow(weight_df)==0){
      stop("All SNPs have weight zero after transfer learning")
    }

    #Gene expression weights
    saveRDS(weight_df,paste0(dir_output,"/",string[1],"_",gene_name,'.rds'))
    #Imputated gene expression
    saveRDS(pred_exp,paste0(dir_output,"/",string[1],"_grex_",gene_name,".rds"))
}