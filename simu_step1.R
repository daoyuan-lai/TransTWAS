#generate the simulated datasets
#The code was adapted from https://github.com/ckhunsr1/PUMICE/blob/c7f006e40309017975ae657c0f92688c74e66949/simulation_multi.R
rm(list = ls())

#automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doParallel",
  "dplyr",
  "tidyr",
  "tidyverse",
  "data.table",
  "glmnet",
  "stringr",
  "corpcor", 
  "MASS", 
  "stats",
  "bigsnpr",
  "RcppArmadillo",
  "Rcpp",
  "liftOver"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
      )
    )
}

#input data
#-----------------------------------------------------------------------------------------------------------
tissue <- "Brain_Frontal_Cortex_BA9"
#gene_name <- "ENSG00000064687.12"                
gene_name<-"ENSG00000224051.6"                   ##Note that we need to ensure that this gene is expressed on the 48 tissues
rho_all <- c(0.7,0.3)                            ##correlation coefficient between the causal and other tissues, {0.3, 0.7}
cor_num_all <- c(47,24,0)                        ##number of correlated tissues, N_corr={0,24,47}
n_causal <- 100                                  ##number of causal snps out of all cis-snps in the region##
# train_size = 200                               ##sample size of training data,1025: we can add a 0.5 later
h2_all<-c(0.5,0.075)                             ##heritibility
sim_num <- 120
ncores<-40
#----------------------------------------------------------------------------------------------------------
tissue_list<-c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary",
"Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia",
"Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus",
"Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c_1","Brain_Substantia_nigra",
"Breast_Mammary_Tissue","Cells_Cultured_fibroblasts" ,"Cells_EBV_transformed_lymphocytes","Colon_Sigmoid", "Colon_Transverse",
"Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle", 
"Liver", "Lung", "Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas", "Pituitary","Prostate",
"Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum","Spleen", "Stomach",
"Testis","Thyroid","Uterus", "Vagina","Whole_Blood")

#----------------------------------------------------------------------------------------------------------
#Pre-processing the genotypes
##Do expression simulations in ext, GTEx, and UKB##

rds_file<-snp_attach(paste0(".../GTEX_cleaned/geno_bigsnpr/",gene_name,".rds"))
geno_gtex<-as.data.frame(rds_file$genotypes[,])

# imputing NA value with mean
geno_gtex<-geno_gtex%>%mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))

## genotype file, centered but not standardized
geno_gtex<-scale(geno_gtex, scale = FALSE)

rownames(geno_gtex)<-rds_file$fam$sample.ID
colnames(geno_gtex)<-paste0("chr",rds_file$map$chromosome[1],"_",rds_file$map$physical.pos)

exp_tis_ind<-list()
for (t in 1:length(tissue_list)){
    exp_tis_ind[[t]]<-fread(file=paste0(".../GTEX_cleaned/expression/",gene_name,"/Y_",tissue_list[t],".txt"),header=FALSE)[,1]
}
names(exp_tis_ind)<-tissue_list
sub_id<-rownames(geno_gtex)

#then produce the simulated expression
for (rho_index in 1:length(rho_all)){
  rho<-rho_all[rho_index]
    for (cor_index in 1:length(cor_num_all)){
        cor_num<-cor_num_all[cor_index]
        for (h2_index in 1:length(h2_all)){
          h2<-h2_all[h2_index]
          cat(paste0("calculation for h2: ",h2,", rho: ",rho,", cor num ",cor_num," start \n"))

          multiexp_sim_list = list()
          exp_sim_list = list()
          #---------------------------------------------------------------------------------------------------------
          #Simulation for GTEx

          # exp_sim_total = as.data.frame(matrix(NA, nrow(geno_gtex), sim_num))

          result_share<-readRDS(file=".../simulation/step1_simu/per_eqtl_shared.rds")%>%as.data.frame()
          rownames(result_share)[18]<-"Brain_Spinal_cord_cervical_c_1"
          rownames(result_share)[21]<-"Cells_EBV_transformed_lymphocytes"
          rownames(result_share)[22]<-"Cells_Cultured_fibroblasts"
          colnames(result_share)[18]<-"Brain_Spinal_cord_cervical_c_1"
          colnames(result_share)[21]<-"Cells_EBV_transformed_lymphocytes"
          colnames(result_share)[22]<-"Cells_Cultured_fibroblasts"

          result_share<-result_share[tissue,, drop = FALSE]

          #Simulate gene expression levels in GTEx
          for (jj in 1:sim_num){
            temp = exp_sim(geno=geno_gtex, result_share=result_share, gene_name=gene_name, cor_num=cor_num,  n_causal=n_causal, n_causal_shared= n_causal_shared, h2=h2, rho=rho,exp_tis_ind=exp_tis_ind)
            rownames(temp)<-rownames(geno_gtex)
            multiexp_sim_list[[jj]] = temp #all tissue simulation. multiexp_sim_list[[i]][[jj]] represents i-th reference panels' jj simulation of expression in all tissues
            rm(temp)
            # print(jj)
          }

          sub_id<-rownames(geno_gtex)

          #filter the training and testing set
          cv_config = cv_helper(length(sub_id), 5)
          cv_perm = cv_config$perm		
          cv_idx = cv_config$idx

          test_index = cv_perm[cv_idx[1,1]:cv_idx[1,2]]   #from now, index stand for number, id stand for "GTEX-..."
          test_id = sub_id[test_index] # sample id in the testing set
          train_id = setdiff(sub_id,test_id)

          #perform the analysis on the training data
          geno_train<-geno_gtex[rownames(geno_gtex) %in% train_id,]
          geno_test<-geno_gtex[rownames(geno_gtex) %in% test_id,]

          exp_train<-list()
          exp_test<-list()

          for (jj in 1:sim_num){
            exp_train[[jj]] = multiexp_sim_list[[jj]][rownames(multiexp_sim_list[[jj]]) %in% train_id,]
            exp_test[[jj]] = multiexp_sim_list[[jj]][rownames(multiexp_sim_list[[jj]]) %in% test_id,]
          }

          rm(multiexp_sim_list)

          timeStart=Sys.time()
          ncore_ut=20
          registerDoParallel(cores=ncore_ut)
          weight_utmost<-foreach(index=1:sim_num,.inorder=FALSE,.combine="cbind")%dopar%{
            res<-UTMOST(geno=geno_train,exp=exp_train[[index]],ncore=ncore_ut)
            cat("UTMOST's weight calculation arrives ",index,"\n")
            return(res)
          }
          stopImplicitCluster()
          timeEnd=Sys.time()
          cat(paste0("UTMOST takes ",difftime(timeEnd, timeStart, units='mins')," mins to finish calculation \n"))

          cat(paste0("calculation for h2: ",h2,", rho: ",rho,", cor num: ",cor_num," underway \n"))
          
          timeStart=Sys.time()
          registerDoParallel(cores=ncores)
          weight_angle<-foreach(index=1:sim_num,.inorder=FALSE,.combine="cbind")%dopar%{
            
            res<-AngleTL(geno=geno_train,exp=exp_train[[index]])

            cat("TransTWAS's weight calculation arrives ",index,"\n")
            return(res)
          }
          timeEnd=Sys.time()
          cat(paste0("AngleTL takes ",difftime(timeEnd, timeStart, units='mins')," mins to finish calculation \n"))
          cat(paste0("calculation for h2: ",h2,", rho: ",rho,", cor num: ",cor_num," underway \n"))
          #warning: the SNPs contain duplication. Thus I can not add rownames as other methods
          #Therefore I added a new column to the weights, which is not the same with previous methods
          #Warning!!!!!!!!
          stopImplicitCluster()
          

          timeStart=Sys.time()
          registerDoParallel(cores=ncores)
          weight_jti<-foreach(index=1:sim_num,.inorder=FALSE,.combine="cbind")%dopar%{

            res<-JTI(geno=geno_train,Y_exp=exp_train[[index]],gene_name=gene_name,tar_tis_name=tissue)

            cat("JTI's weight calculation arrives ",index,"\n")

            return(res)

          }

          timeEnd=Sys.time()
          cat(paste0("JTI takes ",difftime(timeEnd, timeStart, units='mins')," mins to finish calculation \n"))
          stopImplicitCluster()
          cat(paste0("calculation for h2: ",h2,", rho: ",rho,", cor num: ",cor_num," underway \n"))

          registerDoParallel(cores=ncores)
          pred_utmost<-foreach(index=1:sim_num,.inorder=FALSE,.combine="cbind")%dopar%{
            res<-as.matrix(geno_test)%*%as.vector(weight_utmost[,index])
            cat(paste0("UTMOST now at ",index," \n"))
            return(res)
          }
          stopImplicitCluster()

          registerDoParallel(cores=ncores)
          pred_jti<-foreach(index=1:sim_num,.inorder=FALSE,.combine="cbind")%dopar%{
            res<-as.matrix(geno_test)%*%as.vector(weight_jti[,index])
            cat(paste0("JTI now at ",index," \n"))
            return(res)
          }
          stopImplicitCluster()

          registerDoParallel(cores=ncores)
          pred_angle<-foreach(index=1:sim_num,.inorder=FALSE,.combine="cbind")%dopar%{
            res<-as.matrix(geno_test)%*%as.vector(weight_angle[,index])
            cat(paste0("TransTWAS now at ",index," \n"))
            return(res)
          }
          stopImplicitCluster()

          registerDoParallel(cores=ncores)
          r_utmost<-foreach(index=1:sim_num,.inorder=FALSE,.combine="c")%dopar%{
            
            res<-as.numeric(cor.test(exp_test[[index]][,1],as.vector(pred_utmost[,index]))$estimate)
            res<-ifelse(is.na(res)==FALSE,res,0)
            return(res)

          }
          stopImplicitCluster()

          registerDoParallel(cores=ncores)
          r_jti<-foreach(index=1:sim_num,.inorder=FALSE,.combine="c")%dopar%{
            
            res<-as.numeric(cor.test(exp_test[[index]][,1],as.vector(pred_jti[,index]))$estimate)
            res<-ifelse(is.na(res)==FALSE,res,0)
            return(res)

          }
          stopImplicitCluster()

          registerDoParallel(cores=ncores)
          r_angle<-foreach(index=1:sim_num,.inorder=FALSE,.combine="c")%dopar%{
            
            res<-as.numeric(cor.test(exp_test[[index]][,1],as.vector(pred_angle[,index]))$estimate)
            res<-ifelse(is.na(res)==FALSE,res,0)
            return(res)

          }
          stopImplicitCluster()

          #----------------------------------------------------------------------------------------------------
          #Produce the figure as reported in Supplementary Figure 1 of Khunsriraksakul et al.
          length(r_utmost) <- sim_num 
          length(r_jti) <- sim_num                     
          length(r_angle) <- sim_num

          r_pear<-cbind(r_utmost,r_jti,r_angle)
          write.table(r_pear,paste0(".../simulation/step1_simu/data/gtex_",gene_name,"_h2_",h2,"_rho_",rho,"_cor_num_",cor_num,".txt"),row.names=FALSE,col.names=TRUE,sep="\t")
          cat(paste0("calculation for h2: ",h2,", rho: ",rho,", cor num: ",cor_num," finished \n"))
        }
    }
}

#-------------------------------------------------------------------
gene_name<-"ENSG00000224051.6"
rho_all <- c(0.7,0.3)                            ##correlation coefficient between the causal and other tissues, {0.3, 0.7}
cor_num_all <- c(47,24,0)                        ##number of correlated tissues, N_corr={0,24,47}
h2_all<-c(0.3,0.2,0.15,0.1,0.05)
res_ind<-data.frame(matrix(NA,nrow=3,ncol=4))
colnames(res_ind)<-c("h2","JTI","UTMOST","TransTWAS")

for (ii in 1:length(rho_all)){
  rho<-rho_all[ii]
  for (jj in 1:length(cor_num_all)){
    cor_num<-cor_num_all[jj]
    for (kk in 1:length(h2_all)){
      h2<-h2_all[kk]
      
      temp<-read.table(file=paste0(".../simulation/step1_simu/data/gtex_",gene_name,"_h2_",h2,"_rho_",rho,"_cor_num_",cor_num,".txt"),header=TRUE)
      temp<-temp^2
      temp<-temp[,-3]
      
      res_ind[kk,-1]<-as.vector(colMeans(x=temp,na.rm=TRUE))
      res_ind[kk,1]<-h2_all[kk]
      rm(temp)
    }
  
    mm<-reshape2::melt(res_ind, id='h2')
    colnames(mm)<-c("h2","method","r")

    p<-ggplot(data = mm, mapping = aes(x = as.factor(h2), y = r,group=method))+ 
      geom_line(aes(color=method),linewidth = 1)+
      geom_point(aes(color=method),size=3)+
      scale_color_manual("",values = c("UTMOST"="darkseagreen","JTI"="cornflowerblue","TransTWAS"="darkorange"))+
      xlab(expression(h[e]^2)) + 
      ylab(expression(r^2))+
      ggtitle(bquote(N[corr] == .(cor_num)~", "~ rho ==.(rho)))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        legend.key = element_rect(fill = "transparent")
      )


    ggsave(
      filename=paste0("gtex_merged_",gene_name,"_rho_",rho,"_cor_num_",cor_num,".png"),
      plot=p,
      device="png",
      path=paste0(".../simulation/step1_simu/plot/"),
      dpi=1028
    )
  }
}
