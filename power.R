rm(list = ls())

suppressWarnings(library("ggplot2"))
suppressWarnings(library("dplyr"))
suppressWarnings(library("data.table"))
suppressWarnings(library("reshape2"))
suppressWarnings(library("foreach"))
suppressWarnings(library("doParallel"))
suppressWarnings(library("ggpubr"))
suppressWarnings(library("patchwork"))
suppressWarnings(library("gridExtra"))
suppressWarnings(library("stringr"))
suppressWarnings(library("bigsnpr"))
options(scipen=200)
#-----------------------------------------------------------------
#The script is built based on MR-JTI's simulation code
#Special thanks to Dr Dan Zhou
#-----------------------------------------------------------------
# ut: UTMOST
# tt: TransTWAS
# jt: JTI
# sample size
setwd("/put/your/working/directory/here")
methods<-c("TransTWAS","CTIMP","JTI")
tissue_name_all<-c("Brain_Frontal_Cortex_BA9","Cells_EBV_transformed_lymphocytes")
h2_all<-c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)
n_causal_genes_all<-c(40,50,60,70)
n_samples_all<-c(1000,10000,25000,50000,100000,200000,250000,300000,350000,400000,450000,500000)

for(tissue_index in 1:length(tissue_name_all)){
  tissue_name<-tissue_name_all[tissue_index]
  r_threshold<-0.1
  #--------------------------------------------------------------------------------------------------------------
  r_angle<-fread(file=paste0(".../power_data/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"_r_summary_",methods[1],".txt"),header=TRUE)
  r_utmost<-fread(file=paste0(".../simulation/power_data/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"_r_summary_",methods[2],".txt"),header=TRUE)
  r_jti<-fread(file=paste0(".../simulation/power_data/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"_r_summary_",methods[3],".txt"),header=TRUE)
  r_jti<-r_jti[,1:2]
  r_angle<-r_angle[,1:2]
  r_utmost<- r_utmost[,1:2]

  colnames(r_angle)<-c("gene_name_AngleTL","r_angle")
  colnames(r_utmost)<-c("gene_name_UTMOST","r_utmost")
  colnames(r_jti)<-c("gene_name_JTI","r_jti")

  perf_all<-merge(r_utmost,r_angle,by.x="gene_name_UTMOST",by.y="gene_name_AngleTL",all=TRUE)
  perf_all<-merge(perf_all,r_jti,by.x="gene_name_UTMOST",by.y="gene_name_JTI",all=TRUE)
  perf_all[is.na(perf_all)]=0
  colnames(perf_all)<-c("gene_name","r_ut","r_tt","r_jt")

  perf_all$r_max<-apply(perf_all[,c('r_ut','r_tt','r_jt')],1,function(x) max(abs(x)))
  perf_all<-perf_all[perf_all$r_max>=r_threshold,]
  for(h2_index in 1:length(h2_all)){
    h2_all_gene<-h2_all[h2_index]
    for(causal_index in 1:length(n_causal_genes_all)){
      n_causal_genes=n_causal_genes_all[causal_index]
      power_res<-data.frame(TransTWAS=NA,UTMOST=NA,JTI=NA)
      for (sample_index in 1:length(n_samples_all)){
        
        #simulation times per gene
        n_simulations=100
        #n of causal genes
        #heritability for all the causal genes
        n_samples<-n_samples_all[sample_index]

        #generate true expression
        set.seed(2020)
        true_exp_matrix=matrix(rnorm(n_samples*n_causal_genes,0,1),ncol = n_causal_genes)

        #generate effect size of expression on trait
        set.seed(2020)
        beta<-rnorm(n_causal_genes,0,(h2_all_gene/n_causal_genes)^0.5)
        #simulate phenotypes
        set.seed(2021)
        pheno<-true_exp_matrix %*% beta + rnorm(n_samples,mean = 0,sd=(1-h2_all_gene)^0.5)

        #sampling causal genes for simulation
        set.seed(1)

        i=1
        if (nrow(perf_all)>=n_causal_genes){
          perf<-perf_all[sample(seq(1,nrow(perf_all)),n_causal_genes,replace = F),]
          registerDoParallel(cores=50)
          each_gene<-foreach(gene_i=1:n_causal_genes,.inorder=FALSE,.combine="rbind")%dopar%{
            output<-data.frame(sample_size=NA,gene_i=NA,simu_i=NA,p_ut=NA,p_tt=NA,p_jt=NA,ut_sig_bon=NA,tt_sig_bon=NA,jt_sig_bon=NA)
            for (simu_i in 1:n_simulations){
              # #utmost
              set.seed(gene_i+simu_i+2)
              ut_exp=rnorm(n_samples,0,1)
              # angletl
              set.seed(gene_i+simu_i+3)
              tt_exp=rnorm(n_samples,0,1)
              # jti
              set.seed(gene_i+simu_i+4)
              jt_exp=rnorm(n_samples,0,1)
              
              exp<-as.matrix(data.frame(true=true_exp_matrix[,gene_i],ut=ut_exp,tt=tt_exp,jt=jt_exp))
              
              #decorrelate
              c1 <- var(exp)  # find the current correlation matrix
              chol1 <- solve(chol(c1)) # cholesky decomposition to get independence
              exp_decor<-exp %*% chol1 

              perf_matrix <- matrix( 
                c(1, perf$r_ut[gene_i], perf$r_tt[gene_i], perf$r_jt[gene_i],
                  perf$r_ut[gene_i], 1, 0, 0,
                  perf$r_tt[gene_i], 0, 1, 0,
                  perf$r_jt[gene_i], 0, 0, 1), ncol=4)
                            
              chol2 <- try(chol(perf_matrix))
              if('try-error' %in% class(chol2)){next}
              
              #generate simulated expression levels
              exp_simu <- exp_decor %*% chol2 * sd(true_exp_matrix[,gene_i]) + mean(true_exp_matrix[,gene_i])
              colnames(exp_simu)<-c('true','ut','tt','jt')
              
              output[i,'sample_size']<-n_samples
              output[i,'gene_i']<-gene_i
              output[i,'simu_i']<-simu_i

              output[i,'p_ut']<-cor.test(pheno,exp_simu[,'ut'])$p.value
              output[i,'p_tt']<-cor.test(pheno,exp_simu[,'tt'])$p.value
              output[i,'p_jt']<-cor.test(pheno,exp_simu[,'jt'])$p.value

              i=i+1
            }
            return(output)
          }
          stopImplicitCluster()

          each_gene$ut_sig_bon<-ifelse(each_gene$p_ut<=0.05/n_causal_genes,TRUE,FALSE)
          each_gene$tt_sig_bon<-ifelse(each_gene$p_tt<=0.05/n_causal_genes,TRUE,FALSE)
          each_gene$jt_sig_bon<-ifelse(each_gene$p_jt<=0.05/n_causal_genes,TRUE,FALSE)

          #summerize
          utmost<-sum(each_gene$ut_sig_bon, na.rm = TRUE)/(n_simulations*n_causal_genes)
          angle<-sum(each_gene$tt_sig_bon, na.rm = TRUE)/(n_simulations*n_causal_genes)
          jti<-sum(each_gene$jt_sig_bon, na.rm = TRUE)/(n_simulations*n_causal_genes)

          power_res[sample_index,'UTMOST']<-utmost
          power_res[sample_index,'TransTWAS']<-angle
          power_res[sample_index,'JTI']<-jti
   
          cat(paste0("(Tissue is ",tissue_name," h2 is ",h2_all_gene,", causal gene number is ",n_causal_genes,", sample size is ",n_samples,") Power of JTI is ",jti,", UTMOST is ",utmost,", AngleTL is ",angle,"\n"))
          power_res[sample_index,"sample_size"]<-n_samples
        }
      }

      if(all(is.na(power_res))==FALSE){
        power_res<-na.omit(power_res)
        write.table(power_res,paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/h2_",h2_all_gene,"_causal_gene_num_",n_causal_genes,"_power_res.txt"),quote = F,row.names = F,sep = '\t')

        p<-ggplot()+
        geom_line(data = power_res,aes(x = sample_size,y = UTMOST,colour = "UTMOST"),linewidth=1)+
        geom_point(data = power_res,aes(x = sample_size,y = UTMOST,colour = "UTMOST"),size=3)+
        geom_line(data = power_res,aes(x = sample_size,y = TransTWAS,colour = "TransTWAS"),linewidth=1)+
        geom_point(data = power_res,aes(x = sample_size,y = TransTWAS,colour = "TransTWAS"),size=3)+
        geom_line(data = power_res,aes(x = sample_size,y = JTI,colour = "JTI"),linewidth=1)+
        geom_point(data = power_res,aes(x = sample_size,y = JTI,colour = "JTI"),size=3)+
        ylim(0,0.7)+
        scale_colour_manual("",values = c("UTMOST"="darkseagreen","TransTWAS"="darkorange","JTI"="cornflowerblue"))+
        xlab("sample size")+ylab("Power")+
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill='transparent'),
          legend.box.background = element_rect(fill='transparent'),
          legend.key = element_rect(fill = "transparent")
        )+ggtitle(paste0("h2=",h2_all_gene,", causal gene=",n_causal_genes))

        #save the plot as RDS
        saveRDS(p, file = paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/ind_plots/with_JTI_h2_",h2_all_gene,"_causal_gene_number_",n_causal_genes,"_power.rds"))

        ggsave(
          filename=paste0("with_JTI_h2_",h2_all_gene,"_causal_gene_number_",n_causal_genes,"_power.png"),
          plot=p,
          device="png",
          path=paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/ind_plots/"),
          dpi=1028
        )

        rm(power_res,p)
      }
    }
  }
}

library("cowplot")

# after producing figures, we can merge them
for(tissue_index in 1:length(tissue_name_all)){
  tissue_name<-tissue_name_all[tissue_index]
  for (i in 1:length(h2_all)){
    h2_all_gene=h2_all[i]
    p1<-readRDS(paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/ind_plots/with_JTI_h2_",h2_all_gene,"_causal_gene_number_",n_causal_genes_all[1],"_power.rds"))
    p2<-readRDS(paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/ind_plots/with_JTI_h2_",h2_all_gene,"_causal_gene_number_",n_causal_genes_all[2],"_power.rds"))
    p3<-readRDS(paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/ind_plots/with_JTI_h2_",h2_all_gene,"_causal_gene_number_",n_causal_genes_all[3],"_power.rds"))
    p4<-readRDS(paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/ind_plots/with_JTI_h2_",h2_all_gene,"_causal_gene_number_",n_causal_genes_all[4],"_power.rds"))

    figure<-cowplot::plot_grid(
      p4,p3,p2,p1,nrow=4)
    saveRDS(figure, file = paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/sum_plots/h2_",h2_all_gene,".rds"))

    ggsave(
      filename=paste0("h2_",h2_all_gene,".png"),
      plot=figure,
      device="png",
      limitsize=FALSE,
      path=paste0(".../simulation/power/",ifelse(tissue_name=="Cells_EBV_transformed_lymphocytes","GEUV","ROSMAP"),"/sum_plots/")
    )
  }
}
