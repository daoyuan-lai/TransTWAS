#taking the suggestion in MR-JTI
rm(list = ls())
stringsAsFactors=FALSE
suppressWarnings(library(ggplot2))
suppressWarnings(library(dplyr))
suppressWarnings(library(data.table))
suppressWarnings(library(foreach))
suppressWarnings(library(doParallel))

#------------------------------------------------------------------------------------------------------------------------------------
#Step 1: generate the iGene list for each method
gene_list<-read.table("/home/r9user3/UIP-Test/TWAS-TL/real-data-analysis/cleaned_geno_list_MAF_0.05_HWE_0.05_pruned.txt",header=TRUE)
gene_name<-gene_list$ENSG_ID

start=1
end=18380

res_TransTWAS<-data.frame(matrix(NA,nrow=48,ncol=1))
res_utmost<-data.frame(matrix(NA,nrow=48,ncol=1))
res_jti<-data.frame(matrix(NA,nrow=48,ncol=1))

tissue_list<-c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary",
"Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia",
"Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus",
"Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c_1","Brain_Substantia_nigra",
"Breast_Mammary_Tissue","Cells_Cultured_fibroblasts" ,"Cells_EBV_transformed_lymphocytes","Colon_Sigmoid", "Colon_Transverse",
"Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver", "Lung", "Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas", "Pituitary","Prostate",
"Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum","Spleen", 
"Stomach","Testis","Thyroid","Uterus", "Vagina","Whole_Blood")

df_tissue<-data.frame(tissue_list)

res_TransTWAS[,1]<-tissue_list
colnames(res_TransTWAS)<-c("Tissue name")
res_utmost[,1]<-tissue_list
colnames(res_utmost)<-c("Tissue name")
res_jti[,1]<-tissue_list
colnames(res_jti)<-c("Tissue name")

r_res_TransTWAS<-data.frame()
r_res_utmost<-data.frame()
r_res_jti<-data.frame()

setwd("/put/your/working/directory/here")

#-----------------------------------------------------------------------------------------------------------------------------
#Evaluate the type-I error 
tissue_all<-c("Brain_Frontal_Cortex_BA9","Cells_EBV_transformed_lymphocytes")

for(ii in 1:length(tissue_all)){
    
  tissue<-tissue_all[ii]

  main_path<-'/home/r9user3/UIP-Test/TWAS-TL/Results/weight/'
  out_path<-'/home/r9user3/UIP-Test/TWAS-TL/simulation/type1error/'

  simu_times=100 #test 50 100 1000

  #---------------------------------------------------------------------------------------------------------------------------------------------
  #Read the imputable genes of the three methods 
  gene_list_TransTWAS<-read.table(paste0('/home/r9user3/UIP-Test/TWAS-TL/simulation/type1error/data/TransTWAS_',tissue,'.txt'))%>%t()%>%as.vector()#Tranfer learning iGenes
  gene_list_utmost<-read.table(paste0('/home/r9user3/UIP-Test/TWAS-TL/simulation/type1error/data/utmost_',tissue,'.txt'))%>%t()%>%as.vector() #UTMOST iGenes
  gene_list_jti<-read.table(paste0('/home/r9user3/UIP-Test/TWAS-TL/simulation/type1error/data/jti_',tissue,'.txt'))%>%t()%>%as.vector() #JTI iGenes

  #creat dataframes to collect results
  out_c<-out_u<-out_j<-as.data.frame(matrix(data=NA,ncol=3,nrow=0))
  colnames(out_c)<-colnames(out_u)<-colnames(out_j)<-c('gene','i','p')

  #TransTWAS
  out_c_o<-foreach(i=1:length(gene_list_TransTWAS),.inorder=FALSE,.combine="rbind",.errorhandling="remove")%dopar%{
    #for (i in 1:100){
    print(i)
    
    gene<-gene_list_TransTWAS[i]
    #load predicted expression
    d_c<-readRDS(paste0(main_path,'TransTWASTL_var_corrected/',tissue,'/grex_',gene,'.rds'))  #single tissue
    
    for (j in 1:simu_times){   #1000 times
      #single tissue
      y=rnorm(nrow(d_c),mean=0,sd=1) #e~N(0,1^2)
      fit_c<-summary(lm(y~d_c[,1]))
      out_c[j,2]<-j 
      out_c[j,3]<-fit_c$coefficients[2,4]  #p value
      
    }
    
    out_c[,1]<-gene
    return(out_c)
  }

  #utmost
  out_u_o<-foreach(i=1:length(gene_list_utmost),.inorder=FALSE,.combine="rbind")%dopar%{
    #for (i in 1:100){
    print(i)
    
    gene<-gene_list_utmost[i]
    #load predicted expression

    #temp: as we only saved igene weights before, we decide only take exist GreX for now

    dir<-paste0(main_path,'CTIMP/',tissue,'/grex_',gene,'_',tissue,'.rds')
    
    if (file.exists(dir)==TRUE){
      d_u<-readRDS(dir) 

      if(any(d_u==Inf)==FALSE){
        
        for (j in 1:simu_times){   #1000 times
          y=rnorm(nrow(d_u),mean=0,sd=1) #e~N(0,1^2)
          fit_u<-summary(lm(y~d_u[,1]))
          out_u[j,2]<-j 
          out_u[j,3]<-fit_u$coefficients[2,4]  #p value
          
        }
          
        out_u[,1]<-gene
      }

    }
    return(out_u)
  }

  #JTI
  out_j_o<-foreach(i=1:length(gene_list_jti),.inorder=FALSE,.combine="rbind",.errorhandling="remove")%dopar%{
    #for (i in 1:100){
    print(i)
    
    gene<-gene_list_jti[i]
    #load predicted expression

    #temp: as we only saved igene weights before, we decide only take exist GreX for now

    dir<-paste0(main_path,'JTI/',tissue,'/grex_',gene,'_',tissue,'.rds')
    
    if (file.exists(dir)==TRUE){
      d_j<-readRDS(dir) 
      
      for (j in 1:simu_times){   #1000 times
        #print(paste0(i,' ',j))
        # if(var(d_u[,2])==0){next}
        #single tissue
        y=rnorm(nrow(d_j),mean=0,sd=1) #e~N(0,1^2)
        fit_j<-summary(lm(y~d_j[,1]))
        out_j[j,2]<-j 
        out_j[j,3]<-fit_j$coefficients[2,4]  #p value
        
      }
        
      out_j[,1]<-gene

    }
    return(out_j)
  }

  #TransTWAS
  out_c_o<-out_c_o[!is.na(out_c_o$p),]

  #UTMOST
  out_u_o<-out_u_o[!is.na(out_u_o$p),]

  #JTI
  out_j_o<-out_j_o[!is.na(out_j_o$p),]

  gg_qqplot <- function(ps, ci = 0.95,title='QQ-plot',ymax) {
    n  <- length(ps)
    df <- data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(ppoints(n)),
      clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )
    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
      geom_point(aes(expected, observed), shape = 1, size = 1.5) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
      geom_line(aes(expected, cupper), linetype = 2,color='dodgerblue3') +
      geom_line(aes(expected, clower), linetype = 2,color='dodgerblue3') +
      ylim(0,ymax) +
      xlab(log10Pe) +
      ylab(log10Po) +
      ggtitle(title)+
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
  }

  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                      ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
      print(plots[[1]])
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }



  png(paste0(out_path,'plot/',tissue,'_allgenes.png'),width = 2800,height = 1000,res=200,bg = "transparent")
  par(mfcol=c(1,1))

  #UTMOST
  p1=gg_qqplot(out_u_o$p,ymax=6,title=paste0(tissue,' (UTMOST)'))

  ggsave(
    filename=paste0("UTMOST_",tissue,"_allgenes.png"),
    plot=p1,
    device="png",
    path="/home/r9user3/UIP-Test/TWAS-TL/simulation/type1error/plot/",
    dpi=1028
  )

  #TransTWAS
  p2=gg_qqplot(out_c_o$p,ymax=6,title=paste0(tissue,' (TransTWAS)'))
  
  ggsave(
    filename=paste0("TransTWASTL_",tissue,"_allgenes.png"),
    plot=p2,
    device="png",
    path="/home/r9user3/UIP-Test/TWAS-TL/simulation/type1error/plot/",
    dpi=1028
  )

  p3=gg_qqplot(out_j_o$p,ymax=6,title=paste0(tissue,' (JTI)'))
  
  ggsave(
    filename=paste0("JTI_",tissue,"_allgenes.png"),
    plot=p3,
    device="png",
    path="/home/r9user3/UIP-Test/TWAS-TL/simulation/type1error/plot/",
    dpi=1028
  )

  multiplot(p1, p3, p2, cols=3)
  dev.off()

  #--output results--
  write.table(out_c_o,paste0(out_path,'TransTWAS_',tissue,'_allgenes.txt'),sep='\t',quote = F,row.names = F)
  write.table(out_u_o,paste0(out_path,'UTMOST_',tissue,'_allgenes.txt'),sep='\t',quote = F,row.names = F)
  write.table(out_j_o,paste0(out_path,'JTI_',tissue,'_allgenes.txt'),sep='\t',quote = F,row.names = F)

}
