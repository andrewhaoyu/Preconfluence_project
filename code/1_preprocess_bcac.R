library(data.table)
library(dplyr)
setwd('/data/zhangh24/breast_cancer_data_analysis/')
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))

#create snp.id with official name
#SNPs with rs number will be listed as rs SNP
#SNPs without rs number will be listed as chr:pos
temp.str <- strsplit(standard_result$SNP.Onco,":")
n <- nrow(standard_result)
snpid <- rep("c",n)
for(i in 1:n){
  snpid[i] <- temp.str[[i]][1]
}
idx <- which(c(1:n)%in%grep("rs",snpid)!=T)
snpid[idx] <- paste0("chr",standard_result$chr.Onco,":",
                     standard_result$Position.Onco)[idx]
standard_result$snp.id <- snpid

#two data fixed-effect meta-analyses function
Twometa <- function(beta1,var1,beta2,var2){
  var_meta <- 1/(1/var1+1/var2)
  beta_meta <- (var_meta)*(beta1/var1+
                             beta2/var2)
  return(list(beta_meta,sqrt(var_meta)))
}
#generate meta-analysis odds ratio for icog and onco
standard_result = standard_result %>% 
  mutate(BCAC_meta_beta = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[1]],
         BCAC_meta_se = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[2]])
head(standard_result)
standard_result_update = standard_result %>% 
  select(var_name,snp.id, chr.iCOGs, Position.iCOGs, Effect.iCOGs,
         Baseline.iCOGs, EAFcontrols.iCOGs, EAFcases.iCOGs,r2.iCOGs,
         beta.iCOGs,SE.iCOGs,P1df_risk_chi.iCOGs,
         SNP.Onco,chr.Onco,Position.Onco,Effect.Onco,
         Baseline.Onco,EAFcontrols.Onco,EAFcases.Onco,r2.Onco,
         beta.Onco,SE.Onco,P1df_risk_chi.Onco,Effect.Meta,
         BCAC_meta_beta,BCAC_meta_se) %>% 
  mutate(p_meta = 2*pnorm(-abs(BCAC_meta_beta/BCAC_meta_se)),
         #effective sample size is defined as 1/var(G*beta)
         #effective sample size is calculated as 1/{var(beta)*2*f(1-f)}
         N_eff = 1/(SE.iCOGs^2*2*EAFcontrols.iCOGs*(1-EAFcontrols.iCOGs))+
           1/(SE.Onco^2*2*EAFcontrols.Onco*(1-EAFcontrols.Onco))) %>%
     rename(non.effect.iCOGs = Baseline.iCOGs,
            p.iCOGs = P1df_risk_chi.iCOGs,
            #effect.Onco is the effect allele
            #baseline.onco is the non-effect allele
            non.effect.Onco = Baseline.Onco,
            #update the p-value variable name
            p.Onco = P1df_risk_chi.Onco,
            BCAC.meta.beta = BCAC_meta_beta,
            BCAC.meta.se = BCAC_meta_se,
            p.meta = p_meta,
            N.eff = N_eff) %>%
  #change all the . in var name to _ in var name
  rename_at(vars(matches("\\.")), ~gsub("\\.", "_", .))


write.table(standard_result_update, 
            file = )


onco_result <- standard_result %>% mutate(
  or = exp(beta.Onco),
  sample_size = 1/(SE.Onco^2*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>%
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,or,SE.Onco,P1df_risk_chi.Onco,
         r2.Onco,
         MAF,
         sample_size)
colnames(onco_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A1",
                           "A2",
                           "or",
                           "se",
                           "P",
                           "info",
                           "MAF",
                           "N")
write.table(onco_result,file="/data/zhangh24/ldsc/onco_result.txt",col.names = T,quote=F)