
linearRegress<-function(exposure.name,metabolite.name,dataset,covariate.name=NULL,  
                        metabolite.log='F',subset=NULL) 
{ 
  ## exposure.name: the single exposure, e.g. 'BMI' - in this case the SNP or score 
  ## metabolite.name: multiple metabolic outcomes, e.g. c('LDL_C','SERUM_TG') 
  
  tx <- dataset
  if(!is.null(subset)) {ind = with(tx, eval(parse(text = subset))); tx = tx[ind, ]} 
  
  ## linear regression of exposure with metabolites      
  add <- numeric() 
  fom <- formula(paste('met~', paste(c(exposure.name, covariate.name),collapse = '+'))) 
  
  for (j in 1:length(metabolite.name))  
  { 
    met <- tx[[metabolite.name[j]]];    
    fit <- lm(fom, data = tx) 
    temp <- c(summary(fit)$coef[exposure.name, ], confint(fit)[exposure.name, ]); 
    add <- rbind(add,temp);  
  }   
  rownames(add) <- metabolite.name 
  add <- data.frame(add,check.names = F) 
  add$Metabolite <- metabolite.name
  add 
} 