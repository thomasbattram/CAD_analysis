
linearRegress<-function(exposure.name,metabolite.name,dataset,covariate.name=NULL,  
                        metabolite.log='F',subset=NULL) 
{ 
  ## exposure.name: the single exposure, e.g. 'BMI' - in this case the SNP or score 
  ## metabolite.name: multiple metabolic outcomes, e.g. c('LDL_C','SERUM_TG') 
  
  tx=dataset
  if(!is.null(subset)) {ind=with(tx,eval(parse(text=subset))); tx=tx[ind,]} 
  
  ##linear regression of exposure with metabolites      
  add=numeric() 
  #fom produces what will be regressed against what, which is further specified by the function input and the next bit of code
  fom=formula(paste('met~',paste(c(exposure.name, covariate.name),collapse='+'))) 
  #the loop uses fom to regress each metabolite, 1 by 1, against the exposure (scores in this case)
  #then puts the values in "add"
  for (j in 1:length(metabolite.name))  
  { 
    met=tx[[metabolite.name[j]]];    
    fit=lm(fom,data=tx) 
    temp=c(summary(fit)$coef[exposure.name,],confint(fit)[exposure.name,]); 
    add=rbind(add,temp);  
  }   
  #adds row names to "add" and converts it into a data frame
  rownames(add)=metabolite.name 
  add=data.frame(add,check.names = F) 
  add$Metabolite <- metabolite.name
  add 
} 