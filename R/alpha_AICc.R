alpha_aicc<-function(alpha_measure, variables, longitudinal = FALSE, group = NA, time = NA, id = NA) {
  
  #run as lmer(alpha_measure ~ a+group*time+(1|id))
  if(longitudinal == TRUE){
    print("Running lmer model as alpha_measure ~ variable+group*time+(1|id)")
    if(length(group)==1){
      print("ERROR: study groups must be provided if longitudinal=TRUE")
    }
    if(length(time)==1){
      print("ERROR: study time points must be provided if longitudinal=TRUE")
    }
    if(length(id)==1){
      print("ERROR: study id's must be provided if longitudinal=TRUE")
    }
    
    m<-data.frame()
    for(i in 1:ncol(variables)){
      a<-variables[,i] 
      m1 <- lmer(alpha_measure~a+group*time+(1|id))
      temp<-AICc(m1)
      m[i,1]<-colnames(variables[i])
      m[i,2]<-temp
    }
    colnames(m)<-c("variable","AICc")
    
    m1 <- lmer(alpha_measure~group*time+(1|id)) #simple model
    m_simple<-data.frame("simple",AICc(m1))
    colnames(m_simple)<-c("variable","AICc")
    m2 <- lmer(reformulate(response="alpha_measure", termlabels=c(colnames(variables),"group*time","(1|id)")),data=variables) #full model; contains all terms
    m_all<-data.frame("all",AICc(m2))
    colnames(m_all)<-c("variable","AICc")
    m_models<-rbind(m_simple, m, m_all)
  }
  
  
  # run as lm(alpha_measure~a)
  if(longitudinal == FALSE & length(group)==1 & length(time)==1){
    print("Running lm model as alpha_measure ~ variable)")
    if(length(id)>1){
      print("WARNING: study id's were provided but not used. Set longitudinal=TRUE if you want to run a time series analysis")
    }
    
    m<-data.frame()
    for(i in 1:ncol(variables)){
      a<-variables[,i] 
      m1 <- lm(alpha_measure~a)
      temp<-AICc(m1)
      m[i,1]<-colnames(variables[i])
      m[i,2]<-temp
    }
    colnames(m)<-c("variable","AICc")
    
    m2 <- lm(reformulate(response="alpha_measure", termlabels=c(colnames(variables))),data=variables) #full model; contains all terms
    m_all<-data.frame("all",AICc(m2))
    colnames(m_all)<-c("variable","AICc")
    m_models<-rbind(m, m_all)
  }
  
  
  
  #run as lm(alpha_measure ~ a+group)
  if(longitudinal == FALSE & length(group)>1 & length(time)==1){
    print("Running lm model as alpha_measure ~ variable+group")
    if(length(id)>1){
      print("WARNING: study id's were provided but not used. Set longitudinal=TRUE if you want to run a time series analysis")
    }
    
    m<-data.frame()
    for(i in 1:ncol(variables)){
      a<-variables[,i] 
      m1 <- lm(alpha_measure~a+group)
      temp<-AICc(m1)
      m[i,1]<-colnames(variables[i])
      m[i,2]<-temp
    }
    colnames(m)<-c("variable","AICc")
    
    m1 <- lm(alpha_measure~group) #simple model
    m_simple<-data.frame("simple",AICc(m1))
    colnames(m_simple)<-c("variable","AICc")
    m2 <- lm(reformulate(response="alpha_measure", termlabels=c(colnames(variables),"group")),data=variables) #full model; contains all terms
    m_all<-data.frame("all",AICc(m2))
    colnames(m_all)<-c("variable","AICc")
    m_models<-rbind(m_simple, m, m_all)
  }
  
  
  
  #run as lm(alpha_measure ~ a+time)
  if(longitudinal == FALSE & length(group)==1 & length(time)>1){
    print("Running lm model as alpha_measure ~ variable+time")
    if(length(id)>1){
      print("WARNING: study id's were provided but not used. Set longitudinal=TRUE if you want to run a time series analysis")
    }
    
    m<-data.frame()
    for(i in 1:ncol(variables)){
      a<-variables[,i] 
      m1 <- lm(alpha_measure~a+time)
      temp<-AICc(m1)
      m[i,1]<-colnames(variables[i])
      m[i,2]<-temp
    }
    colnames(m)<-c("variable","AICc")
    
    m1 <- lm(alpha_measure~time) #simple model
    m_simple<-data.frame("simple",AICc(m1))
    colnames(m_simple)<-c("variable","AICc")
    m2 <- lm(reformulate(response="alpha_measure", termlabels=c(colnames(variables),"time")),data=variables) #full model; contains all terms
    m_all<-data.frame("all",AICc(m2))
    colnames(m_all)<-c("variable","AICc")
    m_models<-rbind(m_simple, m, m_all)
  }
  
  
  return(m_models)
}


#The code above can run several models depending on the input:
  # a linear mixed effects model is run if longitudinal==TRUE, and you specify groups, timepoints, and id
  # a linear model of just alpha ~ variable
  # a linear model of alpha ~ variable + group/time

#example
m<-alpha_aicc(shan, cho_meta.e,longitudinal = TRUE,group = fam$Group, time=fam$time,id=fam$PTID)