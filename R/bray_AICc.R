bray_AICc<-function(variables, otus, group = NA, time = NA, id = NA){
  n<-data.frame()
  for(i in 1:ncol(variables)){
    a1<-variables[,i]
    df<-vegdist(otus,method="bray")
    perm <- how(nperm = 999)
    setBlocks(perm) <- with(otus, id)
    temp<-adonis2(df~a1+group*time,permutations=perm)
    aic<-AICc.PERMANOVA2(temp)
    n[i,1]<-colnames(variables[i])
    n[i,2]<-aic$AICc
  }
  colnames(n)<-c("model","AICc")
  
  dist<-vegdist(otus,method="bray")
  
  #simple model
  m1<-adonis2(dist~group*time,permutations=perm)
  aic<-AICc.PERMANOVA2(m1)
  aic_simple<-data.frame("simple",aic$AICc)
  colnames(aic_simple)<-c("model","AICc")
  
  #all variables model
  temp<-data.frame(variables, group, time)
  m2<-adonis2(reformulate(response="dist", termlabels=c(colnames(variables),"group*time")),data=temp,permutations=perm)
  aic<-AICc.PERMANOVA2(m2)
  aic_all<-data.frame("all",aic$AICc)
  colnames(aic_all)<-c("model","AICc")
  m<-rbind(aic_simple, n, aic_all)
  
  return(m) 
  
}