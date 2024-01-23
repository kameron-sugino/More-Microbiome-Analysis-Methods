lasso_stats<-function(collect, y_vars, x_vars = NA, groups = NA){
  #by y~x*group
  if(length(x_vars)>1 & length(groups)>1){
    a<-y_vars
    b<-x_vars
    groups<-groups
    overall<-data.frame()
    for(i in 1:nrow(collect)){
      
      y<-a[,colnames(a)%in%collect[i,1]]
      vars<-b[,colnames(b)%in%gsub("\\:.*","",collect[i,2])]
      
      if(grepl("\\:",collect[i,2])){
        m<-summary(aov(y~vars*groups))
        rsq<-summary(lm(y~vars*groups))
        r<-rsq$adj.r.squared
        c<-cor.test(y,vars)
        
        #collect p values
        p<-m[[1]][["Pr(>F)"]][1:3]
        model<-data.frame(cbind(collect[i,1],collect[i,2],r,c$estimate,p[1], p[2], p[3]))
        colnames(model)<-c("measure","xvar","r2","correlation","p-value xvar", "p-value group", "p-value interaction")
        
        overall<-rbind(overall,model)
      }
      if(grepl("\\:",collect[i,2])==FALSE){
        m<-summary(aov(y~vars))
        rsq<-summary(lm(y~vars))
        r<-rsq$adj.r.squared
        c<-cor.test(y,vars)
        
        #collect p values
        p<-m[[1]][["Pr(>F)"]][1]
        model<-data.frame(cbind(collect[i,1],collect[i,2],r,c$estimate,p, NA, NA))
        colnames(model)<-c("measure","xvar","r2","correlation","p-value xvar", "p-value group", "p-value interaction")
        
        overall<-rbind(overall,model)
      }
    }
  }
  
  #by y~x
  if(length(x_vars)>1 & length(groups)==1){
    a<-y_vars
    b<-x_vars
    groups<-groups
    overall<-data.frame()
    for(i in 1:nrow(collect)){
      
      y<-a[,colnames(a)%in%collect[i,1]]
      vars<-b[,colnames(b)%in%gsub("\\:.*","",collect[i,2])]
      
      
      m<-summary(aov(y~vars))
      rsq<-summary(lm(y~vars))
      r<-rsq$adj.r.squared
      c<-cor.test(y,vars)
      
      #collect p values
      p<-m[[1]][["Pr(>F)"]][1]
      model<-data.frame(cbind(collect[i,1],collect[i,2],r,c$estimate,p))
      colnames(model)<-c("measure","xvar","r2","correlation","p-value xvar")
      
      overall<-rbind(overall,model)
    }
  }
  
  
  #by y~group
  if(length(x_vars)==1 & length(groups)>1){
    a<-y_vars
    b<-x_vars
    groups<-groups
    overall<-data.frame()
    for(i in 1:nrow(collect)){
      
      y<-a[,colnames(a)%in%collect[i,1]]
      vars<-groups
      
      m<-summary(aov(y~vars))
      rsq<-summary(lm(y~vars))
      r<-rsq$adj.r.squared
      #c<-cor.test(y,vars)
      
      #collect p values
      p<-m[[1]][["Pr(>F)"]][1]
      model<-data.frame(cbind(collect[i,1],collect[i,2],r,p))
      colnames(model)<-c("measure","group","r2","p-value xvar")
      
      overall<-rbind(overall,model)
    }
  }
  
  return(overall)
}
