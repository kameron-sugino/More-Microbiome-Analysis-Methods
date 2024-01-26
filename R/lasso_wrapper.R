lasso_wrapper<-function(y_vars, x_vars = NA, groups = NA){
  
  if(length(x_vars)>1 & length(groups)>1){
    var.n<-data.frame()
    collect<-data.frame()
    a<-y_vars
    b<-x_vars
    for(i in 1:ncol(a)){
      dat.i.c<-a[,i]
      temp<-data.frame(dat.i.c,b)
      temp.c<-temp[complete.cases(temp),]
      vars<-data.frame(temp.c[,-1])
      temp.groups<-groups[complete.cases(temp)]
      var.n[i,1]<-nrow(vars)
      
      ##Setting up models
      f <- as.formula(y ~ .*temp.groups)
      y <- temp.c$dat.i.c
      # Second step: using model.matrix to take advantage of f
      x <- model.matrix(f, vars)[, -1]
      
      #find best lambda
      cv_model <- cv.glmnet(x, y, alpha = 1, nfolds=50)
      best_lambda <- cv_model$lambda.min
      
      bm<-glmnet(x, y,alpha=1,lambda = best_lambda,standardize = F)
      print(bm)
      coef(bm)
      test<-coef(bm)
      test2<-data.frame(test@Dimnames[1],as.matrix(test))
      test3<-test2[-1,]
      sig<-test3[test3$s0!=0,]
      if(nrow(sig)>0){
        col.b<-data.frame(rep(colnames(a[i]),nrow(sig)),sig)
        collect<-rbind(collect,col.b)
      }
    }
    collect[,2]<-gsub("temp.groups","",collect[,2])
    stat<-summary(var.n$V1)
    print("Running glmnet models as y ~ x*group; make sure your group factor is set up so that the first term is the comparison group")
    print("LASSO requires complete data to estimate parameters, so rows with incomplete data were removed to complete each comparison. The following are summary statistics of the number of rows (e.g., samples, participants, etc.) included for the LASSO estimates:")
    print(stat)
  }
  
  
  if(length(x_vars)>1 & length(groups)==1){
    var.n<-data.frame()
    collect<-data.frame()
    a<-y_vars
    b<-data.frame(x_vars)
    for(i in 1:ncol(a)){
      dat.i.c<-a[,i]
      temp<-data.frame(dat.i.c,b)
      temp.c<-temp[complete.cases(temp),]
      vars<-data.frame(temp.c[,-1])
      var.n[i,1]<-nrow(vars)
      
      ##Setting up models
      f <- as.formula(y ~ .)
      y <- temp.c$dat.i.c
      # Second step: using model.matrix to take advantage of f
      x <- model.matrix(f, vars)[, -1]
      
      #find best lambda
      cv_model <- cv.glmnet(x, y, alpha = 1, nfolds=50)
      best_lambda <- cv_model$lambda.min
      
      bm<-glmnet(x, y,alpha=1,lambda = best_lambda,standardize = F)
      print(bm)
      coef(bm)
      test<-coef(bm)
      test2<-data.frame(test@Dimnames[1],as.matrix(test))
      test3<-test2[-1,]
      sig<-test3[test3$s0!=0,]
      if(nrow(sig)>0){
        col.b<-data.frame(rep(colnames(a[i]),nrow(sig)),sig)
        collect<-rbind(collect,col.b)
      }
    }
    stat<-summary(var.n$V1)
    print("Running glmnet models as y ~ x")
    print("WARNING: LASSO requires complete data to estimate parameters, so rows with incomplete data were removed to complete each comparison. The following are summary statistics of the number of rows (e.g., samples, participants, etc.) included for the LASSO estimates:")
    print(stat)
  }
  
  
  if(length(x_vars)==1 & length(groups)>1){
    var.n<-data.frame()
    collect<-data.frame()
    a<-data.frame(y_vars)
    b<-groups
    for(i in 1:ncol(a)){
      dat.i.c<-a[,i]
      temp<-data.frame(dat.i.c,b)
      temp.c<-temp[complete.cases(temp),]
      temp.groups<-groups[complete.cases(temp)]
      vars<-data.frame(temp.c[,-1])
      var.n[i,1]<-nrow(vars)
      
      ##Setting up models
      f <- as.formula(y ~ .)
      y <- temp.c$dat.i.c
      # Second step: using model.matrix to take advantage of f
      x <- model.matrix(~factor(c(temp.groups)))

      #find best lambda
      cv_model <- cv.glmnet(x, y, alpha = 1, nfolds=50)
      best_lambda <- cv_model$lambda.min
      
      bm<-glmnet(x, y,alpha=1,lambda = best_lambda,standardize = F)
      print(bm)
      coef(bm)
      test<-coef(bm)
      test2<-data.frame(test@Dimnames[1],as.matrix(test))
      test3<-test2[-1,]
      sig<-test3[test3$s0!=0,]
      if(nrow(sig)>0){
        col.b<-data.frame(rep(colnames(a[i]),nrow(sig)),sig)
        collect<-rbind(collect,col.b)
      }
    }
    stat<-summary(var.n$V1)
    print("Running glmnet models as y ~ group")
    print("WARNING: LASSO requires complete data to estimate parameters, so rows with incomplete data were removed to complete each comparison. The following are summary statistics of the number of rows (e.g., samples, participants, etc.) included for the LASSO estimates:")
    print(stat)
  }
  
  
  return(collect)
}
