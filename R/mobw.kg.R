mobw.kg <-
function(par,lenw,lenm,method,span)
  {
   #cat("\n Body length measurements in lenw and lenm must be in the same units and weight data in lenw must be in kg")
   methods.set <- c("BFGS","CG","Nelder-Mead","SANN")
   if(!method%in%methods.set)
     {stop("Use one of the following methods: 'BFGS','CG','Nelder-Mead','SANN' ")}
   if(length(par)!=2)
     {stop("par must be a length two vector of initial values for the two-parameters power model")}
   if(any(is.na(lenw[,1])))
     {
      lenw <- lenw[-which(is.na(lenw[,1])),]
     }
   if(any(is.na(lenw[,2])))
     {
      lenw <- lenw[-which(is.na(lenw[,2])),]
     }
   if(any(is.na(lenm[,1])))
     {
      lenm <- lenm[-which(is.na(lenm[,1])),]
     }
   if(any(is.na(lenm[,2])))
     {
      lenm <- lenm[-which(is.na(lenm[,2])),]
     }
   lwfit <- optim(par=par,fn=.lw,lenw=lenw,method=method,control=list(maxit=50000))
   predw <- lwfit$par[1]*lenm[,2]^lwfit$par[2]
   mobw  <- data.frame(Month=1:12,
                       Obs.mean.kg=NA,
                       Obs.sd.kg=NA,
                       Smooth.mean.kg=NA,
                       Smooth.sd.kg=NA)
   mobw$Obs.mean.kg[match(aggregate(predw,list(lenm[,1]),mean)$Group.1,mobw$Month)]  <- aggregate(predw,list(lenm[,1]),mean)$x
   mobw$Obs.sd.kg[match(aggregate(predw,list(lenm[,1]),sd)$Group.1,mobw$Month)] <- aggregate(predw,list(lenm[,1]),sd)$x
   mobw$Smooth.mean.kg <- predict(loess(predw~lenm[,1],span=span,control=loess.control(surface="direct")),newdata=1:12)
   mobw$Smooth.sd.kg   <- predict(loess(predw~lenm[,1],span=span,control=loess.control(surface="direct")),newdata=1:12,se=TRUE)$se.fit
   return(mobw)
  }
