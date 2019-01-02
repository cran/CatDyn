M.Hoenig <-
function(max.age,time.step)
  {
   time.step.set <- c("day","week","month")
   if(sum(time.step %in% time.step.set) == 0 | sum(time.step %in% time.step.set) > 1)
     {stop("time.step must be either 'day', 'week' or 'month' ")}
   if(!is.numeric(max.age))
     {stop("max.age must be numeric and in years")}
   M.yr       <- exp(1.44-0.982*log(max.age))
   hoenig.cov <- matrix(c(0.123^2,0,0,0.0439^2),2,2)
   M.yr.SE    <-deltamethod(g=list(~exp(x1-x2*log(max.age))),
                            mean=c(1.44,0.982),
                            cov=hoenig.cov,
                            ses=TRUE)
   if(time.step == "day")
     {
      M    <- M.yr/365
      M.SE <- M.yr.SE/365
     }
   else if(time.step == "week")
     {
      M    <- M.yr/53
      M.SE <- M.yr.SE/53
     }
   else if(time.step == "month")
     {
      M <- M.yr/12
      M.SE <- M.yr.SE/12
     }
   M.pred <- data.frame(M,M.SE)
   names(M.pred) <- paste(c("M.pred.","M.pred.SE."),time.step,sep="")
   return(M.pred)
  }
