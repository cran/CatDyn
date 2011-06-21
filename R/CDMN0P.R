CDMN0P <-
function(par,dates,obscat,obseff,obsmbm,M.fixed,M,distr)
  {
                  ts.start  <- head(dates,1);
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff    <- vector("numeric",sealen);
                  effn      <- vector("numeric",sealen);
                  predcat   <- vector("numeric",sealen);
                  res       <- vector("numeric",sealen);
                  likcontr  <- vector("numeric",sealen);
                  catdynmod <- matrix(0,sealen,6);
                  if(M.fixed==TRUE)
                    {
                    M         <- M;
                    logN0     <- par[1];
                    logscale  <- par[2];
                    logalpha  <- par[3];
                    logbeta   <- par[4];
                    mccum[1]  <- 0;
                    nstep[1]  <- exp(logN0)*exp(-M);
                    for(i in 2:sealen)
                       {
                       mccum[i] <- obscat[i-1] + mccum[i-1]*exp(-M);
                       nstep[i] <- exp(logN0)*exp(-M*i) - mccum[i]*exp(-M/2);
                       }
                    effeff     <- obseff^(exp(logalpha));
                    effn       <- nstep^(exp(logbeta));
                    predcat    <- exp(logscale)*(effeff*effn)*exp(-M/2);
                    if(distr=='normal')
                      {
                      res        <- obscat-predcat;
                      likcontr   <- res^2;
                      }
                    else
                      {
                      res        <- ifelse(obscat==0 | predcat==0,0,log(obscat)-log(predcat));
                      likcontr   <- res^2;
                      }
                    }
                  else
                    {
                    logM      <- par[1];
                    logN0     <- par[2];
                    logscale  <- par[3];
                    logalpha  <- par[4];
                    logbeta   <- par[5];
                    mccum[1]  <- 0;
                    nstep[1]  <- exp(logN0)*exp(-exp(logM));
                    for(i in 2:sealen)
                      {
                      mccum[i] <- obscat[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) - mccum[i]*exp(-exp(logM)/2);
                      }
                    effeff     <- obseff^(exp(logalpha));
                    effn       <- nstep^(exp(logbeta));
                    predcat    <- exp(logscale)*(effeff*effn)*exp(-exp(logM)/2);
                    if(distr=='normal')
                      {
                      res        <- obscat-predcat;
                      likcontr   <- res^2;
                      }
                    else
                      {
                      res        <- ifelse(obscat==0 | predcat==0,0,log(obscat)-log(predcat));
                      likcontr   <- res^2;
                      }
                    }
                  biom             <- obsmbm*1e-3*nstep*1e9;
                  catdynmod        <- data.frame(period=ts.start:ts.end,obseff=obseff,obscat=obscat,modcat=predcat,resids=res,npred=nstep,biompred=biom);
                  class(catdynmod) <- "CatDynMod";
                  return(catdynmod);
 }

