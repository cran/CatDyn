CDMN4P <-
function(par,dates,obscat,obseff,obsmbm,M.fixed,M,distr)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff    <- vector("numeric",sealen);
                  effn      <- vector("numeric",sealen);
                  predcat   <- vector("numeric",sealen);
                  res       <- vector("numeric",sealen);
                  likcontr  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  catdynmod <- matrix(0,sealen,6);
                  if(M.fixed==TRUE)
                    {
                    M         <- M;
                    logN0     <- par[1];
                    logP1     <- par[2];
                    logP2     <- par[3];
                    logP3     <- par[4];
                    logP4     <- par[5];
                    logscale  <- par[6];
                    logalpha  <- par[7];
                    logbeta   <- par[8];
                    mccum[1]  <- 0;
                    nstep[1]  <- exp(logN0)*exp(-M);
                    for(i in 2:sealen)
                       {
                       mccum[i] <- obscat[i-1] + mccum[i-1]*exp(-M);
                       nstep[i] <- exp(logN0)*exp(-M*i) +
                                   ind.P1[i]*exp(logP1)*exp(-M*(i-(ts.P1-ts.start)+1)) +
                                   ind.P2[i]*exp(logP2)*exp(-M*(i-(ts.P2-ts.start)+1)) +
                                   ind.P3[i]*exp(logP3)*exp(-M*(i-(ts.P3-ts.start)+1)) +
                                   ind.P4[i]*exp(logP4)*exp(-M*(i-(ts.P4-ts.start)+1)) -
                                   mccum[i]*exp(-M/2);
                       }
                    effeff     <- obseff^(exp(logalpha));
                    effn       <- nstep^(exp(logbeta));
                    predcat    <- exp(logscale)*(effeff*effn)*exp(-M/2);
                    if(distr=='normal')
                      {
                      res        <- obscat - predcat;
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
                   logP1     <- par[3];
                   logP2     <- par[4];
                   logP3     <- par[5];
                   logP4     <- par[6];
                   logscale  <- par[7];
                   logalpha  <- par[8];
                   logbeta   <- par[9];
                   mccum[1]  <- 0;
                   nstep[1]  <- exp(logN0)*exp(-exp(logM));
                   for(i in 2:sealen)
                      {
                      mccum[i] <- obscat[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                  ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                  ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                  ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                      }
                    effeff     <- obseff^(exp(logalpha));
                    effn       <- nstep^(exp(logbeta));
                    predcat    <- exp(logscale)*(effeff*effn)*exp(-exp(logM)/2);
                    if(distr=='normal')
                      {
                      res        <- obscat - predcat;
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

