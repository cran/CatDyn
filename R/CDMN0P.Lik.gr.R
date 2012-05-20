CDMN0P.Lik.gr <-
function(par,dates,obscat,obseff,M.fixed,M,distr)
  {
                  ts.start  <- head(dates,1);
                  ts.end    <- tail(dates,1);
                  period    <- ts.start:ts.end;
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff    <- vector("numeric",sealen);
                  effn      <- vector("numeric",sealen);
                  predcat   <- vector("numeric",sealen);
                  n.res     <- vector("numeric",sealen);
                  ln.res    <- vector("numeric",sealen);
                  likcontr  <- vector("numeric",sealen);
                  grlogMs   <- vector("numeric",sealen);
                  grlogN0s  <- vector("numeric",sealen);
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
                    predcat    <- exp(logscale-M/2)*(effeff*effn);
                    n.res      <- obscat-predcat;
                    ln.res     <- ifelse(obscat==0 | predcat==0,0,log(obscat)-log(predcat));
                    if(distr=='normal')
                       {
                       totlik        <- sum(n.res^2);
                       for(i in 1:sealen)
                           {
                           grlogN0s[i] <- effeff[i]*exp(-i*M-M/2+logN0+logscale+logbeta)*n.res[i]*(nstep[i]^(exp(logbeta)-1));
                           }
                       grad.logN0    <- -(sealen-2)*sum(grlogN0s)/totlik;
                       grad.logscale <- -(sealen-2)*exp(logscale-M/2)*sum(effeff*effn*n.res)/totlik;
                       grad.logalpha <- -(sealen-2)*exp(-M/2+logscale+logalpha)*sum(effeff*effn*log(obseff)*n.res)/totlik;
                       grad.logbeta  <- -(sealen-2)*exp(-M/2+logscale+logbeta)*sum(effeff*effn*log(nstep)*n.res)/totlik;
                       }
                    else
                       {
                       totlik        <- sum(ln.res^2);
                       for(i in 1:sealen)
                           {
                           grlogN0s[i] <- exp(-i*M+logN0+logbeta)*ln.res[i]/nstep[i];
                           }
                       grad.logN0    <- -(sealen-2)*sum(grlogN0s)/totlik;
                       grad.logscale <- -(sealen-2)*sum(ln.res)/totlik;
                       grad.logalpha <- -(sealen-2)*exp(logalpha)*sum(log(obseff)*ln.res)/totlik;
                       grad.logbeta  <- -(sealen-2)*exp(logbeta)*sum(log(nstep)*ln.res)/totlik;
                       }
                    grad       <- c(grad.logN0,grad.logscale,grad.logalpha,grad.logbeta)
                    }
                  else 
                    {
                    logM       <- par[1];
                    logN0      <- par[2];
                    logscale   <- par[3];
                    logalpha   <- par[4];
                    logbeta    <- par[5];
                    mccumgrad1 <- vector("numeric",sealen);
                    mccumgrad2 <- vector("numeric",sealen);
                    mccumgfact <- vector("numeric",sealen);
                    mccum[1]   <- 0;
                    nstep[1]   <- exp(logN0)*exp(-exp(logM));
                    for(i in 2:sealen)
                      {
                      mccum[i] <- obscat[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) - mccum[i]*exp(-exp(logM)/2);
                      }
                    effeff     <- obseff^(exp(logalpha));
                    effn       <- nstep^(exp(logbeta));
                    predcat    <- exp(logscale)*(effeff*effn)*exp(-exp(logM)/2);
                    n.res      <- obscat-predcat;
                    ln.res     <- ifelse(obscat==0 | predcat==0,0,log(obscat)-log(predcat));
                    if(distr=='normal')
                       {
                       totlik        <- sum(n.res^2);
                       mccumgrad1[1] <- 0;
                       mccumgrad2[1] <- 0;
                       mccumgfact[1] <- -exp(logN0-exp(logM)+logM);
                       grlogMs[1]    <- n.res[1]*(exp(logM)*predcat[1]/2-effeff[1]*exp(logbeta+logscale-exp(logM)/2)*nstep[1]^(exp(logbeta)-1)*mccumgfact[1]);
                       for(i in 2:sealen)
                          {
                          mccumgrad1[i]  <- obscat[i-1]+(i-1)*mccumgrad1[i-1]*exp(logM-(i-1)*exp(logM));
                          mccumgrad1[i]  <- mccumgrad1[i]*exp(-exp(logM)/2);
                          mccumgrad2[i]  <- obscat[i-1] + mccumgrad2[i-1]*exp(-exp(logM));
                          mccumgrad2[i]  <- mccumgrad2[i]*(1/2)*exp(logM-exp(logM)/2);
                          mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*exp(logN0-i*exp(logM)+logM);
                          grlogMs[i]     <- n.res[i]*(exp(logM)*predcat[i]/2-effeff[i]*exp(logbeta+logscale-exp(logM)/2)*nstep[i]^(exp(logbeta)-1)*mccumgfact[i]);
                          }
                       grad.logM      <-  (sealen-2)*sum(grlogMs)/totlik;
                       for(i in 1:sealen)
                          {
                          grlogN0s[i] <- effeff[i]*exp(-i*exp(logM)-exp(logM)/2+logN0+logscale+logbeta)*n.res[i]*(nstep[i]^(exp(logbeta)-1));
                          }
                       grad.logN0     <- -(sealen-2)*sum(grlogN0s)/totlik;
                       grad.logscale  <- -(sealen-2)*exp(logscale-exp(logM)/2)*sum(effeff*effn*n.res)/totlik;
                       grad.logalpha  <- -(sealen-2)*exp(-exp(logM)/2+logscale+logalpha)*sum(effeff*effn*log(obseff)*n.res)/totlik;
                       grad.logbeta   <- -(sealen-2)*exp(-exp(logM)/2+logscale+logbeta)*sum(effeff*effn*log(nstep)*n.res)/totlik;
                       }
                    else
                       {
                       totlik        <- sum(ln.res^2);
                       mccumgrad1[1] <- 0;
                       mccumgrad2[1] <- 0;
                       mccumgfact[1] <- -exp(logN0-exp(logM)+logM);
                       grlogMs[1]    <- ln.res[1]*(effeff[1]*exp(-exp(logM)/2+logscale+logbeta)*nstep[1]^(exp(logbeta)-1)*mccumgfact[1]-exp(logM)*predcat[1]/2)/(effeff[1]*nstep[1]);
                       for(i in 2:sealen)
                          {
                          mccumgrad1[i]  <- obscat[i-1]+(i-1)*mccumgrad1[i-1]*exp(logM-(i-1)*exp(logM));
                          mccumgrad1[i]  <- mccumgrad1[i]*exp(-exp(logM)/2);
                          mccumgrad2[i]  <- obscat[i-1] + mccumgrad2[i-1]*exp(-exp(logM));
                          mccumgrad2[i]  <- mccumgrad2[i]*(1/2)*exp(logM-exp(logM)/2);
                          mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*exp(logN0-i*exp(logM)+logM);
                          grlogMs[i]     <- ln.res[i]*(effeff[i]*exp(-exp(logM)/2+logscale+logbeta)*nstep[i]^(exp(logbeta)-1)*mccumgfact[i]-exp(logM)*predcat[i]/2)/(effeff[i]*nstep[i]);
                          }
                       grad.logM      <- -(sealen-2)*exp(-exp(logM)/2-logscale)*sum(grlogMs)/totlik;
                       for(i in 1:sealen)
                          {
                           grlogN0s[i]   <- exp(-i*exp(logM)+logN0+logbeta)*ln.res[i]/nstep[i];
                          }
                       grad.logN0    <- -((sealen-2)*sum(grlogN0s))/totlik;
                       grad.logscale <- -((sealen-2)*sum(ln.res))/totlik;
                       grad.logalpha <- -((sealen-2)*exp(logalpha)*sum(log(obseff)*ln.res))/totlik;
                       grad.logbeta  <- -((sealen-2)*exp(logbeta)*sum(log(nstep)*ln.res))/totlik;
                       }
                    grad       <- c(grad.logM,grad.logN0,grad.logscale,grad.logalpha,grad.logbeta)
                    }
                  grad <-grad;
                  return(grad);
   }
