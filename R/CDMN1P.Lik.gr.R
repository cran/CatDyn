CDMN1P.Lik.gr <-
function(par,dates,obscat,obseff,M.fixed,M,distr)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  period    <- ts.start:ts.end;
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
                  grlogP1s  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  if(M.fixed==TRUE)
                    {
                    M         <- M;
                    logN0     <- par[1];
                    logP1     <- par[2];
                    logscale  <- par[3];
                    logalpha  <- par[4];
                    logbeta   <- par[5];
                    mccum[1]  <- 0;
                    nstep[1]  <- exp(logN0)*exp(-M);
                    for(i in 2:sealen)
                       {
                       mccum[i] <- obscat[i-1] + mccum[i-1]*exp(-M);
                       nstep[i] <- exp(logN0)*exp(-M*i) +
                                   ind.P1[i]*exp(logP1)*exp(-M*(i-(ts.P1-ts.start+1))) -
                                   mccum[i]*exp(-M/2);
                       }
                    effeff     <- obseff^(exp(logalpha));
                    effn       <- nstep^(exp(logbeta));
                    predcat    <- exp(logscale)*(effeff*effn)*exp(-M/2);
                    n.res      <- obscat-predcat;
                    ln.res     <- ifelse(obscat==0 | predcat==0,0,log(obscat)-log(predcat));
                    if(distr=='normal')
                       {
                       totlik        <- sum(n.res^2);
                       for(i in 1:sealen)
                           {
                           grlogN0s[i] <- effeff[i]*exp(-i*M-M/2+logN0+logscale+logbeta)*n.res[i]*(nstep[i]^(exp(logbeta)-1));
                           }
                       grad.logN0    <- -((sealen-2)*sum(grlogN0s))/totlik;
                       for(i in (ts.P1-ts.start+1):sealen)
                           {
                           grlogP1s[i] <- effeff[i]*exp(-M*(i-(ts.P1-ts.start+1))-M/2+logscale+logbeta)*nstep[i]^(exp(logbeta)-1)*n.res[i];
                           }
                       grad.logP1    <- -((sum(ind.P1)-2)*sum(grlogP1s))/sum(n.res[(ts.P1-ts.start+1):sealen]^2);
                       grad.logscale <- -((sealen-2)*exp(logscale-M/2)*sum(effeff*effn*n.res))/totlik;
                       grad.logalpha <- -((sealen-2)*exp(-M/2+logscale+logalpha)*sum(effeff*effn*log(obseff)*n.res))/totlik;
                       grad.logbeta  <- -((sealen-2)*exp(-M/2+logscale+logbeta)*sum(effeff*effn*log(nstep)*n.res))/totlik;
                       }
                    else
                       {
                       totlik        <- sum(ln.res^2);
                       for(i in 1:sealen)
                           {
                           grlogN0s[i] <- exp(-i*M+logN0+logbeta)*ln.res[i]/nstep[i];
                           }
                       grad.logN0    <- -((sealen-2)*sum(grlogN0s))/totlik;
                       for(i in (ts.P1-ts.start+1):sealen)
                           {
                           grlogP1s[i] <- exp(-M*(i-(ts.P1-ts.start+1))+logbeta)*ln.res[i]/nstep[i];
                           }
                       grad.logP1    <- -((sum(ind.P1)-2)*sum(grlogP1s))/sum(ln.res[(ts.P1-ts.start+1):sealen]^2);
                       grad.logscale <- -((sealen-2)*sum(ln.res))/totlik;
                       grad.logalpha <- -((sealen-2)*exp(logalpha)*sum(log(obseff)*ln.res))/totlik;
                       grad.logbeta  <- -((sealen-2)*exp(logbeta)*sum(log(nstep)*ln.res))/totlik;
                       }
                    grad       <- c(grad.logN0,grad.logP1,grad.logscale,grad.logalpha,grad.logbeta)
                    }
                  else
                    {
                    logM      <- par[1];
                    logN0     <- par[2];
                    logP1     <- par[3];
                    logscale  <- par[4];
                    logalpha  <- par[5];
                    logbeta   <- par[6];
                    mccumgrad1 <- vector("numeric",sealen);
                    mccumgrad2 <- vector("numeric",sealen);
                    mccumgfact <- vector("numeric",sealen);
                    mccum[1]   <- 0;
                    nstep[1]   <- exp(logN0)*exp(-exp(logM));
                    for(i in 2:sealen)
                      {
                      mccum[i] <- obscat[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
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
                       mccumgfact[1] <- exp(logN0-exp(logM)+logM);
                       grlogMs[1]    <- n.res[1]*(exp(logM)*predcat[1]/2-effeff[1]*exp(logbeta+logscale-exp(logM)/2)*nstep[1]^(exp(logbeta)-1)*mccumgfact[1]);
                       for(i in 2:sealen)
                          {
                          mccumgrad1[i]  <- obscat[i-1]+(i-1)*mccumgrad1[i-1]*exp(logM-(i-1)*exp(logM));
                          mccumgrad1[i]  <- mccumgrad1[i]*exp(-exp(logM)/2);
                          mccumgrad2[i]  <- obscat[i-1] + mccumgrad2[i-1]*exp(-exp(logM));
                          mccumgrad2[i]  <- mccumgrad2[i]*(1/2)*exp(logM-exp(logM)/2);
                          if(i<=ts.P1)
                             {
                              mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*exp(logN0-i*exp(logM)+logM)
                              }
                          else
                             {
                              mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*exp(logN0-i*exp(logM)+logM) - (i-ts.P1)*exp(logM-(i-ts.P1)*exp(logM))*logP1;
                             }
                          grlogMs[i]     <- n.res[i]*(exp(logM)*predcat[i]/2-effeff[i]*exp(logbeta+logscale-exp(logM)/2)*nstep[i]^(exp(logbeta)-1)*mccumgfact[i]);
                          }
                       grad.logM     <-  (sealen-2)*sum(grlogMs)/totlik;
                       for(i in 1:sealen)
                          {
                          grlogN0s[i] <- effeff[i]*exp(-i*exp(logM)-exp(logM)/2+logN0+logscale+logbeta)*n.res[i]*(nstep[i]^(exp(logbeta)-1));
                          }
                       grad.logN0    <- -(sealen-2)*sum(grlogN0s)/totlik;
                       for(i in (ts.P1-ts.start+1):sealen)
                           {
                           grlogP1s[i] <- effeff[i]*exp(-exp(logM)*(i-(ts.P1-ts.start+1))-exp(logM)/2+logscale+logbeta)*nstep[i]^(exp(logbeta)-1)*n.res[i];
                           }
                       grad.logP1    <- -(sum(ind.P1)-2)*sum(grlogP1s)/sum(n.res[(ts.P1-ts.start+1):sealen]^2);
                       grad.logscale <- -(sealen-2)*exp(logscale-exp(logM)/2)*sum(effeff*effn*n.res)/totlik;
                       grad.logalpha <- -(sealen-2)*exp(-exp(logM)/2+logscale+logalpha)*sum(effeff*effn*log(obseff)*n.res)/totlik;
                       grad.logbeta  <- -(sealen-2)*exp(-exp(logM)/2+logscale+logbeta)*sum(effeff*effn*log(nstep)*n.res)/totlik;
                       }
                    else
                       {
                       totlik        <- sum(ln.res^2);
                       mccumgrad1[1] <- 0;
                       mccumgrad2[1] <- 0;
                       mccumgfact[1] <- exp(logN0-exp(logM)+logM);
                       grlogMs[1]    <- ln.res[1]*(effeff[1]*exp(-exp(logM)/2+logscale+logbeta)*nstep[1]^(exp(logbeta)-1)*mccumgfact[1]-exp(logM)*predcat[1]/2)/(effeff[1]*nstep[1]);
                       for(i in 2:sealen)
                          {
                          mccumgrad1[i]  <- obscat[i-1]+(i-1)*mccumgrad1[i-1]*exp(logM-(i-1)*exp(logM));
                          mccumgrad1[i]  <- mccumgrad1[i]*exp(-exp(logM)/2);
                          mccumgrad2[i]  <- obscat[i-1] + mccumgrad2[i-1]*exp(-exp(logM));
                          mccumgrad2[i]  <- mccumgrad2[i]*(1/2)*exp(logM-exp(logM)/2);
                          if(i<=ts.P1)
                             {
                              mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*exp(logN0-i*exp(logM)+logM);
                             }
                          else
                             {
                              mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*exp(logN0-i*exp(logM)+logM) - (i-ts.P1)*exp(logM-(i-ts.P1)*exp(logM))*logP1;
                             }
                          grlogMs[i]    <- ln.res[i]*(effeff[i]*exp(-exp(logM)/2+logscale+logbeta)*nstep[i]^(exp(logbeta)-1)*mccumgfact[i]-exp(logM)*predcat[i]/2)/(effeff[i]*nstep[i]);
                          }
                       grad.logM     <- -((sealen-2))*exp(exp(logM)/2-logscale)*sum(grlogMs)/totlik;
                       for(i in 1:sealen)
                           {
                           grlogN0s[i] <- exp(-i*exp(logM)+logN0+logbeta)*ln.res[i]/nstep[i];
                           }
                       grad.logN0    <- -((sealen-2)*sum(grlogN0s))/totlik;
                       for(i in (ts.P1-ts.start+1):sealen)
                           {
                           grlogP1s[i] <- exp(-exp(logM)*(i-(ts.P1-ts.start+1))+logbeta)*ln.res[i]/nstep[i];
                           }
                       grad.logP1    <- -((sum(ind.P1)-2)*sum(grlogP1s))/sum(ln.res[(ts.P1-ts.start+1):sealen]^2);
                       grad.logscale <- -((sealen-2)*sum(ln.res))/totlik;
                       grad.logalpha <- -((sealen-2)*exp(logalpha)*sum(log(obseff)*ln.res))/totlik;
                       grad.logbeta  <- -((sealen-2)*exp(logbeta)*sum(log(nstep)*ln.res))/totlik;
                       }
                    grad       <- c(grad.logM,grad.logN0,grad.logP1,grad.logscale,grad.logalpha,grad.logbeta)
                    }
                  grad <- grad;
 }
