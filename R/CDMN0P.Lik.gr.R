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
                  grMs      <- vector("numeric",sealen);
                  grN0s     <- vector("numeric",sealen);
                  if(M.fixed==TRUE)
                    {
                    M         <- M;
                    N0        <- exp(par[1]);
                    Scale     <- exp(par[2]);
                    Alpha     <- exp(par[3]);
                    Beta      <- exp(par[4]);
                    expM      <- exp(-M);
                    expMhalf  <- exp(-M/2);
                    mccum[1]  <- 0;
                    nstep[1]  <- N0*expM;
                    for(i in 2:sealen)
                       {
                       mccum[i] <- obscat[i-1] + mccum[i-1]*expM;
                       nstep[i] <- N0*(expM^i) - mccum[i]*expMhalf;
                       }
                    effeff     <- obseff^Alpha;
                    effn       <- nstep^Beta;
                    predcat    <- Scale*expMhalf*(effeff*effn);
                    n.res      <- obscat-predcat;
                    ln.res     <- ifelse(obscat==0 | predcat==0,0,log(obscat)-log(predcat));
                    if(distr=='normal')
                       {
                       totlik        <- sum(n.res^2);
                       for(i in 1:sealen)
                           {
                           grN0s[i] <- effeff[i]*(expM^i)*expMhalf*N0*Scale*Beta*n.res[i]*(nstep[i]^(Beta-1));
                           }
                       grad.N0    <- -(sealen-2)*sum(grN0s)/totlik;
                       grad.Scale <- -(sealen-2)*Scale*expMhalf*sum(effeff*effn*n.res)/totlik;
                       grad.Alpha <- -(sealen-2)*Scale*expMhalf*Alpha*sum(effeff*effn*ifelse(obseff==0,0,log(obseff))*n.res)/totlik;
                       grad.Beta  <- -(sealen-2)*Scale*expMhalf*Beta*sum(effeff*effn*log(nstep)*n.res)/totlik;
                       }
                    else
                       {
                       totlik        <- sum(ln.res^2);
                       for(i in 1:sealen)
                           {
                           grN0s[i] <- N0*Beta*(expM^i)*ln.res[i]/nstep[i];
                           }
                       grad.N0    <- -(sealen-2)*sum(grN0s)/totlik;
                       grad.Scale <- -(sealen-2)*sum(ln.res)/totlik;
                       grad.Alpha <- -(sealen-2)*Alpha*sum(ifelse(obseff==0,0,log(obseff))*ln.res)/totlik;
                       grad.Beta  <- -(sealen-2)*Beta*sum(log(nstep)*ln.res)/totlik;
                       }
                    grad       <- c(grad.N0,grad.Scale,grad.Alpha,grad.Beta)
                    }
                  else 
                    {
                    M          <- exp(par[1]);
                    N0         <- exp(par[2]);
                    Scale      <- exp(par[3]);
                    Alpha      <- exp(par[4]);
                    Beta       <- exp(par[5]);
                    expM       <- exp(-M);
                    expMhalf   <- exp(-M/2);
                    mccumgrad1 <- vector("numeric",sealen);
                    mccumgrad2 <- vector("numeric",sealen);
                    mccumgfact <- vector("numeric",sealen);
                    mccum[1]   <- 0;
                    nstep[1]   <- N0*expM;
                    for(i in 2:sealen)
                      {
                      mccum[i] <- obscat[i-1] + mccum[i-1]*expM;
                      nstep[i] <- N0*(expM^i) - mccum[i]*expMhalf;
                      }
                    effeff     <- obseff^Alpha;
                    effn       <- nstep^Beta;
                    predcat    <- Scale*(effeff*effn)*expMhalf;
                    n.res      <- obscat-predcat;
                    ln.res     <- ifelse(obscat==0 | predcat==0,0,log(obscat)-log(predcat));
                    if(distr=='normal')
                       {
                       totlik        <- sum(n.res^2);
                       mccumgrad1[1] <- 0;
                       mccumgrad2[1] <- 0;
                       mccumgfact[1] <- -M*N0*expM;
                       grMs[1]       <- n.res[1]*(M*predcat[1]/2-effeff[1]*Beta*Scale*expMhalf*nstep[1]^(Beta-1)*mccumgfact[1]);
                       for(i in 2:sealen)
                          {
                          mccumgrad1[i]  <- obscat[i-1]+(i-1)*mccumgrad1[i-1]*M*expM^(i-1);
                          mccumgrad1[i]  <- mccumgrad1[i]*expMhalf;
                          mccumgrad2[i]  <- obscat[i-1] + mccumgrad2[i-1]*expM;
                          mccumgrad2[i]  <- mccumgrad2[i]*(1/2)*M*expMhalf;
                          mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*N0*M*(expM^i);
                          grMs[i]        <- n.res[i]*(M*predcat[i]/2-effeff[i]*Beta*Scale*expMhalf*(nstep[i]^(Beta-1))*mccumgfact[i]);
                          }
                       grad.M      <-  (sealen-2)*sum(grMs)/totlik;
                       for(i in 1:sealen)
                          {
                          grN0s[i]    <- effeff[i]*(expM^i)*expMhalf*N0*Scale*Beta*n.res[i]*(nstep[i]^(Beta-1));
                          }
                       grad.N0     <- -(sealen-2)*sum(grN0s)/totlik;
                       grad.Scale  <- -(sealen-2)*Scale*expMhalf*sum(effeff*effn*n.res)/totlik;
                       grad.Alpha  <- -(sealen-2)*Scale*expMhalf*Alpha*sum(effeff*effn*ifelse(obseff==0,0,log(obseff))*n.res)/totlik;
                       grad.Beta   <- -(sealen-2)*Scale*expMhalf*Beta*sum(effeff*effn*log(nstep)*n.res)/totlik;
                       }
                    else
                       {
                       totlik        <- sum(ln.res^2);
                       mccumgrad1[1] <- 0;
                       mccumgrad2[1] <- 0;
                       mccumgfact[1] <- -M*N0*expM;
                       grMs[1]       <- ln.res[1]*(effeff[1]*Scale*Beta*expMhalf*nstep[1]^(Beta-1)*mccumgfact[1]-M*predcat[1]/2)/(effeff[1]*effn[1]);
                       for(i in 2:sealen)
                          {
                          mccumgrad1[i]  <- obscat[i-1]+(i-1)*mccumgrad1[i-1]*M*expM^(i-1);
                          mccumgrad1[i]  <- mccumgrad1[i]*expMhalf;
                          mccumgrad2[i]  <- obscat[i-1] + mccumgrad2[i-1]*expM;
                          mccumgrad2[i]  <- mccumgrad2[i]*(1/2)*M*expMhalf;
                          mccumgfact[i]  <- mccumgrad1[i] + mccumgrad2[i] - i*N0*M*(expM^i);
                          grMs[i]        <- ln.res[i]*(effeff[i]*Scale*Beta*expMhalf*(nstep[i]^(Beta-1))*mccumgfact[i]-M*predcat[i]/2)/(effeff[i]*effn[i]);
                          }
                       grad.M        <- -(sealen-2)*(1/(Scale*expMhalf))*sum(grMs,na.rm=TRUE)/totlik;
                       for(i in 1:sealen)
                          {
                           grN0s[i]     <- N0*Beta*(expM^i)*ln.res[i]/nstep[i];
                          }
                       grad.N0       <- -(sealen-2)*sum(grN0s)/totlik;
                       grad.Scale    <- -(sealen-2)*sum(ln.res)/totlik;
                       grad.Alpha    <- -(sealen-2)*Alpha*sum(ifelse(obseff==0,0,log(obseff))*ln.res)/totlik;
                       grad.Beta     <- -(sealen-2)*Beta*sum(log(nstep)*ln.res)/totlik;
                       }
                    grad       <- c(grad.M,grad.N0,grad.Scale,grad.Alpha,grad.Beta)
                    }
                  grad <-grad;
                  return(grad);
   }

