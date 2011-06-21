CDMN1P.Lik <-
function(par,dates,obscat,obseff,M.fixed,M,distr)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
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
                  if(M.fixed==TRUE)
                    {
                    M         <- M;
                    N0        <- exp(par[1]);
                    P1        <- exp(par[2]);
                    Scale     <- exp(par[3]);
                    Alpha     <- exp(par[4]);
                    Beta      <- exp(par[5]);
                    expM      <- exp(-M);
                    expMhalf  <- exp(-M/2);
                    mccum[1]  <- 0;
                    nstep[1]  <- N0*expM;
                    for(i in 2:sealen)
                       {
                       mccum[i] <- obscat[i-1] + mccum[i-1]*expM;
                       nstep[i] <- N0*(expM^i) +
                                   ind.P1[i]*P1*expM^(i-(ts.P1-ts.start)+1) - mccum[i]*expMhalf;
                       }
                    effeff     <- obseff^Alpha;
                    effn       <- nstep^Beta;
                    predcat    <- Scale*(effeff*effn)*expMhalf;
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
                    M         <- exp(par[1]);
                    N0        <- exp(par[2]);
                    P1        <- exp(par[3]);
                    Scale     <- exp(par[4]);
                    Alpha     <- exp(par[5]);
                    Beta      <- exp(par[6]);
                    expM      <- exp(-M);
                    expMhalf  <- exp(-M/2);
                    mccum[1]  <- 0;
                    nstep[1]  <- N0*expM;
                    for(i in 2:sealen)
                      {
                      mccum[i] <- obscat[i-1] + mccum[i-1]*expM;
                      nstep[i] <- N0*(expM^i) +
                                  ind.P1[i]*P1*expM^(i-(ts.P1-ts.start+1)) - mccum[i]*expMhalf;
                      }
                    effeff     <- obseff^Alpha;
                    effn       <- nstep^Beta;
                    predcat    <- Scale*(effeff*effn)*expMhalf;
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
                  negsup <- ((sealen-2)/2)*log(sum(likcontr));
 }

