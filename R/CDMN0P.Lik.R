CDMN0P.Lik <-
function(par,dates,obscat,obseff,M.fixed,M,distr)
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
                      effeff      <- obseff^Alpha;
                      effn        <- nstep^Beta;
                      predcat     <- Scale*(effeff*effn)*expMhalf;
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
                      Scale     <- exp(par[3]);
                      Alpha     <- exp(par[4]);
                      Beta      <- exp(par[5]);
                      mccum[1]  <- 0;
                      expM      <- exp(-M);
                      expMhalf  <- exp(-M/2);
                      nstep[1]  <- N0*expM;
                      for(i in 2:sealen)
                        {
                        mccum[i] <- obscat[i-1] + mccum[i-1]*expM;
                        nstep[i] <- N0*(expM^i) - mccum[i]*expMhalf;
                        }
                      effeff      <- obseff^Alpha;
                      effn        <- nstep^Beta;
                      predcat     <- Scale*(effeff*effn)*expMhalf;
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
                    return(negsup);
 }

