CDMN3P.Lik <-
function(par,dates,obscat,obseff,M.fixed,M,distr)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
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
                  if(M.fixed==TRUE){
                    M         <- M;
                    N0        <- exp(par[1]);
                    P1        <- exp(par[2]);
                    P2        <- exp(par[3]);
                    P3        <- exp(par[4]);
                    Scale     <- exp(par[5]);
                    Alpha     <- exp(par[6]);
                    Beta      <- exp(par[7]);
                    expM      <- exp(-M);
                    expMhalf  <- exp(-M/2);
                    mccum[1]  <- 0;
                    nstep[1]  <- N0*expM;
                    for(i in 2:sealen){
                       mccum[i] <- obscat[i-1] + mccum[i-1]*expM;
                       nstep[i] <- N0*(expM^i) +
                                   ind.P1[i]*P1*expM^(i-(ts.P1-ts.start)+1) +
                                   ind.P2[i]*P2*expM^(i-(ts.P2-ts.start)+1) +
                                   ind.P3[i]*P3*expM^(i-(ts.P3-ts.start)+1) -
                                   mccum[i]*expMhalf;
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
                  else {
                   M         <- exp(par[1]);
                   N0        <- exp(par[2]);
                   P1        <- exp(par[3]);
                   P2        <- exp(par[4]);
                   P3        <- exp(par[5]);
                   Scale     <- exp(par[6]);
                   Alpha     <- exp(par[7]);
                   Beta      <- exp(par[8]);
                   expM      <- exp(-M);
                   expMhalf  <- exp(-M/2);
                   mccum[1]  <- 0;
                   nstep[1]  <- N0*expM;
                   for(i in 2:sealen){
                      mccum[i] <- obscat[i-1] + mccum[i-1]*expM;
                      nstep[i] <- N0*(expM^i) +
                                  ind.P1[i]*P1*expM^(i-(ts.P1-ts.start)+1) +
                                  ind.P2[i]*P2*expM^(i-(ts.P2-ts.start)+1) +
                                  ind.P3[i]*P3*expM^(i-(ts.P3-ts.start)+1) -
                                  mccum[i]*expMhalf;
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

