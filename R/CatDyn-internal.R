.onAttach <-
function (lib, pkg) {
    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="CatDyn"),
                              fields=c("Title","Version","Date")))

   packageStartupMessage(paste(
   "--------------------------------------------------------------\n",
   pkg.info["Title"]),"\n",
   paste(" CatDyn version ", pkg.info["Version"],
   " (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
	 "--------------------------------------------------------------\n"
    )
}
.CDMN0P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=0,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale  <- par[3];
                  logalpha  <- par[4];
                  logbeta   <- par[5];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                    mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                    nstep[i] <- exp(logN0)*exp(-exp(logM)*i) - mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN1P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=1,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logscale  <- par[4];
                  logalpha  <- par[5];
                  logbeta   <- par[6];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                    mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                    nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) -
                                mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN2P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=2,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logscale  <- par[5];
                  logalpha  <- par[6];
                  logbeta   <- par[7];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                    mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                    nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) -
                                mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN3P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=3,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                   ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logscale  <- par[6];
                  logalpha  <- par[7];
                  logbeta   <- par[8];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen){
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN4P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=4,partial)
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
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
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
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN5P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=5,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logscale  <- par[8];
                  logalpha  <- par[9];
                  logbeta   <- par[10];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN6P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=6,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logscale  <- par[9];
                  logalpha  <- par[10];
                  logbeta   <- par[11];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN7P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=7,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logscale  <- par[10];
                  logalpha  <- par[11];
                  logbeta   <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN8P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=8,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logscale  <- par[11];
                  logalpha  <- par[12];
                  logbeta   <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN9P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=9,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logscale  <- par[12];
                  logalpha  <- par[13];
                  logbeta   <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN10P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=10,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logscale  <- par[13];
                  logalpha  <- par[14];
                  logbeta   <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN11P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=11,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logscale  <- par[14];
                  logalpha  <- par[15];
                  logbeta   <- par[16];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN12P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=12,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logscale  <- par[15];
                  logalpha  <- par[16];
                  logbeta   <- par[17];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN13P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=13,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logscale  <- par[16];
                  logalpha  <- par[17];
                  logbeta   <- par[18];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN14P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=14,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logscale  <- par[17];
                  logalpha  <- par[18];
                  logbeta   <- par[19];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN15P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=15,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logscale  <- par[18];
                  logalpha  <- par[19];
                  logbeta   <- par[20];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN16P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=16,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logscale  <- par[19];
                  logalpha  <- par[20];
                  logbeta   <- par[21];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN17P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=17,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logscale  <- par[20];
                  logalpha  <- par[21];
                  logbeta   <- par[22];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN18P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=18,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logscale  <- par[21];
                  logalpha  <- par[22];
                  logbeta   <- par[23];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN19P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=19,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logscale  <- par[22];
                  logalpha  <- par[23];
                  logbeta   <- par[24];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN20P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=20,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.P20    <- dates[21];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logscale  <- par[23];
                  logalpha  <- par[24];
                  logbeta   <- par[25];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN21P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=21,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.P20    <- dates[21];
                  ts.P21    <- dates[22];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logP21    <- par[23];
                  logscale  <- par[24];
                  logalpha  <- par[25];
                  logbeta   <- par[26];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start)+1)) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN22P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=22,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.P20    <- dates[21];
                  ts.P21    <- dates[22];
                  ts.P22    <- dates[23];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logP21    <- par[23];
                  logP22    <- par[24];
                  logscale  <- par[25];
                  logalpha  <- par[26];
                  logbeta   <- par[27];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start)+1)) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start)+1)) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN23P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=23,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.P20    <- dates[21];
                  ts.P21    <- dates[22];
                  ts.P22    <- dates[23];
                  ts.P23    <- dates[24];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start), 0, 1);
                  ind.P23   <- ifelse(1:sealen < (ts.P23-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logP21    <- par[23];
                  logP22    <- par[24];
                  logP23    <- par[25];
                  logscale  <- par[26];
                  logalpha  <- par[27];
                  logbeta   <- par[28];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start)+1)) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start)+1)) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start)+1)) +
                                 ind.P23[i]*exp(logP23)*exp(-exp(logM)*(i-(ts.P23-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN24P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=24,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.P20    <- dates[21];
                  ts.P21    <- dates[22];
                  ts.P22    <- dates[23];
                  ts.P23    <- dates[24];
                  ts.P24    <- dates[25];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start), 0, 1);
                  ind.P23   <- ifelse(1:sealen < (ts.P23-ts.start), 0, 1);
                  ind.P24   <- ifelse(1:sealen < (ts.P24-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logP21    <- par[23];
                  logP22    <- par[24];
                  logP23    <- par[25];
                  logP24    <- par[26];
                  logscale  <- par[27];
                  logalpha  <- par[28];
                  logbeta   <- par[29];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start)+1)) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start)+1)) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start)+1)) +
                                 ind.P23[i]*exp(logP23)*exp(-exp(logM)*(i-(ts.P23-ts.start)+1)) +
                                 ind.P24[i]*exp(logP24)*exp(-exp(logM)*(i-(ts.P24-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN25P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=25,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.P20    <- dates[21];
                  ts.P21    <- dates[22];
                  ts.P22    <- dates[23];
                  ts.P23    <- dates[24];
                  ts.P24    <- dates[25];
                  ts.P25    <- dates[26];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start), 0, 1);
                  ind.P23   <- ifelse(1:sealen < (ts.P23-ts.start), 0, 1);
                  ind.P24   <- ifelse(1:sealen < (ts.P24-ts.start), 0, 1);
                  ind.P25   <- ifelse(1:sealen < (ts.P25-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logP21    <- par[23];
                  logP22    <- par[24];
                  logP23    <- par[25];
                  logP24    <- par[26];
                  logP25    <- par[27];
                  logscale  <- par[28];
                  logalpha  <- par[29];
                  logbeta   <- par[30];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start)+1)) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start)+1)) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start)+1)) +
                                 ind.P23[i]*exp(logP23)*exp(-exp(logM)*(i-(ts.P23-ts.start)+1)) +
                                 ind.P24[i]*exp(logP24)*exp(-exp(logM)*(i-(ts.P24-ts.start)+1)) +
                                 ind.P25[i]*exp(logP25)*exp(-exp(logM)*(i-(ts.P25-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT1P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-1,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logscale  <- par[5];
                     logalpha  <- par[6];
                     logbeta   <- par[7];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logscale  <- par[4];
                     logalpha  <- par[5];
                     logbeta   <- par[6];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT2P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-2,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logscale  <- par[7];
                     logalpha  <- par[8];
                     logbeta   <- par[9];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logscale  <- par[5];
                     logalpha  <- par[6];
                     logbeta   <- par[7];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT3P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-3,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logscale  <- par[9];
                     logalpha  <- par[10];
                     logbeta   <- par[11];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logscale  <- par[6];
                     logalpha  <- par[7];
                     logbeta   <- par[8];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT4P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-4,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[6];
                     logQ3     <- par[7];
                     logP4     <- par[8];
                     logQ4     <- par[9];
                     logscale  <- par[10];
                     logalpha  <- par[11];
                     logbeta   <- par[12];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logscale  <- par[7];
                     logalpha  <- par[8];
                     logbeta   <- par[9];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT5P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-5,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10]
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logscale  <- par[13];
                     logalpha  <- par[14];
                     logbeta   <- par[15];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logscale  <- par[8];
                     logalpha  <- par[9];
                     logbeta   <- par[10];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT6P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-6,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logscale  <- par[15];
                     logalpha  <- par[16];
                     logbeta   <- par[17];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logscale  <- par[9];
                     logalpha  <- par[10];
                     logbeta   <- par[11];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT7P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-7,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logscale  <- par[17];
                     logalpha  <- par[18];
                     logbeta   <- par[19];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logscale  <- par[10];
                     logalpha  <- par[11];
                     logbeta   <- par[12];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT8P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-8,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logscale  <- par[19];
                     logalpha  <- par[20];
                     logbeta   <- par[21];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logscale  <- par[11];
                     logalpha  <- par[12];
                     logbeta   <- par[13];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT9P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-9,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logscale  <- par[21];
                     logalpha  <- par[22];
                     logbeta   <- par[23];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logscale  <- par[12];
                     logalpha  <- par[13];
                     logbeta   <- par[14];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT10P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-10,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logscale  <- par[23];
                     logalpha  <- par[24];
                     logbeta   <- par[25];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logscale  <- par[13];
                     logalpha  <- par[14];
                     logbeta   <- par[15];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT11P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-11,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10]
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logscale  <- par[25];
                     logalpha  <- par[26];
                     logbeta   <- par[27];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logscale  <- par[14];
                     logalpha  <- par[15];
                     logbeta   <- par[16];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP10+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT12P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-12,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logscale  <- par[27];
                     logalpha  <- par[28];
                     logbeta   <- par[29];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logscale  <- par[15];
                     logalpha  <- par[16];
                     logbeta   <- par[17];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT13P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-13,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logscale  <- par[29];
                     logalpha  <- par[30];
                     logbeta   <- par[31];
                     }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logscale  <- par[16];
                     logalpha  <- par[17];
                     logbeta   <- par[18];
                     }
                 mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT14P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-14,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ11    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logscale  <- par[31];
                     logalpha  <- par[32];
                     logbeta   <- par[33];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logscale  <- par[17];
                     logalpha  <- par[18];
                     logbeta   <- par[19];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT15P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-15,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logscale  <- par[33];
                     logalpha  <- par[34];
                     logbeta   <- par[35];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logscale  <- par[18];
                     logalpha  <- par[19];
                     logbeta   <- par[20];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT16P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-16,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logscale  <- par[35];
                     logalpha  <- par[36];
                     logbeta   <- par[37];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logscale  <- par[19];
                     logalpha  <- par[20];
                     logbeta   <- par[21];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) -
                                mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT17P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-17,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logscale  <- par[37];
                     logalpha  <- par[38];
                     logbeta   <- par[39];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logscale  <- par[20];
                     logalpha  <- par[21];
                     logbeta   <- par[22];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT18P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-18,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logscale  <- par[39];
                     logalpha  <- par[40];
                     logbeta   <- par[41];
                      }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logscale  <- par[21];
                     logalpha  <- par[22];
                     logbeta   <- par[23];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT19P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-19,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logP19    <- par[39];
                     logQ19    <- par[40];
                     logscale  <- par[41];
                     logalpha  <- par[42];
                     logbeta   <- par[43];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logP19    <- par[21];
                     logQ19    <- 1;
                     logscale  <- par[22];
                     logalpha  <- par[23];
                     logbeta   <- par[24];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)))
                                          - ind.N19[i]*exp((1-partial)*logP19+(partial*logQ19))*exp(-exp(logM)*(i-((1-partial)*ts.P19+partial*ts.N19-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT20P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-20,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.P20    <- dates[40];
                  ts.N20    <- dates[41];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start+1), 0, 1);
                  ind.N20   <- ifelse(1:sealen < (ts.N20-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logP19    <- par[39];
                     logQ19    <- par[40];
                     logP20    <- par[41];
                     logQ20    <- par[42];
                     logscale  <- par[43];
                     logalpha  <- par[44];
                     logbeta   <- par[45];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logP19    <- par[21];
                     logQ19    <- 1;
                     logP20    <- par[22];
                     logQ20    <- 1;
                     logscale  <- par[23];
                     logalpha  <- par[24];
                     logbeta   <- par[25];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)))
                                          - ind.N19[i]*exp((1-partial)*logP19+(partial*logQ19))*exp(-exp(logM)*(i-((1-partial)*ts.P19+partial*ts.N19-ts.start+1)))
                                          - ind.N20[i]*exp((1-partial)*logP20+(partial*logQ20))*exp(-exp(logM)*(i-((1-partial)*ts.P20+partial*ts.N20-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT21P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-21,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.P20    <- dates[40];
                  ts.N20    <- dates[41];
                  ts.P21    <- dates[42];
                  ts.N21    <- dates[43];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start+1), 0, 1);
                  ind.N20   <- ifelse(1:sealen < (ts.N20-ts.start+1), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start+1), 0, 1);
                  ind.N21   <- ifelse(1:sealen < (ts.N21-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logP19    <- par[39];
                     logQ19    <- par[40];
                     logP20    <- par[41];
                     logQ20    <- par[42];
                     logP21    <- par[43];
                     logQ22    <- par[44];
                     logscale  <- par[45];
                     logalpha  <- par[46];
                     logbeta   <- par[47];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logP19    <- par[21];
                     logQ19    <- 1;
                     logP20    <- par[22];
                     logQ20    <- 1;
                     logP21    <- par[23];
                     logQ21    <- 1;
                     logscale  <- par[24];
                     logalpha  <- par[25];
                     logbeta   <- par[26];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1))) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)))
                                          - ind.N19[i]*exp((1-partial)*logP19+(partial*logQ19))*exp(-exp(logM)*(i-((1-partial)*ts.P19+partial*ts.N19-ts.start+1)))
                                          - ind.N20[i]*exp((1-partial)*logP20+(partial*logQ20))*exp(-exp(logM)*(i-((1-partial)*ts.P20+partial*ts.N20-ts.start+1)))
                                          - ind.N21[i]*exp((1-partial)*logP21+(partial*logQ21))*exp(-exp(logM)*(i-((1-partial)*ts.P21+partial*ts.N21-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT22P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-22,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.P20    <- dates[40];
                  ts.N20    <- dates[41];
                  ts.P21    <- dates[42];
                  ts.N21    <- dates[43];
                  ts.P22    <- dates[44];
                  ts.N22    <- dates[45];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start+1), 0, 1);
                  ind.N20   <- ifelse(1:sealen < (ts.N20-ts.start+1), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start+1), 0, 1);
                  ind.N21   <- ifelse(1:sealen < (ts.N21-ts.start+1), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start+1), 0, 1);
                  ind.N22   <- ifelse(1:sealen < (ts.N22-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ9     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logP19    <- par[39];
                     logQ19    <- par[40];
                     logP20    <- par[41];
                     logQ20    <- par[42];
                     logP21    <- par[43];
                     logQ21    <- par[44];
                     logP22    <- par[45];
                     logQ22    <- par[46];
                     logscale  <- par[47];
                     logalpha  <- par[48];
                     logbeta   <- par[49];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logP19    <- par[21];
                     logQ19    <- 1;
                     logP20    <- par[22];
                     logQ20    <- 1;
                     logP21    <- par[23];
                     logQ21    <- 1;
                     logP22    <- par[24];
                     logQ22    <- 1;
                     logscale  <- par[25];
                     logalpha  <- par[26];
                     logbeta   <- par[27];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1))) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start+1))) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)))
                                          - ind.N19[i]*exp((1-partial)*logP19+(partial*logQ19))*exp(-exp(logM)*(i-((1-partial)*ts.P19+partial*ts.N19-ts.start+1)))
                                          - ind.N20[i]*exp((1-partial)*logP20+(partial*logQ20))*exp(-exp(logM)*(i-((1-partial)*ts.P20+partial*ts.N20-ts.start+1)))
                                          - ind.N21[i]*exp((1-partial)*logP21+(partial*logQ21))*exp(-exp(logM)*(i-((1-partial)*ts.P21+partial*ts.N21-ts.start+1)))
                                          - ind.N22[i]*exp((1-partial)*logP22+(partial*logQ22))*exp(-exp(logM)*(i-((1-partial)*ts.P22+partial*ts.N22-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT23P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-23,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.P20    <- dates[40];
                  ts.N20    <- dates[41];
                  ts.P21    <- dates[42];
                  ts.N21    <- dates[43];
                  ts.P22    <- dates[44];
                  ts.N22    <- dates[45];
                  ts.P23    <- dates[46];
                  ts.N23    <- dates[47];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start+1), 0, 1);
                  ind.N20   <- ifelse(1:sealen < (ts.N20-ts.start+1), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start+1), 0, 1);
                  ind.N21   <- ifelse(1:sealen < (ts.N21-ts.start+1), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start+1), 0, 1);
                  ind.N22   <- ifelse(1:sealen < (ts.N22-ts.start+1), 0, 1);
                  ind.P23   <- ifelse(1:sealen < (ts.P23-ts.start+1), 0, 1);
                  ind.N23   <- ifelse(1:sealen < (ts.N23-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logP19    <- par[39];
                     logQ19    <- par[40];
                     logP20    <- par[41];
                     logQ20    <- par[42];
                     logP21    <- par[43];
                     logQ21    <- par[44];
                     logP22    <- par[45];
                     logQ22    <- par[46];
                     logP23    <- par[47];
                     logQ23    <- par[48];
                     logscale  <- par[49];
                     logalpha  <- par[50];
                     logbeta   <- par[51];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logP19    <- par[21];
                     logQ19    <- 1;
                     logP20    <- par[22];
                     logQ20    <- 1;
                     logP21    <- par[23];
                     logQ21    <- 1;
                     logP22    <- par[24];
                     logQ22    <- 1;
                     logP23    <- par[25];
                     logQ23    <- 1;
                     logscale  <- par[26];
                     logalpha  <- par[27];
                     logbeta   <- par[28];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1))) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start+1))) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start+1))) +
                                 ind.P23[i]*exp(logP23)*exp(-exp(logM)*(i-(ts.P23-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)))
                                          - ind.N19[i]*exp((1-partial)*logP19+(partial*logQ19))*exp(-exp(logM)*(i-((1-partial)*ts.P19+partial*ts.N19-ts.start+1)))
                                          - ind.N20[i]*exp((1-partial)*logP20+(partial*logQ20))*exp(-exp(logM)*(i-((1-partial)*ts.P20+partial*ts.N20-ts.start+1)))
                                          - ind.N21[i]*exp((1-partial)*logP21+(partial*logQ21))*exp(-exp(logM)*(i-((1-partial)*ts.P21+partial*ts.N21-ts.start+1)))
                                          - ind.N22[i]*exp((1-partial)*logP22+(partial*logQ22))*exp(-exp(logM)*(i-((1-partial)*ts.P22+partial*ts.N22-ts.start+1)))
                                          - ind.N23[i]*exp((1-partial)*logP23+(partial*logQ23))*exp(-exp(logM)*(i-((1-partial)*ts.P23+partial*ts.N23-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT24P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-24,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.P20    <- dates[40];
                  ts.N20    <- dates[41];
                  ts.P21    <- dates[42];
                  ts.N21    <- dates[43];
                  ts.P22    <- dates[44];
                  ts.N22    <- dates[45];
                  ts.P23    <- dates[46];
                  ts.N23    <- dates[47];
                  ts.P24    <- dates[48];
                  ts.N24    <- dates[49];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start+1), 0, 1);
                  ind.N20   <- ifelse(1:sealen < (ts.N20-ts.start+1), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start+1), 0, 1);
                  ind.N21   <- ifelse(1:sealen < (ts.N21-ts.start+1), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start+1), 0, 1);
                  ind.N22   <- ifelse(1:sealen < (ts.N22-ts.start+1), 0, 1);
                  ind.P23   <- ifelse(1:sealen < (ts.P23-ts.start+1), 0, 1);
                  ind.N23   <- ifelse(1:sealen < (ts.N23-ts.start+1), 0, 1);
                  ind.P24   <- ifelse(1:sealen < (ts.P24-ts.start+1), 0, 1);
                  ind.N24   <- ifelse(1:sealen < (ts.N24-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logP19    <- par[39];
                     logQ19    <- par[40];
                     logP20    <- par[41];
                     logQ20    <- par[42];
                     logP21    <- par[43];
                     logQ21    <- par[44];
                     logP22    <- par[45];
                     logQ22    <- par[46];
                     logP23    <- par[47];
                     logQ23    <- par[48];
                     logP24    <- par[49];
                     logQ24    <- par[50];
                     logscale  <- par[51];
                     logalpha  <- par[52];
                     logbeta   <- par[53];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logP19    <- par[21];
                     logQ19    <- 1;
                     logP20    <- par[22];
                     logQ20    <- 1;
                     logP21    <- par[23];
                     logQ21    <- 1;
                     logP22    <- par[24];
                     logQ22    <- 1;
                     logP23    <- par[25];
                     logQ23    <- 1;
                     logP24    <- par[26];
                     logQ24    <- 1;
                     logscale  <- par[27];
                     logalpha  <- par[28];
                     logbeta   <- par[29];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1))) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start+1))) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start+1))) +
                                 ind.P23[i]*exp(logP23)*exp(-exp(logM)*(i-(ts.P23-ts.start+1))) +
                                 ind.P24[i]*exp(logP24)*exp(-exp(logM)*(i-(ts.P24-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)))
                                          - ind.N19[i]*exp((1-partial)*logP19+(partial*logQ19))*exp(-exp(logM)*(i-((1-partial)*ts.P19+partial*ts.N19-ts.start+1)))
                                          - ind.N20[i]*exp((1-partial)*logP20+(partial*logQ20))*exp(-exp(logM)*(i-((1-partial)*ts.P20+partial*ts.N20-ts.start+1)))
                                          - ind.N21[i]*exp((1-partial)*logP21+(partial*logQ21))*exp(-exp(logM)*(i-((1-partial)*ts.P21+partial*ts.N21-ts.start+1)))
                                          - ind.N22[i]*exp((1-partial)*logP22+(partial*logQ22))*exp(-exp(logM)*(i-((1-partial)*ts.P22+partial*ts.N22-ts.start+1)))
                                          - ind.N23[i]*exp((1-partial)*logP23+(partial*logQ23))*exp(-exp(logM)*(i-((1-partial)*ts.P23+partial*ts.N23-ts.start+1)))
                                          - ind.N24[i]*exp((1-partial)*logP24+(partial*logQ24))*exp(-exp(logM)*(i-((1-partial)*ts.P24+partial*ts.N24-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT25P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-25,partial)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.P20    <- dates[40];
                  ts.N20    <- dates[41];
                  ts.P21    <- dates[42];
                  ts.N21    <- dates[43];
                  ts.P22    <- dates[44];
                  ts.N22    <- dates[45];
                  ts.P23    <- dates[46];
                  ts.N23    <- dates[47];
                  ts.P24    <- dates[48];
                  ts.N24    <- dates[49];
                  ts.P25    <- dates[50];
                  ts.N25    <- dates[51];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start+1), 0, 1);
                  ind.N20   <- ifelse(1:sealen < (ts.N20-ts.start+1), 0, 1);
                  ind.P21   <- ifelse(1:sealen < (ts.P21-ts.start+1), 0, 1);
                  ind.N21   <- ifelse(1:sealen < (ts.N21-ts.start+1), 0, 1);
                  ind.P22   <- ifelse(1:sealen < (ts.P22-ts.start+1), 0, 1);
                  ind.N22   <- ifelse(1:sealen < (ts.N22-ts.start+1), 0, 1);
                  ind.P23   <- ifelse(1:sealen < (ts.P23-ts.start+1), 0, 1);
                  ind.N23   <- ifelse(1:sealen < (ts.N23-ts.start+1), 0, 1);
                  ind.P24   <- ifelse(1:sealen < (ts.P24-ts.start+1), 0, 1);
                  ind.N24   <- ifelse(1:sealen < (ts.N24-ts.start+1), 0, 1);
                  ind.P25   <- ifelse(1:sealen < (ts.P25-ts.start+1), 0, 1);
                  ind.N25   <- ifelse(1:sealen < (ts.N25-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  if(partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- par[4];
                     logP2     <- par[5];
                     logQ2     <- par[6];
                     logP3     <- par[7];
                     logQ3     <- par[8];
                     logP4     <- par[9];
                     logQ4     <- par[10];
                     logP5     <- par[11];
                     logQ5     <- par[12];
                     logP6     <- par[13];
                     logQ6     <- par[14];
                     logP7     <- par[15];
                     logQ7     <- par[16];
                     logP8     <- par[17];
                     logQ8     <- par[18];
                     logP9     <- par[19];
                     logQ9     <- par[20];
                     logP10    <- par[21];
                     logQ10    <- par[22];
                     logP11    <- par[23];
                     logQ11    <- par[24];
                     logP12    <- par[25];
                     logQ12    <- par[26];
                     logP13    <- par[27];
                     logQ13    <- par[28];
                     logP14    <- par[29];
                     logQ14    <- par[30];
                     logP15    <- par[31];
                     logQ15    <- par[32];
                     logP16    <- par[33];
                     logQ16    <- par[34];
                     logP17    <- par[35];
                     logQ17    <- par[36];
                     logP18    <- par[37];
                     logQ18    <- par[38];
                     logP19    <- par[39];
                     logQ19    <- par[40];
                     logP20    <- par[41];
                     logQ20    <- par[42];
                     logP21    <- par[43];
                     logQ21    <- par[44];
                     logP22    <- par[45];
                     logQ22    <- par[46];
                     logP23    <- par[47];
                     logQ23    <- par[48];
                     logP24    <- par[49];
                     logQ24    <- par[50];
                     logQ25    <- par[51];
                     logscale  <- par[52];
                     logalpha  <- par[53];
                     logbeta   <- par[54];
                    }
                  if(!partial)
                    {
                     logP1     <- par[3];
                     logQ1     <- 1;
                     logP2     <- par[4];
                     logQ2     <- 1;
                     logP3     <- par[5];
                     logQ3     <- 1;
                     logP4     <- par[6];
                     logQ4     <- 1;
                     logP5     <- par[7];
                     logQ5     <- 1;
                     logP6     <- par[8];
                     logQ6     <- 1;
                     logP7     <- par[9];
                     logQ7     <- 1;
                     logP8     <- par[10];
                     logQ8     <- 1;
                     logP9     <- par[11];
                     logQ9     <- 1;
                     logP10    <- par[12];
                     logQ10    <- 1;
                     logP11    <- par[13];
                     logQ11    <- 1;
                     logP12    <- par[14];
                     logQ12    <- 1;
                     logP13    <- par[15];
                     logQ13    <- 1;
                     logP14    <- par[16];
                     logQ14    <- 1;
                     logP15    <- par[17];
                     logQ15    <- 1;
                     logP16    <- par[18];
                     logQ16    <- 1;
                     logP17    <- par[19];
                     logQ17    <- 1;
                     logP18    <- par[20];
                     logQ18    <- 1;
                     logP19    <- par[21];
                     logQ19    <- 1;
                     logP20    <- par[22];
                     logQ20    <- 1;
                     logP21    <- par[23];
                     logQ21    <- 1;
                     logP22    <- par[24];
                     logQ22    <- 1;
                     logP23    <- par[25];
                     logQ23    <- 1;
                     logP24    <- par[26];
                     logQ24    <- 1;
                     logP25    <- par[27];
                     logQ25    <- 1;
                     logscale  <- par[28];
                     logalpha  <- par[29];
                     logbeta   <- par[30];
                    }
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1))) +
                                 ind.P21[i]*exp(logP21)*exp(-exp(logM)*(i-(ts.P21-ts.start+1))) +
                                 ind.P22[i]*exp(logP22)*exp(-exp(logM)*(i-(ts.P22-ts.start+1))) +
                                 ind.P23[i]*exp(logP23)*exp(-exp(logM)*(i-(ts.P23-ts.start+1))) +
                                 ind.P24[i]*exp(logP24)*exp(-exp(logM)*(i-(ts.P24-ts.start+1))) +
                                 ind.P25[i]*exp(logP25)*exp(-exp(logM)*(i-(ts.P25-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp((1-partial)*logP1+(partial*logQ1))*exp(-exp(logM)*(i-((1-partial)*ts.P1+partial*ts.N1-ts.start+1)))
                                          - ind.N2[i]*exp((1-partial)*logP2+(partial*logQ2))*exp(-exp(logM)*(i-((1-partial)*ts.P2+partial*ts.N2-ts.start+1)))
                                          - ind.N3[i]*exp((1-partial)*logP3+(partial*logQ3))*exp(-exp(logM)*(i-((1-partial)*ts.P3+partial*ts.N3-ts.start+1)))
                                          - ind.N4[i]*exp((1-partial)*logP4+(partial*logQ4))*exp(-exp(logM)*(i-((1-partial)*ts.P4+partial*ts.N4-ts.start+1)))
                                          - ind.N5[i]*exp((1-partial)*logP5+(partial*logQ5))*exp(-exp(logM)*(i-((1-partial)*ts.P5+partial*ts.N5-ts.start+1)))
                                          - ind.N6[i]*exp((1-partial)*logP6+(partial*logQ6))*exp(-exp(logM)*(i-((1-partial)*ts.P6+partial*ts.N6-ts.start+1)))
                                          - ind.N7[i]*exp((1-partial)*logP7+(partial*logQ7))*exp(-exp(logM)*(i-((1-partial)*ts.P7+partial*ts.N7-ts.start+1)))
                                          - ind.N8[i]*exp((1-partial)*logP8+(partial*logQ8))*exp(-exp(logM)*(i-((1-partial)*ts.P8+partial*ts.N8-ts.start+1)))
                                          - ind.N9[i]*exp((1-partial)*logP9+(partial*logQ9))*exp(-exp(logM)*(i-((1-partial)*ts.P9+partial*ts.N9-ts.start+1)))
                                          - ind.N10[i]*exp((1-partial)*logP10+(partial*logQ10))*exp(-exp(logM)*(i-((1-partial)*ts.P10+partial*ts.N10-ts.start+1)))
                                          - ind.N11[i]*exp((1-partial)*logP11+(partial*logQ11))*exp(-exp(logM)*(i-((1-partial)*ts.P11+partial*ts.N11-ts.start+1)))
                                          - ind.N12[i]*exp((1-partial)*logP12+(partial*logQ12))*exp(-exp(logM)*(i-((1-partial)*ts.P12+partial*ts.N12-ts.start+1)))
                                          - ind.N13[i]*exp((1-partial)*logP13+(partial*logQ13))*exp(-exp(logM)*(i-((1-partial)*ts.P13+partial*ts.N13-ts.start+1)))
                                          - ind.N14[i]*exp((1-partial)*logP14+(partial*logQ14))*exp(-exp(logM)*(i-((1-partial)*ts.P14+partial*ts.N14-ts.start+1)))
                                          - ind.N15[i]*exp((1-partial)*logP15+(partial*logQ15))*exp(-exp(logM)*(i-((1-partial)*ts.P15+partial*ts.N15-ts.start+1)))
                                          - ind.N16[i]*exp((1-partial)*logP16+(partial*logQ16))*exp(-exp(logM)*(i-((1-partial)*ts.P16+partial*ts.N16-ts.start+1)))
                                          - ind.N17[i]*exp((1-partial)*logP17+(partial*logQ17))*exp(-exp(logM)*(i-((1-partial)*ts.P17+partial*ts.N17-ts.start+1)))
                                          - ind.N18[i]*exp((1-partial)*logP18+(partial*logQ18))*exp(-exp(logM)*(i-((1-partial)*ts.P18+partial*ts.N18-ts.start+1)))
                                          - ind.N19[i]*exp((1-partial)*logP19+(partial*logQ19))*exp(-exp(logM)*(i-((1-partial)*ts.P19+partial*ts.N19-ts.start+1)))
                                          - ind.N20[i]*exp((1-partial)*logP20+(partial*logQ20))*exp(-exp(logM)*(i-((1-partial)*ts.P20+partial*ts.N20-ts.start+1)))
                                          - ind.N21[i]*exp((1-partial)*logP21+(partial*logQ21))*exp(-exp(logM)*(i-((1-partial)*ts.P21+partial*ts.N21-ts.start+1)))
                                          - ind.N22[i]*exp((1-partial)*logP22+(partial*logQ22))*exp(-exp(logM)*(i-((1-partial)*ts.P22+partial*ts.N22-ts.start+1)))
                                          - ind.N23[i]*exp((1-partial)*logP23+(partial*logQ23))*exp(-exp(logM)*(i-((1-partial)*ts.P23+partial*ts.N23-ts.start+1)))
                                          - ind.N24[i]*exp((1-partial)*logP24+(partial*logQ24))*exp(-exp(logM)*(i-((1-partial)*ts.P24+partial*ts.N24-ts.start+1)))
                                          - ind.N25[i]*exp((1-partial)*logP25+(partial*logQ25))*exp(-exp(logM)*(i-((1-partial)*ts.P25+partial*ts.N25-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN0P0P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,0),partial)
  {

                  ts.start  <- head(dates,1);
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logscale2 <- par[6];
                  logalpha2 <- par[7];
                  logbeta2  <- par[8];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] +mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) - mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P1P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,1),partial)
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logscale2 <- par[7];
                  logalpha2 <- par[8];
                  logbeta2  <- par[9];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P2P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,2),partial)
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logscale2 <- par[8];
                  logalpha2 <- par[9];
                  logbeta2  <- par[10];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,3),partial)
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts2.P3    <- dates[4];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logP3F2   <- par[8];
                  logscale2 <- par[9];
                  logalpha2 <- par[10];
                  logbeta2  <- par[11];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,4),partial)
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts2.P3    <- dates[4];
                  ts2.P4    <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logP3F2   <- par[8];
                  logP4F2   <- par[9];
                  logscale2 <- par[10];
                  logalpha2 <- par[11];
                  logbeta2  <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,5),partial)
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts2.P3    <- dates[4];
                  ts2.P4    <- dates[5];
                  ts2.P5    <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logP3F2   <- par[8];
                  logP4F2   <- par[9];
                  logP5F2   <- par[10];
                  logscale2 <- par[11];
                  logalpha2 <- par[12];
                  logbeta2  <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P1P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,1),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start+1), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logscale2 <- par[8];
                  logalpha2 <- par[9];
                  logbeta2  <- par[10];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start+1))) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P2P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,2),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logscale2 <- par[9];
                  logalpha2 <- par[10];
                  logbeta2  <- par[11];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,3),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts2.P3    <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logP3F2   <- par[9];
                  logscale2 <- par[10];
                  logalpha2 <- par[11];
                  logbeta2  <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,4),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts2.P3    <- dates[5];
                  ts2.P4    <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logP3F2   <- par[9];
                  logP4F2   <- par[10];
                  logscale2 <- par[11];
                  logalpha2 <- par[12];
                  logbeta2  <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,5),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts2.P3    <- dates[5];
                  ts2.P4    <- dates[6];
                  ts2.P5    <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logP3F2   <- par[9];
                  logP4F2   <- par[10];
                  logP5F2   <- par[11];
                  logscale2 <- par[12];
                  logalpha2 <- par[13];
                  logbeta2  <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P2P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,2),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logscale2 <- par[10];
                  logalpha2 <- par[11];
                  logbeta2  <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,3),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts2.P3    <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logP3F2   <- par[10];
                  logscale2 <- par[11];
                  logalpha2 <- par[12];
                  logbeta2  <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,4),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts2.P3    <- dates[6];
                  ts2.P4    <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logP3F2   <- par[10];
                  logP4F2   <- par[11];
                  logscale2 <- par[12];
                  logalpha2 <- par[13];
                  logbeta2  <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,5),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts2.P3    <- dates[6];
                  ts2.P4    <- dates[7];
                  ts2.P5    <- dates[8];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logP3F2   <- par[10];
                  logP4F2   <- par[11];
                  logP5F2   <- par[12];
                  logscale2 <- par[13];
                  logalpha2 <- par[14];
                  logbeta2  <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN3P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(3,3),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts2.P1    <- dates[5];
                  ts2.P2    <- dates[6];
                  ts2.P3    <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logscale1 <- par[6];
                  logalpha1 <- par[7];
                  logbeta1  <- par[8];
                  logP1F2   <- par[9];
                  logP2F2   <- par[10];
                  logP3F2   <- par[11];
                  logscale2 <- par[12];
                  logalpha2 <- par[13];
                  logbeta2  <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN3P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(3,4),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts2.P1    <- dates[5];
                  ts2.P2    <- dates[6];
                  ts2.P3    <- dates[7];
                  ts2.P4    <- dates[8];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logscale1 <- par[6];
                  logalpha1 <- par[7];
                  logbeta1  <- par[8];
                  logP1F2   <- par[9];
                  logP2F2   <- par[10];
                  logP3F2   <- par[11];
                  logP4F2   <- par[12];
                  logscale2 <- par[13];
                  logalpha2 <- par[14];
                  logbeta2  <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN3P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(3,5),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts2.P1    <- dates[5];
                  ts2.P2    <- dates[6];
                  ts2.P3    <- dates[7];
                  ts2.P4    <- dates[8];
                  ts2.P5    <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logscale1 <- par[6];
                  logalpha1 <- par[7];
                  logbeta1  <- par[8];
                  logP1F2   <- par[9];
                  logP2F2   <- par[10];
                  logP3F2   <- par[11];
                  logP4F2   <- par[12];
                  logP5F2   <- par[13];
                  logscale2 <- par[14];
                  logalpha2 <- par[15];
                  logbeta2  <- par[16];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN4P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(4,4),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts2.P1    <- dates[6];
                  ts2.P2    <- dates[7];
                  ts2.P3    <- dates[8];
                  ts2.P4    <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logscale1 <- par[7];
                  logalpha1 <- par[8];
                  logbeta1  <- par[9];
                  logP1F2   <- par[10];
                  logP2F2   <- par[11];
                  logP3F2   <- par[12];
                  logP4F2   <- par[13];
                  logscale2 <- par[14];
                  logalpha2 <- par[15];
                  logbeta2  <- par[16];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN4P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(4,5),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts2.P1    <- dates[6];
                  ts2.P2    <- dates[7];
                  ts2.P3    <- dates[8];
                  ts2.P4    <- dates[9];
                  ts2.P5    <- dates[10];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logscale1 <- par[7];
                  logalpha1 <- par[8];
                  logbeta1  <- par[9];
                  logP1F2   <- par[10];
                  logP2F2   <- par[11];
                  logP3F2   <- par[12];
                  logP4F2   <- par[13];
                  logP5F2   <- par[14];
                  logscale2 <- par[15];
                  logalpha2 <- par[16];
                  logbeta2  <- par[17];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN5P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(5,5),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts2.P1    <- dates[7];
                  ts2.P2    <- dates[8];
                  ts2.P3    <- dates[9];
                  ts2.P4    <- dates[10];
                  ts2.P5    <- dates[11];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logscale1 <- par[8];
                  logalpha1 <- par[9];
                  logbeta1  <- par[10];
                  logP1F2   <- par[11];
                  logP2F2   <- par[12];
                  logP3F2   <- par[13];
                  logP4F2   <- par[14];
                  logP5F2   <- par[15];
                  logscale2 <- par[16];
                  logalpha2 <- par[17];
                  logbeta2  <- par[18];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN6P6P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(6,6),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts2.P1    <- dates[8];
                  ts2.P2    <- dates[9];
                  ts2.P3    <- dates[10];
                  ts2.P4    <- dates[11];
                  ts2.P5    <- dates[12];
                  ts2.P6    <- dates[13];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logscale1 <- par[9];
                  logalpha1 <- par[10];
                  logbeta1  <- par[11];
                  logP1F2   <- par[12];
                  logP2F2   <- par[13];
                  logP3F2   <- par[14];
                  logP4F2   <- par[15];
                  logP5F2   <- par[16];
                  logP6F2   <- par[17];
                  logscale2 <- par[18];
                  logalpha2 <- par[19];
                  logbeta2  <- par[20];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                  ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                  ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN7P7P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(7,7),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts2.P1    <- dates[9];
                  ts2.P2    <- dates[10];
                  ts2.P3    <- dates[11];
                  ts2.P4    <- dates[12];
                  ts2.P5    <- dates[13];
                  ts2.P6    <- dates[14];
                  ts2.P7    <- dates[15];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logscale1 <- par[10];
                  logalpha1 <- par[11];
                  logbeta1  <- par[12];
                  logP1F2   <- par[13];
                  logP2F2   <- par[14];
                  logP3F2   <- par[15];
                  logP4F2   <- par[16];
                  logP5F2   <- par[17];
                  logP6F2   <- par[18];
                  logP7F2   <- par[19];
                  logscale2 <- par[20];
                  logalpha2 <- par[21];
                  logbeta2  <- par[22];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN8P8P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(8,8),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts2.P1    <- dates[10];
                  ts2.P2    <- dates[11];
                  ts2.P3    <- dates[12];
                  ts2.P4    <- dates[13];
                  ts2.P5    <- dates[14];
                  ts2.P6    <- dates[15];
                  ts2.P7    <- dates[16];
                  ts2.P8    <- dates[17];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logscale1 <- par[11];
                  logalpha1 <- par[12];
                  logbeta1  <- par[13];
                  logP1F2   <- par[14];
                  logP2F2   <- par[15];
                  logP3F2   <- par[16];
                  logP4F2   <- par[17];
                  logP5F2   <- par[18];
                  logP6F2   <- par[19];
                  logP7F2   <- par[20];
                  logP8F2   <- par[21];
                  logscale2 <- par[22];
                  logalpha2 <- par[23];
                  logbeta2  <- par[24];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN9P9P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(9,9),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts2.P1    <- dates[11];
                  ts2.P2    <- dates[12];
                  ts2.P3    <- dates[13];
                  ts2.P4    <- dates[14];
                  ts2.P5    <- dates[15];
                  ts2.P6    <- dates[16];
                  ts2.P7    <- dates[17];
                  ts2.P8    <- dates[18];
                  ts2.P9    <- dates[19];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logscale1 <- par[12];
                  logalpha1 <- par[13];
                  logbeta1  <- par[14];
                  logP1F2   <- par[15];
                  logP2F2   <- par[16];
                  logP3F2   <- par[17];
                  logP4F2   <- par[18];
                  logP5F2   <- par[19];
                  logP6F2   <- par[20];
                  logP7F2   <- par[21];
                  logP8F2   <- par[22];
                  logP9F2   <- par[23];
                  logscale2 <- par[24];
                  logalpha2 <- par[25];
                  logbeta2  <- par[26];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN10P10P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(10,10),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts2.P1    <- dates[12];
                  ts2.P2    <- dates[13];
                  ts2.P3    <- dates[14];
                  ts2.P4    <- dates[15];
                  ts2.P5    <- dates[16];
                  ts2.P6    <- dates[17];
                  ts2.P7    <- dates[18];
                  ts2.P8    <- dates[19];
                  ts2.P9    <- dates[20];
                  ts2.P10   <- dates[21];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                   ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logscale1 <- par[13];
                  logalpha1 <- par[14];
                  logbeta1  <- par[15];
                  logP1F2   <- par[16];
                  logP2F2   <- par[17];
                  logP3F2   <- par[18];
                  logP4F2   <- par[19];
                  logP5F2   <- par[20];
                  logP6F2   <- par[21];
                  logP7F2   <- par[22];
                  logP8F2   <- par[23];
                  logP9F2   <- par[24];
                  logP10F2  <- par[25];
                  logscale2 <- par[26];
                  logalpha2 <- par[27];
                  logbeta2  <- par[28];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN11P11P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(11,11),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts2.P1    <- dates[13];
                  ts2.P2    <- dates[14];
                  ts2.P3    <- dates[15];
                  ts2.P4    <- dates[16];
                  ts2.P5    <- dates[17];
                  ts2.P6    <- dates[18];
                  ts2.P7    <- dates[19];
                  ts2.P8    <- dates[20];
                  ts2.P9    <- dates[21];
                  ts2.P10   <- dates[22];
                  ts2.P11   <- dates[23];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logscale1 <- par[14];
                  logalpha1 <- par[15];
                  logbeta1  <- par[16];
                  logP1F2   <- par[17];
                  logP2F2   <- par[18];
                  logP3F2   <- par[19];
                  logP4F2   <- par[20];
                  logP5F2   <- par[21];
                  logP6F2   <- par[22];
                  logP7F2   <- par[23];
                  logP8F2   <- par[24];
                  logP9F2   <- par[25];
                  logP10F2  <- par[26];
                  logP11F2  <- par[27];
                  logscale2 <- par[28];
                  logalpha2 <- par[29];
                  logbeta2  <- par[30];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN12P12P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(12,12),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts2.P1    <- dates[14];
                  ts2.P2    <- dates[15];
                  ts2.P3    <- dates[16];
                  ts2.P4    <- dates[17];
                  ts2.P5    <- dates[18];
                  ts2.P6    <- dates[19];
                  ts2.P7    <- dates[20];
                  ts2.P8    <- dates[21];
                  ts2.P9    <- dates[22];
                  ts2.P10   <- dates[23];
                  ts2.P11   <- dates[24];
                  ts2.P12   <- dates[25];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logscale1 <- par[15];
                  logalpha1 <- par[16];
                  logbeta1  <- par[17];
                  logP1F2   <- par[18];
                  logP2F2   <- par[19];
                  logP3F2   <- par[20];
                  logP4F2   <- par[21];
                  logP5F2   <- par[22];
                  logP6F2   <- par[23];
                  logP7F2   <- par[24];
                  logP8F2   <- par[25];
                  logP9F2   <- par[26];
                  logP10F2  <- par[27];
                  logP11F2  <- par[28];
                  logP12F2  <- par[29];
                  logscale2 <- par[30];
                  logalpha2 <- par[31];
                  logbeta2  <- par[32];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN13P13P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(13,13),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts2.P1    <- dates[15];
                  ts2.P2    <- dates[16];
                  ts2.P3    <- dates[17];
                  ts2.P4    <- dates[18];
                  ts2.P5    <- dates[19];
                  ts2.P6    <- dates[20];
                  ts2.P7    <- dates[21];
                  ts2.P8    <- dates[22];
                  ts2.P9    <- dates[23];
                  ts2.P10   <- dates[24];
                  ts2.P11   <- dates[25];
                  ts2.P12   <- dates[26];
                  ts2.P13   <- dates[27];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logscale1 <- par[16];
                  logalpha1 <- par[17];
                  logbeta1  <- par[18];
                  logP1F2   <- par[19];
                  logP2F2   <- par[20];
                  logP3F2   <- par[21];
                  logP4F2   <- par[22];
                  logP5F2   <- par[23];
                  logP6F2   <- par[24];
                  logP7F2   <- par[25];
                  logP8F2   <- par[26];
                  logP9F2   <- par[27];
                  logP10F2  <- par[28];
                  logP11F2  <- par[29];
                  logP12F2  <- par[30];
                  logP13F2  <- par[31];
                  logscale2 <- par[32];
                  logalpha2 <- par[33];
                  logbeta2  <- par[34];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN14P14P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(14,14),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts2.P1    <- dates[16];
                  ts2.P2    <- dates[17];
                  ts2.P3    <- dates[18];
                  ts2.P4    <- dates[19];
                  ts2.P5    <- dates[20];
                  ts2.P6    <- dates[21];
                  ts2.P7    <- dates[22];
                  ts2.P8    <- dates[23];
                  ts2.P9    <- dates[24];
                  ts2.P10   <- dates[25];
                  ts2.P11   <- dates[26];
                  ts2.P12   <- dates[27];
                  ts2.P13   <- dates[28];
                  ts2.P14   <- dates[29];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logscale1 <- par[17];
                  logalpha1 <- par[18];
                  logbeta1  <- par[19];
                  logP1F2   <- par[20];
                  logP2F2   <- par[21];
                  logP3F2   <- par[22];
                  logP4F2   <- par[23];
                  logP5F2   <- par[24];
                  logP6F2   <- par[25];
                  logP7F2   <- par[26];
                  logP8F2   <- par[27];
                  logP9F2   <- par[28];
                  logP10F2  <- par[29];
                  logP11F2  <- par[30];
                  logP12F2  <- par[31];
                  logP13F2  <- par[32];
                  logP14F2  <- par[33];
                  logscale2 <- par[34];
                  logalpha2 <- par[35];
                  logbeta2  <- par[36];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN15P15P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(15,15),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts2.P1    <- dates[17];
                  ts2.P2    <- dates[18];
                  ts2.P3    <- dates[19];
                  ts2.P4    <- dates[20];
                  ts2.P5    <- dates[21];
                  ts2.P6    <- dates[22];
                  ts2.P7    <- dates[23];
                  ts2.P8    <- dates[24];
                  ts2.P9    <- dates[25];
                  ts2.P10   <- dates[26];
                  ts2.P11   <- dates[27];
                  ts2.P12   <- dates[28];
                  ts2.P13   <- dates[29];
                  ts2.P14   <- dates[30];
                  ts2.P15   <- dates[31];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logscale1 <- par[18];
                  logalpha1 <- par[19];
                  logbeta1  <- par[20];
                  logP1F2   <- par[21];
                  logP2F2   <- par[22];
                  logP3F2   <- par[23];
                  logP4F2   <- par[24];
                  logP5F2   <- par[25];
                  logP6F2   <- par[26];
                  logP7F2   <- par[27];
                  logP8F2   <- par[28];
                  logP9F2   <- par[29];
                  logP10F2  <- par[30];
                  logP11F2  <- par[31];
                  logP12F2  <- par[32];
                  logP13F2  <- par[33];
                  logP14F2  <- par[34];
                  logP15F2  <- par[35];
                  logscale2 <- par[36];
                  logalpha2 <- par[37];
                  logbeta2  <- par[38];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN16P16P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(16,16),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts2.P1    <- dates[18];
                  ts2.P2    <- dates[19];
                  ts2.P3    <- dates[20];
                  ts2.P4    <- dates[21];
                  ts2.P5    <- dates[22];
                  ts2.P6    <- dates[23];
                  ts2.P7    <- dates[24];
                  ts2.P8    <- dates[25];
                  ts2.P9    <- dates[26];
                  ts2.P10   <- dates[27];
                  ts2.P11   <- dates[28];
                  ts2.P12   <- dates[29];
                  ts2.P13   <- dates[30];
                  ts2.P14   <- dates[31];
                  ts2.P15   <- dates[32];
                  ts2.P16   <- dates[33];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logscale1 <- par[19];
                  logalpha1 <- par[20];
                  logbeta1  <- par[21];
                  logP1F2   <- par[22];
                  logP2F2   <- par[23];
                  logP3F2   <- par[24];
                  logP4F2   <- par[25];
                  logP5F2   <- par[26];
                  logP6F2   <- par[27];
                  logP7F2   <- par[28];
                  logP8F2   <- par[29];
                  logP9F2   <- par[30];
                  logP10F2  <- par[31];
                  logP11F2  <- par[32];
                  logP12F2  <- par[33];
                  logP13F2  <- par[34];
                  logP14F2  <- par[35];
                  logP15F2  <- par[36];
                  logP16F2  <- par[37];
                  logscale2 <- par[38];
                  logalpha2 <- par[39];
                  logbeta2  <- par[40];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN17P17P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(17,17),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts2.P1    <- dates[19];
                  ts2.P2    <- dates[20];
                  ts2.P3    <- dates[21];
                  ts2.P4    <- dates[22];
                  ts2.P5    <- dates[23];
                  ts2.P6    <- dates[24];
                  ts2.P7    <- dates[25];
                  ts2.P8    <- dates[26];
                  ts2.P9    <- dates[27];
                  ts2.P10   <- dates[28];
                  ts2.P11   <- dates[29];
                  ts2.P12   <- dates[30];
                  ts2.P13   <- dates[31];
                  ts2.P14   <- dates[32];
                  ts2.P15   <- dates[33];
                  ts2.P16   <- dates[34];
                  ts2.P17   <- dates[35];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logscale1 <- par[20];
                  logalpha1 <- par[21];
                  logbeta1  <- par[22];
                  logP1F2   <- par[23];
                  logP2F2   <- par[24];
                  logP3F2   <- par[25];
                  logP4F2   <- par[26];
                  logP5F2   <- par[27];
                  logP6F2   <- par[28];
                  logP7F2   <- par[29];
                  logP8F2   <- par[30];
                  logP9F2   <- par[31];
                  logP10F2  <- par[32];
                  logP11F2  <- par[33];
                  logP12F2  <- par[34];
                  logP13F2  <- par[35];
                  logP14F2  <- par[36];
                  logP15F2  <- par[37];
                  logP16F2  <- par[38];
                  logP17F2  <- par[39];
                  logscale2 <- par[40];
                  logalpha2 <- par[41];
                  logbeta2  <- par[42];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN18P18P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(18,18),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts2.P1    <- dates[20];
                  ts2.P2    <- dates[21];
                  ts2.P3    <- dates[22];
                  ts2.P4    <- dates[23];
                  ts2.P5    <- dates[24];
                  ts2.P6    <- dates[25];
                  ts2.P7    <- dates[26];
                  ts2.P8    <- dates[27];
                  ts2.P9    <- dates[28];
                  ts2.P10   <- dates[29];
                  ts2.P11   <- dates[30];
                  ts2.P12   <- dates[31];
                  ts2.P13   <- dates[32];
                  ts2.P14   <- dates[33];
                  ts2.P15   <- dates[34];
                  ts2.P16   <- dates[35];
                  ts2.P17   <- dates[36];
                  ts2.P18   <- dates[37];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logscale1 <- par[21];
                  logalpha1 <- par[22];
                  logbeta1  <- par[23];
                  logP1F2   <- par[24];
                  logP2F2   <- par[25];
                  logP3F2   <- par[26];
                  logP4F2   <- par[27];
                  logP5F2   <- par[28];
                  logP6F2   <- par[29];
                  logP7F2   <- par[30];
                  logP8F2   <- par[31];
                  logP9F2   <- par[32];
                  logP10F2  <- par[33];
                  logP11F2  <- par[34];
                  logP12F2  <- par[35];
                  logP13F2  <- par[36];
                  logP14F2  <- par[37];
                  logP15F2  <- par[38];
                  logP16F2  <- par[39];
                  logP17F2  <- par[40];
                  logP18F2  <- par[41];
                  logscale2 <- par[42];
                  logalpha2 <- par[43];
                  logbeta2  <- par[44];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN19P19P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(19,19),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts2.P1    <- dates[21];
                  ts2.P2    <- dates[22];
                  ts2.P3    <- dates[23];
                  ts2.P4    <- dates[24];
                  ts2.P5    <- dates[25];
                  ts2.P6    <- dates[26];
                  ts2.P7    <- dates[27];
                  ts2.P8    <- dates[28];
                  ts2.P9    <- dates[29];
                  ts2.P10   <- dates[30];
                  ts2.P11   <- dates[31];
                  ts2.P12   <- dates[32];
                  ts2.P13   <- dates[33];
                  ts2.P14   <- dates[34];
                  ts2.P15   <- dates[35];
                  ts2.P16   <- dates[36];
                  ts2.P17   <- dates[37];
                  ts2.P18   <- dates[38];
                  ts2.P19   <- dates[39];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logscale1 <- par[22];
                  logalpha1 <- par[23];
                  logbeta1  <- par[24];
                  logP1F2   <- par[25];
                  logP2F2   <- par[26];
                  logP3F2   <- par[27];
                  logP4F2   <- par[28];
                  logP5F2   <- par[29];
                  logP6F2   <- par[30];
                  logP7F2   <- par[31];
                  logP8F2   <- par[32];
                  logP9F2   <- par[33];
                  logP10F2  <- par[34];
                  logP11F2  <- par[35];
                  logP12F2  <- par[36];
                  logP13F2  <- par[37];
                  logP14F2  <- par[38];
                  logP15F2  <- par[39];
                  logP16F2  <- par[40];
                  logP17F2  <- par[41];
                  logP18F2  <- par[42];
                  logP19F2  <- par[43];
                  logscale2 <- par[44];
                  logalpha2 <- par[45];
                  logbeta2  <- par[46];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN20P20P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(20,20),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts1.P20   <- dates[21];
                  ts2.P1    <- dates[22];
                  ts2.P2    <- dates[23];
                  ts2.P3    <- dates[24];
                  ts2.P4    <- dates[25];
                  ts2.P5    <- dates[26];
                  ts2.P6    <- dates[27];
                  ts2.P7    <- dates[28];
                  ts2.P8    <- dates[29];
                  ts2.P9    <- dates[30];
                  ts2.P10   <- dates[31];
                  ts2.P11   <- dates[32];
                  ts2.P12   <- dates[33];
                  ts2.P13   <- dates[34];
                  ts2.P14   <- dates[35];
                  ts2.P15   <- dates[36];
                  ts2.P16   <- dates[37];
                  ts2.P17   <- dates[38];
                  ts2.P18   <- dates[39];
                  ts2.P19   <- dates[40];
                  ts2.P20   <- dates[41];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  ind1.P20  <- ifelse(1:sealen < (ts1.P20-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  ind2.P20  <- ifelse(1:sealen < (ts2.P20-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logP20F1  <- par[22];
                  logscale1 <- par[23];
                  logalpha1 <- par[24];
                  logbeta1  <- par[25];
                  logP1F2   <- par[26];
                  logP2F2   <- par[27];
                  logP3F2   <- par[28];
                  logP4F2   <- par[29];
                  logP5F2   <- par[30];
                  logP6F2   <- par[31];
                  logP7F2   <- par[32];
                  logP8F2   <- par[33];
                  logP9F2   <- par[34];
                  logP10F2  <- par[35];
                  logP11F2  <- par[36];
                  logP12F2  <- par[37];
                  logP13F2  <- par[38];
                  logP14F2  <- par[39];
                  logP15F2  <- par[40];
                  logP16F2  <- par[41];
                  logP17F2  <- par[42];
                  logP18F2  <- par[43];
                  logP19F2  <- par[44];
                  logP20F2  <- par[45];
                  logscale2 <- par[46];
                  logalpha2 <- par[47];
                  logbeta2  <- par[48];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind1.P20[i]*exp(logP20F1)*exp(-exp(logM)*(i-(ts1.P20-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) +
                                 ind2.P20[i]*exp(logP20F2)*exp(-exp(logM)*(i-(ts2.P20-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN21P21P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(21,21),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts1.P20   <- dates[21];
                  ts1.P21   <- dates[22];
                  ts2.P1    <- dates[23];
                  ts2.P2    <- dates[24];
                  ts2.P3    <- dates[25];
                  ts2.P4    <- dates[26];
                  ts2.P5    <- dates[27];
                  ts2.P6    <- dates[28];
                  ts2.P7    <- dates[29];
                  ts2.P8    <- dates[30];
                  ts2.P9    <- dates[31];
                  ts2.P10   <- dates[32];
                  ts2.P11   <- dates[33];
                  ts2.P12   <- dates[34];
                  ts2.P13   <- dates[35];
                  ts2.P14   <- dates[36];
                  ts2.P15   <- dates[37];
                  ts2.P16   <- dates[38];
                  ts2.P17   <- dates[39];
                  ts2.P18   <- dates[40];
                  ts2.P19   <- dates[41];
                  ts2.P20   <- dates[42];
                  ts2.P21   <- dates[43];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  ind1.P20  <- ifelse(1:sealen < (ts1.P20-ts.start), 0, 1);
                  ind1.P21  <- ifelse(1:sealen < (ts1.P21-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  ind2.P20  <- ifelse(1:sealen < (ts2.P20-ts.start), 0, 1);
                  ind2.P21  <- ifelse(1:sealen < (ts2.P21-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logP20F1  <- par[22];
                  logP21F1  <- par[23];
                  logscale1 <- par[24];
                  logalpha1 <- par[25];
                  logbeta1  <- par[26];
                  logP1F2   <- par[27];
                  logP2F2   <- par[28];
                  logP3F2   <- par[29];
                  logP4F2   <- par[30];
                  logP5F2   <- par[31];
                  logP6F2   <- par[32];
                  logP7F2   <- par[33];
                  logP8F2   <- par[34];
                  logP9F2   <- par[35];
                  logP10F2  <- par[36];
                  logP11F2  <- par[37];
                  logP12F2  <- par[38];
                  logP13F2  <- par[39];
                  logP14F2  <- par[40];
                  logP15F2  <- par[41];
                  logP16F2  <- par[42];
                  logP17F2  <- par[43];
                  logP18F2  <- par[44];
                  logP19F2  <- par[45];
                  logP20F2  <- par[46];
                  logP21F2  <- par[47];
                  logscale2 <- par[48];
                  logalpha2 <- par[49];
                  logbeta2  <- par[50];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind1.P20[i]*exp(logP20F1)*exp(-exp(logM)*(i-(ts1.P20-ts.start)+1)) +
                                 ind1.P21[i]*exp(logP21F1)*exp(-exp(logM)*(i-(ts1.P21-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) +
                                 ind2.P20[i]*exp(logP20F2)*exp(-exp(logM)*(i-(ts2.P20-ts.start)+1)) +
                                 ind2.P21[i]*exp(logP21F2)*exp(-exp(logM)*(i-(ts2.P21-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN22P22P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(22,22),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts1.P20   <- dates[21];
                  ts1.P21   <- dates[22];
                  ts1.P22   <- dates[23];
                  ts2.P1    <- dates[24];
                  ts2.P2    <- dates[25];
                  ts2.P3    <- dates[26];
                  ts2.P4    <- dates[27];
                  ts2.P5    <- dates[28];
                  ts2.P6    <- dates[29];
                  ts2.P7    <- dates[30];
                  ts2.P8    <- dates[31];
                  ts2.P9    <- dates[32];
                  ts2.P10   <- dates[33];
                  ts2.P11   <- dates[34];
                  ts2.P12   <- dates[35];
                  ts2.P13   <- dates[36];
                  ts2.P14   <- dates[37];
                  ts2.P15   <- dates[38];
                  ts2.P16   <- dates[39];
                  ts2.P17   <- dates[40];
                  ts2.P18   <- dates[41];
                  ts2.P19   <- dates[42];
                  ts2.P20   <- dates[43];
                  ts2.P21   <- dates[44];
                  ts2.P22   <- dates[45];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  ind1.P20  <- ifelse(1:sealen < (ts1.P20-ts.start), 0, 1);
                  ind1.P21  <- ifelse(1:sealen < (ts1.P21-ts.start), 0, 1);
                  ind1.P22  <- ifelse(1:sealen < (ts1.P22-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  ind2.P20  <- ifelse(1:sealen < (ts2.P20-ts.start), 0, 1);
                  ind2.P21  <- ifelse(1:sealen < (ts2.P21-ts.start), 0, 1);
                  ind2.P22  <- ifelse(1:sealen < (ts2.P22-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logP20F1  <- par[22];
                  logP21F1  <- par[23];
                  logP22F1  <- par[24];
                  logscale1 <- par[25];
                  logalpha1 <- par[26];
                  logbeta1  <- par[27];
                  logP1F2   <- par[28];
                  logP2F2   <- par[29];
                  logP3F2   <- par[30];
                  logP4F2   <- par[31];
                  logP5F2   <- par[32];
                  logP6F2   <- par[33];
                  logP7F2   <- par[34];
                  logP8F2   <- par[35];
                  logP9F2   <- par[36];
                  logP10F2  <- par[37];
                  logP11F2  <- par[38];
                  logP12F2  <- par[39];
                  logP13F2  <- par[40];
                  logP14F2  <- par[41];
                  logP15F2  <- par[42];
                  logP16F2  <- par[43];
                  logP17F2  <- par[44];
                  logP18F2  <- par[45];
                  logP19F2  <- par[46];
                  logP20F2  <- par[47];
                  logP21F2  <- par[48];
                  logP22F2  <- par[49];
                  logscale2 <- par[50];
                  logalpha2 <- par[51];
                  logbeta2  <- par[52];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind1.P20[i]*exp(logP20F1)*exp(-exp(logM)*(i-(ts1.P20-ts.start)+1)) +
                                 ind1.P21[i]*exp(logP21F1)*exp(-exp(logM)*(i-(ts1.P21-ts.start)+1)) +
                                 ind1.P22[i]*exp(logP22F1)*exp(-exp(logM)*(i-(ts1.P22-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) +
                                 ind2.P20[i]*exp(logP20F2)*exp(-exp(logM)*(i-(ts2.P20-ts.start)+1)) +
                                 ind2.P21[i]*exp(logP21F2)*exp(-exp(logM)*(i-(ts2.P21-ts.start)+1)) +
                                 ind2.P22[i]*exp(logP22F2)*exp(-exp(logM)*(i-(ts2.P22-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN23P23P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(23,23),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts1.P20   <- dates[21];
                  ts1.P21   <- dates[22];
                  ts1.P22   <- dates[23];
                  ts1.P23   <- dates[24];
                  ts2.P1    <- dates[25];
                  ts2.P2    <- dates[26];
                  ts2.P3    <- dates[27];
                  ts2.P4    <- dates[28];
                  ts2.P5    <- dates[29];
                  ts2.P6    <- dates[30];
                  ts2.P7    <- dates[31];
                  ts2.P8    <- dates[32];
                  ts2.P9    <- dates[33];
                  ts2.P10   <- dates[34];
                  ts2.P11   <- dates[35];
                  ts2.P12   <- dates[36];
                  ts2.P13   <- dates[37];
                  ts2.P14   <- dates[38];
                  ts2.P15   <- dates[39];
                  ts2.P16   <- dates[40];
                  ts2.P17   <- dates[41];
                  ts2.P18   <- dates[42];
                  ts2.P19   <- dates[43];
                  ts2.P20   <- dates[44];
                  ts2.P21   <- dates[45];
                  ts2.P22   <- dates[46];
                  ts2.P23   <- dates[47];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  ind1.P20  <- ifelse(1:sealen < (ts1.P20-ts.start), 0, 1);
                  ind1.P21  <- ifelse(1:sealen < (ts1.P21-ts.start), 0, 1);
                  ind1.P22  <- ifelse(1:sealen < (ts1.P22-ts.start), 0, 1);
                  ind1.P23  <- ifelse(1:sealen < (ts1.P23-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  ind2.P20  <- ifelse(1:sealen < (ts2.P20-ts.start), 0, 1);
                  ind2.P21  <- ifelse(1:sealen < (ts2.P21-ts.start), 0, 1);
                  ind2.P22  <- ifelse(1:sealen < (ts2.P22-ts.start), 0, 1);
                  ind2.P23  <- ifelse(1:sealen < (ts2.P23-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logP20F1  <- par[22];
                  logP21F1  <- par[23];
                  logP22F1  <- par[24];
                  logP23F1  <- par[25];
                  logscale1 <- par[26];
                  logalpha1 <- par[27];
                  logbeta1  <- par[28];
                  logP1F2   <- par[29];
                  logP2F2   <- par[30];
                  logP3F2   <- par[31];
                  logP4F2   <- par[32];
                  logP5F2   <- par[33];
                  logP6F2   <- par[34];
                  logP7F2   <- par[35];
                  logP8F2   <- par[36];
                  logP9F2   <- par[37];
                  logP10F2  <- par[38];
                  logP11F2  <- par[39];
                  logP12F2  <- par[40];
                  logP13F2  <- par[41];
                  logP14F2  <- par[42];
                  logP15F2  <- par[43];
                  logP16F2  <- par[44];
                  logP17F2  <- par[45];
                  logP18F2  <- par[46];
                  logP19F2  <- par[47];
                  logP20F2  <- par[48];
                  logP21F2  <- par[49];
                  logP22F2  <- par[50];
                  logP23F2  <- par[51];
                  logscale2 <- par[52];
                  logalpha2 <- par[53];
                  logbeta2  <- par[54];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind1.P20[i]*exp(logP20F1)*exp(-exp(logM)*(i-(ts1.P20-ts.start)+1)) +
                                 ind1.P21[i]*exp(logP21F1)*exp(-exp(logM)*(i-(ts1.P21-ts.start)+1)) +
                                 ind1.P22[i]*exp(logP22F1)*exp(-exp(logM)*(i-(ts1.P22-ts.start)+1)) +
                                 ind1.P23[i]*exp(logP23F1)*exp(-exp(logM)*(i-(ts1.P23-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) +
                                 ind2.P20[i]*exp(logP20F2)*exp(-exp(logM)*(i-(ts2.P20-ts.start)+1)) +
                                 ind2.P21[i]*exp(logP21F2)*exp(-exp(logM)*(i-(ts2.P21-ts.start)+1)) +
                                 ind2.P22[i]*exp(logP22F2)*exp(-exp(logM)*(i-(ts2.P22-ts.start)+1)) +
                                 ind2.P23[i]*exp(logP23F2)*exp(-exp(logM)*(i-(ts2.P23-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN24P24P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(24,24),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts1.P20   <- dates[21];
                  ts1.P21   <- dates[22];
                  ts1.P22   <- dates[23];
                  ts1.P23   <- dates[24];
                  ts1.P24   <- dates[25];
                  ts2.P1    <- dates[26];
                  ts2.P2    <- dates[27];
                  ts2.P3    <- dates[28];
                  ts2.P4    <- dates[29];
                  ts2.P5    <- dates[30];
                  ts2.P6    <- dates[31];
                  ts2.P7    <- dates[32];
                  ts2.P8    <- dates[33];
                  ts2.P9    <- dates[34];
                  ts2.P10   <- dates[35];
                  ts2.P11   <- dates[36];
                  ts2.P12   <- dates[37];
                  ts2.P13   <- dates[38];
                  ts2.P14   <- dates[39];
                  ts2.P15   <- dates[40];
                  ts2.P16   <- dates[41];
                  ts2.P17   <- dates[42];
                  ts2.P18   <- dates[43];
                  ts2.P19   <- dates[44];
                  ts2.P20   <- dates[45];
                  ts2.P21   <- dates[46];
                  ts2.P22   <- dates[47];
                  ts2.P23   <- dates[48];
                  ts2.P24   <- dates[49];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  ind1.P20  <- ifelse(1:sealen < (ts1.P20-ts.start), 0, 1);
                  ind1.P21  <- ifelse(1:sealen < (ts1.P21-ts.start), 0, 1);
                  ind1.P22  <- ifelse(1:sealen < (ts1.P22-ts.start), 0, 1);
                  ind1.P23  <- ifelse(1:sealen < (ts1.P23-ts.start), 0, 1);
                  ind1.P24  <- ifelse(1:sealen < (ts1.P24-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  ind2.P20  <- ifelse(1:sealen < (ts2.P20-ts.start), 0, 1);
                  ind2.P21  <- ifelse(1:sealen < (ts2.P21-ts.start), 0, 1);
                  ind2.P22  <- ifelse(1:sealen < (ts2.P22-ts.start), 0, 1);
                  ind2.P23  <- ifelse(1:sealen < (ts2.P23-ts.start), 0, 1);
                  ind2.P24  <- ifelse(1:sealen < (ts2.P24-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logP20F1  <- par[22];
                  logP21F1  <- par[23];
                  logP22F1  <- par[24];
                  logP23F1  <- par[25];
                  logP24F1  <- par[26];
                  logscale1 <- par[27];
                  logalpha1 <- par[28];
                  logbeta1  <- par[29];
                  logP1F2   <- par[30];
                  logP2F2   <- par[31];
                  logP3F2   <- par[32];
                  logP4F2   <- par[33];
                  logP5F2   <- par[34];
                  logP6F2   <- par[35];
                  logP7F2   <- par[36];
                  logP8F2   <- par[37];
                  logP9F2   <- par[38];
                  logP10F2  <- par[39];
                  logP11F2  <- par[40];
                  logP12F2  <- par[41];
                  logP13F2  <- par[42];
                  logP14F2  <- par[43];
                  logP15F2  <- par[44];
                  logP16F2  <- par[45];
                  logP17F2  <- par[46];
                  logP18F2  <- par[47];
                  logP19F2  <- par[48];
                  logP20F2  <- par[49];
                  logP21F2  <- par[50];
                  logP22F2  <- par[51];
                  logP23F2  <- par[52];
                  logP24F2  <- par[53];
                  logscale2 <- par[54];
                  logalpha2 <- par[55];
                  logbeta2  <- par[56];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind1.P20[i]*exp(logP20F1)*exp(-exp(logM)*(i-(ts1.P20-ts.start)+1)) +
                                 ind1.P21[i]*exp(logP21F1)*exp(-exp(logM)*(i-(ts1.P21-ts.start)+1)) +
                                 ind1.P22[i]*exp(logP22F1)*exp(-exp(logM)*(i-(ts1.P22-ts.start)+1)) +
                                 ind1.P23[i]*exp(logP23F1)*exp(-exp(logM)*(i-(ts1.P23-ts.start)+1)) +
                                 ind1.P24[i]*exp(logP24F1)*exp(-exp(logM)*(i-(ts1.P24-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) +
                                 ind2.P20[i]*exp(logP20F2)*exp(-exp(logM)*(i-(ts2.P20-ts.start)+1)) +
                                 ind2.P21[i]*exp(logP21F2)*exp(-exp(logM)*(i-(ts2.P21-ts.start)+1)) +
                                 ind2.P22[i]*exp(logP22F2)*exp(-exp(logM)*(i-(ts2.P22-ts.start)+1)) +
                                 ind2.P23[i]*exp(logP23F2)*exp(-exp(logM)*(i-(ts2.P23-ts.start)+1)) +
                                 ind2.P24[i]*exp(logP24F2)*exp(-exp(logM)*(i-(ts2.P24-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN25P25P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(25,25),partial)
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts1.P20   <- dates[21];
                  ts1.P21   <- dates[22];
                  ts1.P22   <- dates[23];
                  ts1.P23   <- dates[24];
                  ts1.P24   <- dates[25];
                  ts1.P25   <- dates[26];
                  ts2.P1    <- dates[27];
                  ts2.P2    <- dates[28];
                  ts2.P3    <- dates[29];
                  ts2.P4    <- dates[30];
                  ts2.P5    <- dates[31];
                  ts2.P6    <- dates[32];
                  ts2.P7    <- dates[33];
                  ts2.P8    <- dates[34];
                  ts2.P9    <- dates[35];
                  ts2.P10   <- dates[36];
                  ts2.P11   <- dates[37];
                  ts2.P12   <- dates[38];
                  ts2.P13   <- dates[39];
                  ts2.P14   <- dates[40];
                  ts2.P15   <- dates[41];
                  ts2.P16   <- dates[42];
                  ts2.P17   <- dates[43];
                  ts2.P18   <- dates[44];
                  ts2.P19   <- dates[45];
                  ts2.P20   <- dates[46];
                  ts2.P21   <- dates[47];
                  ts2.P22   <- dates[48];
                  ts2.P23   <- dates[49];
                  ts2.P24   <- dates[50];
                  ts2.P25   <- dates[51];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  ind1.P20  <- ifelse(1:sealen < (ts1.P20-ts.start), 0, 1);
                  ind1.P21  <- ifelse(1:sealen < (ts1.P21-ts.start), 0, 1);
                  ind1.P22  <- ifelse(1:sealen < (ts1.P22-ts.start), 0, 1);
                  ind1.P23  <- ifelse(1:sealen < (ts1.P23-ts.start), 0, 1);
                  ind1.P24  <- ifelse(1:sealen < (ts1.P24-ts.start), 0, 1);
                  ind1.P25  <- ifelse(1:sealen < (ts1.P25-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  ind2.P20  <- ifelse(1:sealen < (ts2.P20-ts.start), 0, 1);
                  ind2.P21  <- ifelse(1:sealen < (ts2.P21-ts.start), 0, 1);
                  ind2.P22  <- ifelse(1:sealen < (ts2.P22-ts.start), 0, 1);
                  ind2.P23  <- ifelse(1:sealen < (ts2.P23-ts.start), 0, 1);
                  ind2.P24  <- ifelse(1:sealen < (ts2.P24-ts.start), 0, 1);
                  ind2.P25  <- ifelse(1:sealen < (ts2.P25-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logP20F1  <- par[22];
                  logP21F1  <- par[23];
                  logP22F1  <- par[24];
                  logP23F1  <- par[25];
                  logP24F1  <- par[26];
                  logP25F1  <- par[27];
                  logscale1 <- par[28];
                  logalpha1 <- par[29];
                  logbeta1  <- par[30];
                  logP1F2   <- par[31];
                  logP2F2   <- par[32];
                  logP3F2   <- par[33];
                  logP4F2   <- par[34];
                  logP5F2   <- par[35];
                  logP6F2   <- par[36];
                  logP7F2   <- par[37];
                  logP8F2   <- par[38];
                  logP9F2   <- par[39];
                  logP10F2  <- par[40];
                  logP11F2  <- par[41];
                  logP12F2  <- par[42];
                  logP13F2  <- par[43];
                  logP14F2  <- par[44];
                  logP15F2  <- par[45];
                  logP16F2  <- par[46];
                  logP17F2  <- par[47];
                  logP18F2  <- par[48];
                  logP19F2  <- par[49];
                  logP20F2  <- par[50];
                  logP21F2  <- par[51];
                  logP22F2  <- par[52];
                  logP23F2  <- par[53];
                  logP24F2  <- par[54];
                  logP25F2  <- par[55];
                  logscale2 <- par[56];
                  logalpha2 <- par[57];
                  logbeta2  <- par[58];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind1.P20[i]*exp(logP20F1)*exp(-exp(logM)*(i-(ts1.P20-ts.start)+1)) +
                                 ind1.P21[i]*exp(logP21F1)*exp(-exp(logM)*(i-(ts1.P21-ts.start)+1)) +
                                 ind1.P22[i]*exp(logP22F1)*exp(-exp(logM)*(i-(ts1.P22-ts.start)+1)) +
                                 ind1.P23[i]*exp(logP23F1)*exp(-exp(logM)*(i-(ts1.P23-ts.start)+1)) +
                                 ind1.P24[i]*exp(logP24F1)*exp(-exp(logM)*(i-(ts1.P24-ts.start)+1)) +
                                 ind1.P25[i]*exp(logP25F1)*exp(-exp(logM)*(i-(ts1.P25-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) +
                                 ind2.P20[i]*exp(logP20F2)*exp(-exp(logM)*(i-(ts2.P20-ts.start)+1)) +
                                 ind2.P21[i]*exp(logP21F2)*exp(-exp(logM)*(i-(ts2.P21-ts.start)+1)) +
                                 ind2.P22[i]*exp(logP22F2)*exp(-exp(logM)*(i-(ts2.P22-ts.start)+1)) +
                                 ind2.P23[i]*exp(logP23F2)*exp(-exp(logM)*(i-(ts2.P23-ts.start)+1)) +
                                 ind2.P24[i]*exp(logP24F2)*exp(-exp(logM)*(i-(ts2.P24-ts.start)+1)) +
                                 ind2.P25[i]*exp(logP25F2)*exp(-exp(logM)*(i-(ts2.P25-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CatDynLik1F <-
function(obscat1,predcat1,distr,par) # <- put par in the arguments
  {
                  Likel        <- vector("list",4);
                  names(Likel) <- c("Dispersion","Deviance","Likelihood","DevianceResidual");
                  sealen       <- length(obscat1);
                  #1/9
                  if(distr=='poisson')
                    {
                     psi1        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  #2/9
                  else if(distr=='negbin')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  #3/9
                  else if(distr=='normal')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- (obscat1-predcat1)^2;
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  #4/9
                  else if(distr=='apnormal')
                    {
                     psi1        <- NA;
                     dev1        <- NA;
                     likcontr1   <- (obscat1-predcat1)^2;
                     res1        <- obscat1-predcat1;
                    }
                  #5/9
                  else if(distr=='lognormal')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  #6/9
                  else if(distr=='aplnormal')
                    {
                     psi1        <- NA;
                     dev1        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                    }
                  #7/9
                  else if(distr=='gamma')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-(2/psi1)*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  #8/9
                  else if(distr=='roblognormal')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-log(psi1)+log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  #9/9
                  else if(distr=='gumbel')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(-1+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-log(psi1)-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  Likel[["Dispersion"]]       <- psi1
                  Likel[["Deviance"]]         <- dev1
                  Likel[["Likelihood"]]       <- likcontr1
                  Likel[["DevianceResidual"]] <- res1
                  return(Likel);
  }
.CatDynLik2F <-
function(obscat1,predcat1,obscat2,predcat2,distr,par)
  {
                  Likel        <- vector("list",4);
                  names(Likel) <- c("Dispersion","Deviance","Likelihood","DevianceResidual");
                  sealen       <- length(obscat1);
                  #1/81
                  if(sum(distr==c('poisson','poisson')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #2/81
                  else if(sum(distr==c('poisson','negbin')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #3/81
                  else if(sum(distr==c('poisson','normal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- (obscat2-predcat2)^2;
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #4/81
                  else if(sum(distr==c('poisson','apnormal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- NA;
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #5/81
                  else if(sum(distr==c('poisson','lognormal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #6/81
                  else if(sum(distr==c('poisson','aplnormal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- NA;
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #7/81
                  else if(sum(distr==c('poisson','gamma')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #8/81
                  else if(sum(distr==c('poisson','roblognormal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #9/81
                  else if(sum(distr==c('poisson','gumbel')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #10/81
                  else if(sum(distr==c('negbin','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #11/81
                  else if(sum(distr==c('negbin','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #12/81
                  else if(sum(distr==c('negbin', 'normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #13/81
                  else if(sum(distr==c('negbin', 'apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- obscat2-predcat2;
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- dev2^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- dev2;
                    }
                  #14/81
                  else if(sum(distr==c('negbin', 'lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #15/81
                  else if(sum(distr==c('negbin', 'aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- NA;
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- log(obscat2)-log(predcat2);
                    }
                  #16/81
                  else if(sum(distr==c('negbin', 'gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #17/81
                  else if(sum(distr==c('negbin', 'roblognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #18/81
                  else if(sum(distr==c('negbin', 'gumbel')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #19/81
                  else if(sum(distr==c('normal','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ((obscat1-predcat1)^2)/sealen;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #20/81
                  else if(sum(distr==c('normal', 'negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #21/81
                  else if(sum(distr==c('normal','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #22/81
                  else if(sum(distr==c('normal','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- NA;
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #23/81
                  else if(sum(distr==c('normal','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #24/81
                  else if(sum(distr==c('normal','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- NA;
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #25/81
                  else if(sum(distr==c('normal','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #26/81
                  else if(sum(distr==c('normal','roblognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #27/81
                  else if(sum(distr==c('normal','gumbel')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #28/81
                  else if(sum(distr==c('apnormal','poisson')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- 1;
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #29/81
                  else if(sum(distr==c('apnormal','negbin')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat1/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #30/81
                  else if(sum(distr==c('apnormal','normal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ((obscat2-predcat2)^2)/sealen;
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #31/81
                  else if(sum(distr==c('apnormal','apnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- obscat1-predcat1;
                     res2        <- obscat2-predcat2;
                    }
                  #32/81
                  else if(sum(distr==c('apnormal','lognormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #33/81
                  else if(sum(distr==c('apnormal','aplnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- obscat1-predcat1;
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #34/81
                  else if(sum(distr==c('apnormal','gamma')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #35/81
                  else if(sum(distr==c('apnormal','roblognormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #36/81
                  else if(sum(distr==c('apnormal','gumbel')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #37/81
                  else if(sum(distr==c('lognormal', 'poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2)/sealen;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #38/81
                  else if(sum(distr==c('lognormal','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #39/81
                  else if(sum(distr==c('lognormal','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ((obscat2-predcat2)^2)/sealen;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #40/81
                  else if(sum(distr==c('lognormal','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(sealen/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #41/81
                  else if(sum(distr==c('lognormal','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #42/81
                  else if(sum(distr==c('lognormal','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #43/81
                  else if(sum(distr==c('lognormal','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #44/81
                  else if(sum(distr==c('lognormal','roblognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #45/81
                  else if(sum(distr==c('lognormal','gumbel')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #46/81
                  else if(sum(distr==c('aplnormal','poisson')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- 1;
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #47/81
                  else if(sum(distr==c('aplnormal','negbin')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #48/81
                  else if(sum(distr==c('aplnormal','normal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #49/81
                  else if(sum(distr==c('aplnormal','apnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- obscat2-predcat2;
                    }
                  #50/81
                  else if(sum(distr==c('aplnormal','lognormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #51/81
                  else if(sum(distr==c('aplnormal','aplnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #52/81
                  else if(sum(distr==c('aplnormal','gamma')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #53/81
                  else if(sum(distr==c('aplnormal','roblognormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #54/81
                  else if(sum(distr==c('aplnormal','gumbel')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #55/81
                  else if(sum(distr==c('gamma','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #56/81
                  else if(sum(distr==c('gamma','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #57/81
                  else if(sum(distr==c('gamma','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #58/81
                  else if(sum(distr==c('gamma','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #59/81
                  else if(sum( distr==c('gamma','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #60/81
                  else if(sum(distr==c('gamma','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #61/81
                  else if(sum(distr==c('gamma','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #62/81
                  else if(sum( distr==c('gamma','roblognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #63/81
                  else if(sum( distr==c('gamma','gumbel')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #64/81
                  else if(sum(distr==c('roblognormal','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #65/81
                  else if(sum(distr==c('roblognormal','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #66/81
                  else if(sum(distr==c('roblognormal','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #67/81
                  else if(sum(distr==c('roblognormal','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #68/81
                  else if(sum( distr==c('roblognormal','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #69/81
                  else if(sum(distr==c('roblognormal','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #70/81
                  else if(sum(distr==c('roblognormal','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #71/81
                  else if(sum( distr==c('roblognormal','roblognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #72/81
                  else if(sum( distr==c('roblognormal','gumbel')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)-log(exp(-0.5*(log(obscat1/predcat1)/psi1+0.5*psi1)^2)+0.01));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #73/81
                  else if(sum(distr==c('gumbel','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #74/81
                  else if(sum(distr==c('gumbel','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #75/81
                  else if(sum(distr==c('gumbel','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #76/81
                  else if(sum(distr==c('gumbel','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #77/81
                  else if(sum( distr==c('gumbel','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #78/81
                  else if(sum(distr==c('gumbel','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #79/81
                  else if(sum(distr==c('gumbel','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #80/81
                  else if(sum( distr==c('gumbel','roblognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)-log(exp(-0.5*(log(obscat2/predcat2)/psi2+0.5*psi2)^2)+0.01));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #81/81
                  else if(sum( distr==c('gumbel','gumbel')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(1-(obscat1-predcat1)/psi1-exp(-(obscat1-predcat1)/psi1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(1-(obscat2-predcat2)/psi2-exp(-(obscat2-predcat2)/psi2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,log(psi1)+(obscat1-predcat1)/psi1+exp(-(obscat1-predcat1)/psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,log(psi2)+(obscat2-predcat2)/psi2+exp(-(obscat2-predcat2)/psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  Likel[[1]] <- c(psi1,psi2)
                  Likel[[2]] <- data.frame(fleet1=dev1,fleet2=dev2)
                  Likel[[3]] <- data.frame(fleet1=likcontr1,fleet2=likcontr2)
                  Likel[[4]] <- data.frame(fleet1=res1,fleet2=res2)
                  return(Likel);
  }
.CatDynExp1F.Res <-
function(properties,nstep,obsmbm1,ppp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,partial)
  {
                  partial   <- partial
                  sdistr    <- c("poisson","apnormal","aplnormal")
                  biom      <- vector("numeric",length(nstep))
                  obscat.w  <- vector("numeric",length(nstep))
                  predcat.w <- vector("numeric",length(nstep))
                  if(properties$Units[2] == "ind")
                    {
                     biom      <- NA
                     obscat.w  <- NA
                     predcat.w <- NA
                    }
                  if(properties$Units[2] != "ind")
                    {
                     #Set biom in tonnes for both units of individual body weight
                     #Set obscat.w, predcat.w in kg for both units of individual body weight
                     biom <- obsmbm1*(1e-6*(properties$Units[3] == "g")+
                                      1e-3*(properties$Units[3] == "kg"))*
                             nstep*(1e9*(properties$Units[4] == "bill")+
                                    1e6*(properties$Units[4] == "mill")+
                                    1e3*(properties$Units[4] == "thou")+
                                    1e2*(properties$Units[4] == "hund"))
                     obscat.w <- obsmbm1*(1e-3*(properties$Units[3] == "g")+
                                          1e0*(properties$Units[3] == "kg"))*
                                 obscat1*(1e9*(properties$Units[4] == "bill")+
                                          1e6*(properties$Units[4] == "mill")+
                                          1e3*(properties$Units[4] == "thou")+
                                          1e2*(properties$Units[4] == "hund"))
                     predcat.w <- obsmbm1*(1e-3*(properties$Units[3] == "g")+
                                           1e0*(properties$Units[3] == "kg"))*
                                  predcat1*(1e9*(properties$Units[4] == "bill")+
                                            1e6*(properties$Units[4] == "mill")+
                                            1e3*(properties$Units[4] == "thou")+
                                            1e2*(properties$Units[4] == "hund"))
                    }
                  catdynres                         <- vector("list",2);
                  names(catdynres)                  <- c("Properties","Model");
                  catdynres$Properties              <- properties;
                  catdynres$Model                   <- vector("list",5);
                  names(catdynres$Model)            <- c("Type","Dates","Distribution","Parameters","Results");
                  catdynres$Model$Type              <- ppp;
                  catdynres$Model$Dates             <- dates;
                  catdynres$Model$Distribution      <- distr;
                  catdynres$Model$Parameters        <- exp(par);
                  if(ppp==0)
                    {
                     if(distr %in% sdistr)
                       {
                        names(catdynres$Model$Parameters) <- c("M","N0","k","alpha","beta")
                       }
                     if(!distr %in% sdistr)
                       {
                        names(catdynres$Model$Parameters) <- c("M","N0","k","alpha","beta","psi")
                       }
                    }
                  if(ppp>0)
                    {
                     if(distr %in% sdistr)
                       {
                        names(catdynres$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta")
                       }
                     if(!distr %in% sdistr)
                       {
                        names(catdynres$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta","psi")
                       }
                    }
                  if(ppp<0)
                    {
                     ppp <- abs(ppp);
                     if(distr %in% sdistr)
                       {
                        if(!partial)
                          {
                           names(catdynres$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta")
                          }
                        if(partial)
                          {
                           names(catdynres$Model$Parameters) <- c("M","N0",paste(c("P","Q"),sort(rep(1:ppp,2)),sep=""),"k","alpha","beta")
                          }
                       }
                     if(!distr %in% sdistr)
                       {
                        if(!partial)
                          {
                           names(catdynres$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta","psi")
                          }
                        if(partial)
                          {
                           names(catdynres$Model$Parameters) <- c("M","N0",paste(c("P","Q"),sort(rep(1:ppp,2)),sep=""),"k","alpha","beta","psi")
                          }
                       }
                    }
                  catdynres$Model$Results           <- data.frame(period=ts.start:ts.end,
                                                                  obseff1=obseff1,
                                                                  obscat.n=obscat1,
                                                                  predcat.n=predcat1,
                                                                  obscat.w=obscat.w,
                                                                  predcat.w=predcat.w,
                                                                  deviance=Likel[[2]],
                                                                  likcontr=Likel[[3]],
                                                                  devresid=Likel[[4]],
                                                                  npred=nstep,
                                                                  biompred=biom,
                                                                  obsexplot=NA,
                                                                  predexplot=NA,
                                                                  F.obsexplot=NA,
                                                                  F.predexplot=NA);
                  if(properties$Units[2] == "ind")
                    {
                     catdynres$Model$Results$obsexplot  <- obscat1/nstep
                     catdynres$Model$Results$predexplot <- predcat1/nstep
                     obsexplot.n                        <- obscat1/nstep
                     predexplot.n                       <- predcat1/nstep
                    }
                  if(properties$Units[2] != "ind")
                    {
                     catdynres$Model$Results$obsexplot   <- 1e-3*obscat.w/biom
                     catdynres$Model$Results$predexplot  <- 1e-3*predcat.w/biom
                     obsexplot.n                         <- obscat1/nstep
                     predexplot.n                        <- predcat1/nstep
                    }
                  for(t in 1:(ts.end-ts.start+1))
                    {
                     #if(biom[t]==0)
                     #  {
                     #   catdynres$Model$Results$biompred[t]  <- NA
                     #   catdynres$Model$Results$obsexplot[t]  <- NA
                     #   catdynres$Model$Results$predexplot[t] <- NA
                     #  }
                     catdynres$Model$Results$F.obsexplot[t]  <- ifelse(obsexplot.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),obsexplot.n[t],extendInt="yes")$root)
                     catdynres$Model$Results$F.predexplot[t] <- ifelse(predexplot.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),predexplot.n[t],extendInt="yes")$root)
                    }
                  names(catdynres$Model$Results)    <- c(paste("Period",properties$Units[1],sep="."),
                                                         paste("Effort",properties$Fleets[1,2],sep="."),
                                                         paste("Observed.Catch",properties$Units[4],sep="."),
                                                         paste("Predicted.Catch",properties$Units[4],sep="."),
                                                         "Observed.Catch.kg",
                                                         "Predicted.Catch.kg",
                                                         "Deviance",
                                                         "Likelihood",
                                                         "Deviance.Residuals",
                                                         paste("Predicted.Abundance",properties$Units[4],sep="."),
                                                         "Predicted.Biomass.tonnes",
                                                         "Observed.Explotrate",
                                                         "Predicted.Explotrate",
                                                         paste("Observed.F.1/",properties$Units[1],sep=""),
                                                         paste("Predicted.F.1/",properties$Units[1],sep=""))
                  return(catdynres);
  }
.CatDynExp2F.Res <-
function(properties,nstep,obsmbm1,obsmbm2,ppp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
  {
                  sdistr.set <- c("poisson","apnormal","aplnormal");
                  fleet.name <- properties$Fleets$Fleet
                  if(properties$Units[2] == "ind")
                    {
                     biom <- NA
                    }
                  else
                    {
                     #Set biom in tonnes for both units of individual body weight
                     #Set obscat.w, predcat.w in kg for both units of individual body weight
                     biom <- mean(c(obsmbm1*(1e-6*(properties$Units[3] == "g")+
                                             1e-3*(properties$Units[3] == "kg")),
                                    obsmbm2*(1e-6*(properties$Units[3] == "g")+
                                             1e-3*(properties$Units[3] == "kg"))),na.rm=TRUE)*
                             nstep*(1e9*(properties$Units[4] == "bill")+
                                    1e6*(properties$Units[4] == "mill")+
                                    1e3*(properties$Units[4] == "thou")+
                                    1e2*(properties$Units[4] == "hund"))
                     obscat1.w <- obsmbm1*(1e-3*(properties$Units[3] == "g")+
                                           1e0*(properties$Units[3] == "kg"))*
                                  obscat1*(1e9*(properties$Units[4] == "bill")+
                                           1e6*(properties$Units[4] == "mill")+
                                           1e3*(properties$Units[4] == "thou")+
                                           1e2*(properties$Units[4] == "hund"))
                     predcat1.w <- obsmbm1*(1e-3*(properties$Units[3] == "g")+
                                            1e0*(properties$Units[3] == "kg"))*
                                   predcat1*(1e9*(properties$Units[4] == "bill")+
                                             1e6*(properties$Units[4] == "mill")+
                                             1e3*(properties$Units[4] == "thou")+
                                             1e2*(properties$Units[4] == "hund"))
                     obscat2.w <- obsmbm2*(1e-3*(properties$Units[3] == "g")+
                                           1e0*(properties$Units[3] == "kg"))*
                                  obscat2*(1e9*(properties$Units[4] == "bill")+
                                           1e6*(properties$Units[4] == "mill")+
                                           1e3*(properties$Units[4] == "thou")+
                                           1e2*(properties$Units[4] == "hund"))
                     predcat2.w <- obsmbm2*(1e-3*(properties$Units[3] == "g")+
                                            1e0*(properties$Units[3] == "kg"))*
                                   predcat2*(1e9*(properties$Units[4] == "bill")+
                                             1e6*(properties$Units[4] == "mill")+
                                             1e3*(properties$Units[4] == "thou")+
                                             1e2*(properties$Units[4] == "hund"))
                    }
                  catdynres                         <- vector("list",2);
                  names(catdynres)                  <- c("Properties","Model");
                  catdynres$Properties              <- properties;
                  catdynres$Model                   <- vector("list",5);
                  names(catdynres$Model)            <- c("Type","Dates","Distribution","Parameters","Results");
                  catdynres$Model$Type              <- ppp;
                  catdynres$Model$Dates             <- dates;
                  catdynres$Model$Distribution      <- distr;
                  catdynres$Model$Parameters        <- exp(par);
                  if(ppp == c(0,0))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,1))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,2))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,3))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,4))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,5))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[2],sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                   paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[2],sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                   paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[2],sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[2],sep=""),
                                                paste("psi.",fleet.name[1],sep=""),paste("psi.",fleet.name[2],sep=""));
                       }
                    }
                  names(catdynres$Model$Parameters) <- par.names;
                  catdynres$Model$Results           <- data.frame(period=ts.start:ts.end,
                                                                  obseff1=obseff1,
                                                                  obscat1.n=obscat1,
                                                                  predcat1.n=predcat1,
                                                                  obscat1.w=obscat1.w,
                                                                  predcat1.w=predcat1.w,
                                                                  deviance1=Likel[[2]][1],
                                                                  likcontr1=Likel[[3]][1],
                                                                  devresid1=Likel[[4]][1],
                                                                  obseff2=obseff2,
                                                                  obscat2.n=obscat2,
                                                                  predcat2.n=predcat2,
                                                                  obscat2.w=obscat2.w,
                                                                  predcat2.w=predcat2.w,
                                                                  deviance2=Likel[[2]][2],
                                                                  likcontr2=Likel[[3]][2],
                                                                  devresid2=Likel[[4]][2],
                                                                  npred=nstep,
                                                                  biompred=biom,
                                                                  obsexplot1=NA,
                                                                  predexplot1=NA,
                                                                  obsexplot2=NA,
                                                                  predexplot2=NA,
                                                                  obsexplot=NA,
                                                                  predexplot=NA,
                                                                  F.obsexplot1=NA,
                                                                  F.predexplot1=NA,
                                                                  F.obsexplot2=NA,
                                                                  F.predexplot2=NA,
                                                                  F.obsexplot=NA,
                                                                  F.predexplot=NA);
                  if(properties$Units[2] == "ind")
                    {
                     catdynres$Model$Results$obsexplot1  <- obscat1/nstep
                     catdynres$Model$Results$predexplot1 <- predcat1/nstep
                     catdynres$Model$Results$obsexplot2  <- obscat2/nstep
                     catdynres$Model$Results$predexplot2 <- predcat2/nstep
                     catdynres$Model$Results$obsexplot   <- (obscat1+obscat2)/nstep
                     catdynres$Model$Results$predexplot  <- (predcat1+predcat2)/nstep
                    }
                  if(properties$Units[2] != "ind")
                    {
                     catdynres$Model$Results$obsexplot1  <- 1e-3*obscat1/biom
                     catdynres$Model$Results$predexplot1 <- 1e-3*predcat1/biom
                     catdynres$Model$Results$obsexplot2  <- 1e-3*obscat2/biom
                     catdynres$Model$Results$predexplot2 <- 1e-3*predcat2/biom
                     catdynres$Model$Results$obsexplot   <- 1e-3*(obscat1+obscat2)/biom
                     catdynres$Model$Results$predexplot  <- 1e-3*(predcat1+predcat2)/biom
                     obsexplot1.n                        <- obscat1/nstep
                     predexplot1.n                       <- predcat1/nstep
                     obsexplot2.n                        <- obscat2/nstep
                     predexplot2.n                       <- predcat2/nstep
                     obsexplot.n                         <- (obscat1+obscat2)/nstep
                     predexplot.n                        <- (predcat1+predcat2)/nstep
                    }
                  for(t in 1:(ts.end-ts.start+1))
                    {
                     if(biom[t]==0)
                       {
                        catdynres$Model$Results$biompred[t]    <- NA
                        catdynres$Model$Results$obsexplot1[t]  <- NA
                        catdynres$Model$Results$predexplot1[t] <- NA
                        catdynres$Model$Results$obsexplot2[t]  <- NA
                        catdynres$Model$Results$predexplot2[t] <- NA
                        catdynres$Model$Results$obsexplot[t]   <- NA
                        catdynres$Model$Results$predexplot[t]  <- NA
                       }
                     catdynres$Model$Results$F.obsexplot1[t]  <- ifelse(obsexplot1.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),obsexplot1.n[t], extendInt = "yes")$root)
                     catdynres$Model$Results$F.predexplot1[t] <- ifelse(predexplot1.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),predexplot1.n[t], extendInt = "yes")$root)
                     catdynres$Model$Results$F.obsexplot2[t]  <- ifelse(obsexplot2.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),obsexplot2.n[t], extendInt = "yes")$root)
                     catdynres$Model$Results$F.predexplot2[t] <- ifelse(predexplot2.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),predexplot2.n[t], extendInt = "yes")$root)
                     catdynres$Model$Results$F.obsexplot[t]   <- ifelse(obsexplot.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),obsexplot.n[t], extendInt = "yes")$root)
                     catdynres$Model$Results$F.predexplot[t]  <- ifelse(predexplot.n[t]==0,0,uniroot(.F.explotrate,c(1e-9,1e1),exp(par[1]),predexplot.n[t], extendInt = "yes")$root)
                    }
                  names(catdynres$Model$Results) <- c(paste("Period.",properties$Units[1],sep=""),
                                                      paste("Effort",properties$Fleets[1,1],properties$Fleets[1,2],sep="."),
                                                      paste("Observed.Catch",properties$Fleets[1,1],properties$Units[4],sep="."),
                                                      paste("Predicted.Catch",properties$Fleets[1,1],properties$Units[4],sep="."),
                                                      paste("Observed.Catch",properties$Fleets[1,1],"kg",sep="."),
                                                      paste("Predicted.Catch",properties$Fleets[1,1],"kg",sep="."),
                                                      paste("Deviance",properties$Fleets[1,1],properties$Units[4],sep="."),
                                                      paste("Likelihood",properties$Fleets[1,1],properties$Units[4],sep="."),
                                                      paste("Deviance.Residuals",properties$Fleets[1,1],properties$Units[4],sep="."),
                                                      paste("Effort",properties$Fleets[2,1],properties$Fleets[2,2],sep="."),
                                                      paste("Observed.Catch",properties$Fleets[2,1],properties$Units[4],sep="."),
                                                      paste("Predicted.Catch",properties$Fleets[2,1],properties$Units[4],sep="."),
                                                      paste("Observed.Catch",properties$Fleets[2,1],"kg",sep="."),
                                                      paste("Predicted.Catch",properties$Fleets[2,1],"kg",sep="."),
                                                      paste("Deviance",properties$Fleets[2,1],properties$Units[4],sep="."),
                                                      paste("Likelihood",properties$Fleets[2,1],properties$Units[4],sep="."),
                                                      paste("Deviance.Residuals",properties$Fleets[2,1],properties$Units[4],sep="."),
                                                      paste("Predicted.Abundance",properties$Units[4],sep="."),
                                                      "Predicted.Biomass.tonnes",
                                                      paste("Observed.Explotrate",properties$Fleets[1,1],sep="."),
                                                      paste("Predicted.Explotrate",properties$Fleets[1,1],sep="."),
                                                      paste("Observed.Explotrate",properties$Fleets[2,1],sep="."),
                                                      paste("Predicted.Explotrate",properties$Fleets[2,1],sep="."),
                                                      "Observed.Explotrate.Total",
                                                      "Predicted.Explotrate.Total",
                                                      paste("Observed.F",properties$Fleets[1,1],"1/",properties$Units[1],sep="."),
                                                      paste("Predicted.F",properties$Fleets[1,1],"1/",properties$Units[1],sep="."),
                                                      paste("Observed.F",properties$Fleets[2,1],"1/",properties$Units[1],sep="."),
                                                      paste("Predicted.F",properties$Fleets[2,1],"1/",properties$Units[1],sep="."),
                                                      paste("Observed.F.Total","1/",properties$Units[1],sep="."),
                                                      paste("Predicted.F.Total","1/",properties$Units[1],sep="."))
                  return(catdynres);
  }
.lw <-
function(par,lenw,method)
                       {
                        len   <- lenw[,1]
                        wgt   <- lenw[,2]
                        res <- sum((par[1]*len^par[2]-wgt)^2)
                        return(res)
                       }
.F.explotrate <-
function(F.explot,M,explot)
  {
   (F.explot/(F.explot+M))*(1-exp(-F.explot-M)) - explot
  }
