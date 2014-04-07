CatDynFit <-
function(x, p, par, dates, distr, method, control=list(), hessian=TRUE, itnmax)
          {
          options(warn=-1)
          fleet.name <- x$Properties$Fleets$Fleet
          if(class(x) != "CatDynData") {stop("Pass an object of class CatDynData as first argument - see the help for as.CatDynData")}
          if(sum(sort(p) == p)<length(p)) {stop("Number of perturbations per fleet must not be arranged in descending order")}
          if(any(is.na(p))) {stop("NAs are not allowed in the perturbations per fleet vector 'p'")}
          if(length(p) < 1 || length(p) > 3) {stop("The integer vector 'p' determines the number of perturbations per fleet; its length must be 1 for single fleet, 2 for two fleets, or 3 for three fleets")}
          if(any(p < 0) || any(p > 20)) {stop("The number of perturbations per fleet shall not be less than 0 nor higher than 20")}
          if(class(p) != "numeric") {stop("'p' must be a numeric vector")}
          if(any(is.na(par))) {stop("NAs are not allowed in the parameter vector 'par'")}
          if(length(par) < 4 || length(par) > 48) {stop("For any of the model versions the number of parameters is > 4 and < 49")}
          if(class(par) != "numeric") {stop("'par' must be a numeric vector")}
          if(any(is.na(dates))){stop("NAs are not allowed in the dates vector")}
          if(length(dates) != (sum(p)+2)) {stop("The dates vector must contain initial time step, time steps of all perturbations, and final time step")}
          #if(class(dates) != "integer") {stop("'dates' must be a vector of integers")}
          if(dates[1] > dates[2]) {stop("Initial date must be less than next date")}
          if(tail(dates,1) < tail(dates,2)[1]) {stop("Final date must not be less than previous date")}
          if(distr != "normal" && distr != "lognormal" ) {stop("'distr' must be 'normal' for the additive model, or 'lognormal' for the multiplicative model")}
          if(length(fleet.name) == 1)
            {
             if(sum(sort(dates) == dates)<length(dates)){stop("Dates must be arranged in ascending order")} 
             if(length(p) != 1)
               {
                stop("For a 1-fleet model 'p' must be of length 1")
               }
             #pure depletion
             if(p == 0)
                {
                 if(length(dates) != 2) {stop("For a 1-fleet 0P (pure depletion) model 'dates' must be a vector of length = 2")}
                 if(length(par) != 5) {stop("For a 1-fleet 0P (pure depletion) model 'par' must be a vector of lengt = 5")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.end");
                    par.names                      <- c("M","N0",paste(c("k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev         <- NA;
                       results2[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                      }
                    else
                      {
                       results2[[i]]$converg          <- results1[i,length(par)+5];
                       results2[[i]]$kkt              <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC              <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par           <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads        <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #1P
             else if(p == 1)
                {
                 if(length(dates) != 3) {stop("For a 1-fleet 1P model 'dates' must be a vector of length = 3")}
                 if(length(par) != 6) {stop("For a 1-fleet 1P model 'par' must be a vector of lengt = 6")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #2P
             else if(p == 2)
                {
                 if(length(dates) != 4) {stop("For a 1-fleet 2P model 'dates' must be a vector of length = 4")}
                 if(length(par) != 7) {stop("For a 1-fleet 2P model 'par' must be a vector of lengt = 7")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7)),
                                                                            mean=as.numeric(results1[i,1:length(par)]),
                                                                            cov=try(solve(temp[i,]$nhatend)),
                                                                            ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));                                                                         
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }                                                   
                    }
                rm(temp);
                }
             #3P
             else if(p == 3)
                {
                 if(length(dates) != 5) {stop("For a 1-fleet 3P model 'dates' must be a vector of length = 5")}
                 if(length(par) != 8) {stop("For a 1-fleet 3P model 'par' must be a vector of lengt = 8")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #4P
             else if(p == 4)
                {
                 if(length(dates) != 6) {stop("For a 1-fleet 4P model 'dates' must be a vector of length = 6")}
                 if(length(par) != 9) {stop("For a 1-fleet 4P model 'par' must be a vector of lengt = 9")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #5P
             else if(p == 5)
                {
                 if(length(dates) != 7) {stop("For a 1-fleet 5P model 'dates' must be a vector of length = 7")}
                 if(length(par) != 10) {stop("For a 1-fleet 5P model 'par' must be a vector of lengt = 10")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev         <- rep(NA, length(par.names))
                       results2[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev        <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                     ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10)),
                                                                              mean=as.numeric(results1[i,1:length(par)]),
                                                                              cov=try(solve(temp[i,]$nhatend)),
                                                                              ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #6P
             else if(p == 6)
                {
                 if(length(dates) != 8) {stop("For a 1-fleet 6P model 'dates' must be a vector of length = 8")}
                 if(length(par) != 11) {stop("For a 1-fleet 6P model 'par' must be a vector of lengt = 11")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #7P
             else if(p == 7)
                {
                 if(length(dates) != 9) {stop("For a 1-fleet 7P model 'dates' must be a vector of length = 9")}
                 if(length(par) != 12) {stop("For a 1-fleet 7P model 'par' must be a vector of lengt = 12")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6","ts.P7","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #8P
             else if(p == 8)
                {
                 if(length(dates) != 10) {stop("For a 1-fleet 8P model 'dates' must be a vector of length = 10")}
                 if(length(par) != 13) {stop("For a 1-fleet 8P model 'par' must be a vector of lengt = 13")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #9P
             else if(p == 9)
                {
                 if(length(dates) != 11) {stop("For a 1-fleet 9P model 'dates' must be a vector of length = 11")}
                 if(length(par) != 14) {stop("For a 1-fleet 9P model 'par' must be a vector of lengt = 14")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp)
                }
             #10P
             else if(p == 10)
                {
                 if(length(dates) != 12) {stop("For a 1-fleet 10P model 'dates' must be a vector of length = 12")}
                 if(length(par) != 15) {stop("For a 1-fleet 10P model 'par' must be a vector of lengt = 15")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #11P
             else if(p == 11)
                {
                 if(length(dates) != 13) {stop("For a 1-fleet 11P model 'dates' must be a vector of length = 13")}
                 if(length(par) != 16) {stop("For a 1-fleet 11P model 'par' must be a vector of lengt = 16")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #12P
             else if(p == 12)
                {
                 if(length(dates) != 14) {stop("For a 1-fleet 12P model 'dates' must be a vector of length = 14")}
                 if(length(par) != 17) {stop("For a 1-fleet 12P model 'par' must be a vector of lengt = 17")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #13P
             else if(p == 13)
                {
                 if(length(dates) != 15) {stop("For a 1-fleet 13P model 'dates' must be a vector of length = 15")}
                 if(length(par) != 18) {stop("For a 1-fleet 13P model 'par' must be a vector of lengt = 18")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #14P
             else if(p == 14)
                {
                 if(length(dates) != 16) {stop("For a 1-fleet 14P model 'dates' must be a vector of length = 16")}
                 if(length(par) != 19) {stop("For a 1-fleet 14P model 'par' must be a vector of lengt = 19")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13",
                                                        "ts.P14","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18),~exp(x19)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #15P
             else if(p == 15)
                {
                 if(length(dates) != 17) {stop("For a 1-fleet 15P model 'dates' must be a vector of length = 17")}
                 if(length(par) != 20) {stop("For a 1-fleet 15P model 'par' must be a vector of lengt = 20")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13",
                                                        "ts.P14","ts.P15","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #16P
             else if(p == 16)
                {
                 if(length(dates) != 18) {stop("For a 1-fleet 16P model 'dates' must be a vector of length = 18")}
                 if(length(par) != 21) {stop("For a 1-fleet 16P model 'par' must be a vector of lengt = 21")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13",
                                                        "ts.P14","ts.P15","ts.P16","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                              ~exp(x21)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #17P
             else if(p == 17)
                {
                 if(length(dates) != 19) {stop("For a 1-fleet 17P model 'dates' must be a vector of length = 19")}
                 if(length(par) != 22) {stop("For a 1-fleet 17P model 'par' must be a vector of lengt = 22")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13",
                                                        "ts.P14","ts.P15","ts.P16","ts.P17","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                              ~exp(x21),~exp(x22)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #18P
             else if(p == 18)
                {
                 if(length(dates) != 20) {stop("For a 1-fleet 18P model 'dates' must be a vector of length = 20")}
                 if(length(par) != 23) {stop("For a 1-fleet 18P model 'par' must be a vector of lengt = 23")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,3], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13",
                                                        "ts.P14","ts.P15","ts.P16","ts.P17","ts.P18","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "P18.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                              ~exp(x21),~exp(x22),~exp(x23)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #19P
             else if(p == 19)
                {
                 if(length(dates) != 21) {stop("For a 1-fleet 19P model 'dates' must be a vector of length = 21")}
                 if(length(par) != 24) {stop("For a 1-fleet 19P model 'par' must be a vector of lengt = 24")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13",
                                                        "ts.P14","ts.P15","ts.P16","ts.P17","ts.P18","ts.P19","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "P18.","P19.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                              ~exp(x21),~exp(x22),~exp(x23),~exp(x24)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #20P
             else if(p == 20)
                {
                 if(length(dates) != 22) {stop("For a 1-fleet 20P model 'dates' must be a vector of length = 22")}
                 if(length(par) != 25) {stop("For a 1-fleet 20P model 'par' must be a vector of lengt = 25")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name]][,2], 
                                                         obscat1=x$Data[[fleet.name]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian, 
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1","ts.P2","ts.P3","ts.P4","ts.P5","ts.P6",
                                                        "ts.P7","ts.P8","ts.P9","ts.P10","ts.P11","ts.P12","ts.P13",
                                                        "ts.P14","ts.P15","ts.P16","ts.P17","ts.P18","ts.P19","ts.P20","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "P18.","P19.","P20.","k.","alpha.","beta."),sort(rep(fleet.name,3+p)),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                              ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                              ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                              ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                              ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25)),
                                                                       mean=as.numeric(results1[i,1:length(par)]),
                                                                       cov=try(solve(temp[i,]$nhatend)),
                                                                       ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             else 
                {
                 stop("Parameter p -number of perturbations- does not correspond with any of the allowed integers for one fleet, 0<=p<=20")
                }
            } #Fleet 1
          else if(length(fleet.name) == 2)
            {
             if(sum(sort(dates[1:(p[1]+1)]) == dates[1:(p[1]+1)]) < length(dates[1:(p[1]+1)]) || sum(sort(dates[(p[1]+2):(p[1]+p[2]+2)]) == dates[(p[1]+2):(p[1]+p[2]+2)]) < length(dates[(p[1]+2):(p[1]+p[2]+2)]))
               {stop("Perturbation dates for each fleet should be arranged in ascending order")}
#             if(length(unique(dates[1:(p[1]+1)])) != length(dates[1:(p[1]+1)]) || length(unique(dates[(p[1]+2):(p[1]+p[2]+2)])) != length(dates[(p[1]+2):(p[1]+p[2]+2)]))   
#               {stop("Perturbation dates inside each fleet should be all distinct and distinct from initial and final date")}
             if(length(p) != 2) {stop("For a 2-fleet model 'p' must be of length 2")}
             #0P0P
             if(p[1] == 0 & p[2] == 0)
                {
                 if(length(par) != 8) {stop("For a 2-Fleet 0P 0P (pure depletion) model 'par' must be a vector of length = 8")}
                 if(length(dates) != 2) {stop("For a 2-fleet 0P 0P (pure depletion) model 'dates' must be a vector of length = 2")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.end");
                    par.names                      <- c("M","N0",paste(rep(c("k.","alpha.","beta."),2),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),
                                                                                ~exp(x7),~exp(x8)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #0P1P
             else if(p[1] == 0 & p[2] == 1)  
                {
                 if(length(par) != 9) {stop("For a 2-Fleet 0P 1P (pure depletion) model 'par' must be a vector of length = 9")}
                 if(length(dates) != 3) {stop("For a 2-fleet 0P 1P (pure depletion) model 'dates' must be a vector of length = 3")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F2","ts.end");
                    par.names                      <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #0P2P 
             else if(p[1] == 0 & p[2] == 2)
                {
                 if(length(par) != 10) {stop("For a 2-Fleet 0P 2P model 'par' must be a vector of length = 10")}
                 if(length(dates) != 4) {stop("For a 2-Fleet 0P 2P model 'dates' must be a vector of length = 4")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F2","ts.P2F2","ts.end");
                    par.names                      <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #0P3P
             else if(p[1] == 0 & p[2] == 3)
                {
                 if(length(par) != 11) {stop("For a 2-Fleet 0P 3P model M 'par' must be a vector of length = 11")}
                 if(length(dates) != 5) {stop("For a 2-Fleet 0P 3P model 'dates' must be a vector of length = 5")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F2","ts.P2F2","ts.P3F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #0P4P
             else if (p[1] == 0 & p[2] == 4)
                {
                 if(length(par) != 12) {stop("For a 2-Fleet 0P 4P model 'par' must be a vector of length = 12")}
                 if(length(dates) != 6) {stop("For a 2-Fleet 0P 4P model 'dates' must be a vector of length = 6")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #0P5P
             else if (p[1] == 0 & p[2] == 5)
                {
                 if(length(par) != 13) {stop("For a 2-Fleet 0P 5P model 'par' must be a vector of length = 13")}
                 if(length(dates) != 7) {stop("For a 2-Fleet 0P 5P model 'dates' must be a vector of length = 7")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
              #1P1P
              else if (p[1] == 1 & p[2] == 1)
                {
                 if(length(par) != 10) {stop("For a 2-Fleet 1P 1P model 'par' must be a vector of length = 10")}
                 if(length(dates) != 4) {stop("For a 2-Fleet 1P 1P model 'dates' must be a vector of length = 4")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P1F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","k.","alpha.","beta.","P1.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #1P2P
             else if (p[1] == 1 & p[2] == 2)
                {
                 if(length(par) != 11) {stop("For a 2-Fleet 1P 2P model 'par' must be a vector of length = 11")}
                 if(length(dates) != 5) {stop("For a 2-Fleet 1P 2P model 'dates' must be a vector of length = 5")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P1F2","ts.P2F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #1P3P
             else if (p[1] == 1 & p[2] == 3)
                {
                 if(length(par) != 12) {stop("For a 2-Fleet 1P 3P model 'par' must be a vector of length = 12")}
                 if(length(dates) != 6) {stop("For a 2-Fleet 1P 3P model 'dates' must be a vector of length = 6")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P1F2","ts.P2F2","ts.P3F2","ts.end");
                    par.names                      <- c("M","N0",paste(c("P1.","k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                  }
                 rm(temp);
                }
             #1P4P
             else if (p[1] == 1 & p[2] == 4)
                {
                 if(length(par) != 13) {stop("For a 2-Fleet 1P 4P model 'par' must be a vector of length = 13")}
                 if(length(dates) != 7) {stop("For a 2-Fleet 1P 4P model 'dates' must be a vector of length = 7")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                  }
                 rm(temp);
                }
            #2P2P
             else if (p[1] == 2 & p[2] == 2)
                {
                 if(length(par) != 12) {stop("For a 2-Fleet 2P 2P model 'par' must be a vector of length = 12")}
                 if(length(dates) != 6) {stop("For a 2-Fleet 2P 2P model 'dates' must be a vector of length = 6")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P1F2","ts.P2F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
            #2P3P
             else if (p[1] == 2 & p[2] == 3)
                {
                 if(length(par) != 13) {stop("For a 2-Fleet 2P 3P model 'par' must be a vector of length = 13")}
                 if(length(dates) != 7) {stop("For a 2-Fleet 2P 3P model 'dates' must be a vector of length = 7")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P1F2","ts.P2F2","ts.P3F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                        c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #2P4P
             else if (p[1] == 2 & p[2] == 4)
                {
                 if(length(par) != 14) {stop("For a 2-Fleet 2P 4P model 'par' must be a vector of length = 14")}
                 if(length(dates) != 8) {stop("For a 2-Fleet 2P 4P model 'dates' must be a vector of length = 8")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.",
                                                        "beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #2P5P
             else if (p[1] == 2 & p[2] == 5)
                {
                 if(length(par) != 15) {stop("For a 2-Fleet 2P 5P model 'par' must be a vector of length = 15")}
                 if(length(dates) != 9) {stop("For a 2-Fleet 2P 5P model 'dates' must be a vector of length = 9")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P1F2","ts.P2F2","ts.P3F2",
                                                        "ts.P4F2","ts.P5F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","k.","alpha.","beta.","P1.","P2.",
                                                                         "P3.","P4.","P5.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                 }
             #3P3P
             else if (p[1] == 3 & p[2] == 3)
                {
                 if(length(par) != 14) {stop("For a 2-Fleet 3P 3P model 'par' must be a vector of length = 14")}
                 if(length(dates) != 8) {stop("For a 2-Fleet 3P 3P model 'dates' must be a vector of length = 8")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P1F2","ts.P2F2","ts.P3F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #3P4P
             else if (p[1] == 3 & p[2] == 4)
                {
                 if(length(par) != 15) {stop("For a 2-Fleet 3P 4P model 'par' must be a vector of length = 15")}
                 if(length(dates) != 9) {stop("For a 2-Fleet 3P 4P model 'dates' must be a vector of length = 9")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #3P5P
             else if (p[1] == 3 & p[2] == 5)
                {
                 if(length(par) != 16) {stop("For a 2-Fleet 3P 5P model 'par' must be a vector of length = 16")}
                 if(length(dates) != 10) {stop("For a 2-Fleet 3P 5P model 'dates' must be a vector of length = 10")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P1F2","ts.P2F2",
                                                        "ts.P3F2","ts.P4F2","ts.P5F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","k.","alpha.","beta.","P1.","P2.",
                                                        "P3.","P4.","P5.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #4P4P
             else if (p[1] == 4 & p[2] == 4)
                {
                 if(length(par) != 16) {stop("For a 2-Fleet 4P 4P model 'par' must be a vector of length = 16")}
                 if(length(dates) != 10) {stop("For a 2-Fleet 4P 4P model 'dates' must be a vector of length = 10")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P1F2",
                                                        "ts.P2F2","ts.P3F2","ts.P4F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #4P5P
             else if (p[1] == 4 & p[2] == 5)                                                              
                {
                 if(length(par) != 17) {stop("For a 2-Fleet 4P 5P model 'par' must be a vector of length = 17")}
                 if(length(dates) != 11) {stop("For a 2-Fleet 4P 5P model 'dates' must be a vector of length = 11")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P1F2",
                                                        "ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),
                                                                                      ~exp(x17)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #5P5P
             else if (p[1] == 5 & p[2] == 5)                                                              
                {
                 if(length(par) != 18) {stop("For a 2-Fleet 5P 5P model 'par' must be a vector of length = 18")}
                 if(length(dates) != 12) {stop("For a 2-Fleet 5P 5P model 'dates' must be a vector of length = 12")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #6P6P
             else if (p[1] == 6 & p[2] == 6)                                                              
                {
                 if(length(par) != 20) {stop("For a 2-Fleet 6P 6P model 'par' must be a vector of length = 20")}
                 if(length(dates) != 14) {stop("For a 2-Fleet 6P 6P model 'dates' must be a vector of length = 14")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1","ts.P6F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #7P7P
             else if (p[1] == 7 & p[2] == 7)                                                              
                {
                 if(length(par) != 22) {stop("For a 2-Fleet 7P 7P model 'par' must be a vector of length = 22")}
                 if(length(dates) != 16) {stop("For a 2-Fleet 7P 7P model 'dates' must be a vector of length = 16")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #8P8P
             else if (p[1] == 8 & p[2] == 8)                                                              
                {
                 if(length(par) != 24) {stop("For a 2-Fleet 8P 8P model 'par' must be a vector of length = 24")}
                 if(length(dates) != 18) {stop("For a 2-Fleet 8P 8P model 'dates' must be a vector of length = 18")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #9P9P
             else if (p[1] == 9 & p[2] == 9)                                                              
                {
                 if(length(par) != 26) {stop("For a 2-Fleet 9P 9P model 'par' must be a vector of length = 26")}
                 if(length(dates) != 20) {stop("For a 2-Fleet 9P 9P model 'dates' must be a vector of length = 20")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #10P10P
             else if (p[1] == 10 & p[2] == 10)                                                              
                {
                 if(length(par) != 28) {stop("For a 2-Fleet 10P 10P model 'par' must be a vector of length = 28")}
                 if(length(dates) != 22) {stop("For a 2-Fleet 10P 10P model 'dates' must be a vector of length = 22")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #11P11P
             else if (p[1] == 11 & p[2] == 11)                                                              
                {
                 if(length(par) != 30) {stop("For a 2-Fleet 11P 11P model 'par' must be a vector of length = 30")}
                 if(length(dates) != 24) {stop("For a 2-Fleet 11P 11P model 'dates' must be a vector of length = 24")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #12P12P
             else if (p[1] == 12 & p[2] == 12)                                                              
                {
                 if(length(par) != 32) {stop("For a 2-Fleet 12P 12P model 'par' must be a vector of length = 32")}
                 if(length(dates) != 26) {stop("For a 2-Fleet 12P 12P model 'dates' must be a vector of length = 26")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #13P13P
             else if (p[1] == 13 & p[2] == 13)
                {
                 if(length(par) != 34) {stop("For a 2-Fleet 13P 13P model 'par' must be a vector of length = 34")}
                 if(length(dates) != 28) {stop("For a 2-Fleet 13P 13P model 'dates' must be a vector of length = 28")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             else if (p[1] == 14 & p[2] == 14)
                {
                 if(length(par) != 36) {stop("For a 2-Fleet 14P 14P model 'par' must be a vector of length = 36")}
                 if(length(dates) != 30) {stop("For a 2-Fleet 14P 14P model 'dates' must be a vector of length = 30")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1","ts.P14F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.P14F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.","P14.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                      ~exp(x36)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                ~exp(x36)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #15P15P
             else if (p[1] == 15 & p[2] == 15)
                {
                 if(length(par) != 38) {stop("For a 2-Fleet 15P 15P model 'par' must be a vector of length = 38")}
                 if(length(dates) != 32) {stop("For a 2-Fleet 15P 15P model 'dates' must be a vector of length = 32")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1","ts.P14F1","ts.P15F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.P14F2","ts.P15F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.","P14.","P15.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                      ~exp(x36),~exp(x37),~exp(x38)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                ~exp(x36),~exp(x37),~exp(x38)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #16P16P
             else if (p[1] == 16 & p[2] == 16)
                {
                 if(length(par) != 40) {stop("For a 2-Fleet 16P 16P model 'par' must be a vector of length = 40")}
                 if(length(dates) != 34) {stop("For a 2-Fleet 16P 16P model 'dates' must be a vector of length = 34")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1","ts.P14F1","ts.P15F1","ts.P16F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.P14F2","ts.P15F2","ts.P16F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.","P14.","P15.","P16.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                      ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #17P17P
             else if (p[1] == 17 & p[2] == 17)
                {
                 if(length(par) != 42) {stop("For a 2-Fleet 17P 17P model 'par' must be a vector of length = 42")}
                 if(length(dates) != 36) {stop("For a 2-Fleet 17P 17P model 'dates' must be a vector of length = 36")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1","ts.P14F1","ts.P15F1","ts.P16F1",
                                                        "ts.P17F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.P14F2","ts.P15F2","ts.P16F2","ts.P17F2",
                                                        "ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                      ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                      ~exp(x41),~exp(x42)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                ~exp(x41),~exp(x42)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
             #18P18P
             else if (p[1] == 18 & p[2] == 18)
                {
                 if(length(par) != 44) {stop("For a 2-Fleet 18P 18P model 'par' must be a vector of length = 44")}
                 if(length(dates) != 38) {stop("For a 2-Fleet 18P 18P model 'dates' must be a vector of length = 38")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1","ts.P14F1","ts.P15F1","ts.P16F1",
                                                        "ts.P17F1","ts.P18F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.P14F2","ts.P15F2","ts.P16F2","ts.P17F2",
                                                        "ts.P18F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "P18.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.","P14.","P15.","P16.","P17.","P18.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                      ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                      ~exp(x41),~exp(x42),~exp(x43),~exp(x44)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                ~exp(x41),~exp(x42),~exp(x43),~exp(x44)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
            #19P19P
            else if (p[1] == 19 & p[2] == 19)
                {
                 if(length(par) != 46) {stop("For a 2-Fleet 19P 19P model 'par' must be a vector of length = 46")}
                 if(length(dates) != 40) {stop("For a 2-Fleet 19P 19P model 'dates' must be a vector of length = 40")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1","ts.P14F1","ts.P15F1","ts.P16F1",
                                                        "ts.P17F1","ts.P18F1","ts.P19F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.P14F2","ts.P15F2","ts.P16F2","ts.P17F2",
                                                        "ts.P18F2","ts.P19F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "P18.","P19.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.","P14.","P15.","P16.","P17.","P18.","P19.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                      ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                      ~exp(x41),~exp(x42),~exp(x43),~exp(x44),~exp(x45),
                                                                                      ~exp(x46)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                ~exp(x41),~exp(x42),~exp(x43),~exp(x44),~exp(x45),
                                                                                ~exp(x46)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
            #20P20P
            else if (p[1] == 20 & p[2] == 20)
                {
                 if(length(par) != 48) {stop("For a 2-Fleet 20P 20P model 'par' must be a vector of length = 48")}
                 if(length(dates) != 42) {stop("For a 2-Fleet 20P 20P model 'dates' must be a vector of length = 42")}
                 results1        <- do.call(optimx, list(par=par, 
                                                         fn=eval(as.name(paste('.CDMN',as.character(p[1]),'P',as.character(p[2]),'P.Lik',sep=""))), 
                                                         gr=NULL, 
                                                         dates=dates, 
                                                         obseff1=x$Data[[fleet.name[1]]][,2], 
                                                         obscat1=x$Data[[fleet.name[1]]][,5], 
                                                         obseff2=x$Data[[fleet.name[2]]][,2], 
                                                         obscat2=x$Data[[fleet.name[2]]][,5], 
                                                         distr=distr, 
                                                         method=method, 
                                                         lower = -Inf, 
                                                         upper = Inf, 
                                                         control=list(), 
                                                         hessian=hessian,
                                                         itnmax=itnmax));
                 results2        <- vector("list",length(method));
                 names(results2) <- method;
                 temp            <- attr(results1, "details")
                 for(i in 1:length(method))
                   {
                    results2[[i]]$Type       <- p;
                    results2[[i]]$Dates      <- dates;
                    results2[[i]]$Distr      <- distr;
                    names(results2[[i]]$Dates)     <- c("ts.start","ts.P1F1","ts.P2F1","ts.P3F1","ts.P4F1","ts.P5F1",
                                                        "ts.P6F1","ts.P7F1","ts.P8F1","ts.P9F1","ts.P10F1","ts.P11F1",
                                                        "ts.P12F1","ts.P13F1","ts.P14F1","ts.P15F1","ts.P16F1",
                                                        "ts.P17F1","ts.P18F1","ts.P19F1","ts.P20F1",
                                                        "ts.P1F2","ts.P2F2","ts.P3F2","ts.P4F2","ts.P5F2","ts.P6F2",
                                                        "ts.P7F2","ts.P8F2","ts.P9F2","ts.P10F2","ts.P11F2","ts.P12F2",
                                                        "ts.P13F2","ts.P14F2","ts.P15F2","ts.P16F2","ts.P17F2",
                                                        "ts.P18F2","ts.P19F2","ts.P20F2","ts.end"); 
                    par.names                      <- c("M","N0",paste(c("P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.",
                                                        "P9.","P10.","P11.","P12.","P13.","P14.","P15.","P16.","P17.",
                                                        "P18.","P19.","P20.",
                                                        "k.","alpha.","beta.",
                                                        "P1.","P2.","P3.","P4.","P5.","P6.","P7.","P8.","P9.","P10.",
                                                        "P11.","P12.","P13.","P14.","P15.","P16.","P17.","P18.","P19.","P20.",
                                                        "k.","alpha.","beta."),c(rep(fleet.name[1],3+p[1]),rep(fleet.name[2],3+p[2])),sep=""));
                    if(length(temp[i,]$ngatend) < length(par.names) || any(is.na(temp[i,]$nhatend)) || 1/kappa(temp[i,]$nhatend) < 1e-15)
                      {
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev <- rep(NA, length(par.names))
                       results2[[i]]$Cor      <- matrix(NA,length(par.names),length(par.names))
                      }
                    else
                      {
                       results2[[i]]$converg    <- results1[i,length(par)+5];
                       results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                       results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                       results2[[i]]$bt.par     <- exp(results1[i,1:length(par)])
                       results2[[i]]$num.grads  <- temp[i,]$ngatend;
                       results2[[i]]$bt.stdev         <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                      ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                      ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                      ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                      ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                      ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                      ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                      ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                      ~exp(x41),~exp(x42),~exp(x43),~exp(x44),~exp(x45),
                                                                                      ~exp(x46),~exp(x47),~exp(x48)),
                                                                               mean=as.numeric(results1[i,1:length(par)]),
                                                                               cov=try(solve(temp[i,]$nhatend)),
                                                                               ses=FALSE)));
                       results2[[i]]$Cor              <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),
                                                                                ~exp(x6),~exp(x7),~exp(x8),~exp(x9),~exp(x10),
                                                                                ~exp(x11),~exp(x12),~exp(x13),~exp(x14),~exp(x15),
                                                                                ~exp(x16),~exp(x17),~exp(x18),~exp(x19),~exp(x20),
                                                                                ~exp(x21),~exp(x22),~exp(x23),~exp(x24),~exp(x25),
                                                                                ~exp(x26),~exp(x27),~exp(x28),~exp(x29),~exp(x30),
                                                                                ~exp(x31),~exp(x32),~exp(x33),~exp(x34),~exp(x35),
                                                                                ~exp(x36),~exp(x37),~exp(x38),~exp(x39),~exp(x40),
                                                                                ~exp(x41),~exp(x42),~exp(x43),~exp(x44),~exp(x45),
                                                                                ~exp(x46),~exp(x47),~exp(x48)),
                                                                         mean=as.numeric(results1[i,1:length(par)]),
                                                                         cov=try(solve(temp[i,]$nhatend)),
                                                                         ses=FALSE));
                       names(results2[[i]]$num.grads) <- par.names;
                       names(results2[[i]]$bt.par)    <- par.names;
                       names(results2[[i]]$bt.stdev)  <- par.names;
                       colnames(results2[[i]]$Cor)    <- par.names;
                       rownames(results2[[i]]$Cor)    <- par.names;
                      }
                   }
                 rm(temp);
                }
            #stop
            else 
                {
                 stop("The perturbations parameter vector p does not correspond with any of the allowed models for 2 fleets")
                }
             } #Fleet 2
            results                     <- vector("list",3);
            names(results)              <- c("Data","Initial","Model")
            results$Data                <- x;
            results$Initial             <- par;
            names(results$Initial)      <- par.names;
            results$Model               <- results2;
            class(results)              <- "catdyn";
            return(results);
           }
