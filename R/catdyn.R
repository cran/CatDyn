catdyn <-
function(p, par, dates, obseff, obscat, M.fixed, M, distr, method, control=list(), hessian=TRUE, itnmax)
          {
          require(optimx);
          require(msm);
               if(M.fixed == TRUE)
                  {
                   if(p == 0)
                      {
                       if(length(par) != 4 | length(dates) != 2)
                          {
                           stop('For depletion model with fixed M par must be a vector of length = 4 and dates must be a vector of length = 2')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 1)
                      {
                       if(length(par) != 5 | length(dates) != 3)
                          {
                           stop('For 1-perturbation model with fixed M par must be a vector of length = 5 and dates must be a vector of length = 3')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 2)
                      {
                       if(length(par) != 6 | length(dates) != 4)
                          {
                           stop('For 2-perturbations model with fixed M par must be a vector of length = 6 and dates must be a vector of length = 4')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 3)
                      {
                       if(length(par) != 7 | length(dates) != 5)
                          {
                           stop('For 3-perturbations model with fixed M par must be a vector of length = 7 and dates must be a vector of length = 5')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 4)
                      {
                       if(length(par) != 8 | length(dates) != 6)
                          {
                           stop('For 4-perturbations model with fixed M par must be a vector of length = 8 and dates must be a vector of length = 6')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                       }
                   else if(p>4 | p<0)
                      {
                      stop('p is the number of perturbations, either 0, 1, 2, 3, or 4')
                      }
                   results2 <- vector("list",length(method));
                       for(i in 1:length(method))
                           {
                           results2[[i]]$Model     <- paste('CDMN',as.character(p),'P',sep="");
                           results2[[i]]$Distr     <- distr;
                           results2[[i]]$method    <- attr(results1, "details")[[i]]$method;
                           results2[[i]]$M         <- M;
                           results2[[i]]$converg   <- attr(results1, "details")[[i]]$convergence;
                           results2[[i]]$kkt       <- c(attr(results1, "details")[[i]]$kkt1,attr(results1, "details")[[i]]$kkt2);
                           results2[[i]]$AIC       <- 2*length(attr(results1, "details")[[i]]$par)-2*(-attr(results1, "details")[[i]]$value);
                           results2[[i]]$bt.par    <- exp(attr(results1, "details")[[i]]$par);
                           results2[[i]]$num.grads <- attr(results1, "details")[[i]]$ngatend;
                           if(p == 0)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN0P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 1)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN1P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 2)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN2P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 3)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN3P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 4)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN4P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           }

                   }
               else
                  {
                   if(p == 0)
                      {
                       if(length(par) != 5 | length(dates) != 2)
                          {
                           stop('For depletion model with fixed M par must be a vector of length = 5 and dates must be a vector of length = 2')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 1)
                      {
                       if(length(par) != 6 | length(dates) != 3)
                          {
                           stop('For 1-perturbation model with free M par must be a vector of length = 6 and dates must be a vector of length = 3')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 2)
                      {
                       if(length(par) != 7 | length(dates) != 4)
                          {
                           stop('For 2-perturbations model with free M par must be a vector of length = 7 and dates must be a vector of length = 4')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 3)
                      {
                       if(length(par) != 8 | length(dates) != 5)
                          {
                           stop('For 3-perturbations model with free M par must be a vector of length = 8 and dates must be a vector of length = 5')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                      }
                   else if(p == 4)
                      {
                       if(length(par) != 9 | length(dates) != 6)
                          {
                           stop('For 4-perturbations model with free M par must be a vector of length = 9 and dates must be a vector of length = 6')
                          }
                       else
                          {
                           results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                          }
                       }
                   else if(p>4 | p<0)
                      {
                      stop('p is the number of perturbations, either 0, 1, 2, 3, or 4')
                      }
                   results2 <- vector("list",length(method));
                       for(i in 1:length(method))
                           {
                           results2[[i]]$Model     <- paste('CDMN',as.character(p),'P',sep="");
                           results2[[i]]$Distr     <- distr;
                           results2[[i]]$method    <- attr(results1, "details")[[i]]$method;
                           results2[[i]]$M         <- "Free (bt.par[1])"
                           results2[[i]]$converg   <- attr(results1, "details")[[i]]$convergence;
                           results2[[i]]$kkt       <- c(attr(results1, "details")[[i]]$kkt1,attr(results1, "details")[[i]]$kkt2);
                           results2[[i]]$AIC       <- 2*length(attr(results1, "details")[[i]]$par)-2*(-attr(results1, "details")[[i]]$value);
                           results2[[i]]$bt.par    <- exp(attr(results1, "details")[[i]]$par);
                           results2[[i]]$num.grads <- attr(results1, "details")[[i]]$ngatend;
                           if(p == 0)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN0P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 1)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN1P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 2)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN2P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 3)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN3P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           else if(p == 4)
                              {
                              results2[[i]]$ana.grads <- do.call(CDMN4P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, M=M, distr=distr));
                              results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                           mean=attr(results1, "details")[[i]]$par,
                                                           cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                           ses=FALSE)));
                              results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                          mean=attr(results1, "details")[[i]]$par,
                                                          cov=solve(attr(results1, "details")[[i]]$nhatend),
                                                          ses=FALSE));
                              }
                           }
                  }
               return(results2);
          }

