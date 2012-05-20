catdyn <-
function(p, par, dates, obseff, obscat, M.fixed, M, distr, num.wrapper, 
                   method, control=list(), hessian, itnmax)
          {
          require(optimx);
          require(msm);
               if(p == 0)
                  {
                   if(M.fixed == FALSE)
                      {
                       if(length(par) != 5)
                          {
                           stop('For depletion model with free M par must be a vector of length = 5')
                          }
                      }
                   else
                      {
                       if(length(par) != 4)
                          {
                           stop('For depletion model with fixed M par must be a vector of length = 4')
                          }
                      }
                   if(length(dates) != 2)
                      {
                       stop('For depletion model dates must be a vector of length = 2')
                      }
                  }
               else if(p == 1)
                  {
                   if(M.fixed == FALSE)
                      {
                       if(length(par) != 6)
                          {
                           stop('For 1-Perturbation model with free M par must be a vector of length = 6')
                          }
                      }
                   else
                      {
                       if(length(par) != 5)
                          {
                           stop('For 1-Perturbation model with fixed M par must be a vector of length = 5')
                          }
                      }
                   if(length(dates) != 3)
                      {
                       stop('For 1-Perturbation model dates must be a vector of length = 3')
                      }
                   }
               else if(p == 2)
                  {
                   if(M.fixed == FALSE)
                      {
                       if(length(par) != 7)
                          {
                           stop('For 2-Perturbations model with free M par must be a vector of length = 7')
                          }
                      }
                   else
                      {
                       if(length(par) != 6)
                          {
                           stop('For 2-Perturbations model with fixed M par must be a vector of length = 6')
                          }
                      }
                   if(length(dates) != 4)
                      {
                       stop('For 2-Perturbations model dates must be a vector of length = 4')
                      }
                   }
               else if(p == 3)
                  {
                   if(M.fixed == FALSE)
                      {
                       if(length(par) != 8)
                          {
                           stop('For 3-Perturbations model with free M par must be a vector of length = 8')
                          }
                      }
                   else
                      {
                       if(length(par) != 7)
                          {
                           stop('For 3-Perturbations model with fixed M par must be a vector of length = 7')
                          }
                      }
                   if(length(dates) != 5)
                      {
                       stop('For 3-Perturbations model dates must be a vector of length = 5')
                      }
                   }
               else if (p == 4)
                  {
                   if(M.fixed == FALSE)
                      {
                       if(length(par) != 9)
                          {
                           stop('For 4-Perturbations model with free M par must be a vector of length = 9')
                          }
                      }
                   else
                      {
                       if(length(par) != 8)
                          {
                           stop('For 4-Perturbations model with fixed M par must be a vector of length = 8')
                          }
                      }
                   if(length(dates) != 6)
                      {
                       stop('For 3-Perturbations model dates must be a vector of length = 6')
                      }
                   }
               if(length(obseff) != length(obscat))
                  {
                   stop('The number of time steps with observations of catch and nominal effort must be equal')
                  }
               else if(any(obseff < 0))
                  {
                   stop('Effort observations must be non-negative')
                  }
               else if(any(obscat < 0))
                  {
                   stop('Catch observations must be non-negative')
                  }
               else
                  {
                   if(num.wrapper == "optimx")
                      {
                      results1 <- do.call(optimx, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=NULL, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian, itnmax=itnmax));
                      results2 <- vector("list",length(method));
                      for(i in 1:length(method))
                          {
                          results2[[i]]$Model     <- paste('CDMN',as.character(p),'P',sep="");
                          results2[[i]]$Distr     <- distr;
                          results2[[i]]$method    <- attr(results1, "details")[[i]]$method;
                          results2[[i]]$converg   <- attr(results1, "details")[[i]]$convergence;
                          results2[[i]]$kkt       <- c(attr(results1, "details")[[i]]$kkt1,attr(results1, "details")[[i]]$kkt2);
                          results2[[i]]$AIC       <- 2*length(attr(results1, "details")[[i]]$par)-2*(-attr(results1, "details")[[i]]$value);
                          results2[[i]]$bt.par    <- exp(attr(results1, "details")[[i]]$par);
                          results2[[i]]$num.grads <- attr(results1, "details")[[i]]$ngatend;
                          if(p == 0)
                             {
                             results2[[i]]$ana.grads <- do.call(CDMN0P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                             results2[[i]]$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=attr(results1, "details")[[i]]$par,
                                                                              cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                              ses=FALSE)));
                             results2[[i]]$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                        mean=attr(results1, "details")[[i]]$par,
                                                                        cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                        ses=FALSE));
                             }
                          else if(p == 1)
                             {
                             results2[[i]]$ana.grads <- do.call(CDMN1P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                             results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                             mean=attr(results1, "details")[[i]]$par,
                                                                             cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                             ses=FALSE)));
                             results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                       mean=attr(results1, "details")[[i]]$par,
                                                                       cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                       ses=FALSE));
                             }
                          else if(p == 2)
                             {
                             results2[[i]]$ana.grads <- do.call(CDMN2P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                             results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                                             mean=attr(results1, "details")[[i]]$par,
                                                                             cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                             ses=FALSE)));
                             results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                                       mean=attr(results1, "details")[[i]]$par,
                                                                       cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                       ses=FALSE));
                             }
                          else if(p == 3)
                             {
                             results2[[i]]$ana.grads <- do.call(CDMN3P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                             results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                                             mean=attr(results1, "details")[[i]]$par,
                                                                             cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                             ses=FALSE)));
                             results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                                       mean=attr(results1, "details")[[i]]$par,
                                                                       cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                       ses=FALSE));
                             }
                          else if(p == 4)
                             {
                             results2[[i]]$ana.grads <- do.call(CDMN4P.Lik.gr, list(par=attr(results1, "details")[[i]]$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                             results2[[i]]$bt.stdev <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                                             mean=attr(results1, "details")[[i]]$par,
                                                                             cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                             ses=FALSE)));
                             results2[[i]]$Cor      <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                                       mean=attr(results1, "details")[[i]]$par,
                                                                       cov=try(solve(attr(results1, "details")[[i]]$nhatend)),
                                                                       ses=FALSE));
                             }
                          }
                      }
                   else if(num.wrapper == "optim")
                      {
                      results1 <- do.call(optim, list(par=par, fn=eval(as.name(paste('CDMN',as.character(p),'P.Lik',sep=""))), gr=eval(as.name(paste('CDMN',as.character(p),'P.Lik.gr',sep=""))), dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr, method=method, lower = -Inf, upper = Inf, control=list(), hessian=hessian));                      
                      results2 <- vector("list",1);
                      results2$Model     <- paste('CDMN',as.character(p),'P',sep="");
                      results2$Distr     <- distr;
                      results2$method    <- method;
                      results2$converg   <- results1$convergence;
                      results2$AIC       <- 2*length(results1$par)-2*(-(results1$value));
                      results2$bt.par    <- exp(results1$par);
                      if(p == 0)
                         {
                         results2$ana.grads <- do.call(CDMN0P.Lik.gr, list(par=results1$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                         results2$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                     mean=results1$par,
                                                                     cov=try(solve(results1$hessian)),
                                                                     ses=FALSE)));
                         results2$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                               mean=results1$par,
                                                               cov=try(solve(results1$hessian)),
                                                               ses=FALSE));
                         }
                      else if(p == 1)
                         {
                         results2$ana.grads <- do.call(CDMN1P.Lik.gr, list(par=results1$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                         results2$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                     mean=results1$par,
                                                                     cov=try(solve(results1$hessian)),
                                                                     ses=FALSE)));
                         results2$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                               mean=results1$par,
                                                               cov=try(solve(results1$hessian)),
                                                               ses=FALSE));
                         }
                      else if(p == 2)
                         {
                         results2$ana.grads <- do.call(CDMN2P.Lik.gr, list(par=results1$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                         results2$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                                     mean=results1$par,
                                                                     cov=try(solve(results1$hessian)),
                                                                     ses=FALSE)));
                         results2$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),
                                                               mean=results1$par,
                                                               cov=try(solve(results1$hessian)),
                                                               ses=FALSE));
                         }
                      else if(p == 3)
                         {
                         results2$ana.grads <- do.call(CDMN3P.Lik.gr, list(par=results1$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                         results2$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                                     mean=results1$par,
                                                                     cov=try(solve(results1$hessian)),
                                                                     ses=FALSE)));
                         results2$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8)),
                                                               mean=results1$par,
                                                               cov=try(solve(results1$hessian)),
                                                               ses=FALSE));
                         }
                      else if(p == 4)
                         {
                         results2$ana.grads <- do.call(CDMN4P.Lik.gr, list(par=results1$par, dates=dates, obseff=obseff, obscat=obscat, M.fixed=M.fixed, distr=distr));
                         results2$bt.stdev  <- sqrt(diag(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                                     mean=results1$par,
                                                                     cov=try(solve(results1$nhatend)),
                                                                     ses=FALSE)));
                         results2$Cor       <- cor(deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8),~exp(x9)),
                                                     mean=results1$par,
                                                     cov=try(solve(results1$hessian)),
                                                     ses=FALSE));
                         }
                      }
                   else
                      {
                       stop('Numerical engine of choice (optim or optimx) must be specified')
                      }
               return(results2);
                  }
          }
