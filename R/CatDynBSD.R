CatDynBSD <-
function(x,method,multi,mbw.sd) #mbw.sd must be given in kg
  {
   #time step month
   if(multi)
     {
      if(class(x)!="catdyn")
        {
         stop("In multi-annual models 'x' must be a single object of class 'catdyn' run at monthly time steps")
        }
      if(x$Data$Properties$Units[3] == "ind")
        {
        stop("This function is used to calculate standard deviation of annual biomass when the catch is recorded in weight")
        }
      if(x$Data$Properties$Units["Time Step"]!="month")
        {
         stop("In multi-annual models 'x' must be a single object of class 'catdyn' run at monthly time steps")
        }
      Thou.scaler <- 1e6*(x$Data$Properties$Units[4]=="bill") +
                     1e3*(x$Data$Properties$Units[4]=="mill") +
                     1e0*(x$Data$Properties$Units[4]=="thou") +
                     1e-1*(x$Data$Properties$Units[4]=="hund")
      PopDyn <- data.frame(M=x$Model[[method]]$bt.par$M,
                           SE.M=x$Model[[method]]$bt.stdev["M"],
                           N0=Thou.scaler*x$Model[[method]]$bt.par$N0,
                           SE.N0=Thou.scaler*x$Model[[method]]$bt.stdev["N0"])
      #Replacing missing standard errors, if any
      if(is.na(PopDyn[2]))
        {
         PopDyn[2] <- PopDyn[1]*mean(unlist(x$Model[[method]]$bt.stdev[which(!is.na(x$Model[[method]]$bt.stdev))])/
                                     unlist(x$Model[[method]]$bt.par[which(!is.na(x$Model[[method]]$bt.stdev))]))
        }
      if(is.na(PopDyn[4]))
        {
         PopDyn[4] <- PopDyn[3]*mean(unlist(x$Model[[method]]$bt.stdev[which(!is.na(x$Model[[method]]$bt.stdev))])/
                                     unlist(x$Model[[method]]$bt.par[which(!is.na(x$Model[[method]]$bt.stdev))]))
        }
      Perts  <- data.frame(Pest=unlist(x$Model[[method]]$bt.par[grep("P",names(x$Model[[method]]$bt.par))])*Thou.scaler,
                           SE.Pest=unlist(x$Model[[method]]$bt.stdev[grep("P.",names(x$Model[[method]]$bt.par))])*Thou.scaler,
                           tsteps=x$Model[[method]]$Dates[grep("ts.P",names(x$Model[[method]]$Dates))])
      #Replacing missing standard errors
      if(any(is.na(Perts$SE.Pest)))
        {
         Perts$SE.Pest[which(is.na(Perts$SE.Pest))] <- Perts$Pest[which(is.na(Perts$SE.Pest))]*mean(Perts$SE.Pest[which(!is.na(Perts$SE.Pest))]/Perts$Pest[which(!is.na(Perts$SE.Pest))])
        }
      #One fleet
      if(length(x$Data$Properties$Fleets$Fleet)==1)
       {
         mt <- x$Model[[method]]$Type
         Timing <- matrix(0,12*mt,1)
         Timing[1:12]    <- ifelse(row(Timing)[1:12] >= Perts$tsteps[1],1,0)
         Timing[13:24]   <- ifelse(row(Timing)[13:24] >= Perts$tsteps[2],1,0)
         Timing[25:36]   <- ifelse(row(Timing)[25:36] >= Perts$tsteps[3],1,0)
         Timing[37:48]   <- ifelse(row(Timing)[37:48] >= Perts$tsteps[4],1,0)
         Timing[49:60]   <- ifelse(row(Timing)[49:60] >= Perts$tsteps[5],1,0)
         Timing[61:72]   <- ifelse(row(Timing)[61:72] >= Perts$tsteps[6],1,0)
         Timing[73:84]   <- ifelse(row(Timing)[73:84] >= Perts$tsteps[7],1,0)
         Timing[85:96]   <- ifelse(row(Timing)[85:96] >= Perts$tsteps[8],1,0)
         Timing[97:108]  <- ifelse(row(Timing)[97:108] >= Perts$tsteps[9],1,0)
         Timing[109:120] <- ifelse(row(Timing)[109:120] >= Perts$tsteps[10],1,0)
         Timing[121:132] <- ifelse(row(Timing)[121:132] >= Perts$tsteps[11],1,0)
         Timing[133:144] <- ifelse(row(Timing)[133:144] >= Perts$tsteps[12],1,0)
         Timing[145:156] <- ifelse(row(Timing)[145:156] >= Perts$tsteps[13],1,0)
         Timing[157:168] <- ifelse(row(Timing)[157:168] >= Perts$tsteps[14],1,0)
         Timing[169:180] <- ifelse(row(Timing)[169:180] >= Perts$tsteps[15],1,0)
         fleet1 <- x$Data$Properties$Fleets[1,1]
         if(mt <= 14)
           {
            stop("This function is intended to be used to calculate the standard deviation of annual biomass \n to fit a biomass dynamic model from the output of a multi-annual generalized depletion model (MAGD) \n with the catch recorded in biomass. At least 15 years of data must have been used in fitting \n the MAGD for its outputs to be used in this manner.")
           }
         Cov.Mat   <- cor2cov(cor.mat=x$Model[[method]]$Cor[c(1:(mt+2)),
                                                            c(1:(mt+2))],
                              sd=c(PopDyn$SE.M,PopDyn$SE.N0,Perts$SE.Pest))
         if(length(mbw.sd) != 12 & length(mbw.sd) != 12*mt)
           {stop("mbw.sd must be a vector of length 12 (monthly mean weight) or 12*number of years (in kg)")}
         yr1 <- as.numeric(format(as.Date(x$Data$Properties$Dates[1]),"%Y"))
         yr2 <- as.numeric(format(as.Date(x$Data$Properties$Dates[2]),"%Y"))
         z   <- CatDynPred(x,method)
         PredStock <- data.frame(Year=sort(rep(yr1:yr2,12)),
                                 Month=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                                 TimeStep=1:((yr2-yr1+1)*12),
                                 Mmw.kg=x$Data$Data[[fleet1]]$obsmbw.kg,
                                 SDmw.kg=mbw.sd,
                                 N.thou=z$Model$Results[,10],
                                 N.thou.SE=0,
                                 B.ton=z$Model$Results[,11],
                                 B.ton.SE=0)
         #Year 1
         for(m in 1:12)
          {
           PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                 mean=c(PopDyn$M,PopDyn$N0,Timing[m]*Perts$Pest[1]),
                                                 cov=Cov.Mat[c(1:3),c(1:3)])
           PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                          (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
          }
         #Year 2
         for(m in 13:24)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[12],Timing[m]*Perts$Pest[2]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,4],
                                                               0,PredStock$N.thou.SE[12]^2,0,
                                                               Cov.Mat[4,1],0,Cov.Mat[4,4]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 3
         for(m in 25:36)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[24],Timing[m]*Perts$Pest[3]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,5],
                                                               0,PredStock$N.thou.SE[24]^2,0,
                                                               Cov.Mat[5,1],0,Cov.Mat[5,5]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 4
         for(m in 37:48)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[36],Timing[m]*Perts$Pest[4]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,6],
                                                               0,PredStock$N.thou.SE[36]^2,0,
                                                               Cov.Mat[6,1],0,Cov.Mat[6,6]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 5
         for(m in 49:60)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[48],Timing[m]*Perts$Pest[5]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,7],
                                                               0,PredStock$N.thou.SE[48]^2,0,
                                                               Cov.Mat[7,1],0,Cov.Mat[7,7]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 6
         for(m in 61:72)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[60],Timing[m]*Perts$Pest[6]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,8],
                                                               0,PredStock$N.thou.SE[60]^2,0,
                                                               Cov.Mat[8,1],0,Cov.Mat[8,8]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 7
         for(m in 73:84)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[72],Timing[m]*Perts$Pest[7]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,9],
                                                               0,PredStock$N.thou.SE[72]^2,0,
                                                               Cov.Mat[9,1],0,Cov.Mat[9,9]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 8
         for(m in 85:96)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[84],Timing[m]*Perts$Pest[8]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,10],
                                                               0,PredStock$N.thou.SE[84]^2,0,
                                                               Cov.Mat[10,1],0,Cov.Mat[10,10]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 9
         for(m in 97:108)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[96],Timing[m]*Perts$Pest[9]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,11],
                                                               0,PredStock$N.thou.SE[96]^2,0,
                                                               Cov.Mat[11,1],0,Cov.Mat[11,11]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 10
         for(m in 109:120)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[108],Timing[m]*Perts$Pest[10]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,12],
                                                               0,PredStock$N.thou.SE[108]^2,0,
                                                               Cov.Mat[12,1],0,Cov.Mat[12,12]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 11
         for(m in 121:132)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[120],Timing[m,1]*Perts$Pest[11]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,13],
                                                               0,PredStock$N.thou.SE[120]^2,0,
                                                               Cov.Mat[13,1],0,Cov.Mat[13,13]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 12
         for(m in 133:144)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[132],Timing[m]*Perts$Pest[12]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,14],
                                                               0,PredStock$N.thou.SE[132]^2,0,
                                                               Cov.Mat[14,1],0,Cov.Mat[14,14]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 13
         for(m in 145:156)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[144],Timing[m]*Perts$Pest[13]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,15],
                                                               0,PredStock$N.thou.SE[144]^2,0,
                                                               Cov.Mat[15,1],0,Cov.Mat[15,15]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 14
         for(m in 157:168)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[156],Timing[m]*Perts$Pest[14]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,16],
                                                               0,PredStock$N.thou.SE[156]^2,0,
                                                               Cov.Mat[16,1],0,Cov.Mat[16,16]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 15
         for(m in 169:180)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[168],Timing[m]*Perts$Pest[15]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,17],
                                                               0,PredStock$N.thou.SE[168]^2,0,
                                                               Cov.Mat[17,1],0,Cov.Mat[17,17]),3,3))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 16
         if(mt == 16)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 17)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 18)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 19)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228] <- ifelse(row(Timing)[217:228] >= Perts$tsteps[19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m]*Perts$Pest[19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],
                                                                  0,PredStock$N.thou.SE[216]^2,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 20)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228] <- ifelse(row(Timing)[217:228] >= Perts$tsteps[19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m]*Perts$Pest[19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],
                                                                  0,PredStock$N.thou.SE[216]^2,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240] <- ifelse(row(Timing)[229:240] >= Perts$tsteps[20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m]*Perts$Pest[20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],
                                                                  0,PredStock$N.thou.SE[228]^2,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 21)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228] <- ifelse(row(Timing)[217:228] >= Perts$tsteps[19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m]*Perts$Pest[19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],
                                                                  0,PredStock$N.thou.SE[216]^2,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240] <- ifelse(row(Timing)[229:240] >= Perts$tsteps[20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m]*Perts$Pest[20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],
                                                                  0,PredStock$N.thou.SE[228]^2,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252] <- ifelse(row(Timing)[241:252] >= Perts$tsteps[21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m]*Perts$Pest[21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],
                                                                  0,PredStock$N.thou.SE[240]^2,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 22)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228] <- ifelse(row(Timing)[217:228] >= Perts$tsteps[19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m]*Perts$Pest[19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],
                                                                  0,PredStock$N.thou.SE[216]^2,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240] <- ifelse(row(Timing)[229:240] >= Perts$tsteps[20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m]*Perts$Pest[20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],
                                                                  0,PredStock$N.thou.SE[228]^2,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252] <- ifelse(row(Timing)[241:252] >= Perts$tsteps[21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m]*Perts$Pest[21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],
                                                                  0,PredStock$N.thou.SE[240]^2,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264] <- ifelse(row(Timing)[253:264] >= Perts$tsteps[22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m]*Perts$Pest[22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],
                                                                  0,PredStock$N.thou.SE[252]^2,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 23)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228] <- ifelse(row(Timing)[217:228] >= Perts$tsteps[19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m]*Perts$Pest[19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],
                                                                  0,PredStock$N.thou.SE[216]^2,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240] <- ifelse(row(Timing)[229:240] >= Perts$tsteps[20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m]*Perts$Pest[20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],
                                                                  0,PredStock$N.thou.SE[228]^2,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252] <- ifelse(row(Timing)[241:252] >= Perts$tsteps[21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m]*Perts$Pest[21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],
                                                                  0,PredStock$N.thou.SE[240]^2,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264] <- ifelse(row(Timing)[253:264] >= Perts$tsteps[22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m]*Perts$Pest[22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],
                                                                  0,PredStock$N.thou.SE[252]^2,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[265:276] <- ifelse(row(Timing)[265:276] >= Perts$tsteps[23],1,0)
            for(m in 265:276)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[264],Timing[m]*Perts$Pest[23]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,25],
                                                                  0,PredStock$N.thou.SE[264]^2,0,
                                                                  Cov.Mat[25,1],0,Cov.Mat[25,25]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 24)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228] <- ifelse(row(Timing)[217:228] >= Perts$tsteps[19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m]*Perts$Pest[19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],
                                                                  0,PredStock$N.thou.SE[216]^2,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240] <- ifelse(row(Timing)[229:240] >= Perts$tsteps[20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m]*Perts$Pest[20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],
                                                                  0,PredStock$N.thou.SE[228]^2,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252] <- ifelse(row(Timing)[241:252] >= Perts$tsteps[21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m]*Perts$Pest[21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],
                                                                  0,PredStock$N.thou.SE[240]^2,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264] <- ifelse(row(Timing)[253:264] >= Perts$tsteps[22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m]*Perts$Pest[22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],
                                                                  0,PredStock$N.thou.SE[252]^2,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[265:276] <- ifelse(row(Timing)[265:276] >= Perts$tsteps[23],1,0)
            for(m in 265:276)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[264],Timing[m]*Perts$Pest[23]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,25],
                                                                  0,PredStock$N.thou.SE[264]^2,0,
                                                                  Cov.Mat[25,1],0,Cov.Mat[25,25]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[277:288] <- ifelse(row(Timing)[277:288] >= Perts$tsteps[24],1,0)
            for(m in 277:288)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[276],Timing[m]*Perts$Pest[24]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,26],
                                                                  0,PredStock$N.thou.SE[276]^2,0,
                                                                  Cov.Mat[26,1],0,Cov.Mat[26,26]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 25)
           {
            Timing[181:192] <- ifelse(row(Timing)[181:192] >= Perts$tsteps[16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m]*Perts$Pest[16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],
                                                                  0,PredStock$N.thou.SE[180]^2,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204] <- ifelse(row(Timing)[193:204] >= Perts$tsteps[17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m]*Perts$Pest[17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],
                                                                  0,PredStock$N.thou.SE[192]^2,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216] <- ifelse(row(Timing)[205:216] >= Perts$tsteps[18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m]*Perts$Pest[18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],
                                                                  0,PredStock$N.thou.SE[204]^2,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228] <- ifelse(row(Timing)[217:228] >= Perts$tsteps[19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m]*Perts$Pest[19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],
                                                                  0,PredStock$N.thou.SE[216]^2,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240] <- ifelse(row(Timing)[229:240] >= Perts$tsteps[20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m]*Perts$Pest[20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],
                                                                  0,PredStock$N.thou.SE[228]^2,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252] <- ifelse(row(Timing)[241:252] >= Perts$tsteps[21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m]*Perts$Pest[21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],
                                                                  0,PredStock$N.thou.SE[240]^2,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264] <- ifelse(row(Timing)[253:264] >= Perts$tsteps[22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m]*Perts$Pest[22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],
                                                                  0,PredStock$N.thou.SE[252]^2,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[265:276] <- ifelse(row(Timing)[265:276] >= Perts$tsteps[23],1,0)
            for(m in 265:276)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[264],Timing[m]*Perts$Pest[23]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,25],
                                                                  0,PredStock$N.thou.SE[264]^2,0,
                                                                  Cov.Mat[25,1],0,Cov.Mat[25,25]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[277:288] <- ifelse(row(Timing)[277:288] >= Perts$tsteps[24],1,0)
            for(m in 277:288)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[276],Timing[m]*Perts$Pest[24]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,26],
                                                                  0,PredStock$N.thou.SE[276]^2,0,
                                                                  Cov.Mat[26,1],0,Cov.Mat[26,26]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[289:300] <- ifelse(row(Timing)[289:300] >= Perts$tsteps[25],1,0)
            for(m in 289:300)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[288],Timing[m]*Perts$Pest[25]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,27],
                                                                  0,PredStock$N.thou.SE[288]^2,0,
                                                                  Cov.Mat[27,1],0,Cov.Mat[27,27]),3,3))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
       #End of one fleet
       }
      if(length(x$Data$Properties$Fleets$Fleet)==2)
       {
         mt <- x$Model[[method]]$Type[1]
         Timing <- matrix(0,12*mt,2)
         Timing[1:12,1]    <- ifelse(row(Timing)[1:12,1] >= Perts$tsteps[1],1,0)
         Timing[1:12,2]    <- ifelse(row(Timing)[1:12,2] >= Perts$tsteps[mt+1],1,0)
         Timing[13:24,1]   <- ifelse(row(Timing)[13:24,1] >= Perts$tsteps[2],1,0)
         Timing[13:24,2]   <- ifelse(row(Timing)[13:24,2] >= Perts$tsteps[mt+2],1,0)
         Timing[25:36,1]   <- ifelse(row(Timing)[25:36,1] >= Perts$tsteps[3],1,0)
         Timing[25:36,2]   <- ifelse(row(Timing)[25:36,2] >= Perts$tsteps[mt+3],1,0)
         Timing[37:48,1]   <- ifelse(row(Timing)[37:48,1] >= Perts$tsteps[4],1,0)
         Timing[37:48,2]   <- ifelse(row(Timing)[37:48,2] >= Perts$tsteps[mt+4],1,0)
         Timing[49:60,1]   <- ifelse(row(Timing)[49:60,1] >= Perts$tsteps[5],1,0)
         Timing[49:60,2]   <- ifelse(row(Timing)[49:60,2] >= Perts$tsteps[mt+5],1,0)
         Timing[61:72,1]   <- ifelse(row(Timing)[61:72,1] >= Perts$tsteps[6],1,0)
         Timing[61:72,2]   <- ifelse(row(Timing)[61:72,2] >= Perts$tsteps[mt+6],1,0)
         Timing[73:84,1]   <- ifelse(row(Timing)[73:84,1] >= Perts$tsteps[7],1,0)
         Timing[73:84,2]   <- ifelse(row(Timing)[73:84,2] >= Perts$tsteps[mt+7],1,0)
         Timing[85:96,1]   <- ifelse(row(Timing)[85:96,1] >= Perts$tsteps[8],1,0)
         Timing[85:96,2]   <- ifelse(row(Timing)[85:96,2] >= Perts$tsteps[mt+8],1,0)
         Timing[97:108,1]  <- ifelse(row(Timing)[97:108,1] >= Perts$tsteps[9],1,0)
         Timing[97:108,2]  <- ifelse(row(Timing)[97:108,2] >= Perts$tsteps[mt+9],1,0)
         Timing[109:120,1] <- ifelse(row(Timing)[109:120,1] >= Perts$tsteps[10],1,0)
         Timing[109:120,2] <- ifelse(row(Timing)[109:120,2] >= Perts$tsteps[mt+10],1,0)
         Timing[121:132,1] <- ifelse(row(Timing)[121:132,1] >= Perts$tsteps[11],1,0)
         Timing[121:132,2] <- ifelse(row(Timing)[121:132,2] >= Perts$tsteps[mt+11],1,0)
         Timing[133:144,1] <- ifelse(row(Timing)[133:144,1] >= Perts$tsteps[12],1,0)
         Timing[133:144,2] <- ifelse(row(Timing)[133:144,2] >= Perts$tsteps[mt+12],1,0)
         Timing[145:156,1] <- ifelse(row(Timing)[145:156,1] >= Perts$tsteps[13],1,0)
         Timing[145:156,2] <- ifelse(row(Timing)[145:156,2] >= Perts$tsteps[mt+13],1,0)
         Timing[157:168,1] <- ifelse(row(Timing)[157:168,1] >= Perts$tsteps[14],1,0)
         Timing[157:168,2] <- ifelse(row(Timing)[157:168,2] >= Perts$tsteps[mt+14],1,0)
         Timing[169:180,1] <- ifelse(row(Timing)[169:180,1] >= Perts$tsteps[15],1,0)
         Timing[169:180,2] <- ifelse(row(Timing)[169:180,2] >= Perts$tsteps[mt+15],1,0)
         fleet1 <- x$Data$Properties$Fleets[1,1]
         fleet2 <- x$Data$Properties$Fleets[2,1]
         if(mt <= 14)
           {
            stop("This function is intended to be used to calculate the standard deviation of annual biomass \n to fit a biomass dynamic model from the output of a multi-annual generalized depletion model (MAGD) \n with the catch recorded in biomass. At least 15 years of data must have been used in fitting \n the MAGD for its outputs to be used in this manner.")
           }
         Cov.Mat   <- cor2cov(cor.mat=x$Model[[method]]$Cor[c(1:(mt+2),(mt+6):(2*(mt)+5)),
                                                            c(1:(mt+2),(mt+6):(2*(mt)+5))],
                              sd=c(PopDyn$SE.M,PopDyn$SE.N0,Perts$SE.Pest))
         if(length(mbw.sd) != 12 & length(mbw.sd) != 12*mt)
           {stop("mbw.sd must be a vector of length 12 (monthly mean weight) or 12*number of years, all in kg")}
         yr1 <- as.numeric(format(as.Date(x$Data$Properties$Dates[1]),"%Y"))
         yr2 <- as.numeric(format(as.Date(x$Data$Properties$Dates[2]),"%Y"))
         z   <- CatDynPred(x,method)
         PredStock <- data.frame(Year=sort(rep(yr1:yr2,12)),
                                 Month=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                                 TimeStep=1:((yr2-yr1+1)*12),
                                 Mmw.kg=(x$Data$Data[[fleet1]]$obsmbw.kg+x$Data$Data[[fleet2]]$obsmbw.kg)/2,
                                 SDmw.kg=mbw.sd,
                                 N.thou=z$Model$Results[,18],
                                 N.thou.SE=0,
                                 B.ton=z$Model$Results[,19],
                                 B.ton.SE=0)
         #Year 1
         for(m in 1:12)
          {
           PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                 mean=c(PopDyn$M,PopDyn$N0,Timing[m,1]*Perts$Pest[1],Timing[m,2]*Perts$Pest[mt+1]),
                                                 cov=Cov.Mat[c(1:3,23),c(1:3,23)])
           PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                          (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
          }
         #Year 2
         for(m in 13:24)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[12],Timing[m,1]*Perts$Pest[2],Timing[m,2]*Perts$Pest[mt+2]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,4],Cov.Mat[1,24],
                                                               0,PredStock$N.thou.SE[12]^2,0,0,
                                                               Cov.Mat[4,1],0,Cov.Mat[4,4],Cov.Mat[4,24],
                                                               Cov.Mat[24,1],0,Cov.Mat[24,4],Cov.Mat[24,24]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 3
         for(m in 25:36)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[24],Timing[m,1]*Perts$Pest[3],Timing[m,2]*Perts$Pest[mt+3]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,5],Cov.Mat[1,25],
                                                               0,PredStock$N.thou.SE[24]^2,0,0,
                                                               Cov.Mat[5,1],0,Cov.Mat[5,5],Cov.Mat[5,25],
                                                               Cov.Mat[25,1],0,Cov.Mat[25,5],Cov.Mat[25,25]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 4
         for(m in 37:48)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[36],Timing[m,1]*Perts$Pest[4],Timing[m,2]*Perts$Pest[mt+4]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,6],Cov.Mat[1,26],
                                                               0,PredStock$N.thou.SE[36]^2,0,0,
                                                               Cov.Mat[6,1],0,Cov.Mat[6,6],Cov.Mat[6,26],
                                                               Cov.Mat[26,1],0,Cov.Mat[26,6],Cov.Mat[26,26]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 5
         for(m in 49:60)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[48],Timing[m,1]*Perts$Pest[5],Timing[m,2]*Perts$Pest[mt+5]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,7],Cov.Mat[1,27],
                                                               0,PredStock$N.thou.SE[48]^2,0,0,
                                                               Cov.Mat[7,1],0,Cov.Mat[7,7],Cov.Mat[7,27],
                                                               Cov.Mat[27,1],0,Cov.Mat[27,7],Cov.Mat[27,27]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 6
         for(m in 61:72)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[60],Timing[m,1]*Perts$Pest[6],Timing[m,2]*Perts$Pest[mt+6]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,8],Cov.Mat[1,28],
                                                               0,PredStock$N.thou.SE[60]^2,0,0,
                                                               Cov.Mat[8,1],0,Cov.Mat[8,8],Cov.Mat[8,28],
                                                               Cov.Mat[28,1],0,Cov.Mat[28,8],Cov.Mat[28,28]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 7
         for(m in 73:84)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[72],Timing[m,1]*Perts$Pest[7],Timing[m,2]*Perts$Pest[mt+7]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,9],Cov.Mat[1,29],
                                                               0,PredStock$N.thou.SE[72]^2,0,0,
                                                               Cov.Mat[9,1],0,Cov.Mat[9,9],Cov.Mat[9,29],
                                                               Cov.Mat[29,1],0,Cov.Mat[29,9],Cov.Mat[29,29]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 8
         for(m in 85:96)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[84],Timing[m,1]*Perts$Pest[8],Timing[m,2]*Perts$Pest[mt+8]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,10],Cov.Mat[1,30],
                                                               0,PredStock$N.thou.SE[84]^2,0,0,
                                                               Cov.Mat[10,1],0,Cov.Mat[10,10],Cov.Mat[10,30],
                                                               Cov.Mat[30,1],0,Cov.Mat[30,10],Cov.Mat[30,30]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 9
         for(m in 97:108)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[96],Timing[m,1]*Perts$Pest[9],Timing[m,2]*Perts$Pest[mt+9]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,11],Cov.Mat[1,31],
                                                               0,PredStock$N.thou.SE[96]^2,0,0,
                                                               Cov.Mat[11,1],0,Cov.Mat[11,11],Cov.Mat[11,31],
                                                               Cov.Mat[31,1],0,Cov.Mat[31,11],Cov.Mat[31,31]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 10
         for(m in 109:120)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[108],Timing[m,1]*Perts$Pest[10],Timing[m,2]*Perts$Pest[mt+10]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,12],Cov.Mat[1,32],
                                                               0,PredStock$N.thou.SE[108]^2,0,0,
                                                               Cov.Mat[12,1],0,Cov.Mat[12,12],Cov.Mat[12,32],
                                                               Cov.Mat[32,1],0,Cov.Mat[32,12],Cov.Mat[32,32]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 11
         for(m in 121:132)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[120],Timing[m,1]*Perts$Pest[11],Timing[m,2]*Perts$Pest[mt+11]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,13],Cov.Mat[1,33],
                                                               0,PredStock$N.thou.SE[120]^2,0,0,
                                                               Cov.Mat[13,1],0,Cov.Mat[13,13],Cov.Mat[13,33],
                                                               Cov.Mat[33,1],0,Cov.Mat[33,13],Cov.Mat[33,33]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 12
         for(m in 133:144)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[132],Timing[m,1]*Perts$Pest[12],Timing[m,2]*Perts$Pest[mt+12]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,14],Cov.Mat[1,34],
                                                               0,PredStock$N.thou.SE[132]^2,0,0,
                                                               Cov.Mat[14,1],0,Cov.Mat[14,14],Cov.Mat[14,34],
                                                               Cov.Mat[34,1],0,Cov.Mat[34,14],Cov.Mat[34,34]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 13
         for(m in 145:156)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[144],Timing[m,1]*Perts$Pest[13],Timing[m,2]*Perts$Pest[mt+13]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,15],Cov.Mat[1,35],
                                                               0,PredStock$N.thou.SE[144]^2,0,0,
                                                               Cov.Mat[15,1],0,Cov.Mat[15,15],Cov.Mat[15,35],
                                                               Cov.Mat[35,1],0,Cov.Mat[35,15],Cov.Mat[35,35]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 14
         for(m in 157:168)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[156],Timing[m,1]*Perts$Pest[14],Timing[m,2]*Perts$Pest[mt+14]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,16],Cov.Mat[1,36],
                                                               0,PredStock$N.thou.SE[156]^2,0,0,
                                                               Cov.Mat[16,1],0,Cov.Mat[16,16],Cov.Mat[16,36],
                                                               Cov.Mat[36,1],0,Cov.Mat[36,16],Cov.Mat[36,36]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 15
         for(m in 169:180)
           {
            PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                  mean=c(PopDyn$M,PredStock$N.thou[168],Timing[m,1]*Perts$Pest[15],Timing[m,2]*Perts$Pest[mt+15]),
                                                  cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,17],Cov.Mat[1,37],
                                                               0,PredStock$N.thou.SE[168]^2,0,0,
                                                               Cov.Mat[17,1],0,Cov.Mat[17,17],Cov.Mat[17,37],
                                                               Cov.Mat[37,1],0,Cov.Mat[37,17],Cov.Mat[37,37]),4,4))
            PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                           (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
           }
         #Year 16
         if(mt == 16)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 17)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 18)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 19)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228,1] <- ifelse(row(Timing)[217:228,1] >= Perts$tsteps[19],1,0)
            Timing[217:228,2] <- ifelse(row(Timing)[217:228,2] >= Perts$tsteps[mt+19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m,1]*Perts$Pest[19],Timing[m,2]*Perts$Pest[mt+19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],Cov.Mat[1,41],
                                                                  0,PredStock$N.thou.SE[216]^2,0,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21],Cov.Mat[21,41],
                                                                  Cov.Mat[41,1],0,Cov.Mat[41,21],Cov.Mat[41,41]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 20)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228,1] <- ifelse(row(Timing)[217:228,1] >= Perts$tsteps[19],1,0)
            Timing[217:228,2] <- ifelse(row(Timing)[217:228,2] >= Perts$tsteps[mt+19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m,1]*Perts$Pest[19],Timing[m,2]*Perts$Pest[mt+19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],Cov.Mat[1,41],
                                                                  0,PredStock$N.thou.SE[216]^2,0,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21],Cov.Mat[21,41],
                                                                  Cov.Mat[41,1],0,Cov.Mat[41,21],Cov.Mat[41,41]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240,1] <- ifelse(row(Timing)[229:240,1] >= Perts$tsteps[20],1,0)
            Timing[229:240,2] <- ifelse(row(Timing)[229:240,2] >= Perts$tsteps[mt+20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m,1]*Perts$Pest[20],Timing[m,2]*Perts$Pest[mt+20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],Cov.Mat[1,42],
                                                                  0,PredStock$N.thou.SE[228]^2,0,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22],Cov.Mat[22,42],
                                                                  Cov.Mat[42,1],0,Cov.Mat[42,22],Cov.Mat[42,42]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 21)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228,1] <- ifelse(row(Timing)[217:228,1] >= Perts$tsteps[19],1,0)
            Timing[217:228,2] <- ifelse(row(Timing)[217:228,2] >= Perts$tsteps[mt+19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m,1]*Perts$Pest[19],Timing[m,2]*Perts$Pest[mt+19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],Cov.Mat[1,41],
                                                                  0,PredStock$N.thou.SE[216]^2,0,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21],Cov.Mat[21,41],
                                                                  Cov.Mat[41,1],0,Cov.Mat[41,21],Cov.Mat[41,41]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240,1] <- ifelse(row(Timing)[229:240,1] >= Perts$tsteps[20],1,0)
            Timing[229:240,2] <- ifelse(row(Timing)[229:240,2] >= Perts$tsteps[mt+20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m,1]*Perts$Pest[20],Timing[m,2]*Perts$Pest[mt+20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],Cov.Mat[1,42],
                                                                  0,PredStock$N.thou.SE[228]^2,0,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22],Cov.Mat[22,42],
                                                                  Cov.Mat[42,1],0,Cov.Mat[42,22],Cov.Mat[42,42]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252,1] <- ifelse(row(Timing)[241:252,1] >= Perts$tsteps[21],1,0)
            Timing[241:252,2] <- ifelse(row(Timing)[241:252,2] >= Perts$tsteps[mt+21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m,1]*Perts$Pest[21],Timing[m,2]*Perts$Pest[mt+21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],Cov.Mat[1,43],
                                                                  0,PredStock$N.thou.SE[240]^2,0,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23],Cov.Mat[23,43],
                                                                  Cov.Mat[43,1],0,Cov.Mat[43,23],Cov.Mat[43,43]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 22)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228,1] <- ifelse(row(Timing)[217:228,1] >= Perts$tsteps[19],1,0)
            Timing[217:228,2] <- ifelse(row(Timing)[217:228,2] >= Perts$tsteps[mt+19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m,1]*Perts$Pest[19],Timing[m,2]*Perts$Pest[mt+19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],Cov.Mat[1,41],
                                                                  0,PredStock$N.thou.SE[216]^2,0,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21],Cov.Mat[21,41],
                                                                  Cov.Mat[41,1],0,Cov.Mat[41,21],Cov.Mat[41,41]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240,1] <- ifelse(row(Timing)[229:240,1] >= Perts$tsteps[20],1,0)
            Timing[229:240,2] <- ifelse(row(Timing)[229:240,2] >= Perts$tsteps[mt+20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m,1]*Perts$Pest[20],Timing[m,2]*Perts$Pest[mt+20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],Cov.Mat[1,42],
                                                                  0,PredStock$N.thou.SE[228]^2,0,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22],Cov.Mat[22,42],
                                                                  Cov.Mat[42,1],0,Cov.Mat[42,22],Cov.Mat[42,42]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252,1] <- ifelse(row(Timing)[241:252,1] >= Perts$tsteps[21],1,0)
            Timing[241:252,2] <- ifelse(row(Timing)[241:252,2] >= Perts$tsteps[mt+21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m,1]*Perts$Pest[21],Timing[m,2]*Perts$Pest[mt+21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],Cov.Mat[1,43],
                                                                  0,PredStock$N.thou.SE[240]^2,0,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23],Cov.Mat[23,43],
                                                                  Cov.Mat[43,1],0,Cov.Mat[43,23],Cov.Mat[43,43]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264,1] <- ifelse(row(Timing)[253:264,1] >= Perts$tsteps[22],1,0)
            Timing[253:264,2] <- ifelse(row(Timing)[253:264,2] >= Perts$tsteps[mt+22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m,1]*Perts$Pest[22],Timing[m,2]*Perts$Pest[mt+22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],Cov.Mat[1,44],
                                                                  0,PredStock$N.thou.SE[252]^2,0,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24],Cov.Mat[24,44],
                                                                  Cov.Mat[44,1],0,Cov.Mat[44,24],Cov.Mat[44,44]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 23)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228,1] <- ifelse(row(Timing)[217:228,1] >= Perts$tsteps[19],1,0)
            Timing[217:228,2] <- ifelse(row(Timing)[217:228,2] >= Perts$tsteps[mt+19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m,1]*Perts$Pest[19],Timing[m,2]*Perts$Pest[mt+19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],Cov.Mat[1,41],
                                                                  0,PredStock$N.thou.SE[216]^2,0,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21],Cov.Mat[21,41],
                                                                  Cov.Mat[41,1],0,Cov.Mat[41,21],Cov.Mat[41,41]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240,1] <- ifelse(row(Timing)[229:240,1] >= Perts$tsteps[20],1,0)
            Timing[229:240,2] <- ifelse(row(Timing)[229:240,2] >= Perts$tsteps[mt+20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m,1]*Perts$Pest[20],Timing[m,2]*Perts$Pest[mt+20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],Cov.Mat[1,42],
                                                                  0,PredStock$N.thou.SE[228]^2,0,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22],Cov.Mat[22,42],
                                                                  Cov.Mat[42,1],0,Cov.Mat[42,22],Cov.Mat[42,42]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252,1] <- ifelse(row(Timing)[241:252,1] >= Perts$tsteps[21],1,0)
            Timing[241:252,2] <- ifelse(row(Timing)[241:252,2] >= Perts$tsteps[mt+21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m,1]*Perts$Pest[21],Timing[m,2]*Perts$Pest[mt+21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],Cov.Mat[1,43],
                                                                  0,PredStock$N.thou.SE[240]^2,0,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23],Cov.Mat[23,43],
                                                                  Cov.Mat[43,1],0,Cov.Mat[43,23],Cov.Mat[43,43]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264,1] <- ifelse(row(Timing)[253:264,1] >= Perts$tsteps[22],1,0)
            Timing[253:264,2] <- ifelse(row(Timing)[253:264,2] >= Perts$tsteps[mt+22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m,1]*Perts$Pest[22],Timing[m,2]*Perts$Pest[mt+22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],Cov.Mat[1,44],
                                                                  0,PredStock$N.thou.SE[252]^2,0,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24],Cov.Mat[24,44],
                                                                  Cov.Mat[44,1],0,Cov.Mat[44,24],Cov.Mat[44,44]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[265:276,1] <- ifelse(row(Timing)[265:276,1] >= Perts$tsteps[23],1,0)
            Timing[265:276,2] <- ifelse(row(Timing)[265:276,2] >= Perts$tsteps[mt+23],1,0)
            for(m in 265:276)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[264],Timing[m,1]*Perts$Pest[23],Timing[m,2]*Perts$Pest[mt+23]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,25],Cov.Mat[1,45],
                                                                  0,PredStock$N.thou.SE[264]^2,0,0,
                                                                  Cov.Mat[25,1],0,Cov.Mat[25,25],Cov.Mat[25,45],
                                                                  Cov.Mat[45,1],0,Cov.Mat[45,25],Cov.Mat[45,45]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 24)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228,1] <- ifelse(row(Timing)[217:228,1] >= Perts$tsteps[19],1,0)
            Timing[217:228,2] <- ifelse(row(Timing)[217:228,2] >= Perts$tsteps[mt+19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m,1]*Perts$Pest[19],Timing[m,2]*Perts$Pest[mt+19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],Cov.Mat[1,41],
                                                                  0,PredStock$N.thou.SE[216]^2,0,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21],Cov.Mat[21,41],
                                                                  Cov.Mat[41,1],0,Cov.Mat[41,21],Cov.Mat[41,41]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240,1] <- ifelse(row(Timing)[229:240,1] >= Perts$tsteps[20],1,0)
            Timing[229:240,2] <- ifelse(row(Timing)[229:240,2] >= Perts$tsteps[mt+20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m,1]*Perts$Pest[20],Timing[m,2]*Perts$Pest[mt+20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],Cov.Mat[1,42],
                                                                  0,PredStock$N.thou.SE[228]^2,0,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22],Cov.Mat[22,42],
                                                                  Cov.Mat[42,1],0,Cov.Mat[42,22],Cov.Mat[42,42]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252,1] <- ifelse(row(Timing)[241:252,1] >= Perts$tsteps[21],1,0)
            Timing[241:252,2] <- ifelse(row(Timing)[241:252,2] >= Perts$tsteps[mt+21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m,1]*Perts$Pest[21],Timing[m,2]*Perts$Pest[mt+21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],Cov.Mat[1,43],
                                                                  0,PredStock$N.thou.SE[240]^2,0,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23],Cov.Mat[23,43],
                                                                  Cov.Mat[43,1],0,Cov.Mat[43,23],Cov.Mat[43,43]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264,1] <- ifelse(row(Timing)[253:264,1] >= Perts$tsteps[22],1,0)
            Timing[253:264,2] <- ifelse(row(Timing)[253:264,2] >= Perts$tsteps[mt+22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m,1]*Perts$Pest[22],Timing[m,2]*Perts$Pest[mt+22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],Cov.Mat[1,44],
                                                                  0,PredStock$N.thou.SE[252]^2,0,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24],Cov.Mat[24,44],
                                                                  Cov.Mat[44,1],0,Cov.Mat[44,24],Cov.Mat[44,44]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[265:276,1] <- ifelse(row(Timing)[265:276,1] >= Perts$tsteps[23],1,0)
            Timing[265:276,2] <- ifelse(row(Timing)[265:276,2] >= Perts$tsteps[mt+23],1,0)
            for(m in 265:276)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[264],Timing[m,1]*Perts$Pest[23],Timing[m,2]*Perts$Pest[mt+23]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,25],Cov.Mat[1,45],
                                                                  0,PredStock$N.thou.SE[264]^2,0,0,
                                                                  Cov.Mat[25,1],0,Cov.Mat[25,25],Cov.Mat[25,45],
                                                                  Cov.Mat[45,1],0,Cov.Mat[45,25],Cov.Mat[45,45]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[277:288,1] <- ifelse(row(Timing)[277:288,1] >= Perts$tsteps[24],1,0)
            Timing[277:288,2] <- ifelse(row(Timing)[277:288,2] >= Perts$tsteps[mt+24],1,0)
            for(m in 277:288)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[276],Timing[m,1]*Perts$Pest[24],Timing[m,2]*Perts$Pest[mt+24]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,26],Cov.Mat[1,46],
                                                                  0,PredStock$N.thou.SE[276]^2,0,0,
                                                                  Cov.Mat[26,1],0,Cov.Mat[26,26],Cov.Mat[26,46],
                                                                  Cov.Mat[46,1],0,Cov.Mat[46,26],Cov.Mat[46,46]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
         if(mt == 25)
           {
            Timing[181:192,1] <- ifelse(row(Timing)[181:192,1] >= Perts$tsteps[16],1,0)
            Timing[181:192,2] <- ifelse(row(Timing)[181:192,2] >= Perts$tsteps[mt+16],1,0)
            for(m in 181:192)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[180],Timing[m,1]*Perts$Pest[16],Timing[m,2]*Perts$Pest[mt+16]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,18],Cov.Mat[1,38],
                                                                  0,PredStock$N.thou.SE[180]^2,0,0,
                                                                  Cov.Mat[18,1],0,Cov.Mat[18,18],Cov.Mat[18,38],
                                                                  Cov.Mat[38,1],0,Cov.Mat[38,18],Cov.Mat[38,38]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[193:204,1] <- ifelse(row(Timing)[193:204,1] >= Perts$tsteps[17],1,0)
            Timing[193:204,2] <- ifelse(row(Timing)[193:204,2] >= Perts$tsteps[mt+17],1,0)
            for(m in 193:204)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[192],Timing[m,1]*Perts$Pest[17],Timing[m,2]*Perts$Pest[mt+17]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,19],Cov.Mat[1,39],
                                                                  0,PredStock$N.thou.SE[192]^2,0,0,
                                                                  Cov.Mat[19,1],0,Cov.Mat[19,19],Cov.Mat[19,39],
                                                                  Cov.Mat[39,1],0,Cov.Mat[39,19],Cov.Mat[39,39]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[205:216,1] <- ifelse(row(Timing)[205:216,1] >= Perts$tsteps[18],1,0)
            Timing[205:216,2] <- ifelse(row(Timing)[205:216,2] >= Perts$tsteps[mt+18],1,0)
            for(m in 205:216)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[204],Timing[m,1]*Perts$Pest[18],Timing[m,2]*Perts$Pest[mt+18]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,20],Cov.Mat[1,40],
                                                                  0,PredStock$N.thou.SE[204]^2,0,0,
                                                                  Cov.Mat[20,1],0,Cov.Mat[20,20],Cov.Mat[20,40],
                                                                  Cov.Mat[40,1],0,Cov.Mat[40,20],Cov.Mat[40,40]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[217:228,1] <- ifelse(row(Timing)[217:228,1] >= Perts$tsteps[19],1,0)
            Timing[217:228,2] <- ifelse(row(Timing)[217:228,2] >= Perts$tsteps[mt+19],1,0)
            for(m in 217:228)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[216],Timing[m,1]*Perts$Pest[19],Timing[m,2]*Perts$Pest[mt+19]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,21],Cov.Mat[1,41],
                                                                  0,PredStock$N.thou.SE[216]^2,0,0,
                                                                  Cov.Mat[21,1],0,Cov.Mat[21,21],Cov.Mat[21,41],
                                                                  Cov.Mat[41,1],0,Cov.Mat[41,21],Cov.Mat[41,41]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[229:240,1] <- ifelse(row(Timing)[229:240,1] >= Perts$tsteps[20],1,0)
            Timing[229:240,2] <- ifelse(row(Timing)[229:240,2] >= Perts$tsteps[mt+20],1,0)
            for(m in 229:240)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[228],Timing[m,1]*Perts$Pest[20],Timing[m,2]*Perts$Pest[mt+20]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,22],Cov.Mat[1,42],
                                                                  0,PredStock$N.thou.SE[228]^2,0,0,
                                                                  Cov.Mat[22,1],0,Cov.Mat[22,22],Cov.Mat[22,42],
                                                                  Cov.Mat[42,1],0,Cov.Mat[42,22],Cov.Mat[42,42]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[241:252,1] <- ifelse(row(Timing)[241:252,1] >= Perts$tsteps[21],1,0)
            Timing[241:252,2] <- ifelse(row(Timing)[241:252,2] >= Perts$tsteps[mt+21],1,0)
            for(m in 241:252)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[240],Timing[m,1]*Perts$Pest[21],Timing[m,2]*Perts$Pest[mt+21]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,23],Cov.Mat[1,43],
                                                                  0,PredStock$N.thou.SE[240]^2,0,0,
                                                                  Cov.Mat[23,1],0,Cov.Mat[23,23],Cov.Mat[23,43],
                                                                  Cov.Mat[43,1],0,Cov.Mat[43,23],Cov.Mat[43,43]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[253:264,1] <- ifelse(row(Timing)[253:264,1] >= Perts$tsteps[22],1,0)
            Timing[253:264,2] <- ifelse(row(Timing)[253:264,2] >= Perts$tsteps[mt+22],1,0)
            for(m in 253:264)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[252],Timing[m,1]*Perts$Pest[22],Timing[m,2]*Perts$Pest[mt+22]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,24],Cov.Mat[1,44],
                                                                  0,PredStock$N.thou.SE[252]^2,0,0,
                                                                  Cov.Mat[24,1],0,Cov.Mat[24,24],Cov.Mat[24,44],
                                                                  Cov.Mat[44,1],0,Cov.Mat[44,24],Cov.Mat[44,44]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[265:276,1] <- ifelse(row(Timing)[265:276,1] >= Perts$tsteps[23],1,0)
            Timing[265:276,2] <- ifelse(row(Timing)[265:276,2] >= Perts$tsteps[mt+23],1,0)
            for(m in 265:276)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[264],Timing[m,1]*Perts$Pest[23],Timing[m,2]*Perts$Pest[mt+23]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,25],Cov.Mat[1,45],
                                                                  0,PredStock$N.thou.SE[264]^2,0,0,
                                                                  Cov.Mat[25,1],0,Cov.Mat[25,25],Cov.Mat[25,45],
                                                                  Cov.Mat[45,1],0,Cov.Mat[45,25],Cov.Mat[45,45]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[277:288,1] <- ifelse(row(Timing)[277:288,1] >= Perts$tsteps[24],1,0)
            Timing[277:288,2] <- ifelse(row(Timing)[277:288,2] >= Perts$tsteps[mt+24],1,0)
            for(m in 277:288)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[276],Timing[m,1]*Perts$Pest[24],Timing[m,2]*Perts$Pest[mt+24]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,26],Cov.Mat[1,46],
                                                                  0,PredStock$N.thou.SE[276]^2,0,0,
                                                                  Cov.Mat[26,1],0,Cov.Mat[26,26],Cov.Mat[26,46],
                                                                  Cov.Mat[46,1],0,Cov.Mat[46,26],Cov.Mat[46,46]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
            Timing[289:300,1] <- ifelse(row(Timing)[289:300,1] >= Perts$tsteps[25],1,0)
            Timing[289:300,2] <- ifelse(row(Timing)[289:300,2] >= Perts$tsteps[mt+25],1,0)
            for(m in 289:300)
              {
               PredStock$N.thou.SE[m] <- deltamethod(g=list(~x2*exp(-x1)+x3*exp(-x1)+x4*exp(-x1)),
                                                     mean=c(PopDyn$M,PredStock$N.thou[288],Timing[m,1]*Perts$Pest[25],Timing[m,2]*Perts$Pest[mt+25]),
                                                     cov=matrix(c(Cov.Mat[1,1],0,Cov.Mat[1,27],Cov.Mat[1,47],
                                                                  0,PredStock$N.thou.SE[288]^2,0,0,
                                                                  Cov.Mat[27,1],0,Cov.Mat[27,27],Cov.Mat[27,47],
                                                                  Cov.Mat[47,1],0,Cov.Mat[47,27],Cov.Mat[47,47]),4,4))
               PredStock$B.ton.SE[m]  <- sqrt((1e3*PredStock$N.thou.SE[m])^2*(PredStock$Mmw.kg[m]*1e-3)^2 +
                                              (1e3*PredStock$N.thou[m])^2*(PredStock$SDmw.kg[m]*1e-3)^2)
              }
           }
       #End of two fleets
       }
    #End of month
    }
  #time step day or week
  if(!multi)
    {
     if(length(unique(sapply(1:length(x), function(u) length(x[[u]]$Data$Properties$Fleets$Fleet)))) > 1)
       {stop("All catdyn objects in the list 'x' must be either 1-fleet or 2-fleets, not some 1-fleet and some 2-fleets")}
     #One fleet
     if(length(unique(as.vector(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet)))) == 1)
       {
        if(length(unique(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet))) > 1)
          {stop("All catdyn objects in the list 'x' must have the same name of fleet")}
        if(class(x) != "list" | length(x) < 15)
          {stop("For intra-annual models (time step is daily or weekly) x must be a list of objects of class 'catdyn' \n from succesful fit of models using CatDynFit, and the number of objects in the list (intra-annual fits) must be 15 consecutive years or more")}
        ny <- length(x)
        if(sum(dim(mbw.sd) != c(ny,3)) != 0)
          {stop(" 'mbw.sd' must be a data.frame with as many rows as the length of 'x' and three columns: \n year, mean weight, and standard deviation of mean weight")}
        if(length(method) != ny)
          {stop("One numerical method must be supplied for each catdyn object in 'x' ")}
        if(any(sapply(1:length(x), function(u) x[[u]]$Model[[method[u]]]$Type)>5))
          {stop("The maximum number of perturbations in any of the elements of 'x' must not be higher than 5")}
        PredStock <- data.frame(Year=as.numeric(format(as.Date(x[[1]]$Data$Properties$Dates[1]),"%Y")):as.numeric(format(as.Date(x[[ny]]$Data$Properties$Dates[1]),"%Y")),
                                Mw.kg=mbw.sd[,2],
                                SDmw.kg=mbw.sd[,3],
                                N0Tot.thou=0,
                                N0Tot.thou.SE=0,
                                B0Tot.ton=0,
                                B0Tot.ton.SE=0)
        for(i in 1:ny)
          {
           #Replacing missing standard errors
           if(any(is.na(x[[i]]$Model[[method[i]]]$bt.stdev)))
             {
              x[[i]]$Model[[method[i]]]$bt.stdev[which(is.na(x[[i]]$Model[[method[i]]]$bt.stdev))] <-
              x[[i]]$Model[[method[i]]]$bt.par[which(is.na(x[[i]]$Model[[method[i]]]$bt.stdev))]*
              mean(unlist(x[[i]]$Model[[method[i]]]$bt.stdev[which(!is.na(x[[i]]$Model[[method[i]]]$bt.stdev))])/
              unlist(x[[i]]$Model[[method[i]]]$bt.par[which(!is.na(x[[i]]$Model[[method[i]]]$bt.stdev))]))
             }
           mt             <- x[[i]]$Model[[method[i]]]$Type
           Thou.scaler    <- 1e6*(x[[i]]$Data$Properties$Units[4]=="bill") +
                             1e3*(x[[i]]$Data$Properties$Units[4]=="mill") +
                             1e0*(x[[i]]$Data$Properties$Units[4]=="thou") +
                             1e-1*(x[[i]]$Data$Properties$Units[4]=="hund")
           M     <- unlist(x[[i]]$Model[[method[i]]]$bt.par["M"])
           SE.M  <- unlist(x[[i]]$Model[[method[i]]]$bt.stdev["M"])
           N0    <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par["N0"])
           SE.N0 <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev["N0"])
           if(mt == 0)
             {
              PredStock$N0Tot.thou[i]    <- N0
              PredStock$N0Tor.thou.SE[i] <- SE.N0
              PredStock$B0Tot.ton[i]     <- N0*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]  <- sqrt((SE.N0)^2*(mbw.sd[i,2])^2+(N0)^2*(mbw.sd[i,3])^2)
             }
           if(mt == 1)
             {
              P1                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1.back                     <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[1:3,1:3],sd=c(SE.M,SE.N0,P1.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1*exp(M*P1.back)
              form                        <- sprintf("~x2 + x3*exp(x1*%i)",
                                                     P1.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
              }
           if(mt == 2)
             {
              P1                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1.back                     <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2.back                     <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[1:4,1:4],sd=c(SE.M,SE.N0,P1.SE,P2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1*exp(M*P1.back) +
                                                  P2*exp(M*P2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)",
                                                     P1.back,P2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1,P2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(mt == 3)
             {
              P1                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1.back                     <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2.back                     <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3.back                     <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[1:5,1:5],sd=c(SE.M,SE.N0,P1.SE,P2.SE,P3.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1*exp(M*P1.back) +
                                                  P2*exp(M*P2.back) +
                                                  P3*exp(M*P3.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)",
                                                     P1.back,P2.back,P3.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1,P2,P3),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(mt == 4)
             {
              P1                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1.back                     <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2.back                     <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3.back                     <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P4.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P4.back                     <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[1:6,1:6],sd=c(SE.M,SE.N0,P1.SE,P2.SE,P3.SE,P4.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1*exp(M*P1.back) +
                                                  P2*exp(M*P2.back) +
                                                  P3*exp(M*P3.back) +
                                                  P4*exp(M*P4.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)",
                                                     P1.back,P2.back,P3.back,P4.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1,P2,P3,P4),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(mt == 5)
             {
              P1                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1.back                     <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2.back                     <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3.back                     <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P4.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P4.back                     <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P5                          <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P5.SE                       <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P5.back                     <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[1:7,1:7],sd=c(SE.M,SE.N0,P1.SE,P2.SE,P3.SE,P4.SE,P5.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1*exp(M*P1.back) +
                                                  P2*exp(M*P2.back) +
                                                  P3*exp(M*P3.back) +
                                                  P4*exp(M*P4.back) +
                                                  P5*exp(M*P5.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)",
                                                     P1.back,P2.back,P3.back,P4.back,P5.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1,P2,P3,P4,P5),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
          }
      #End of one fleet
      }
    #Two fleets
    if(length(unique(as.vector(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet)))) == 2)
      {
        if(length(unique(as.vector(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet)))) != 2)
          {stop("All catdyn objects in the list 'x' must have the same two fleets")}
        if(class(x) != "list" | length(x) < 15)
          {stop("For intra-annual models (time step is daily or weekly) x must be a list of objects of class 'catdyn' \n from succesful fit of models using CatDynFit, and the number of objects in the list (intra-annual fits) must be 15 consecutive years or more")}
        ny <- length(x)
        if(sum(dim(mbw.sd) != c(ny,3)) != 0)
          {stop(" 'mbw.sd' must be a data.frame with as many rows as the length of 'x' and three columns: \n year, mean weight, and standard deviation of mean length ")}
        if(length(method) != ny)
          {stop("One numerical method must be supplied for each catdyn object in 'x' ")}
        if(max(as.vector(sapply(1:length(x), function(u) x[[u]]$Model[[method[u]]]$Type)))>5)
          {stop("The maximum number of perturbations in any of the elements and fleets of 'x' must not be higher than 5; \n any of the models in 'x' can go from a minimum of c(0,0) (pure depletion both fleets) to a maximum of c(5,5)")}
        PredStock <- data.frame(Year=as.numeric(format(as.Date(x[[1]]$Data$Properties$Dates[1]),"%Y")):as.numeric(format(as.Date(x[[ny]]$Data$Properties$Dates[1]),"%Y")),
                                Mw.kg=mbw.sd[,2],
                                SDmw.kg=mbw.sd[,3],
                                N0Tot.thou=0,
                                N0Tot.thou.SE=0,
                                B0Tot.ton=0,
                                B0Tot.ton.SE=0)
        for(i in 1:ny)
          {
           #Replacing missing standard errors
           if(any(is.na(x[[i]]$Model[[method[i]]]$bt.stdev)))
             {
              x[[i]]$Model[[method[i]]]$bt.stdev[which(is.na(x[[i]]$Model[[method[i]]]$bt.stdev))] <-
              x[[i]]$Model[[method[i]]]$bt.par[which(is.na(x[[i]]$Model[[method[i]]]$bt.stdev))]*
              mean(unlist(x[[i]]$Model[[method[i]]]$bt.stdev[which(!is.na(x[[i]]$Model[[method[i]]]$bt.stdev))])/
              unlist(x[[i]]$Model[[method[i]]]$bt.par[which(!is.na(x[[i]]$Model[[method[i]]]$bt.stdev))]))
             }
           mt             <- x[[i]]$Model[[method[i]]]$Type
           Thou.scaler    <- 1e6*(x[[i]]$Data$Properties$Units[4]=="bill") +
                             1e3*(x[[i]]$Data$Properties$Units[4]=="mill") +
                             1e0*(x[[i]]$Data$Properties$Units[4]=="thou") +
                             1e-1*(x[[i]]$Data$Properties$Units[4]=="hund")
           M     <- unlist(x[[i]]$Model[[method[i]]]$bt.par["M"])
           SE.M  <- unlist(x[[i]]$Model[[method[i]]]$bt.stdev["M"])
           N0    <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par["N0"])
           SE.N0 <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev["N0"])
           if(sum(mt == c(0,0))==2)
             {
              PredStock$N0Tot.thou[i]    <- N0
              PredStock$N0Tor.thou.SE[i] <- SE.N0
              PredStock$B0Tot.ton[i]     <- N0*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]  <- sqrt((SE.N0)^2*(mbw.sd[i,2])^2+(N0)^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(0,1))==2)
             {
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              ts1                         <- x[[i]]$Model[[method[i]]]$Dates[1]
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,6),c(1,2,6)],sd=c(SE.M,SE.N0,P1F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F2*exp(M*P1F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)",
                                                     P1F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(0,2))==2)
             {
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,6,7),c(1,2,6,7)],sd=c(SE.M,SE.N0,P1F2.SE,P2F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)",
                                                     P1F2.back,P2F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F2,P2F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(0,3))==2)
             {
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,6,7,8),c(1,2,6,7,8)],sd=c(SE.M,SE.N0,P1F2.SE,P2F2.SE,P3F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)",
                                                     P1F2.back,P2F2.back,P3F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F2,P2F2,P3F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(0,4))==2)
             {
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              ts1                         <- x[[i]]$Model[[method[i]]]$Dates[1]
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,6:9),c(1,2,6:9)],sd=c(SE.M,SE.N0,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)",
                                                     P1F2.back,P2F2.back,P3F2.back,P4F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F2,P2F2,P3F2,P4F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(0,5))==2)
             {
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P5F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P5F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P5F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,6:10),c(1,2,6:10)],sd=c(SE.M,SE.N0,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE,P5F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back) +
                                                  P5F2*exp(M*P5F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)",
                                                     P1F2.back,P2F2.back,P3F2.back,P4F2.back,P5F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F2,P2F2,P3F2,P4F2,P5F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(1,1))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,7),c(1,2,3,7)],sd=c(SE.M,SE.N0,P1F1.SE,P1F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P1F2*exp(M*P1F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)",
                                                     P1F1.back,P1F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P1F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(1,2))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                     <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P1F2.back                     <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P2F2.back                     <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,7,8),c(1,2,3,7,8)],sd=c(SE.M,SE.N0,P1F1.SE,P1F2.SE,P2F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)",
                                                     P1F1.back,P1F2.back,P2F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P1F2,P2F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(1,3))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,7:9),c(1,2,3,7:9)],sd=c(SE.M,SE.N0,P1F1.SE,P1F2.SE,P2F2.SE,P3F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back)
              form                       <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)",
                                                     P1F1.back,P1F2.back,P2F2.back,P3F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P1F2,P2F2,P3F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(1,4))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,7:10),c(1,2,3,7:10)],sd=c(SE.M,SE.N0,P1F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)",
                                                     P1F1.back,P1F2.back,P2F2.back,P3F2.back,P4F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P1F2,P2F2,P3F2,P4F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(1,5))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P5F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P5F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P5F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[7]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              ts1                         <-
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,7:11),c(1,2,3,7:11)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE,P5F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back) +
                                                  P5F2*exp(M*P5F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)+x8*exp(x1*%i)",
                                                     P1F1.back,P1F2.back,P2F2.back,P3F2.back,P4F2.back,P5F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P1F2,P2F2,P3F2,P4F2,P5F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(2,2))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,4,8,9),c(1,2,3,4,8,9)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P1F2.SE,P2F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P2F1*exp(M*P1F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)",
                                                     P1F1.back,P2F1.back,P1F2.back,P2F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P2F1,P1F2,P2F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(2,3))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,4,8:10),c(1,2,3,4,8:10)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P1F2.SE,P2F2.SE,P3F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P2F1*exp(M*P2F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)",
                                                     P1F1.back,P2F1.back,P1F2.back,P2F2.back,P3F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P2F1,P1F2,P2F2,P3F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(2,4))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[7]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,4,8:11),c(1,2,3,4,8:11)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P2F1*exp(M*P2F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)+x8*exp(x1*%i)",
                                                     P1F1.back,P2F1.back,P1F2.back,P2F2.back,P3F2.back,P4F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P2F1,P1F2,P2F2,P3F2,P4F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(2,5))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[8])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[8])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[7]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P5F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[12])
              P5F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[12])
              P5F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[8]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1,2,3,4,8:12),c(1,2,3,4,8:12)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE,P5F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P2F1*exp(M*P2F2.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back) +
                                                  P5F2*exp(M*P5F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)+x8*exp(x1*%i)+x9*exp(x1*%i)",
                                                     P1F1.back,P2F1.back,P1F2.back,P2F2.back,P3F2.back,P4F2.back,P5F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P2F1,P1F2,P2F2,P3F2,P4F2,P5F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(3,3))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[7]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1:5,9:11),c(1:5,9:11)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P3F1.SE,P1F2.SE,P2F2.SE,P3F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P2F1*exp(M*P2F1.back) +
                                                  P3F1*exp(M*P3F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)+x8*exp(x1*%i)",
                                                     P1F1.back,P2F1.back,P3F1.back,P1F2.back,P2F2.back,P3F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P2F1,P3F1,P1F2,P2F2,P3F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(3,4))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[7]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[12])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[12])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[8]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1:5,9:12),c(1:5,9:12)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P3F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P2F1*exp(M*P2F1.back) +
                                                  P3F1*exp(M*P3F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)+x8*exp(x1*%i)+x9*exp(x1*%i)",
                                                     P1F1.back,P2F1.back,P3F1.back,P1F2.back,P2F2.back,P3F2.back,P4F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P2F1,P3F1,P1F2,P2F2,P3F2,P4F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(3,5))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[2]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[3]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3F1.back                   <- x[[i]]$Model[[method[i]]]$Dates[4]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[9])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[9])
              P1F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[5]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P2F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[6]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P3F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[7]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[12])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[12])
              P4F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[8]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              P5F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[13])
              P5F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[13])
              P5F2.back                   <- x[[i]]$Model[[method[i]]]$Dates[9]-x[[i]]$Model[[method[i]]]$Dates[1]+1
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1:5,9:13),c(1:5,9:13)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P3F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE,P5F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*P1F1.back) +
                                                  P2F1*exp(M*P2F1.back) +
                                                  P3F1*exp(M*P3F1.back) +
                                                  P1F2*exp(M*P1F2.back) +
                                                  P2F2*exp(M*P2F2.back) +
                                                  P3F2*exp(M*P3F2.back) +
                                                  P4F2*exp(M*P4F2.back) +
                                                  P5F2*exp(M*P5F2.back)
              form                        <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)+x7*exp(x1*%i)+x8*exp(x1*%i)+x9*exp(x1*%i)+x10*exp(x1*%i)",
                                                     P1F1.back,P2F1.back,P3F1.back,P1F2.back,P2F2.back,P3F2.back,P4F2.back,P5F2.back)
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=as.formula(form),
                                                         mean=c(M,N0,P1F1,P2F1,P3F1,P1F2,P2F2,P3F2,P4F2,P5F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(4,4))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[2]
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[3]
              P3F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[4]
              P4F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P4F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P4F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[5]
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P1F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[6]
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P2F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[7]
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[12])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[12])
              P3F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[8]
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[13])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[13])
              P4F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[9]
              ts1                         <- x[[i]]$Model[[method[i]]]$Dates[1]
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1:6,10:13),c(1:6,10:13)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P3F1.SE,P4F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*(P1F1.ts-ts1+1)) + P2F1*exp(M*(P2F1.ts-ts1+1)) + P3F1*exp(M*(P3F1.ts-ts1+1)) + P4F1*exp(M*(P4F1.ts-ts1+1)) +
                                                  P1F2*exp(M*(P1F2.ts-ts1+1)) + P2F2*exp(M*(P2F2.ts-ts1+1)) + P3F2*exp(M*(P3F2.ts-ts1+1)) + P4F2*exp(M*(P4F2.ts-ts1+1))
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=list(~x2 + x3*exp(x1*(P1F1.ts-ts1+1)) + x4*exp(x1*(P2F1.ts-ts1+1)) + x5*exp(x1*(P3F1.ts-ts1+1)) + x6*exp(x1*(P4F1.ts-ts1+1)) +
                                                                      x7*exp(x1*(P1F2.ts-ts1+1)) + x7*exp(x1*(P2F2.ts-ts1+1)) + x9*exp(x1*(P3F2.ts-ts1+1)) + x10*exp(x1*(P4F2.ts-ts1+1))),
                                                         mean=c(M,N0,P1F1,P2F1,P3F1,P4F1,P1F2,P2F2,P3F2,P4F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(4,5))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[2]
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[3]
              P3F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[4]
              P4F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P4F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P4F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[5]
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[10])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[10])
              P1F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[6]
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P2F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[7]
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[12])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[12])
              P3F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[8]
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[13])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[13])
              P4F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[9]
              P5F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[14])
              P5F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[14])
              P5F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[10]
              ts1                         <- x[[i]]$Model[[method[i]]]$Dates[1]
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1:6,10:14),c(1:6,10:14)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P3F1.SE,P4F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE,P5F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*(P1F1.ts-ts1+1)) + P2F1*exp(M*(P2F1.ts-ts1+1)) + P3F1*exp(M*(P3F1.ts-ts1+1)) + P4F1*exp(M*(P4F1.ts-ts1+1)) +
                                                  P1F2*exp(M*(P1F2.ts-ts1+1)) + P2F2*exp(M*(P2F2.ts-ts1+1)) + P3F2*exp(M*(P3F2.ts-ts1+1)) + P4F2*exp(M*(P4F2.ts-ts1+1)) + P5F2*exp(M*(P5F2.ts-ts1+1))
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=list(~x2 + x3*exp(x1*(P1F1.ts-ts1+1)) + x4*exp(x1*(P2F1.ts-ts1+1)) + x5*exp(x1*(P3F1.ts-ts1+1)) + x6*exp(x1*(P4F1.ts-ts1+1)) +
                                                                      x7*exp(x1*(P1F2.ts-ts1+1)) + x7*exp(x1*(P2F2.ts-ts1+1)) + x9*exp(x1*(P3F2.ts-ts1+1)) + x10*exp(x1*(P4F2.ts-ts1+1)) + x11*exp(x1*(P5F2.ts-ts1+1))),
                                                         mean=c(M,N0,P1F1,P2F1,P3F1,P4F1,P1F2,P2F2,P3F2,P4F2,P5F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
           if(sum(mt == c(5,5))==2)
             {
              P1F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
              P1F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
              P1F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[2]
              P2F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
              P2F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
              P2F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[3]
              P3F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
              P3F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
              P3F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[4]
              P4F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
              P4F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
              P4F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[5]
              P5F1                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[7])
              P5F1.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[7])
              P5F1.ts                     <- x[[i]]$Model[[method[i]]]$Dates[6]
              P1F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[11])
              P1F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[11])
              P1F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[7]
              P2F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[12])
              P2F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[12])
              P2F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[8]
              P3F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[13])
              P3F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[13])
              P3F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[9]
              P4F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[14])
              P4F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[14])
              P4F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[10]
              P5F2                        <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.par[15])
              P5F2.SE                     <- Thou.scaler*unlist(x[[i]]$Model[[method[i]]]$bt.stdev[15])
              P5F2.ts                     <- x[[i]]$Model[[method[i]]]$Dates[11]
              ts1                         <- x[[i]]$Model[[method[i]]]$Dates[1]
              cov                         <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[c(1:7,11:15),c(1:7,11:15)],
                                                     sd=c(SE.M,SE.N0,P1F1.SE,P2F1.SE,P3F1.SE,P4F1.SE,P5F1.SE,P1F2.SE,P2F2.SE,P3F2.SE,P4F2.SE,P5F2.SE))
              PredStock$N0Tot.thou[i]     <- N0 + P1F1*exp(M*(P1F1.ts-ts1+1)) + P2F1*exp(M*(P2F1.ts-ts1+1)) + P3F1*exp(M*(P3F1.ts-ts1+1)) + P4F1*exp(M*(P4F1.ts-ts1+1)) + P5F1*exp(M*(P5F1.ts-ts1+1)) +
                                                  P1F2*exp(M*(P1F2.ts-ts1+1)) + P2F2*exp(M*(P2F2.ts-ts1+1)) + P3F2*exp(M*(P3F2.ts-ts1+1)) + P4F2*exp(M*(P4F2.ts-ts1+1)) + P5F2*exp(M*(P5F2.ts-ts1+1))
              PredStock$N0Tot.thou.SE[i]  <- deltamethod(g=list(~x2 + x3*exp(x1*(P1F1.ts-ts1+1)) + x4*exp(x1*(P2F1.ts-ts1+1)) + x5*exp(x1*(P3F1.ts-ts1+1)) + x6*exp(x1*(P4F1.ts-ts1+1)) + x7*exp(x1*(P5F1.ts-ts1+1)) +
                                                                      x8*exp(x1*(P1F2.ts-ts1+1)) + x9*exp(x1*(P2F2.ts-ts1+1)) + x10*exp(x1*(P3F2.ts-ts1+1)) + x11*exp(x1*(P4F2.ts-ts1+1)) + x12*exp(x1*(P5F2.ts-ts1+1))),
                                                         mean=c(M,N0,P1F1,P2F1,P3F1,P4F1,P1F2,P2F2,P3F2,P4F2,P5F2),
                                                         cov=cov)
              PredStock$B0Tot.ton[i]      <- PredStock$N0Tot.thou[i]*mbw.sd[i,2]
              PredStock$B0Tot.ton.SE[i]   <- sqrt((PredStock$N0Tot.thou.SE[i])^2*(mbw.sd[i,2])^2+(PredStock$N0Tot.thou[i])^2*(mbw.sd[i,3])^2)
             }
          }
      #End of two fleets
      }
    #End week or day
    }
    return(PredStock)
  #End of function
  }
