CatDynPar <-
function(x,method,partial=TRUE)
   {
   if(class(x) != "catdyn")
     {stop("'x' must be an object of class 'catdyn' (created by function CatDynFit)")}
    sdistr.set <- c("poisson","apnormal","aplnormal")
    if(x$Data$Properties$Units["Time Step"]=="month")
      {
       #Month 1 Fleet
       if(length(x$Data$Properties$Fleets$Fleet)==1)
         {
          month.F1 <- x$Model[[method]]$Dates[2:(x$Model[[method]]$Type+1)]
          years.F1 <- as.numeric(format(as.Date(x$Data$Properties$Dates["StartDate"]),"%Y"))+floor(month.F1/12)
          month.F1 <- as.integer((month.F1/12-floor(month.F1/12))*12)
          #Month 1 Fleet Likelihood without dispersion
          if(x$Model[[method]]$Distr%in%sdistr.set)
            {
             partable <- data.frame(Parameter=c("M.1/month",
                                                paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".",years.F1,sep=""),
                                                paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                "alpha","beta"),
                                    Timing=c("","",paste(years.F1,month.F1,sep="-"),"","",""),
                                    Estimates=unlist(x$Model[[method]]$bt.par),
                                    CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
            }
          #Month 1 Fleet Likelihood with dispersion
          if(!x$Model[[method]]$Distr%in%sdistr.set)
            {
             partable <- data.frame(Parameter=c("M.1/month",
                                                paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".",years.F1,sep=""),
                                                paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                "alpha","beta",
                                                paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                    Timing=c("","",paste(years.F1,month.F1,sep="-"),"","","",""),
                                    Estimates=unlist(x$Model[[method]]$bt.par),
                                    CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
            }
         #end of Month 1 Fleet
         }
       #Month 2 Fleets
       if(length(x$Data$Properties$Fleets$Fleet)==2)
         {
          month.F1 <- x$Model[[method]]$Dates[2:(x$Model[[method]]$Type[1]+1)]
          month.F2 <- x$Model[[method]]$Dates[(x$Model[[method]]$Type[1]+2):(x$Model[[method]]$Type[1]+1+x$Model[[method]]$Type[2])]
          years.F1 <- as.numeric(format(as.Date(x$Data$Properties$Dates["StartDate"]),"%Y"))+floor(month.F1/12)
          years.F2 <- as.numeric(format(as.Date(x$Data$Properties$Dates["StartDate"]),"%Y"))+floor(month.F2/12)
          month.F1 <- as.integer((month.F1/12-floor(month.F1/12))*12)
          month.F2 <- as.integer((month.F2/12-floor(month.F2/12))*12)
          if(any(month.F1==0)){month.F1[which(month.F1==0)] <- 12}
          if(any(month.F2==0)){month.F2[which(month.F2==0)] <- 12}
          #Month 2 Fleets Both fleets with likelihood without dispersion (sdristr.set)
          if(sum(x$Model[[method]]$Distr%in%sdistr.set)==2)
            {
             partable <- data.frame(Parameter=c("M.1/month",
                                                paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".",years.F1,sep=""),
                                                paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                "alpha","beta"),
                                    Timing.F1=c("","",paste(years.F1,month.F1,sep="-"),"","",""),
                                    Estimates.F1=unlist(x$Model[[method]]$bt.par[1:(2+length(years.F1)+3)]),
                                    CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[1:(2+length(years.F1)+3)])/
                                                         unlist(x$Model[[method]]$bt.par[1:(2+length(years.F1)+3)]),1),
                                    Timing.F2=c("","",paste(years.F2,month.F2,sep="-"),"","",""),
                                    Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,(2+length(years.F2)+3+1):length(x$Model[[method]]$bt.par))]),
                                    CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+length(years.F2)+4):(length(x$Model[[method]]$bt.stdev)))])/
                                                         unlist(x$Model[[method]]$bt.par[(c(1,2,(2+length(years.F2)+4):length(x$Model[[method]]$bt.par)))]),1),row.names=NULL)
            names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
            names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
            }
          #Month 2 Fleets First fleet with likelihood without dispersion
          if(diff(x$Model[[method]]$Distr%in%sdistr.set)==-1)
            {
             partable <- data.frame(Parameter=c("M.1/month",
                                                paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".",years.F1,sep=""),
                                                paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                "alpha","beta",
                                                paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                    Timing.F1=c("","",paste(years.F1,month.F1,sep="-"),"","","",""),
                                    Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1:(2+length(years.F1)+3))]),NA),
                                    CVpCent.F1=round(c(100*unlist(x$Model[[method]]$bt.stdev[c(1:(2+length(years.F1)+3))])/
                                                           unlist(x$Model[[method]]$bt.par[c(1:(2+length(years.F1)+3))]),NA),1),
                                    Timing.F2=c("","",paste(years.F2,month.F2,sep="-"),"","","",""),
                                    Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,(2+length(years.F2)+3+1):length(x$Model[[method]]$bt.par))]),
                                    CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+length(years.F2)+3+1):length(x$Model[[method]]$bt.stdev))])/
                                                         unlist(x$Model[[method]]$bt.par[c(1,2,(2+length(years.F2)+3+1):length(x$Model[[method]]$bt.par))]),1),row.names=NULL)
            names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
            names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
            }
          #Month 2 Fleets Second fleet with likelihood without dispersion
          if(diff(x$Model[[method]]$Distr%in%sdistr.set)==1)
            {
             partable <- data.frame(Parameter=c("M.1/month",
                                                paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".",years.F1,sep=""),
                                                paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                "alpha","beta",
                                                paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                    Timing.F1=c("","",paste(years.F1,month.F1,sep="-"),"","","",""),
                                    Estimates.F1=unlist(x$Model[[method]]$bt.par[c(1:(2+length(years.F1)+3),length(unlist(x$Model[[method]]$bt.par)))]),
                                    CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[c(1:(2+length(years.F1)+3),length(unlist(x$Model[[method]]$bt.stdev)))])/
                                                         unlist(x$Model[[method]]$bt.par[c(1:(2+length(years.F1)+3),length(unlist(x$Model[[method]]$bt.par)))]),1),
                                    Timing.F2=c("","",paste(years.F2,month.F2,sep="-"),"","","",""),
                                    Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+length(years.F2)+3+1):(length(x$Model[[method]]$bt.par)-1))]),NA),
                                    CVpCent.F2=c(round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+length(years.F2)+3+1):(length(x$Model[[method]]$bt.stdev)-1))])/
                                                           unlist(x$Model[[method]]$bt.par[c(1,2,(2+length(years.F2)+3+1):(length(x$Model[[method]]$bt.par)-1))]),1),NA),row.names=NULL)
            names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
            names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
            }
          #Month 2 Fleets Both fleets with likelihood with dispersion
          if(sum(x$Model[[method]]$Distr%in%sdistr.set)==0)
            {
             partable <- data.frame(Parameter=c("M.1/month",
                                                paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".",years.F1,sep=""),
                                                paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                "alpha","beta",
                                                paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                    Timing.F1=c("","",paste(years.F1,month.F1,sep="-"),"","","",""),
                                    Estimates.F1=unlist(x$Model[[method]]$bt.par[c(1:(2+length(years.F1)+3),length(unlist(x$Model[[method]]$bt.par))-1)]),
                                    CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[c(1:(2+length(years.F1)+3),length(unlist(x$Model[[method]]$bt.stdev))-1)])/
                                                         unlist(x$Model[[method]]$bt.par[c(1:(2+length(years.F1)+3),length(unlist(x$Model[[method]]$bt.par))-1)]),1),
                                    Timing.F2=c("","",paste(years.F2,month.F2,sep="-"),"","","",""),
                                    Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,(2+length(years.F2)+3+1):(length(x$Model[[method]]$bt.par)-2),(length(x$Model[[method]]$bt.par)))]),
                                    CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+length(years.F2)+3+1):(length(x$Model[[method]]$bt.stdev)-2),(length(x$Model[[method]]$bt.stdev)))])/
                                                         unlist(x$Model[[method]]$bt.par[c(1,2,(2+length(years.F2)+3+1):(length(x$Model[[method]]$bt.par)-2),(length(x$Model[[method]]$bt.par)))]),1),row.names=NULL)
            names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
            names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
            }
         #End of Month 2 Fleets
         }
      #End of month
      }
    if(x$Data$Properties$Units["Time Step"]=="day" | x$Data$Properties$Units["Time Step"]=="week")
      {
       #Week/Day 1 Fleet
       if(length(x$Data$Properties$Fleets$Fleet)==1)
         {
          #Week/Day 1 fleet Pure depletion
          if(x$Model[[method]]$Type==0)
            {
             #Week/Day 1 fleet Pure depletion Likelihood without dispersion
             if(x$Model[[method]]$Distr%in%sdistr.set)
               {
                if(x$Data$Properties$Units["Time Step"]=="week")
                  {
                   partable <- data.frame(Parameter=c("M.1/week",
                                                      paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                      paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                      "alpha","beta"),
                                          Estimates=unlist(x$Model[[method]]$bt.par),
                                          CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                  }
                if(x$Data$Properties$Units["Time Step"]=="day")
                  {
                   partable <- data.frame(Parameter=c("M.1/day",
                                                      paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                      paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                      "alpha","beta"),
                                          Estimates=unlist(x$Model[[method]]$bt.par),
                                          CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                  }
               }
             #Week/Day 1 fleet Pure depletion Likelihood with dispersion
             if(!x$Model[[method]]$Distr%in%sdistr.set)
               {
                if(x$Data$Properties$Units["Time Step"]=="week")
                  {
                   partable <- data.frame(Parameter=c("M.1/week",
                                                      paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                      paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                      "alpha","beta",
                                                      paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                          Estimates=unlist(x$Model[[method]]$bt.par),
                                          CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                  }
                if(x$Data$Properties$Units["Time Step"]=="day")
                  {
                   partable <- data.frame(Parameter=c("M.1/day",
                                                      paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                      paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                      "alpha","beta",
                                                      paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                          Estimates=unlist(x$Model[[method]]$bt.par),
                                          CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                  }
               }
            #End of Week/Day 1 fleet Pure depletion
            }
          if(x$Model[[method]]$Type!=0)
            #Week/Day 1 fleet With perturbations
            {
             #Week/Day 1 fleet With perturbations Likelihood without dispersion
             waves <- abs(x$Model[[method]]$Type)
             if(x$Model[[method]]$Distr %in% sdistr.set)
               {
                if(x$Data$Properties$Units["Time Step"]=="week")
                  {
                   year1           <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   finalweek.year1 <- as.numeric(format(as.Date(paste(year1,"-12-31",sep="")), "%W"))
                   year2           <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   #weeks    <- x$Model[[method]]$Dates[2:(waves+1)]
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   #dates.ranges <- vector("character",waves)
                   if(partial & x$Model[[method]]$Type < 0)
                     {
                      dates.ranges <- vector("character",2*waves)
                      if(year2 > year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(2*waves+1)] - finalweek.year1
                         for(w in 1:2*waves)
                           {
                            dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year2, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      if(year2 == year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(2*waves+1)]
                         for(w in 1:(2*waves))
                          {
                           dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year1, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                          }
                       }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Recruitment.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("Spawning.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta"),
                                             Timing=c("","",dates.ranges,"","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(partial & x$Model[[method]]$Type > 0)
                     {
                      dates.ranges <- vector("character",waves)
                      if(year2 > year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(2*waves+1)] - finalweek.year1
                         for(w in 1:waves)
                           {
                            dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year2, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      if(year2 == year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(waves+1)]
                         for(w in 1:(waves))
                          {
                           dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year1, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                          }
                       }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Recruitment.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta"),
                                             Timing=c("","",dates.ranges,"","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(!partial)
                     {
                      weeks        <- x$Model[[method]]$Dates[2:(waves+1)]
                      dates.ranges <- vector("character",waves)
                      for(w in 1:waves)
                        {
                         dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year1, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                        }
                         partable <- data.frame(Parameter=c("M.1/week",
                                                            paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                            paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                            paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                            "alpha","beta"),
                                                Timing=c("","",dates.ranges,"","",""),
                                                Estimates=unlist(x$Model[[method]]$bt.par),
                                                CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                  }
                if(x$Data$Properties$Units["Time Step"]=="day")
                  {
                   if(partial & x$Model[[method]]$Type < 0)
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("Spawning.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta"),
                                             Timing=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves+2)]),
                                                      "","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(partial & x$Model[[method]]$Type > 0)
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta"),
                                             Timing=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves+1)]),
                                                      "","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(!partial)
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta"),
                                             Timing=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves+1)]),"","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                 }
               }
             #Week/Day 1 fleet With perturbations Likelihood with dispersion
             if(!x$Model[[method]]$Distr %in% sdistr.set)
               {
                if(x$Data$Properties$Units["Time Step"]=="week")
                  {
                   year1           <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   finalweek.year1 <- as.numeric(format(as.Date(paste(year1,"-12-31",sep="")), "%W"))
                   year2           <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days            <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(partial & x$Model[[method]]$Type < 0)
                     {
                      dates.ranges <- vector("character",2*waves)
                      if(year2 > year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(2*waves+1)] - finalweek.year1
                         for(w in 1:2*waves)
                           {
                            dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year2, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      if(year2 == year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(2*waves+1)]
                         for(w in 1:2*waves)
                          {
                           dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year1, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                          }
                       }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Recruitment.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("Spawning.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing=c("","",dates.ranges,"","","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(partial & x$Model[[method]]$Type > 0)
                     {
                      dates.ranges <- vector("character",waves)
                      if(year2 > year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(waves+1)] - finalweek.year1
                         for(w in 1:waves)
                           {
                            dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year2, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      if(year2 == year1)
                        {
                         weeks <- x$Model[[method]]$Dates[2:(waves+1)]
                         for(w in 1:waves)
                          {
                           dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year1, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                          }
                       }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Recruitment.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing=c("","",dates.ranges,"","","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(!partial)
                     {
                      weeks        <- x$Model[[method]]$Dates[2:(waves+1)]
                      dates.ranges <- vector("character",waves)
                      for(w in 1:waves)
                        {
                         dates.ranges[w] <- paste(range(days[sprintf("%d %02d", year1, weeks[w]) == format(days, "%Y %U")]),collapse=' ')
                        }
                         partable <- data.frame(Parameter=c("M.1/week",
                                                            paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                            paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                            paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                            "alpha","beta",
                                                            paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                                Timing=c("","",dates.ranges,"","","",""),
                                                Estimates=unlist(x$Model[[method]]$bt.par),
                                                CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                  }
                if(x$Data$Properties$Units["Time Step"]=="day")
                  {
                   if(partial & x$Model[[method]]$Type < 0)
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("Spawning.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves+2)]),
                                                      "","","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(partial & x$Model[[method]]$Type > 0)
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves+1)]),
                                                      "","","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                   if(!partial)
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units,sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves+1)]),"","","",""),
                                             Estimates=unlist(x$Model[[method]]$bt.par),
                                             CVpCent=round(100*unlist(x$Model[[method]]$bt.stdev)/unlist(x$Model[[method]]$bt.par),1),row.names=NULL)
                     }
                  }
               }
            #End of Week/Day 1 fleet With perturbations
            }
         #End of Week/Day one fleet
         }
       #Week/Day 2 Fleets
       if(length(x$Data$Properties$Fleets$Fleet)==2)
         {
          #Week/Day 2 Fleets Pure depletion
          if(sum(x$Model[[method]]$Type)==0)
            {
            #Week/Day 2 Fleets Pure depletion Both fleets with likelihood without dispersion
            if(sum(x$Model[[method]]$Distr%in%sdistr.set)==2)
              {
               if(x$Data$Properties$Units["Time Step"]=="week")
                 {
                  partable <- data.frame(Parameter=c("M.1/week",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta"),
                                         Estimates.F1=unlist(x$Model[[method]]$bt.par[1:5]),
                                         CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[1:5])/
                                                              unlist(x$Model[[method]]$bt.par[1:5]),1),
                                         Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),
                                         CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:8)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),1),row.names=NULL)
                 names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                 names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                 }
               if(x$Data$Properties$Units["Time Step"]=="day")
                 {
                  partable <- data.frame(Parameter=c("M.1/day",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta"),
                                         Estimates.F1=unlist(x$Model[[method]]$bt.par[1:5]),
                                         CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[1:5])/
                                                              unlist(x$Model[[method]]$bt.par[1:5]),1),
                                         Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),
                                         CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:8)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),1),row.names=NULL)
                  names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                  names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                  }
               }
            #Week/Day 2 Fleets Pure depletion First fleet with likelihood without dispersion
            if(diff(x$Model[[method]]$Distr%in%sdistr.set)==-1)
              {
               if(x$Data$Properties$Units["Time Step"]=="week")
                 {
                  partable <- data.frame(Parameter=c("M.1/week",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta",
                                                     paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                         Estimates.F1=c(unlist(x$Model[[method]]$bt.par[1:5]),NA),
                                         CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[1:5]),NA)/
                                                              c(unlist(x$Model[[method]]$bt.par[1:5]),NA),1),
                                         Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:9)]),
                                         CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:9)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1,2,6:9)]),1),row.names=NULL)
                 names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                 names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                 }
               if(x$Data$Properties$Units["Time Step"]=="day")
                 {
                  partable <- data.frame(Parameter=c("M.1/day",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta",
                                                     paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                         Estimates.F1=c(unlist(x$Model[[method]]$bt.par[1:5]),NA),
                                         CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[1:5]),NA)/
                                                              c(unlist(x$Model[[method]]$bt.par[1:5]),NA),1),
                                         Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:9)]),
                                         CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:9)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1,2,6:9)]),1),row.names=NULL)
                  names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                  names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                  }
              }
            #Week/Day 2 Fleets Pure depletion Second fleet with likelihood without dispersion
            if(diff(x$Model[[method]]$Distr%in%sdistr.set)==1)
              {
               if(x$Data$Properties$Units["Time Step"]=="week")
                 {
                  partable <- data.frame(Parameter=c("M.1/week",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta",
                                                     paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                         Estimates.F1=unlist(x$Model[[method]]$bt.par[c(1:5,9)]),
                                         CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[c(1:5,9)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1:5,9)]),1),
                                         Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),NA),
                                         CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:8)]),NA)/
                                                              c(unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),NA),1),row.names=NULL)
                 names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                 names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                 }
               if(x$Data$Properties$Units["Time Step"]=="day")
                 {
                  partable <- data.frame(Parameter=c("M.1/day",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta",
                                                     paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                         Estimates.F1=unlist(x$Model[[method]]$bt.par[c(1:5,9)]),
                                         CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[c(1:5,9)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1:5,9)]),1),
                                         Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),NA),
                                         CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:8)]),NA)/
                                                              c(unlist(x$Model[[method]]$bt.par[c(1,2,6:8)]),NA),1),row.names=NULL)
                  names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                  names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                  }
              }
            #Week/Day 2 Fleets Pure depletion Both fleets with likelihood with dispersion
            if(sum(x$Model[[method]]$Distr%in%sdistr.set)==0)
              {
               if(x$Data$Properties$Units["Time Step"]=="week")
                 {
                  partable <- data.frame(Parameter=c("M.1/week",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta",
                                                     paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                         Estimates.F1=unlist(x$Model[[method]]$bt.par[c(1:5,9)]),
                                         CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[c(1:5,9)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1:5,9)]),1),
                                         Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:8,10)]),
                                         CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:8,10)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1,2,6:8,10)]),1),row.names=NULL)
                 names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                 names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                 }
               if(x$Data$Properties$Units["Time Step"]=="day")
                 {
                  partable <- data.frame(Parameter=c("M.1/day",
                                                     paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                     paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                     "alpha","beta",
                                                     paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                         Estimates.F1=unlist(x$Model[[method]]$bt.par[c(1:5,9)]),
                                         CVpCent.F1=round(100*unlist(x$Model[[method]]$bt.stdev[c(1:5,9)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1:5,9)]),1),
                                         Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:8,10)]),
                                         CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:8,10)])/
                                                              unlist(x$Model[[method]]$bt.par[c(1,2,6:8,10)]),1),row.names=NULL)
                  names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                  names(partable)[4:5] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                  }
             }
            #End of Week/Day 2 Fleets Pure depletion
            }
          if(sum(x$Model[[method]]$Type)!=0)
            #Week/Day 2 Fleets With perturbations
            {
             year     <- as.numeric(substring(x$Data$Properties$Dates[["StartDate"]],1,4))
             waves.F1 <- x$Model[[method]]$Type[1]
             waves.F2 <- x$Model[[method]]$Type[2]
             #Week/Day 2 Fleets With perturbations Both fleets with likelihood without dispersion (sdristr.set)
             if(sum(x$Model[[method]]$Distr%in%sdistr.set)==2)
               {
                #Only second fleet with perturbations
                if(waves.F1==0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks2    <- x$Model[[method]]$Dates[2:(waves.F2+1)]
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta"),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5])),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5])),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","",""),
                                             Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),
                                             CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))])/
                                                                  unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta"),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5])),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5])),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F2+1)]),"","",""),
                                             Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),
                                             CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))])/
                                                                  unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
                #Both fleets with perturbations
                if(waves.F1!=0 & waves.F2!=0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks1    <- x$Model[[method]]$Dates[2:(waves.F1+1)]
                      weeks2    <- x$Model[[method]]$Dates[(waves.F1+2):(waves.F1+waves.F2+1)]
                      dates.ranges.F1 <- vector("character",waves.F1)
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F1)
                        {
                         if(weeks1[w] <= 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year1, weeks1[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks1[w] > 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year2, weeks1[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta"),
                                             Timing.F1=c("","",dates.ranges.F1,rep("",3+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)])),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)])),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","",""),
                                             Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+4):(2+waves.F1+3+waves.F2+3))]),
                                             CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+waves.F1+4):(2+waves.F1+3+waves.F2+3))])/
                                                                  unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+4):(2+waves.F1+3+waves.F2+3))]),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta"),
                                             Timing.F1=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F1+1)]),rep("",3+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)])),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)])),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[(waves.F1+2):(waves.F1+1+waves.F2)]),"","",""),
                                             Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),
                                             CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(waves.F1+6):(waves.F1+5+waves.F2+3))])/
                                                                  unlist(x$Model[[method]]$bt.par[c(1,2,(waves.F1+6):(waves.F1+5+waves.F2+3))]),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
               }
             #First fleet with likelihood without dispersion
             if(diff(x$Model[[method]]$Distr%in%sdistr.set)==-1)
               {
                #Only second fleet with perturbations
                if(waves.F1==0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks2   <- x$Model[[method]]$Dates[2:(waves.F2+1)]
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),NA),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),NA),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),NA),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),NA),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F2+1)]),"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
                #Both fleets with perturbations
                if(waves.F1!=0 & waves.F2!=0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks1    <- x$Model[[method]]$Dates[2:(waves.F1+1)]
                      weeks2    <- x$Model[[method]]$Dates[(waves.F1+2):(waves.F1+waves.F2+1)]
                      dates.ranges.F1 <- vector("character",waves.F1)
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F1)
                        {
                         if(weeks1[w] <= 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year1, weeks1[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks1[w] > 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year2, weeks1[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing.F1=c("","",dates.ranges.F1,rep("",4+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),NA),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),NA),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","","",""),
                                             Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+4):(2+waves.F1+3+waves.F2+4))]),
                                             CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+waves.F1+4):(2+waves.F1+3+waves.F2+4))])/
                                                                  unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+4):(2+waves.F1+3+waves.F2+4))]),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing.F1=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F1+1)]),rep("",4+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),NA),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),NA),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[(waves.F1+2):(waves.F1+1+waves.F2)]),"","","",""),
                                             Estimates.F2=unlist(x$Model[[method]]$bt.par[c(1,2,(waves.F1+6):((waves.F1+6+waves.F2+3)))]),
                                             CVpCent.F2=round(100*unlist(x$Model[[method]]$bt.stdev[c(1,2,(waves.F1+6):((waves.F1+6+waves.F2+3)))])/
                                                                  unlist(x$Model[[method]]$bt.par[c(1,2,(waves.F1+6):((waves.F1+6+waves.F2+3)))]),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
               }
             #Second fleet with likelihood without dispersion
             if(diff(x$Model[[method]]$Distr%in%sdistr.set)==1)
               {
                #Only second fleet with perturbations
                if(waves.F1==0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks2   <- x$Model[[method]]$Dates[2:(waves.F2+1)]
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),1)),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+length(waves.F2)+3))]),NA),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),NA),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),1)),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F2+1)]),"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),NA),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),NA),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
                #Both fleets with perturbations
                if(waves.F1!=0 & waves.F2!=0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks1    <- x$Model[[method]]$Dates[2:(waves.F1+1)]
                      weeks2    <- x$Model[[method]]$Dates[(waves.F1+2):(waves.F1+waves.F2+1)]
                      dates.ranges.F1 <- vector("character",waves.F1)
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F1)
                        {
                         if(weeks1[w] <= 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year1, weeks1[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks1[w] > 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year2, weeks1[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing.F1=c("","",dates.ranges.F1,rep("",4+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),1)),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),NA),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),NA),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing.F1=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F1+1)]),rep("",4+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),1)),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[(waves.F1+2):(waves.F1+1+waves.F2)]),"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),NA),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),NA)/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),NA),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
               }
             #Both fleets with likelihood with dispersion
             else if(sum(x$Model[[method]]$Distr%in%sdistr.set)==0)
               {
                #Only second fleet with perturbations
                if(waves.F1==0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks2    <- x$Model[[method]]$Dates[2:(waves.F2+1)]
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]),tail(unlist(x$Model[[method]]$bt.stdev),2)[1])/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2)]),rep(0,waves.F2),unlist(x$Model[[method]]$bt.stdev[3:5]),tail(unlist(x$Model[[method]]$bt.stdev),2)[1])/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2)]),rep(1,waves.F2),unlist(x$Model[[method]]$bt.par[3:5]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F2+1)]),"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,6:(5+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),1),row.names=NULL)
                     names(partable)[2:3] <- paste(c("Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[4:6] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
                #Both fleets with perturbations
                else if(waves.F1!=0 & waves.F2!=0)
                  {
                   year1    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["StartDate"]]), "%Y"))
                   year2    <- as.numeric(format(as.Date(x$Data$Properties$Dates[["EndDate"]]), "%Y"))
                   days     <- as.Date(paste(year1, 1, 1, sep = "-"))+0:365*(year2-year1+1)
                   if(x$Data$Properties$Units["Time Step"]=="week")
                     {
                      weeks1    <- x$Model[[method]]$Dates[2:(waves.F1+1)]
                      weeks2    <- x$Model[[method]]$Dates[2:(waves.F2+1)]
                      dates.ranges.F1 <- vector("character",waves.F1)
                      dates.ranges.F2 <- vector("character",waves.F2)
                      for(w in 1:waves.F1)
                        {
                         if(weeks1[w] <= 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year1, weeks1[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks1[w] > 53)
                           {
                            dates.ranges.F1[w] <- paste(range(days[sprintf("%d %02d", year2, weeks1[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      for(w in 1:waves.F2)
                        {
                         if(weeks2[w] <= 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year1, weeks2[w]) == format(days, "%Y %U")]),collapse=' ')
                           }
                         if(weeks2[w] > 53)
                           {
                            dates.ranges.F2[w] <- paste(range(days[sprintf("%d %02d", year2, weeks2[w]-53) == format(days, "%Y %U")]),collapse=' ')
                           }
                        }
                      partable <- data.frame(Parameter=c("M.1/week",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing.F1=c("","",dates.ranges.F1,rep("",4+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.stdev),2)[1])/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),1),
                                             Timing.F2=c("","",dates.ranges.F2,"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                   if(x$Data$Properties$Units["Time Step"]=="day")
                     {
                      partable <- data.frame(Parameter=c("M.1/day",
                                                         paste("N0.",x$Data$Properties$Units["NumbersMultiplier"],sep=""),
                                                         paste("Rec.",x$Data$Properties$Units["NumbersMultiplier"],".Wave",1:waves.F2,sep=""),
                                                         paste("k.1/",x$Data$Properties$Fleets$Units[1],sep=""),
                                                         "alpha","beta",
                                                         paste("psi.",x$Data$Properties$Units["NumbersMultiplier"],".squared",sep="")),
                                             Timing.F1=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[2:(waves.F1+1)]),rep("",4+waves.F2-waves.F1)),
                                             Estimates.F1=c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),
                                             CVpCent.F1=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.stdev[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.stdev),2)[1])/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,3:(waves.F1+2))]),rep(0,waves.F2-waves.F1),unlist(x$Model[[method]]$bt.par[(waves.F1+3):(waves.F1+5)]),tail(unlist(x$Model[[method]]$bt.par),2)[1]),1),
                                             Timing.F2=c("","",as.character(as.Date(x$Data$Properties$Dates[["StartDate"]])+x$Model[[method]]$Dates[(waves.F1+2):(waves.F1+1+waves.F2)]),"","","",""),
                                             Estimates.F2=c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),
                                             CVpCent.F2=round(100*c(unlist(x$Model[[method]]$bt.stdev[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.stdev),1))/
                                                                  c(unlist(x$Model[[method]]$bt.par[c(1,2,(2+waves.F1+3+1):(2+waves.F1+3+waves.F2+3))]),tail(unlist(x$Model[[method]]$bt.par),1)),1),row.names=NULL)
                     names(partable)[2:4] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[1],sep="")
                     names(partable)[5:7] <- paste(c("Timing.","Estimates.","CVpCent."),x$Data$Properties$Fleets$Fleet[2],sep="")
                     }
                  }
               }
            }
         #End of two fleets
         }
      #End of week/day time step
      }
    return(partable)
   #End of function
   }
