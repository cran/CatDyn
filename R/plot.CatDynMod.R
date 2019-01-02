plot.CatDynMod <-
function(x,leg.pos,Biom.tstep,Cat.tstep,Biom.xpos,Biom.ypos,Cat.xpos,Cat.ypos,diagnostics.panels=2,top.lab=TRUE,...)
    {
     if(class(x) != "CatDynMod")
       {stop("Pass an object of class 'CatDynMod' to plot.CatDynMod")}
     if(!is.numeric(Biom.tstep) | !is.numeric(Biom.xpos) | !is.numeric(Biom.ypos) | !is.numeric(Cat.xpos) | !is.numeric(Cat.ypos))
       {stop("Biom.tstep must be integer, Biom.xpos, Biom.ypos, Cat.xpos, Cat.ypos must be numeric")}
     if(Biom.tstep > tail(x$Model$Results[,1],1)) {stop("Biommas cannot be averaqed over a number of time steps as high as Biom.tstep")}
     if(Cat.tstep > tail(x$Model$Results[,1],1)) {stop("Catch cannot be summed over a number of time steps as high as Cat.tstep")}
     fleet.name <- x$Properties$Fleets$Fleet
     period     <- x$Model$Results$Period
     tstep      <- x$Properties$Units[1]
     catunits   <- paste("Catch (",x$Properties$Units[4],")",sep="")
     index      <- 1:length(x$Model$Results[,dim(x$Model$Results)[2]])
     Biom       <- round(mean(tail(x$Model$Results[,"Predicted.Biomass.tonnes"],Biom.tstep),na.rm=TRUE))
     Cat        <- ifelse(length(fleet.name) == 1,
                          round(1e-3*sum(tail(x$Model$Results[,5],Cat.tstep),na.rm=TRUE)),
                          round(1e-3*sum(tail(x$Model$Results[,5],Cat.tstep)+
                                         tail(x$Model$Results[,13],Cat.tstep),na.rm=TRUE)))
     options(warn=-1)
     for(i in 1:length(fleet.name))
       {
       obscat <- x$Model$Results[,i*3+5*(i-1)]   #1Fleet 3 2Fleet 3 11
       modcat <- x$Model$Results[,i*4+4*(i-1)]   #1Fleet 4 2Fleet 4 12
       resids <- x$Model$Results[,i*9-1*(i-1)]   #1Fleet 9 2Fleet 9 17
       if(diagnostics.panels==0)
         {
          par(mfrow=c(1,1))
          plot(x=period,y=obscat,pch=19, xlab=gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", tstep, perl=TRUE),ylab=catunits,ylim=c(0,max(obscat,modcat)),main="",cex=0.75)
          lines(x=period,y=modcat,col="black",lwd=2)
          legend(leg.pos,c("Observed Catch", "Predicted Catch"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
          if(x$Properties$Units[3] != "ind")
            {
             text(x=Biom.xpos*max(period), y=Biom.ypos*max(obscat), adj=0,lab=list(bquote("Biomass_(tons)" ==.(Biom))))
             text(x=Cat.xpos*max(period), y=Cat.ypos*max(obscat), adj=0,lab=list(bquote("Catch.ton" ==.(Cat))))
            }
          if(x$Model$Type[i] != 0)
            {
             if(i == 1)
               {
                if(x$Model$Type[i] > 0)
                  {
                   points(x=x$Model$Dates[2:(x$Model$Type[i]+1)],
                          y=rep(0,x$Model$Type[i]),
                          pch=10,cex=3)
                  }
                else
                  {
                   points(x=x$Model$Dates[seq(2,(2*abs(x$Model$Type[i])+1),2)],
                          y=rep(0,abs(x$Model$Type[i])),
                          pch=10,cex=3,col="red")
                   points(x=x$Model$Dates[seq(3,(2*abs(x$Model$Type[i])+1),2)],
                          y=rep(0,abs(x$Model$Type[i])),
                          pch=10,cex=3,col="blue")
                  }
               }
             else
               {
                points(x=x$Model$Dates[(2+x$Model$Type[i-1]):(1+x$Model$Type[i-1]+x$Model$Type[i])],
                       y=rep(0,x$Model$Type[i]),
                       pch=10,cex=3)
               }
            }
          if(top.lab==TRUE)
            {
             mtext(side=3,outer=TRUE,text=paste("Fleet = ",fleet.name[i],", Perturbations = ",x$Model$Type[i],", Distribution = ",gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", x$Model$Distribution[i], perl=TRUE),", Numerical algorithm = ", x$Model$Method, sep=""))
            }
          options(warn=0)
          devAskNewPage(ask=TRUE)
         }
       if(diagnostics.panels==1)
         {
          par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(4,4,2,2))
          plot(x=period,y=obscat,pch=19, xlab=gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", tstep, perl=TRUE),ylab=catunits,ylim=c(0,max(obscat,modcat)),main="",cex=0.75)
          lines(x=period,y=modcat,col="black",lwd=2)
          legend(leg.pos,c("Observed Catch", "Predicted Catch"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
          if(x$Properties$Units[3] != "ind")
            {
             text(x=Biom.xpos*max(period), y=Biom.ypos*max(obscat), adj=0,lab=list(bquote("Biomass_(tons)" ==.(Biom))))
             text(x=Cat.xpos*max(period), y=Cat.ypos*max(obscat), adj=0,lab=list(bquote("Catch.ton" ==.(Cat))))
            }
          if(x$Model$Type[i] != 0)
            {
             if(i == 1)
               {
                if(x$Model$Type[i] > 0)
                  {
                   points(x=x$Model$Dates[2:(x$Model$Type[i]+1)],
                          y=rep(0,x$Model$Type[i]),
                          pch=10,cex=3)
                  }
                else
                  {
                   points(x=x$Model$Dates[seq(2,(2*abs(x$Model$Type[i])+1),2)],
                          y=rep(0,abs(x$Model$Type[i])),
                          pch=10,cex=3,col="red")
                   points(x=x$Model$Dates[seq(3,(2*abs(x$Model$Type[i])+1),2)],
                          y=rep(0,abs(x$Model$Type[i])),
                          pch=10,cex=3,col="blue")
                  }
               }
             else
               {
                points(x=x$Model$Dates[(2+x$Model$Type[i-1]):(1+x$Model$Type[i-1]+x$Model$Type[i])],
                       y=rep(0,x$Model$Type[i]),
                       pch=10,cex=3)
               }
            }
          hist(x=resids,main="",xlab="Deviance Residuals",ylab="Frequency")
          plot(x=period,y=resids,xlab=gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", tstep, perl=TRUE), ylab="Deviance Residuals",pch=1,type="n")
          abline(h=0,lwd=2)
          text(x=period,y=resids,lab=format(period),cex=0.75)
          qqplot(x=obscat,y=modcat,xlab=paste("Observed Catch (",x$Properties$Units[4],")", sep=""),ylab=paste("Predicted Catch (",x$Properties$Units[4],")", sep=""),pch=1)
          abline(a=0,b=1,lwd=2)
          if(top.lab==TRUE)
            {
             mtext(side=3,outer=TRUE,text=paste("Fleet = ",fleet.name[i],", Perturbations = ",x$Model$Type[i],", Distribution = ",gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", x$Model$Distribution[i], perl=TRUE),", Numerical algorithm = ", x$Model$Method, sep=""))
            }
          options(warn=0)
          devAskNewPage(ask=TRUE)
         }
       if(diagnostics.panels==2)
         {
          layout(matrix(c(1,1,1,1,1,1,1,1,1,2,3,4,2,3,4),5,3, byrow = TRUE))
          par(oma=c(2,2,2,1), mar=c(4,4.5,2,2))
          plot(x=period,y=obscat,pch=19, xlab=gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", tstep, perl=TRUE),ylab=catunits,ylim=c(0,max(obscat,modcat)),main="",cex.lab=1.4,cex.axis=1.5)
          lines(x=period,y=modcat,col="black",lwd=2)
          legend(leg.pos,c("Observed Catch", "Predicted Catch"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"),cex=1.5)
          if(x$Properties$Units[3] != "ind")
            {
             text(x=Biom.xpos*max(period), y=Biom.ypos*max(obscat), adj=0,lab=list(bquote("Biomass_(tons)" ==.(Biom))))
             text(x=Cat.xpos*max(period), y=Cat.ypos*max(obscat), adj=0,lab=list(bquote("Catch.ton" ==.(Cat))))
            }
          if(x$Model$Type[i] != 0)
            {
             if(i == 1)
               {
                if(x$Model$Type[i] > 0)
                  {
                   points(x=x$Model$Dates[2:(x$Model$Type[i]+1)],
                          y=rep(0,x$Model$Type[i]),
                          pch=10,cex=3)
                  }
                else
                  {
                   points(x=x$Model$Dates[seq(2,(2*abs(x$Model$Type[i])+1),2)],
                          y=rep(0,abs(x$Model$Type[i])),
                          pch=10,cex=3,col="red")
                   points(x=x$Model$Dates[seq(3,(2*abs(x$Model$Type[i])+1),2)],
                          y=rep(0,abs(x$Model$Type[i])),
                          pch=10,cex=3,col="blue")
                  }
               }
             else
               {
                points(x=x$Model$Dates[(2+x$Model$Type[i-1]):(1+x$Model$Type[i-1]+x$Model$Type[i])],
                       y=rep(0,x$Model$Type[i]),
                       pch=10,cex=3)
               }
            }
          hist(x=resids,main="",xlab="Deviance Residuals",ylab="Frequency",cex.lab=1.4,cex.axis=1.5)
          plot(x=period,y=resids,xlab=gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", tstep, perl=TRUE),ylab="Deviance Residuals",pch=1,type="n",cex.lab=1.4,cex.axis=1.5)
          abline(h=0,lwd=2)
          text(x=period,y=resids,lab=format(period),cex=0.75)
          qqplot(x=obscat,y=modcat,xlab=paste("Observed Catch (",x$Properties$Units[4],")", sep=""),ylab=paste("Predicted Catch (",x$Properties$Units[4],")", sep=""),pch=1,cex.lab=1.4,cex.axis=1.5)
          abline(a=0,b=1,lwd=2)
          if(top.lab==TRUE)
            {
             mtext(side=3,outer=TRUE,text=paste("Fleet = ",fleet.name[i],", Perturbations = ",x$Model$Type[i],", Distribution = ",gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", x$Model$Distribution[i], perl=TRUE),", Numerical algorithm = ", x$Model$Method, sep=""))
            }
          options(warn=0)
          devAskNewPage(ask=TRUE)
         }
       }
     devAskNewPage(ask=FALSE)
    }
