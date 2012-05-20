plot.CatDynMod <-
function(x,tstep,mult,Biom,AIC,top.text,leg.pos,AIC.xpos,AIC.ypos,Biom.tstep,Biom.xpos,Biom.ypos,p.dates,...)
    {
     options(warn=-1)
     if(any(is.nan(x$modcat))) {stop("Change initial values to increase predicted catch")}
     par(mfrow=c(2,2),oma=c(0,0,2,0),mar=c(4,4,2,2))
     plot(x=x$period,y=x$obscat,pch=19, xlab=tstep, ylab=paste("Catch (",mult,")",sep=""),ylim=c(0,max(x$obscat,x$modcat)),main="",cex=0.75)
     lines(x=x$period,y=x$modcat,col="black",lwd=2)
     legend(leg.pos,c("Observed Catch", "Predicted Catch"), bty="n", pch=c(19, NA), lty=c(0, 1), lwd=c(0,2), col=c("black", "black"))
     if(Biom.tstep=="Fin"  |  Biom.tstep=="fin")
        {
        text(x=Biom.xpos*max(x$period), y=Biom.ypos*max(x$obscat), adj=0,lab=ifelse(sum(x$modcat)>Biom | sum(x$obscat)>Biom,"Collapse",list(bquote("Biom.End_(tons)" ==.(Biom)))))
        }
     else
        {
        if(Biom.tstep=="Mid"  |  Biom.tstep=="mid")
           {
           text(x=Biom.xpos*max(x$period), y=Biom.ypos*max(x$obscat), adj=0,lab=ifelse(sum(x$modcat)>Biom | sum(x$obscat)>Biom,"Collapse",list(bquote("Biom.Mid_(tons)" ==.(Biom)))))
           }
        else
           {
           text(x=Biom.xpos*max(x$period), y=Biom.ypos*max(x$obscat), adj=0,lab=ifelse(sum(x$modcat)>Biom | sum(x$obscat)>Biom,"Collapse",list(bquote("Biom.Ini_(tons)" ==.(Biom)))))
           }
        }
     text(x=AIC.xpos*max(x$period), y=AIC.ypos*max(x$obscat), adj=0, lab=bquote("AIC" == .(AIC)))
     if(p.dates!=0)
        {
        points(x=p.dates,y=x$obscat[p.dates-x$period[1]+1],pch=10,cex=3)
        }
     hist(x=x$resids,main="",xlab="Residuals",ylab="Relative Frequency",prob=TRUE)
     plot(x=x$period,y=x$resids,xlab=tstep, ylab="Residuals",pch=1,type="n")
     abline(h=0,lwd=2)
     text(x=x$period,y=x$resids,lab=format(x$period),cex=0.75)
     qqplot(x=x$obscat,y=x$modcat,xlab=paste("Observed Catch (",mult,")", sep=""),ylab=paste("Predicted Catch (",mult,")", sep=""),pch=1)
     abline(a=0,b=1,lwd=2)
     mtext(side=3,outer=TRUE,text=top.text)
     options(warn=0)
    }
