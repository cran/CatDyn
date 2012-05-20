BioFishpDay.Fk <-
function(BioData,CatEffData,Year,Season,StartDate,EndDate,FGTable,span,sds,plot)
                   {
                   CatEffData.Year.Season <- CatEffData[CatEffData$Year==Year & CatEffData$Season==Season,]
                   SeasonSeq   <- seq(StartDate,EndDate,by='day');
                   u <- aggregate(BioData$Mlen.cm,list(BioData$Date),mean)
                   names(u) <- c("Date","Mlen.cm")
                   u <- u[u$Date>=StartDate & u$Date<=EndDate,]
                   u.date.series <- data.frame(SeasonSeq, NA)
                   names(u.date.series) <- c("Date","Mlen.cm")
                   u.date.series[match(unique(u$Date),unique(u.date.series$Date)),2] <- u$Mlen.cm
                   Mlen.cm.Trend <- loess(u.date.series$Mlen.cm~julian(u.date.series$Date,origin=as.Date(paste(Year,"-01-01",sep=""))),span=span,control=loess.control(surface = "direct"))
                   Pred.Mlen.cm  <- predict(Mlen.cm.Trend,data.frame(julian(SeasonSeq,origin=as.Date(paste(Year,"-01-01",sep="")))),se=TRUE)
                   u.date.series$Mlen.cm[which(is.na(u.date.series$Mlen.cm))] <- Pred.Mlen.cm$fit[which(is.na(u.date.series$Mlen.cm))]
                   #Standard Deviation of Mean Length (cm)
                   v <- aggregate(BioData$Slen.cm,list(BioData$Date),mean)
                   names(v) <- c("Date","Slen.cm")
                   v <- v[v$Date>=StartDate & v$Date<=EndDate,]
                   v.date.series <- data.frame(SeasonSeq, NA)
                   names(v.date.series) <- c("Date","Slen.cm")
                   v.date.series[match(unique(v$Date),unique(v.date.series$Date)),2] <- v$Slen.cm
                   Slen.cm.Trend <- loess(v.date.series$Slen.cm~julian(v.date.series$Date,origin=as.Date(paste(Year,"-01-01",sep=""))),span=span,control=loess.control(surface = "direct"))
                   Pred.Slen.cm  <- predict(Slen.cm.Trend,data.frame(julian(SeasonSeq,origin=as.Date(paste(Year,"-01-01",sep="")))),se=TRUE)
                   v.date.series$Slen.cm[which(is.na(v.date.series$Slen.cm))] <- Pred.Slen.cm$fit[which(is.na(v.date.series$Slen.cm))]
                   #Mean Body Mass (kg)
                   x <- aggregate(BioData$Mbm.kg,list(BioData$Date),mean)
                   names(x) <- c("Date","Mbm.kg")
                   x <- x[x$Date>=StartDate & x$Date<=EndDate,]
                   x.date.series <- data.frame(SeasonSeq, NA)
                   names(x.date.series) <- c("Date","Mbm.kg")
                   x.date.series[match(unique(x$Date),unique(x.date.series$Date)),2] <- x$Mbm.kg
                   Mbm.kg.Trend <- loess(x.date.series$Mbm.kg~julian(x.date.series$Date,origin=as.Date(paste(Year,"-01-01",sep=""))),span=span,control=loess.control(surface="direct"))
                   Pred.Mbm.kg  <- predict(Mbm.kg.Trend,data.frame(julian(SeasonSeq,origin=as.Date(paste(Year,"-01-01",sep="")))),se=TRUE)
                   x.date.series$Mbm.kg[which(is.na(x.date.series$Mbm.kg))] <- Pred.Mbm.kg$fit[which(is.na(x.date.series$Mbm.kg))]
                   y <- aggregate(BioData$Sbm.kg,list(BioData$Date),mean)
                   names(y) <- c("Date","Sbm.kg")
                   y <- y[y$Date>=StartDate & y$Date<=EndDate,]
                   y.date.series <- data.frame(SeasonSeq, NA)
                   names(y.date.series) <- c("Date","Sbm.kg")
                   y.date.series[match(unique(y$Date),unique(y.date.series$Date)),2] <- y$Sbm.kg
                   Sbm.kg.Trend <- loess(y.date.series$Sbm.kg~julian(y.date.series$Date,origin=as.Date(paste(Year,"-01-01",sep=""))),span=span,,control=loess.control(surface="direct"))
                   Pred.Sbm.kg  <- predict(Sbm.kg.Trend,data.frame(julian(SeasonSeq,origin=as.Date(paste(Year,"-01-01",sep="")))),se=TRUE)
                   y.date.series$Sbm.kg[which(is.na(y.date.series$Sbm.kg))] <- Pred.Sbm.kg$fit[which(is.na(y.date.series$Sbm.kg))]
                   if(plot==TRUE)
                      {
                      par(mfrow=c(4,2),oma=c(2,2,0,0),mar=c(2,2,0,0))
                      plot(u.date.series$Date,u.date.series$Mlen.cm, ylim=c(min(Pred.Mlen.cm$fit-sds*Pred.Mlen.cm$se.fit),max(Pred.Mlen.cm$fit+sds*Pred.Mlen.cm$se.fit)))
                      lines(u.date.series$Date,Pred.Mlen.cm$fit,col="blue",lwd=2)
                      lines(u.date.series$Date,Pred.Mlen.cm$fit+sds*Pred.Mlen.cm$se.fit,lwd=2,col="red")
                      lines(u.date.series$Date,Pred.Mlen.cm$fit-sds*Pred.Mlen.cm$se.fit,lwd=2,col="red")
                      legend("topleft",legend="Mean Mantle Length (cm)",cex=0.75)
                      plot(v.date.series$Date,v.date.series$Slen.cm, ylim=c(min(Pred.Slen.cm$fit-sds*Pred.Slen.cm$se.fit),max(Pred.Slen.cm$fit+sds*Pred.Slen.cm$se.fit)))
                      lines(v.date.series$Date,Pred.Slen.cm$fit,col="blue",lwd=2)
                      lines(v.date.series$Date,Pred.Slen.cm$fit+sds*Pred.Slen.cm$se.fit,lwd=2,col="red")
                      lines(v.date.series$Date,Pred.Slen.cm$fit-sds*Pred.Slen.cm$se.fit,lwd=2,col="red")
                      legend("topleft",legend="SD of Mean Mantle Length (cm)",cex=0.75)
                      plot(x.date.series$Date,x.date.series$Mbm.kg, ylim=c(min(Pred.Mbm.kg$fit-sds*Pred.Mbm.kg$se.fit),max(Pred.Mbm.kg$fit+sds*Pred.Mbm.kg$se.fit)))
                      lines(x.date.series$Date,Pred.Mbm.kg$fit,col="blue",lwd=2)
                      lines(x.date.series$Date,Pred.Mbm.kg$fit+sds*Pred.Mbm.kg$se.fit,lwd=2,col="red")
                      lines(x.date.series$Date,Pred.Mbm.kg$fit-sds*Pred.Mbm.kg$se.fit,lwd=2,col="red")
                      legend("topleft",legend="Mean Body Mass (g)",cex=0.75)
                      plot(y.date.series$Date,y.date.series$Sbm.kg, ylim=c(min(Pred.Sbm.kg$fit-sds*Pred.Sbm.kg$se.fit),max(Pred.Sbm.kg$fit+sds*Pred.Sbm.kg$se.fit)))
                      lines(y.date.series$Date,Pred.Sbm.kg$fit,col="blue",lwd=2)
                      lines(y.date.series$Date,Pred.Sbm.kg$fit+sds*Pred.Sbm.kg$se.fit,lwd=2,col="red")
                      lines(y.date.series$Date,Pred.Sbm.kg$fit-sds*Pred.Sbm.kg$se.fit,lwd=2,col="red")
                      legend("topleft",legend="SD of Mean Body Mass (g)",cex=0.75)
                      }
                   #Replacing outliers
                   #Standard Deviation of Mean Body Mass
                   #Mean length (cm)
                   u.date.series$Mlen.cm[which(u.date.series$Mlen.cm>Pred.Mlen.cm$fit+sds*Pred.Mlen.cm$se.fit)] <- Pred.Mlen.cm$fit[which(u.date.series$Mlen.cm>Pred.Mlen.cm$fit+sds*Pred.Mlen.cm$se.fit)]
                   u.date.series$Mlen.cm[which(u.date.series$Mlen.cm<Pred.Mlen.cm$fit-sds*Pred.Mlen.cm$se.fit)] <- Pred.Mlen.cm$fit[which(u.date.series$Mlen.cm<Pred.Mlen.cm$fit-sds*Pred.Mlen.cm$se.fit)]
                   #Standard Deviation of Mean Length (cm)
                   v.date.series$Slen.cm[which(v.date.series$Slen.cm>Pred.Slen.cm$fit+sds*Pred.Slen.cm$se.fit)] <- Pred.Slen.cm$fit[which(v.date.series$Slen.cm>Pred.Slen.cm$fit+sds*Pred.Slen.cm$se.fit)]
                   v.date.series$Slen.cm[which(v.date.series$Slen.cm<Pred.Slen.cm$fit-sds*Pred.Slen.cm$se.fit)] <- Pred.Slen.cm$fit[which(v.date.series$Slen.cm<Pred.Slen.cm$fit-sds*Pred.Slen.cm$se.fit)]
                   #Mean Body Mass (kg)
                   x.date.series$Mbm.kg[which(x.date.series$Mbm.kg>Pred.Mbm.kg$fit+sds*Pred.Mbm.kg$se.fit)] <- Pred.Mbm.kg$fit[which(x.date.series$Mbm.kg>Pred.Mbm.kg$fit+sds*Pred.Mbm.kg$se.fit)]
                   x.date.series$Mbm.kg[which(x.date.series$Mbm.kg<Pred.Mbm.kg$fit-sds*Pred.Mbm.kg$se.fit)] <- Pred.Mbm.kg$fit[which(x.date.series$Mbm.kg<Pred.Mbm.kg$fit-sds*Pred.Mbm.kg$se.fit)]
                   #Standard Deviation of Mean Body Mass (kg)
                   y.date.series$Sbm.kg[which(y.date.series$Sbm.kg>Pred.Sbm.kg$fit+sds*Pred.Sbm.kg$se.fit)] <- Pred.Sbm.kg$fit[which(y.date.series$Sbm.kg>Pred.Sbm.kg$fit+sds*Pred.Sbm.kg$se.fit)]
                   y.date.series$Sbm.kg[which(y.date.series$Sbm.kg<Pred.Sbm.kg$fit-sds*Pred.Sbm.kg$se.fit)] <- Pred.Sbm.kg$fit[which(y.date.series$Sbm.kg<Pred.Sbm.kg$fit-sds*Pred.Sbm.kg$se.fit)]
                   if(plot==TRUE)
                      {
                      plot(u.date.series$Date,u.date.series$Mlen.cm, ylim=c(min(Pred.Mlen.cm$fit-sds*Pred.Mlen.cm$se.fit),max(Pred.Mlen.cm$fit+sds*Pred.Mlen.cm$se.fit)))
                      lines(u.date.series$Date,Pred.Mlen.cm$fit,col="blue",lwd=2)
                      lines(u.date.series$Date,Pred.Mlen.cm$fit+sds*Pred.Mlen.cm$se.fit,lwd=2,col="red")
                      lines(u.date.series$Date,Pred.Mlen.cm$fit-sds*Pred.Mlen.cm$se.fit,lwd=2,col="red")
                      legend("topleft",legend="Mean Mantle Length (cm)",cex=0.75)
                      plot(v.date.series$Date,v.date.series$Slen.cm, ylim=c(min(Pred.Slen.cm$fit-sds*Pred.Slen.cm$se.fit),max(Pred.Slen.cm$fit+sds*Pred.Slen.cm$se.fit)))
                      lines(v.date.series$Date,Pred.Slen.cm$fit,col="blue",lwd=2)
                      lines(v.date.series$Date,Pred.Slen.cm$fit+sds*Pred.Slen.cm$se.fit,lwd=2,col="red")
                      lines(v.date.series$Date,Pred.Slen.cm$fit-sds*Pred.Slen.cm$se.fit,lwd=2,col="red")
                      legend("topleft",legend="SD of Mean Mantle Length (cm)",cex=0.75)
                      plot(x.date.series$Date,x.date.series$Mbm.kg, ylim=c(min(Pred.Mbm.kg$fit-sds*Pred.Mbm.kg$se.fit),max(Pred.Mbm.kg$fit+sds*Pred.Mbm.kg$se.fit)))
                      lines(x.date.series$Date,Pred.Mbm.kg$fit,col="blue",lwd=2)
                      lines(x.date.series$Date,Pred.Mbm.kg$fit+sds*Pred.Mbm.kg$se.fit,lwd=2,col="red")
                      lines(x.date.series$Date,Pred.Mbm.kg$fit-sds*Pred.Mbm.kg$se.fit,lwd=2,col="red")
                      legend("topleft",legend="Mean Body Mass (g)",cex=0.75)
                      plot(y.date.series$Date,y.date.series$Sbm.kg, ylim=c(min(Pred.Sbm.kg$fit-sds*Pred.Sbm.kg$se.fit),max(Pred.Sbm.kg$fit+sds*Pred.Sbm.kg$se.fit)))
                      lines(y.date.series$Date,Pred.Sbm.kg$fit,col="blue",lwd=2)
                      lines(y.date.series$Date,Pred.Sbm.kg$fit+sds*Pred.Sbm.kg$se.fit,lwd=2,col="red")
                      lines(y.date.series$Date,Pred.Sbm.kg$fit-sds*Pred.Sbm.kg$se.fit,lwd=2,col="red")
                      legend("topleft",legend="SD of Mean Body Mass (g)",cex=0.75)
                      }
                   #Process catch effort data
                   #Beauchene
                   l <- if(FGTable["Beauchene"]==0) data.frame(SeasonSeq,0) else
                        aggregate(CatEffData.Year.Season$Catch.kg[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Beauchene"],
                        list(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Beauchene"]),
                        sum)
                   names(l) <- c("Date","B.Catch.kg")
                   l.date.series <- data.frame(SeasonSeq, NA)
                   names(l.date.series) <- c("Date","B.Catch.kg")
                   l.date.series[match(unique(l$Date),unique(l.date.series$Date)),2] <- l$B.Catch.kg
                   l.date.series$B.Catch.kg[which(is.na(l.date.series$B.Catch.kg))] <- 0
                   #Trawl.Time.h
                   m <- if(FGTable["Beauchene"]==0) data.frame(SeasonSeq,0) else
                        aggregate(CatEffData.Year.Season$Trawl.Time.h[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Beauchene"],
                        list(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Beauchene"]),
                        sum)
                   names(m) <- c("Date","B.Trawl.Time.h")
                   m.date.series <- data.frame(SeasonSeq, NA)
                   names(m.date.series) <- c("Date","B.Trawl.Time.h")
                   m.date.series[match(unique(m$Date),unique(m.date.series$Date)),2] <- m$B.Trawl.Time.h
                   m.date.series$B.Trawl.Time.h[which(is.na(m.date.series$B.Trawl.Time.h))] <- 0
                   #VesselDays
                   n <- if(FGTable["Beauchene"]==0) data.frame(SeasonSeq,0) else
                        data.frame(table(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Beauchene"]))
                   names(n) <- c("Date","B.N.Vessels")
                   n$Date   <- as.Date(n$Date)
                   n.date.series <- data.frame(SeasonSeq, NA)
                   names(n.date.series) <- c("Date","B.N.Vessels")
                   n.date.series[match(unique(n$Date),unique(n.date.series$Date)),2] <- n$B.N.Vessels
                   n.date.series$B.N.Vessels[which(is.na(n.date.series$B.N.Vessels))] <- 0
                   #Central
                   #Catch.kg
                   o <- if(FGTable["Central"]==0) data.frame(SeasonSeq,0) else
                        aggregate(CatEffData.Year.Season$Catch.kg[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Central"],
                        list(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Central"]),
                        sum)
                   names(o) <- c("Date","C.Catch.kg")
                   o.date.series <- data.frame(SeasonSeq, NA)
                   names(o.date.series) <- c("Date","C.Catch.kg")
                   o.date.series[match(unique(o$Date),unique(o.date.series$Date)),2] <- o$C.Catch.kg
                   o.date.series$C.Catch.kg[which(is.na(o.date.series$C.Catch.kg))] <- 0
                   #Trawl.Time.h
                   p <- if(FGTable["Central"]==0) data.frame(SeasonSeq,0) else
                        aggregate(CatEffData.Year.Season$Trawl.Time.h[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Central"],
                        list(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Central"]),
                        sum)
                   names(p) <- c("Date","C.Trawl.Time.h")
                   p.date.series <- data.frame(SeasonSeq, NA)
                   names(p.date.series) <- c("Date","C.Trawl.Time.h")
                   p.date.series[match(unique(p$Date),unique(p.date.series$Date)),2] <- p$C.Trawl.Time.h
                   p.date.series$C.Trawl.Time.h[which(is.na(p.date.series$C.Trawl.Time.h))] <- 0
                   #VesselDays
                   q <- if(FGTable["Central"]==0) data.frame(SeasonSeq,0) else
                        data.frame(table(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="Central"]))
                   names(q) <- c("Date","C.N.Vessels")
                   q$Date   <- as.Date(q$Date)
                   q.date.series <- data.frame(SeasonSeq, NA)
                   names(q.date.series) <- c("Date","C.N.Vessels")
                   q.date.series[match(unique(q$Date),unique(q.date.series$Date)),2] <- q$C.N.Vessels
                   q.date.series$C.N.Vessels[which(is.na(q.date.series$C.N.Vessels))] <- 0
                   #North
                   #Catch.kg
                   r <- if(FGTable["North"]==0) data.frame(SeasonSeq,0) else
                        aggregate(CatEffData.Year.Season$Catch.kg[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="North"],
                        list(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="North"]),
                        sum)
                   names(r) <- c("Date","N.Catch.kg")
                   r.date.series <- data.frame(SeasonSeq, NA)
                   names(r.date.series) <- c("Date","N.Catch.kg")
                   r.date.series[match(unique(r$Date),unique(r.date.series$Date)),2] <- r$N.Catch.kg
                   r.date.series$N.Catch.kg[which(is.na(r.date.series$N.Catch.kg))] <- 0
                   #Trawl.Time.h
                   s <- if(FGTable["North"]==0) data.frame(SeasonSeq,0) else
                        aggregate(CatEffData.Year.Season$Trawl.Time.h[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="North"],
                        list(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="North"]),
                        sum)
                   names(s) <- c("Date","N.Trawl.Time.h")
                   s.date.series <- data.frame(SeasonSeq, NA)
                   names(s.date.series) <- c("Date","N.Trawl.Time.h")
                   s.date.series[match(unique(s$Date),unique(s.date.series$Date)),2] <- s$N.Trawl.Time.h
                   s.date.series$N.Trawl.Time.h[which(is.na(s.date.series$N.Trawl.Time.h))] <- 0
                   #VesselDays
                   t <- if(FGTable["Central"]==0) data.frame(SeasonSeq,0) else
                        data.frame(table(CatEffData.Year.Season$Date[CatEffData.Year.Season$Date <= EndDate & CatEffData.Year.Season$Date >= StartDate & CatEffData.Year.Season$Fishing.ground=="North"]))
                   names(t) <- c("Date","N.N.Vessels")
                   t$Date   <- as.Date(t$Date)
                   t.date.series <- data.frame(SeasonSeq, NA)
                   names(t.date.series) <- c("Date","N.N.Vessels")
                   t.date.series[match(unique(t$Date),unique(t.date.series$Date)),2] <- t$N.N.Vessels
                   t.date.series$N.N.Vessels[which(is.na(t.date.series$N.N.Vessels))] <- 0
                   SeasonData.Fk <- data.frame(Date=SeasonSeq,
                                               Day=julian(SeasonSeq,origin=as.Date(paste(Year,"-01-01",sep=""))),
                                               Mlen.cm=u.date.series$Mlen.cm,
                                               Slen.cm=v.date.series$Slen.cm,
                                               Mbm.kg=x.date.series$Mbm.kg,
                                               Sbm.kg=y.date.series$Sbm.kg,
                                               B.Catch.kg=l.date.series$B.Catch.kg,
                                               B.Catch.bill=l.date.series$B.Catch.kg/x.date.series$Mbm.kg/1e9,
                                               B.Trawl.Time.h=m.date.series$B.Trawl.Time.h,
                                               B.N.Vessels=n.date.series$B.N.Vessels,
                                               C.Catch.kg=o.date.series$C.Catch.kg,
                                               C.Catch.bill=o.date.series$C.Catch.kg/x.date.series$Mbm.kg/1e9,
                                               C.Trawl.Time.h=p.date.series$C.Trawl.Time.h,
                                               C.N.Vessels=q.date.series$C.N.Vessels,
                                               N.Catch.kg=r.date.series$N.Catch.kg,
                                               N.Catch.bill=r.date.series$N.Catch.kg/x.date.series$Mbm.kg/1e9,
                                               N.Trawl.Time.h=s.date.series$N.Trawl.Time.h,
                                               N.N.Vessels=t.date.series$N.N.Vessels)
                   SeasonData.Fk.B <- data.frame(period=SeasonData.Fk$Day,
                                                 obseff1=SeasonData.Fk$B.Trawl.Time.h,
                                                 obseff2=SeasonData.Fk$B.N.Vessels,
                                                 obscat=SeasonData.Fk$B.Catch.bill,
                                                 mn.wt=SeasonData.Fk$Mbm.kg)
                   SeasonData.Fk.C <- data.frame(period=SeasonData.Fk$Day,
                                                 obseff1=SeasonData.Fk$C.Trawl.Time.h,
                                                 obseff2=SeasonData.Fk$C.N.Vessels,
                                                 obscat=SeasonData.Fk$C.Catch.bill,
                                                 mn.wt=SeasonData.Fk$Mbm.kg)
                   SeasonData.Fk.N <- data.frame(period=SeasonData.Fk$Day,
                                                 obseff1=SeasonData.Fk$N.Trawl.Time.h,
                                                 obseff2=SeasonData.Fk$N.N.Vessels,
                                                 obscat=SeasonData.Fk$N.Catch.bill,
                                                 mn.wt=SeasonData.Fk$Mbm.kg)

                   class(SeasonData.Fk.B) <- "CatDynData"
                   class(SeasonData.Fk.C) <- "CatDynData"
                   class(SeasonData.Fk.N) <- "CatDynData"
                   SeasonData.Fk <- list(FullSeason=SeasonData.Fk,
                                         Beauchene=SeasonData.Fk.B,
                                         Central=SeasonData.Fk.C,
                                         North=SeasonData.Fk.N)
                   SeasonData.Fk
                   }
