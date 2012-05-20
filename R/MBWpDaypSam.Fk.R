MBWpDaypSam.Fk <-
function(BioData,Year,Season,col1,col2,len1,len2,w1,sd.w1,w2,sd.w2,cor.w1.w2)
                   {
                    n  <- 0;
                    x1 <- 0;
                    x2 <- 0;
                    x3 <- 0;
                    cov.w1.w2 <- cor.w1.w2*sqrt(sd.w1^2*sd.w2^2);
                    BioData.Year.Season <- BioData[BioData$Year==Year & BioData$Season==Season,]
                    for(i in 1:dim(BioData.Year.Season)[1])
                        {
                        n[i] <- sum(BioData.Year.Season[i,col1:col2]);
                        BioData.Year.Season$Mlen.cm[i]<-sum(seq(len1,len2,0.5)*BioData.Year.Season[i,col1:col2])/n[i]; #MEAN LENGTH cm
                        BioData.Year.Season$Slen.cm[i]<-sqrt((sum(seq(len1,len2,0.5)^2*BioData.Year.Season[i,col1:col2])-(sum(seq(len1,len2,0.5)*BioData.Year.Season[i,col1:col2])^2/n[i]))/(n[i]-1)); #SD Mlen cm
                        BioData.Year.Season$Mbm.kg[i]<-sum(w1*seq(len1,len2,0.5)^w2*BioData.Year.Season[i,col1:col2])/n[i]/1000; #MLE MEAN BODY WEIGHT
                        x1[i]<-sd.w1*(sum(w1*seq(len1,len2,0.5)^w2*BioData.Year.Season[i,col1:col2])/w1)^2; #SD MLE MEAN BODY WEIGHT
                        x2[i]<-2*w1*cov.w1.w2*(sum(w1*seq(len1,len2,0.5)^w2*BioData.Year.Season[i,col1:col2])/w1)*sum(log(seq(len1,len2,0.5))*(w1*seq(len1,len2,0.5)^w2*BioData.Year.Season[i,col1:col2]))/w1; #SD MLE MEAN BODY WEIGHT
                        x3[i]<-w1*sd.w2*sum(log(seq(len1,len2,0.5))*((w1*seq(len1,len2,0.5)^w2*BioData.Year.Season[i,col1:col2])/w1))^2; #SD MLE MEAN BODY WEIGHT
                        BioData.Year.Season$Sbm.kg[i]<-sqrt((1/n[i])^2*(x1[i]+x2[i]+x3[i]))/1000; #SD MLE MEAN BODY WEIGHT
                        }
                    rm(n,x1,x2,x3)
                    return(BioData.Year.Season)
                   }
