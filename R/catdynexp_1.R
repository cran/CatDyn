catdynexp <-
function(x, p, par, dates, distr)
    {
     fleet.name <- x$Properties$Fleets$Fleet
     distr.set  <- c("normal","apnormal","lognormal","aplnormal","gamma","poisson","negbin")
     if(class(x) != "CatDynData") {stop("Pass an object of class CatDynData as first argument - see the help for as.CatDynData")}
     if(any(is.na(p))) {stop("NAs are not allowed in the perturbations per fleet vector 'p'")}
     if(length(p) < 1 || length(p) > 2) {stop("The integer vector 'p' determines the number of perturbations per fleet; its length must be 1 for single fleet, 2 for two fleets")}
     if(any(p < -20) || any(p > 100)) {stop("The number of perturbations per fleet shall not be more than 100, nor more than 20 for transit fisheries")}
     if(class(p) != "numeric") {stop("'p' must be a numeric vector or a scalar")}
     if(any(is.na(par))) {stop("NAs are not allowed in the parameter vector 'par'")}
     if(length(par) < 4 || length(par) > 208) {stop("For any of the model versions the number of parameters is > 3 and < 208")}
     if(class(par) != "numeric") {stop("'par' must be a numeric vector")}
     if(any(is.na(dates))){stop("NAs are not allowed in the dates vector")}
     if(dates[1] > dates[2]) {stop("Initial date must not be less than next date")}
     if(tail(dates,1) < tail(dates,2)[1]) {stop("Final date must not be less than previous date")}
     options(warn=-1)
     if(length(fleet.name) == 1)
       {
        if(p >= 0)
          {
           if(sum(sort(p) == p)<length(p)) {stop("Number of perturbations per fleet must not be arranged in descending order")}
           if(length(dates) != (sum(p)+2)) {stop("The dates vector must contain initial time step, time steps of all perturbations, and final time step")}
           if(length(fleet.name) == 1)
             {
              if(sum(sort(dates) == dates)<length(dates)){stop("Dates must be arranged in ascending order")} 
             }
           else if(length(fleet.name) == 2)
             {
              if(sum(sort(dates[1:(p[1]+1)]) == dates[1:(p[1]+1)]) < length(dates[1:(p[1]+1)]) || sum(sort(dates[(p[1]+2):(p[1]+p[2]+2)]) == dates[(p[1]+2):(p[1]+p[2]+2)]) < length(dates[(p[1]+2):(p[1]+p[2]+2)]))
                {stop("Perturbation dates for each fleet must be arranged in ascending order")}
              #if(length(unique(dates[1:(p[1]+1)])) != length(dates[1:(p[1]+1)]) || length(unique(dates[(p[1]+2):(p[1]+p[2]+2)])) != length(dates[(p[1]+2):(p[1]+p[2]+2)]))   
              #  {stop("Perturbation dates inside each fleet should be all distinct and distinct from initial and final date")}
             }
          }
        else
          {
          if(length(dates) != (2*abs(p)+2)) {stop("The dates vector must contain initial time step, time step of entry of each perturbation immediately followed by time step of its exit, for all perturbations ordered by time of entry, and final time step")}
          }           
        #if(class(dates) != "integer") {stop("'dates' must be a vector of integers")}
        if(!distr%in%distr.set)
          {stop("distr must be 'normal','apnormal','lognormal','aplnormal','gamma', 'poisson' or 'negbin', see help pages for CatDynFit")}
        parlist <- list(par=par, 
                        dates=dates, 
                        obseff1=x$Data[[fleet.name]][,2], 
                        obscat1=x$Data[[fleet.name]][,5], 
                        obsmbm1=x$Data[[fleet.name]][,4], 
                        distr=distr, 
                        properties=x$Properties,
                        output="predict");
        if(p==0)
          {
          if(length(dates) != 2)
            {
             stop("For a 1-fleet 0 perturbation (simple depletion model) dates must be a vector with the following time step marks: initial, final")
            }
          if(length(par) != 5)
            {
             stop("For a 1-fleet 0 perturbation (simple depletion model) par must be a vector of length 5")
            }
          results <- do.call(.CDMN0P, parlist);
          }
        else if(p==1)
          {
          if(length(dates) != 3)
            {
             stop("For a 1-fleet 1-perturbation model dates must be a vector with the following time step marks: initial, perturbation 1, final")
            }
          if(length(par) != 6)
            {
             stop("For a 1-fleet 1-perturbation model par must be a vector of length 6")
            }
          results <- do.call(.CDMN1P, parlist);
          }
        else if(p==2)
          {
          if(length(dates) != 4)
            {
             stop("For a 1-fleet 2-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, perturbation 2, final")
            }
          if(length(par) != 7)
            {
             stop("For a 1-fleet 2-perturbation model par must be a vector of length 7")
            }
          results <- do.call(.CDMN2P, parlist);
          }
        else if(p==3)
          {
          if(length(dates) != 5)
            {
             stop("For a 1-fleet 3-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, perturbation 2, perturbation 3, final")
            }
          if(length(par) != 8)
            {
             stop("For a 1-fleet 3-perturbation model par must be a vector of length 8")
            }
          results <- do.call(.CDMN3P, parlist);
          }
        else if(p==4)
          {
          if(length(dates) != 6)
            {
             stop("For a 1-fleet 4-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 4, final")
            }
          if(length(par) != 9)
            {
             stop("For a 1-fleet 4-perturbation model par must be a vector of length 9")
            }
          results <- do.call(.CDMN4P, parlist);
          }
        else if(p==5)
          {
          if(length(dates) != 7)
            {
             stop("For a 1-fleet 5-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 5, final")
            }
          if(length(par) != 10)
            {
             stop("For a 1-fleet 5-perturbation model par must be a vector of length 10")
            }
          results <- do.call(.CDMN5P, parlist);
          }
        else if(p==6)
          {
          if(length(dates) != 8)
            {
             stop("For a 1-fleet 6-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 6, final")
            }
          if(length(par) != 11)
            {
             stop("For a 1-fleet 6-perturbation model par must be a vector of length 11")
            }
          results <- do.call(.CDMN6P, parlist);
          }
        else if(p==7)
          {
          if(length(dates) != 9)
            {
             stop("For a 1-fleet 7-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 7, final")
            }
          if(length(par) != 12)
            {
             stop("For a 1-fleet 7-perturbation model par must be a vector of length 12")
            }
          results <- do.call(.CDMN7P, parlist);
          }
        else if(p==8)
          {
          if(length(dates) != 10)
            {
             stop("For a 1-fleet 8-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 8, final")
            }
          if(length(par) != 13)
            {
             stop("For a 1-fleet 8-perturbation model par must be a vector of length 13")
            }
          results <- do.call(.CDMN8P, parlist);
          }
        else if(p==9)
          {
          if(length(dates) != 11)
            {
             stop("For a 1-fleet 9-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 9, final")
            }
          if(length(par) != 14)
            {
             stop("For a 1-fleet 9-perturbation model par must be a vector of length 14")
            }
          results <- do.call(.CDMN9P, parlist);
          }
        else if(p==10)
          {
          if(length(dates) != 12)
            {
             stop("For a 1-fleet 10-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 10, final")
            }
          if(length(par) != 15)
            {
             stop("For a 1-fleet 10-perturbation model par must be a vector of length 15")
            }
          results <- do.call(.CDMN10P, parlist);
          }
        else if(p==11)
          {
          if(length(dates) != 13)
            {
             stop("For a 1-fleet 11-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 11, final")
            }
          if(length(par) != 16)
            {
             stop("For a 1-fleet 11-perturbation model par must be a vector of length 16")
            }
          results <- do.call(.CDMN11P, parlist);
          }
        else if(p==12)
          {
          if(length(dates) != 14)
            {
             stop("For a 1-fleet 12-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 12, final")
            }
          if(length(par) != 17)
            {
             stop("For a 1-fleet 12-perturbation model par must be a vector of length 17")
            }
          results <- do.call(.CDMN12P, parlist);
          }
        else if(p==13)
          {
          if(length(dates) != 15)
            {
             stop("For a 1-fleet 13-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 13, final")
            }
          if(length(par) != 18)
            {
             stop("For a 1-fleet 13-perturbation model par must be a vector of length 18")
            }
          results <- do.call(.CDMN13P, parlist);
          }
        else if(p==14)
          {
          if(length(dates) != 16)
            {
             stop("For a 1-fleet 14-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 14, final")
            }
          if(length(par) != 19)
            {
             stop("For a 1-fleet 14-perturbation model par must be a vector of length 19")
            }
          results <- do.call(.CDMN14P, parlist);
          }
        else if(p==15)
          {
          if(length(dates) != 17)
            {
             stop("For a 1-fleet 15-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 15, final")
            }
          if(length(par) != 20)
            {
             stop("For a 1-fleet 15-perturbation model par must be a vector of length 20")
            }
          results <- do.call(.CDMN15P, parlist);
          }
        else if(p==16)
          {
          if(length(dates) != 18)
            {
             stop("For a 1-fleet 16-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 16, final")
            }
          if(length(par) != 21)
            {
             stop("For a 1-fleet 16-perturbation model par must be a vector of length 21")
            }
          results <- do.call(.CDMN16P, parlist);
          }
        else if(p==17)
          {
          if(length(dates) != 19)
            {
             stop("For a 1-fleet 17-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 17, final")
            }
          if(length(par) != 22)
            {
             stop("For a 1-fleet 17-perturbation model par must be a vector of length 22")
            }
          results <- do.call(.CDMN17P, parlist);
          }
        else if(p==18)
          {
          if(length(dates) != 20)
            {
             stop("For a 1-fleet 18-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 18, final")
            }
          if(length(par) != 23)
            {
             stop("For a 1-fleet 18-perturbation model par must be a vector of length 23")
            }
          results <- do.call(.CDMN18P, parlist);
          }
        else if(p==19)
          {
          if(length(dates) != 21)
            {
             stop("For a 1-fleet 19-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 19, final")
            }
          if(length(par) != 24)
            {
             stop("For a 1-fleet 19-perturbation model par must be a vector of length 24")
            }
          results <- do.call(.CDMN19P, parlist);
          }
        else if(p==20)
          {
          if(length(dates) != 22)
            {
             stop("For a 1-fleet 20-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 20, final")
            }
          if(length(par) != 25)
            {
             stop("For a 1-fleet 20-perturbation model par must be a vector of length 25")
            }
          results <- do.call(.CDMN20P, parlist);
          }
        else if(p==21)
          {
          if(length(dates) != 23)
            {
             stop("For a 1-fleet 23-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 21, final")
            }
          if(length(par) != 26)
            {
             stop("For a 1-fleet 21-perturbation model par must be a vector of length 26")
            }
          results <- do.call(.CDMN21P, parlist);
          }
        else if(p==22)
          {
          if(length(dates) != 24)
            {
             stop("For a 1-fleet 22-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 22, final")
            }
          if(length(par) != 27)
            {
             stop("For a 1-fleet 22-perturbation model par must be a vector of length 27")
            }
          results <- do.call(.CDMN22P, parlist);
          }
        else if(p==23)
          {
          if(length(dates) != 25)
            {
             stop("For a 1-fleet 23-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 23, final")
            }
          if(length(par) != 28)
            {
             stop("For a 1-fleet 23-perturbation model par must be a vector of length 28")
            }
          results <- do.call(.CDMN23P, parlist);
          }
        else if(p==24)
          {
          if(length(dates) != 26)
            {
             stop("For a 1-fleet 24-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 24, final")
            }
          if(length(par) != 29)
            {
             stop("For a 1-fleet 24-perturbation model par must be a vector of length 29")
            }
          results <- do.call(.CDMN24P, parlist);
          }
        else if(p==25)
          {
          if(length(dates) != 27)
            {
             stop("For a 1-fleet 25-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 25, final")
            }
          if(length(par) != 30)
            {
             stop("For a 1-fleet 25-perturbation model par must be a vector of length 30")
            }
          results <- do.call(.CDMN25P, parlist);
          }
        else if(p==26)
          {
          if(length(dates) != 28)
            {
             stop("For a 1-fleet 26-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 26, final")
            }
          if(length(par) != 31)
            {
             stop("For a 1-fleet 26-perturbation model par must be a vector of length 31")
            }
          results <- do.call(.CDMN26P, parlist);
          }
        else if(p==27)
          {
          if(length(dates) != 29)
            {
             stop("For a 1-fleet 27-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 27, final")
            }
          if(length(par) != 32)
            {
             stop("For a 1-fleet 27-perturbation model par must be a vector of length 32")
            }
          results <- do.call(.CDMN27P, parlist);
          }
        else if(p==28)
          {
          if(length(dates) != 30)
            {
             stop("For a 1-fleet 28-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 28, final")
            }
          if(length(par) != 33)
            {
             stop("For a 1-fleet 28-perturbation model par must be a vector of length 33")
            }
          results <- do.call(.CDMN28P, parlist);
          }
        else if(p==29)
          {
          if(length(dates) != 31)
            {
             stop("For a 1-fleet 28-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 29, final")
            }
          if(length(par) != 34)
            {
             stop("For a 1-fleet 29-perturbation model par must be a vector of length 34")
            }
          results <- do.call(.CDMN29P, parlist);
          }
        else if(p==30)
          {
          if(length(dates) != 32)
            {
             stop("For a 1-fleet 30-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 30, final")
            }
          if(length(par) != 35)
            {
             stop("For a 1-fleet 30-perturbation model par must be a vector of length 35")
            }
          results <- do.call(.CDMN30P, parlist);
          }
        else if(p==31)
          {
          if(length(dates) != 33)
            {
             stop("For a 1-fleet 31-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 31, final")
            }
          if(length(par) != 36)
            {
             stop("For a 1-fleet 31-perturbation model par must be a vector of length 36")
            }
          results <- do.call(.CDMN31P, parlist);
          }
        else if(p==32)
          {
          if(length(dates) != 34)
            {
             stop("For a 1-fleet 32-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 32, final")
            }
          if(length(par) != 37)
            {
             stop("For a 1-fleet 32-perturbation model par must be a vector of length 37")
            }
          results <- do.call(.CDMN32P, parlist);
          }
        else if(p==33)
          {
          if(length(dates) != 35)
            {
             stop("For a 1-fleet 33-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 33, final")
            }
          if(length(par) != 38)
            {
             stop("For a 1-fleet 33-perturbation model par must be a vector of length 38")
            }
          results <- do.call(.CDMN33P, parlist);
          }
        else if(p==34)
          {
          if(length(dates) != 36)
            {
             stop("For a 1-fleet 34-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 34, final")
            }
          if(length(par) != 39)
            {
             stop("For a 1-fleet 34-perturbation model par must be a vector of length 39")
            }
          results <- do.call(.CDMN34P, parlist);
          }
        else if(p==35)
          {
          if(length(dates) != 37)
            {
             stop("For a 1-fleet 35-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 35, final")
            }
          if(length(par) != 40)
            {
             stop("For a 1-fleet 35-perturbation model par must be a vector of length 40")
            }
          results <- do.call(.CDMN35P, parlist);
          }
        else if(p==36)
          {
          if(length(dates) != 38)
            {
             stop("For a 1-fleet 36-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 36, final")
            }
          if(length(par) != 41)
            {
             stop("For a 1-fleet 36-perturbation model par must be a vector of length 41")
            }
          results <- do.call(.CDMN36P, parlist);
          }
        else if(p==37)
          {
          if(length(dates) != 39)
            {
             stop("For a 1-fleet 37-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 37, final")
            }
          if(length(par) != 42)
            {
             stop("For a 1-fleet 37-perturbation model par must be a vector of length 42")
            }
          results <- do.call(.CDMN37P, parlist);
          }
        else if(p==38)
          {
          if(length(dates) != 40)
            {
             stop("For a 1-fleet 38-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 38, final")
            }
          if(length(par) != 43)
            {
             stop("For a 1-fleet 38-perturbation model par must be a vector of length 43")
            }
          results <- do.call(.CDMN38P, parlist);
          }
        else if(p==39)
          {
          if(length(dates) != 41)
            {
             stop("For a 1-fleet 39-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 39, final")
            }
          if(length(par) != 44)
            {
             stop("For a 1-fleet 39-perturbation model par must be a vector of length 44")
            }
          results <- do.call(.CDMN39P, parlist);
          }
        else if(p==40)
          {
          if(length(dates) != 42)
            {
             stop("For a 1-fleet 40-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 40, final")
            }
          if(length(par) != 45)
            {
             stop("For a 1-fleet 40-perturbation model par must be a vector of length 45")
            }
          results <- do.call(.CDMN40P, parlist);
          }
        else if(p==41)
          {
          if(length(dates) != 43)
            {
             stop("For a 1-fleet 41-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 41, final")
            }
          if(length(par) != 46)
            {
             stop("For a 1-fleet 41-perturbation model par must be a vector of length 46")
            }
          results <- do.call(.CDMN41P, parlist);
          }
        else if(p==42)
          {
          if(length(dates) != 44)
            {
             stop("For a 1-fleet 42-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 42, final")
            }
          if(length(par) != 47)
            {
             stop("For a 1-fleet 42-perturbation model par must be a vector of length 47")
            }
          results <- do.call(.CDMN42P, parlist);
          }
        else if(p==43)
          {
          if(length(dates) != 45)
            {
             stop("For a 1-fleet 43-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 43, final")
            }
          if(length(par) != 48)
            {
             stop("For a 1-fleet 43-perturbation model par must be a vector of length 48")
            }
          results <- do.call(.CDMN43P, parlist);
          }
        else if(p==44)
          {
          if(length(dates) != 46)
            {
             stop("For a 1-fleet 44-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 44, final")
            }
          if(length(par) != 49)
            {
             stop("For a 1-fleet 44-perturbation model par must be a vector of length 49")
            }
          results <- do.call(.CDMN44P, parlist);
          }
        else if(p==45)
          {
          if(length(dates) != 47)
            {
             stop("For a 1-fleet 45-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 45, final")
            }
          if(length(par) != 50)
            {
             stop("For a 1-fleet 45-perturbation model par must be a vector of length 50")
            }
          results <- do.call(.CDMN45P, parlist);
          }
        else if(p==46)
          {
          if(length(dates) != 48)
            {
             stop("For a 1-fleet 45-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 46, final")
            }
          if(length(par) != 51)
            {
             stop("For a 1-fleet 46-perturbation model par must be a vector of length 51")
            }
          results <- do.call(.CDMN46P, parlist);
          }
        else if(p==47)
          {
          if(length(dates) != 49)
            {
             stop("For a 1-fleet 47-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 47, final")
            }
          if(length(par) != 52)
            {
             stop("For a 1-fleet 47-perturbation model par must be a vector of length 52")
            }
          results <- do.call(.CDMN47P, parlist);
          }
        else if(p==48)
          {
          if(length(dates) != 50)
            {
             stop("For a 1-fleet 48-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 48, final")
            }
          if(length(par) != 53)
            {
             stop("For a 1-fleet 48-perturbation model par must be a vector of length 53")
            }
          results <- do.call(.CDMN48P, parlist);
          }
        else if(p==49)
          {
          if(length(dates) != 51)
            {
             stop("For a 1-fleet 49-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 49, final")
            }
          if(length(par) != 54)
            {
             stop("For a 1-fleet 49-perturbation model par must be a vector of length 54")
            }
          results <- do.call(.CDMN49P, parlist);
          }
        else if(p==50)
          {
          if(length(dates) != 52)
            {
             stop("For a 1-fleet 50-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 50, final")
            }
          if(length(par) != 55)
            {
             stop("For a 1-fleet 50-perturbation model par must be a vector of length 55")
            }
          results <- do.call(.CDMN50P, parlist);
          }
        else if(p==51)
          {
          if(length(dates) != 53)
            {
             stop("For a 1-fleet 51-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 51, final")
            }
          if(length(par) != 56)
            {
             stop("For a 1-fleet 51-perturbation model par must be a vector of length 56")
            }
          results <- do.call(.CDMN51P, parlist);
          }
        else if(p==52)
          {
          if(length(dates) != 54)
            {
             stop("For a 1-fleet 52-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 52, final")
            }
          if(length(par) != 57)
            {
             stop("For a 1-fleet 52-perturbation model par must be a vector of length 57")
            }
          results <- do.call(.CDMN52P, parlist);
          }
        else if(p==53)
          {
          if(length(dates) != 55)
            {
             stop("For a 1-fleet 53-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 53, final")
            }
          if(length(par) != 58)
            {
             stop("For a 1-fleet 53-perturbation model par must be a vector of length 58")
            }
          results <- do.call(.CDMN53P, parlist);
          }
        else if(p==54)
          {
          if(length(dates) != 56)
            {
             stop("For a 1-fleet 54-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 54, final")
            }
          if(length(par) != 59)
            {
             stop("For a 1-fleet 54-perturbation model par must be a vector of length 59")
            }
          results <- do.call(.CDMN54P, parlist);
          }
        else if(p==55)
          {
          if(length(dates) != 57)
            {
             stop("For a 1-fleet 55-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 55, final")
            }
          if(length(par) != 60)
            {
             stop("For a 1-fleet 55-perturbation model par must be a vector of length 60")
            }
          results <- do.call(.CDMN55P, parlist);
          }
        else if(p==56)
          {
          if(length(dates) != 58)
            {
             stop("For a 1-fleet 56-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 56, final")
            }
          if(length(par) != 61)
            {
             stop("For a 1-fleet 56-perturbation model par must be a vector of length 61")
            }
          results <- do.call(.CDMN56P, parlist);
          }
        else if(p==57)
          {
          if(length(dates) != 59)
            {
             stop("For a 1-fleet 57-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 57, final")
            }
          if(length(par) != 62)
            {
             stop("For a 1-fleet 57-perturbation model par must be a vector of length 62")
            }
          results <- do.call(.CDMN57P, parlist);
          }
        else if(p==58)
          {
          if(length(dates) != 60)
            {
             stop("For a 1-fleet 58-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 58, final")
            }
          if(length(par) != 63)
            {
             stop("For a 1-fleet 58-perturbation model par must be a vector of length 68")
            }
          results <- do.call(.CDMN58P, parlist);
          }
        else if(p==59)
          {
          if(length(dates) != 61)
            {
             stop("For a 1-fleet 59-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 59, final")
            }
          if(length(par) != 64)
            {
             stop("For a 1-fleet 59-perturbation model par must be a vector of length 64")
            }
          results <- do.call(.CDMN59P, parlist);
          }
        else if(p==60)
          {
          if(length(dates) != 62)
            {
             stop("For a 1-fleet 60-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 60, final")
            }
          if(length(par) != 65)
            {
             stop("For a 1-fleet 60-perturbation model par must be a vector of length 65")
            }
          results <- do.call(.CDMN60P, parlist);
          }
        else if(p==61)
          {
          if(length(dates) != 63)
            {
             stop("For a 1-fleet 61-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 61, final")
            }
          if(length(par) != 66)
            {
             stop("For a 1-fleet 61-perturbation model par must be a vector of length 66")
            }
          results <- do.call(.CDMN61P, parlist);
          }
        else if(p==62)
          {
          if(length(dates) != 64)
            {
             stop("For a 1-fleet 62-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 62, final")
            }
          if(length(par) != 67)
            {
             stop("For a 1-fleet 62-perturbation model par must be a vector of length 67")
            }
          results <- do.call(.CDMN62P, parlist);
          }
        else if(p==63)
          {
          if(length(dates) != 65)
            {
             stop("For a 1-fleet 63-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 63, final")
            }
          if(length(par) != 68)
            {
             stop("For a 1-fleet 63-perturbation model par must be a vector of length 68")
            }
          results <- do.call(.CDMN63P, parlist);
          }
        else if(p==64)
          {
          if(length(dates) != 66)
            {
             stop("For a 1-fleet 64-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 64, final")
            }
          if(length(par) != 69)
            {
             stop("For a 1-fleet 64-perturbation model par must be a vector of length 69")
            }
          results <- do.call(.CDMN64P, parlist);
          }
        else if(p==65)
          {
          if(length(dates) != 67)
            {
             stop("For a 1-fleet 65-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 65, final")
            }
          if(length(par) != 70)
            {
             stop("For a 1-fleet 65-perturbation model par must be a vector of length 70")
            }
          results <- do.call(.CDMN65P, parlist);
          }
        else if(p==66)
          {
          if(length(dates) != 68)
            {
             stop("For a 1-fleet 66-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 66, final")
            }
          if(length(par) != 71)
            {
             stop("For a 1-fleet 66-perturbation model par must be a vector of length 71")
            }
          results <- do.call(.CDMN66P, parlist);
          }
        else if(p==67)
          {
          if(length(dates) != 69)
            {
             stop("For a 1-fleet 67-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 67, final")
            }
          if(length(par) != 72)
            {
             stop("For a 1-fleet 67-perturbation model par must be a vector of length 72")
            }
          results <- do.call(.CDMN67P, parlist);
          }
        else if(p==68)
          {
          if(length(dates) != 70)
            {
             stop("For a 1-fleet 68-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 68, final")
            }
          if(length(par) != 73)
            {
             stop("For a 1-fleet 68-perturbation model par must be a vector of length 73")
            }
          results <- do.call(.CDMN68P, parlist);
          }
        else if(p==69)
          {
          if(length(dates) != 71)
            {
             stop("For a 1-fleet 69-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 69, final")
            }
          if(length(par) != 74)
            {
             stop("For a 1-fleet 69-perturbation model par must be a vector of length 74")
            }
          results <- do.call(.CDMN69P, parlist);
          }
        else if(p==70)
          {
          if(length(dates) != 72)
            {
             stop("For a 1-fleet 70-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 70, final")
            }
          if(length(par) != 75)
            {
             stop("For a 1-fleet 70-perturbation model par must be a vector of length 75")
            }
          results <- do.call(.CDMN70P, parlist);
          }
        else if(p==71)
          {
          if(length(dates) != 73)
            {
             stop("For a 1-fleet 71-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 71, final")
            }
          if(length(par) != 76)
            {
             stop("For a 1-fleet 71-perturbation model par must be a vector of length 76")
            }
          results <- do.call(.CDMN71P, parlist);
          }
        else if(p==72)
          {
          if(length(dates) != 74)
            {
             stop("For a 1-fleet 72-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 72, final")
            }
          if(length(par) != 77)
            {
             stop("For a 1-fleet 72-perturbation model par must be a vector of length 77")
            }
          results <- do.call(.CDMN72P, parlist);
          }
        else if(p==73)
          {
          if(length(dates) != 75)
            {
             stop("For a 1-fleet 73-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 73, final")
            }
          if(length(par) != 78)
            {
             stop("For a 1-fleet 73-perturbation model par must be a vector of length 78")
            }
          results <- do.call(.CDMN73P, parlist);
          }
        else if(p==74)
          {
          if(length(dates) != 76)
            {
             stop("For a 1-fleet 74-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 74, final")
            }
          if(length(par) != 79)
            {
             stop("For a 1-fleet 74-perturbation model par must be a vector of length 79")
            }
          results <- do.call(.CDMN74P, parlist);
          }
        else if(p==75)
          {
          if(length(dates) != 77)
            {
             stop("For a 1-fleet 75-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 75, final")
            }
          if(length(par) != 80)
            {
             stop("For a 1-fleet 75-perturbation model par must be a vector of length 80")
            }
          results <- do.call(.CDMN75P, parlist);
          }
        else if(p==76)
          {
          if(length(dates) != 78)
            {
             stop("For a 1-fleet 76-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 76, final")
            }
          if(length(par) != 81)
            {
             stop("For a 1-fleet 76-perturbation model par must be a vector of length 81")
            }
          results <- do.call(.CDMN76P, parlist);
          }
        else if(p==77)
          {
          if(length(dates) != 79)
            {
             stop("For a 1-fleet 77-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 77, final")
            }
          if(length(par) != 82)
            {
             stop("For a 1-fleet 77-perturbation model par must be a vector of length 82")
            }
          results <- do.call(.CDMN77P, parlist);
          }
        else if(p==78)
          {
          if(length(dates) != 80)
            {
             stop("For a 1-fleet 78-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 78, final")
            }
          if(length(par) != 83)
            {
             stop("For a 1-fleet 78-perturbation model par must be a vector of length 83")
            }
          results <- do.call(.CDMN78P, parlist);
          }
        else if(p==79)
          {
          if(length(dates) != 81)
            {
             stop("For a 1-fleet 79-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 79, final")
            }
          if(length(par) != 84)
            {
             stop("For a 1-fleet 79-perturbation model par must be a vector of length 84")
            }
          results <- do.call(.CDMN79P, parlist);
          }
        else if(p==80)
          {
          if(length(dates) != 82)
            {
             stop("For a 1-fleet 80-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 80, final")
            }
          if(length(par) != 85)
            {
             stop("For a 1-fleet 80-perturbation model par must be a vector of length 85")
            }
          results <- do.call(.CDMN80P, parlist);
          }
        else if(p==81)
          {
          if(length(dates) != 83)
            {
             stop("For a 1-fleet 81-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 81, final")
            }
          if(length(par) != 86)
            {
             stop("For a 1-fleet 81-perturbation model par must be a vector of length 86")
            }
          results <- do.call(.CDMN81P, parlist);
          }
        else if(p==82)
          {
          if(length(dates) != 84)
            {
             stop("For a 1-fleet 82-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 82, final")
            }
          if(length(par) != 87)
            {
             stop("For a 1-fleet 82-perturbation model par must be a vector of length 87")
            }
          results <- do.call(.CDMN82P, parlist);
          }
        else if(p==83)
          {
          if(length(dates) != 85)
            {
             stop("For a 1-fleet 83-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 83, final")
            }
          if(length(par) != 88)
            {
             stop("For a 1-fleet 83-perturbation model par must be a vector of length 88")
            }
          results <- do.call(.CDMN83P, parlist);
          }
        else if(p==84)
          {
          if(length(dates) != 86)
            {
             stop("For a 1-fleet 84-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 84, final")
            }
          if(length(par) != 89)
            {
             stop("For a 1-fleet 84-perturbation model par must be a vector of length 89")
            }
          results <- do.call(.CDMN84P, parlist);
          }
        else if(p==85)
          {
          if(length(dates) != 87)
            {
             stop("For a 1-fleet 85-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 85, final")
            }
          if(length(par) != 90)
            {
             stop("For a 1-fleet 85-perturbation model par must be a vector of length 90")
            }
          results <- do.call(.CDMN85P, parlist);
          }
        else if(p==86)
          {
          if(length(dates) != 88)
            {
             stop("For a 1-fleet 86-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 86, final")
            }
          if(length(par) != 91)
            {
             stop("For a 1-fleet 86-perturbation model par must be a vector of length 91")
            }
          results <- do.call(.CDMN86P, parlist);
          }
        else if(p==87)
          {
          if(length(dates) != 89)
            {
             stop("For a 1-fleet 87-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 87, final")
            }
          if(length(par) != 92)
            {
             stop("For a 1-fleet 87-perturbation model par must be a vector of length 92")
            }
          results <- do.call(.CDMN87P, parlist);
          }
        else if(p==88)
          {
          if(length(dates) != 90)
            {
             stop("For a 1-fleet 88-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 88, final")
            }
          if(length(par) != 93)
            {
             stop("For a 1-fleet 88-perturbation model par must be a vector of length 93")
            }
          results <- do.call(.CDMN88P, parlist);
          }
        else if(p==89)
          {
          if(length(dates) != 91)
            {
             stop("For a 1-fleet 89-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 89, final")
            }
          if(length(par) != 94)
            {
             stop("For a 1-fleet 89-perturbation model par must be a vector of length 94")
            }
          results <- do.call(.CDMN89P, parlist);
          }
        else if(p==90)
          {
          if(length(dates) != 92)
            {
             stop("For a 1-fleet 90-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 90, final")
            }
          if(length(par) != 95)
            {
             stop("For a 1-fleet 90-perturbation model par must be a vector of length 95")
            }
          results <- do.call(.CDMN90P, parlist);
          }
        else if(p==91)
          {
          if(length(dates) != 93)
            {
             stop("For a 1-fleet 91-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 91, final")
            }
          if(length(par) != 96)
            {
             stop("For a 1-fleet 91-perturbation model par must be a vector of length 96")
            }
          results <- do.call(.CDMN91P, parlist);
          }
        else if(p==92)
          {
          if(length(dates) != 94)
            {
             stop("For a 1-fleet 92-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 92, final")
            }
          if(length(par) != 97)
            {
             stop("For a 1-fleet 92-perturbation model par must be a vector of length 97")
            }
          results <- do.call(.CDMN92P, parlist);
          }
        else if(p==93)
          {
          if(length(dates) != 95)
            {
             stop("For a 1-fleet 93-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 93, final")
            }
          if(length(par) != 98)
            {
             stop("For a 1-fleet 93-perturbation model par must be a vector of length 98")
            }
          results <- do.call(.CDMN93P, parlist);
          }
        else if(p==94)
          {
          if(length(dates) != 96)
            {
             stop("For a 1-fleet 94-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 94, final")
            }
          if(length(par) != 99)
            {
             stop("For a 1-fleet 94-perturbation model par must be a vector of length 99")
            }
          results <- do.call(.CDMN94P, parlist);
          }
        else if(p==95)
          {
          if(length(dates) != 97)
            {
             stop("For a 1-fleet 95-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 95, final")
            }
          if(length(par) != 100)
            {
             stop("For a 1-fleet 95-perturbation model par must be a vector of length 100")
            }
          results <- do.call(.CDMN95P, parlist);
          }
        else if(p==96)
          {
          if(length(dates) != 98)
            {
             stop("For a 1-fleet 96-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 96, final")
            }
          if(length(par) != 101)
            {
             stop("For a 1-fleet 96-perturbation model par must be a vector of length 101")
            }
          results <- do.call(.CDMN96P, parlist);
          }
        else if(p==97)
          {
          if(length(dates) != 99)
            {
             stop("For a 1-fleet 97-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 97, final")
            }
          if(length(par) != 102)
            {
             stop("For a 1-fleet 97-perturbation model par must be a vector of length 102")
            }
          results <- do.call(.CDMN97P, parlist);
          }
        else if(p==98)
          {
          if(length(dates) != 100)
            {
             stop("For a 1-fleet 98-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 98, final")
            }
          if(length(par) != 103)
            {
             stop("For a 1-fleet 98-perturbation model par must be a vector of length 103")
            }
          results <- do.call(.CDMN98P, parlist);
          }
        else if(p==99)
          {
          if(length(dates) != 101)
            {
             stop("For a 1-fleet 99-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 99, final")
            }
          if(length(par) != 104)
            {
             stop("For a 1-fleet 99-perturbation model par must be a vector of length 104")
            }
          results <- do.call(.CDMN99P, parlist);
          }
        else if(p==100)
          {
          if(length(dates) != 102)
            {
             stop("For a 1-fleet 100-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 100, final")
            }
          if(length(par) != 105)
            {
             stop("For a 1-fleet 100-perturbation model par must be a vector of length 105")
            }
          results <- do.call(.CDMN100P, parlist);
          }
        else if(p==-1)
          {
          if(length(dates) != 4)
            {
             stop("For a 1-fleet 1-perturbation transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, final")
            }
          if(length(par) != 6)
            {
             stop("For a 1-fleet 1-perturbation transit model par must be a vector of length 6")
            }
          results <- do.call(.CDMNT1P, parlist);
          }
        else if(p==-2)
          {
          if(length(dates) != 6)
            {
             stop("For a 1-fleet 2-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, perturbation 2, exit 2, final")
            }
          if(length(par) != 7)
            {
             stop("For a 1-fleet 2-perturbation transit model par must be a vector of length 7")
            }
          results <- do.call(.CDMNT2P, parlist);
          }
        else if(p==-3)
          {
          if(length(dates) != 8)
            {
             stop("For a 1-fleet 3-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, perturbation 2, exit 2, perturbation 3, exit 3, final")
            }
          if(length(par) != 8)
            {
             stop("For a 1-fleet 3-perturbation transit model par must be a vector of length 8")
            }
          results <- do.call(.CDMNT3P, parlist);
          }
        else if(p==-4)
          {
          if(length(dates) != 10)
            {
             stop("For a 1-fleet 4-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 4, exit 4, final")
            }
          if(length(par) != 9)
            {
             stop("For a 1-fleet 4-perturbation transit model par must be a vector of length 9")
            }
          results <- do.call(.CDMNT4P, parlist);
          }
        else if(p==-5)
          {
          if(length(dates) != 12)
            {
             stop("For a 1-fleet 5-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 5, exit 6, final")
            }
          if(length(par) != 10)
            {
             stop("For a 1-fleet 5-perturbation transit model par must be a vector of length 10")
            }
          results <- do.call(.CDMNT5P, parlist);
          }
        else if(p==-6)
          {
          if(length(dates) != 14)
            {
             stop("For a 1-fleet 6-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 6, exit 6, final")
            }
          if(length(par) != 11)
            {
             stop("For a 1-fleet 6-perturbation transit model par must be a vector of length 11")
            }
          results <- do.call(.CDMNT6P, parlist);
          }
        else if(p==-7)
          {
          if(length(dates) != 16)
            {
             stop("For a 1-fleet 7-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 7, exit 7, final")
            }
          if(length(par) != 12)
            {
             stop("For a 1-fleet 7-perturbation transit model par must be a vector of length 12")
            }
          results <- do.call(.CDMNT7P, parlist);
          }
        else if(p==-8)
          {
          if(length(dates) != 18)
            {
             stop("For a 1-fleet 8-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 8, exit 8, final")
            }
          if(length(par) != 13)
            {
             stop("For a 1-fleet 8-perturbation transit model par must be a vector of length 13")
            }
          results <- do.call(.CDMNT8P, parlist);
          }
        else if(p==-9)
          {
          if(length(dates) != 20)
            {
             stop("For a 1-fleet 9-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 9, exit 9, final")
            }
          if(length(par) != 14)
            {
             stop("For a 1-fleet 9-perturbation transit model par must be a vector of length 14")
            }
          results <- do.call(.CDMNT9P, parlist);
          }
        else if(p==-10)
          {
          if(length(dates) != 22)
            {
             stop("For a 1-fleet 10-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 10, exit 10, final")
            }
          if(length(par) != 15)
            {
             stop("For a 1-fleet 10-perturbation transit model par must be a vector of length 15")
            }
          results <- do.call(.CDMNT10P, parlist);
          }
        else if(p==-11)
          {
          if(length(dates) != 24)
            {
             stop("For a 1-fleet 11-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 11, exit 11, final")
            }
          if(length(par) != 16)
            {
             stop("For a 1-fleet 11-perturbation transit model par must be a vector of length 16")
            }
          results <- do.call(.CDMNT11P, parlist);
          }
        else if(p==-12)
          {
          if(length(dates) != 26)
            {
             stop("For a 1-fleet 12-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 12, exit 12, final")
            }
          if(length(par) != 17)
            {
             stop("For a 1-fleet 12-perturbation transit model par must be a vector of length 17")
            }
          results <- do.call(.CDMNT12P, parlist);
          }
        else if(p==-13)
          {
          if(length(dates) != 28)
            {
             stop("For a 1-fleet 13-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 13, exit 13, final")
            }
          if(length(par) != 18)
            {
             stop("For a 1-fleet 13-perturbation transit model par must be a vector of length 18")
            }
          results <- do.call(.CDMNT13P, parlist);
          }
        else if(p==-14)
          {
          if(length(dates) != 30)
            {
             stop("For a 1-fleet 14-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 14, exit 14, final")
            }
          if(length(par) != 19)
            {
             stop("For a 1-fleet 14-perturbation transit model par must be a vector of length 19")
            }
          results <- do.call(.CDMNT14P, parlist);
          }
        else if(p==-15)
          {
          if(length(dates) != 32)
            {
             stop("For a 1-fleet 15-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 15, exit 15, final")
            }
          if(length(par) != 20)
            {
             stop("For a 1-fleet 15-perturbation transit model par must be a vector of length 20")
            }
          results <- do.call(.CDMNT15P, parlist);
          }
        else if(p==-16)
          {
          if(length(dates) != 34)
            {
             stop("For a 1-fleet 16-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 16, exit 16, final")
            }
          if(length(par) != 21)
            {
             stop("For a 1-fleet 16-perturbation transit model par must be a vector of length 21")
            }
          results <- do.call(.CDMNT16P, parlist);
          }
        else if(p==-17)
          {
          if(length(dates) != 36)
            {
             stop("For a 1-fleet 17-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 17, exit 17, final")
            }
          if(length(par) != 22)
            {
             stop("For a 1-fleet 17-perturbation transit model par must be a vector of length 22")
            }
          results <- do.call(.CDMNT17P, parlist);
          }
        else if(p==-18)
          {
          if(length(dates) != 38)
            {
             stop("For a 1-fleet 18-perturbations model transit dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 18, exit 18, final")
            }
          if(length(par) != 23)
            {
             stop("For a 1-fleet 18-perturbation transit model par must be a vector of length 23")
            }
          results <- do.call(.CDMNT18P, parlist);
          }
        else if(p==-19)
          {
          if(length(dates) != 40)
            {
             stop("For a 1-fleet 19-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 19, exit 18, final")
            }
          if(length(par) != 24)
            {
             stop("For a 1-fleet 19-perturbation transit model par must be a vector of length 24")
            }
          results <- do.call(.CDMNT19P, parlist);
          }
        else if(p==-20)
          {
          if(length(dates) != 42)
            {
             stop("For a 1-fleet 20-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 20, exit 18, final")
            }
          if(length(par) != 25)
            {
             stop("For a 1-fleet 20-perturbation transit model par must be a vector of length 25")
            }
          results <- do.call(.CDMNT20P, parlist);
          }
       }
     else 
       {
        if(!distr%in%distr.set)
          {stop("distr must be 'normal','apnormal','lognormal','aplnormal','gamma', 'poisson' or 'negbin', see help pages for CatDynFit")}
        parlist <- list(par=par,
                        dates=dates, 
                        obseff1=x$Data[[fleet.name[1]]][,2], 
                        obscat1=x$Data[[fleet.name[1]]][,5], 
                        obsmbm1=x$Data[[fleet.name[1]]][,4], 
                        obseff2=x$Data[[fleet.name[2]]][,2], 
                        obscat2=x$Data[[fleet.name[2]]][,5], 
                        obsmbm2=x$Data[[fleet.name[2]]][,4],
                        distr=distr, 
                        properties=x$Properties,
                        output="predict");
        if(sum(p == c(0,0)) == length(p))
          {
          if(length(dates) != 2)
            {
             stop("For a 2-fleet 0-perturbation 0-perturbation model dates must be a vector with the following time step marks: initial, and final time steps in the season")
            }
          if(length(par) != 8)
            {
             stop("For a 0-fleet 0-perturbation 0-perturbation model par must be a vector of length 8")
            }
          results <- do.call(.CDMN0P0P, parlist);
          }
        else if(sum(p == c(0,1)) == length(p))
          {
          if(length(dates) != 3)
            {
             stop("For a 2-fleet 0-perturbation 1-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, final")
            }
          if(length(par) != 9)
            {
             stop("For a 2-fleet 0-perturbation 1-perturbation model par must be a vector of length 9")
            }
          results <- do.call(.CDMN0P1P, parlist);
          }
        else if(sum(p == c(0,2)) == length(p))
          {
          if(length(dates) != 4)
            {
             stop("For a 2-fleet 0-perturbation 2-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, fleet 2 perturbation 2, final")
            }
          if(length(par) != 10)
            {
             stop("For a 2-fleet 0-perturbation 2-perturbation model par must be a vector of length 10")
            }
          results <- do.call(.CDMN0P2P, parlist);
          }
        else if(sum(p == c(0,3)) == length(p))
          {
          if(length(dates) != 5)
            {
             stop("For a 2-fleet 0-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, fleet 2 perturbation 2, fleet 2 perturbation 3, final")
            }
          if(length(par) != 11)
            {
             stop("For a 2-fleet 0-perturbation 3-perturbation model par must be a vector of length 11")
            }
          results <- do.call(.CDMN0P3P, parlist);
          }
        else if(sum(p == c(0,4)) == length(p))
          {
          if(length(dates) != 6)
            {
             stop("For a 2-fleet 0-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(length(par) != 12)
            {
             stop("For a 2-fleet 0-perturbation 4-perturbation model par must be a vector of length 12")
            }
          results <- do.call(.CDMN0P4P, parlist);
          }
        else if(sum(p == c(0,5)) == length(p))
          {
          if(length(dates) != 7)
            {
             stop("For a 2-fleet 0-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(length(par) != 13)
            {
             stop("For a 2-fleet 0-perturbation 5-perturbation model par must be a vector of length 13")
            }
          results <- do.call(.CDMN0P5P, parlist);
          }
        else if(sum(p == c(1,1)) == length(p))
          {
          if(length(dates) != 4)
            {
             stop("For a 2-fleet 1-perturbation 1-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, final")
            }
          if(length(par) != 10)
            {
             stop("For a 2-fleet 1-perturbation 1-perturbation model par must be a vector of length 10")
            }
          results <- do.call(.CDMN1P1P, parlist);
          }                                                                                                                                            
        else if(sum(p == c(1,2)) == length(p))
          {
          if(length(dates) != 5)
            {
             stop("For a 2-fleet 1-perturbation 2-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, fleet 2 perturbation 2, final")
            }
          if(length(par) != 11)
            {
             stop("For a 2-fleet 1-perturbation 2-perturbation model par must be a vector of length 11")
            }
          results <- do.call(.CDMN1P2P, parlist);
          }
        else if(sum(p == c(1,3)) == length(p))
          {
          if(length(dates) != 6)
            {
             stop("For a 2-fleet 1-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, ..., fleet 2 perturbation 3, final")
            }
          if(length(par) != 12)
            {
             stop("For a 2-fleet 1-perturbation 3-perturbation model par must be a vector of length 12")
            }
          results <- do.call(.CDMN1P3P, parlist);
          }
        else if(sum(p == c(1,4)) == length(p))
          {
          if(length(dates) != 7)
            {
             stop("For a 2-fleet 1-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(length(par) != 13)
            {
             stop("For a 2-fleet 1-perturbation 4-perturbation model par must be a vector of length 13")
            }
          results <- do.call(.CDMN1P4P, parlist);
          }
        else if(sum(p == c(1,5)) == length(p))
          {
          if(length(dates) != 8)
            {
             stop("For a 2-fleet 1-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(length(par) != 14)
            {
             stop("For a 2-fleet 1-perturbation 5-perturbation model par must be a vector of length 14")
            }
          results <- do.call(.CDMN1P5P, parlist);
          }
        else if(sum(p == c(2,2)) == length(p))
          {
          if(length(dates) != 6)
            {
             stop("For a 2-fleet 2-perturbation 2-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, fleet 2 perturbation 2, final")
            }
          if(length(par) != 12)
            {
             stop("For a 2-fleet 2-perturbation 2-perturbation model par must be a vector of length 12")
            }
          results <- do.call(.CDMN2P2P, parlist);
          }
        else if(sum(p == c(2,3)) == length(p))
          {
          if(length(dates) != 7)
            {
             stop("For a 2-fleet 2-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, ..., fleet 2 perturbation 3, final")
            }
          if(length(par) != 13)
            {
             stop("For a 2-fleet 2-perturbation 3-perturbation model par must be a vector of length 13")
            }
          results <- do.call(.CDMN2P3P, parlist);
          }
        else if(sum(p == c(2,4)) == length(p))
          {
          if(length(dates) != 8)
            {
             stop("For a 2-fleet 2-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(length(par) != 14)
            {
             stop("For a 2-fleet 2-perturbation 4-perturbation model par must be a vector of length 14")
            }
          results <- do.call(.CDMN2P4P, parlist);
          }
        else if(sum(p == c(2,5)) == length(p))
          {
          if(length(dates) != 9)
            {
             stop("For a 2-fleet 2-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(length(par) != 15)
            {
             stop("For a 2-fleet 2-perturbation 5-perturbation model par must be a vector of length 15")
            }
          results <- do.call(.CDMN2P5P, parlist);
          }
        else if(sum(p == c(3,3)) == length(p))
          {
          if(length(dates) != 8)
            {
             stop("For a 2-fleet 3-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, ..., fleet 2 perturbation 2, fleet 2 perturbation 3, final")
            }
          if(length(par) != 14)
            {
             stop("For a 2-fleet 3-perturbation 3-perturbation model par must be a vector of length 14")
            }
          results <- do.call(.CDMN3P3P, parlist);
          }
        else if(sum(p == c(3,4)) == length(p))
          {
          if(length(dates) != 9)
            {
             stop("For a 2-fleet 3-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(length(par) != 15)
            {
             stop("For a 2-fleet 3-perturbation 4-perturbation model par must be a vector of length 15")
            }
          results <- do.call(.CDMN3P4P, parlist);
          }
        else if(sum(p == c(3,5)) == length(p))
          {
          if(length(dates) != 10)
            {
             stop("For a 2-fleet 3-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(length(par) != 16)
            {
             stop("For a 2-fleet 3-perturbation 5-perturbation model par must be a vector of length 16")
            }
          results <- do.call(.CDMN3P5P, parlist);
          }
        else if(sum(p == c(4,4)) == length(p))
          {
          if(length(dates) != 10)
            {
             stop("For a 2-fleet 4-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, fleet 1 perturbation 4, fleet 2 perturbation 1, fleet 2 perturbation 2, fleet 2 perturbation 3, fleet 2 perturbation 4, final")
            }
          if(length(par) != 16)
            {
             stop("For a 2-fleet 4-perturbation 4-perturbation model par must be a vector of length 16")
            }
          results <- do.call(.CDMN4P4P, parlist);
          }
        else if(sum(p == c(4,5)) == length(p))
          {
          if(length(dates) != 11)
            {
             stop("For a 2-fleet 4-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 4, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(length(par) != 17)
            {
             stop("For a 2-fleet 4-perturbation 5-perturbation model par must be a vector of length 17")
            }
          results <- do.call(.CDMN4P5P, parlist);
          }
        else if(sum(p == c(5,5)) == length(p))
          {
          if(length(dates) != 12)
            {
             stop("For a 2-fleet 5-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 5, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(length(par) != 18)
            {
             stop("For a 2-fleet 5-perturbation 5-perturbation model par must be a vector of length 18")
            }
          results <- do.call(.CDMN5P5P, parlist);
          }
        else if(sum(p == c(6,6)) == length(p))
          {
          if(length(dates) != 14)
            {
             stop("For a 2-fleet 6-perturbation 6-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 6, fleet 2 perturbation 1, ..., fleet 2 perturbation 6, final")
            }
          if(length(par) != 20)
            {
             stop("For a 2-fleet 6-perturbation 6-perturbation model par must be a vector of length 20")
            }
          results <- do.call(.CDMN6P6P, parlist);
          }
        else if(sum(p == c(7,7)) == length(p))
          {
          if(length(dates) != 16)
            {
             stop("For a 2-fleet 7-perturbation 7-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 7, fleet 2 perturbation 1, ..., fleet 2 perturbation 7, final")
            }
          if(length(par) != 22)
            {
             stop("For a 2-fleet 7-perturbation 7-perturbation model par must be a vector of length 22")
            }
          results <- do.call(.CDMN7P7P, parlist);
          }
        else if(sum(p == c(8,8)) == length(p))
          {
          if(length(dates) != 18)
            {
             stop("For a 2-fleet 8-perturbation 8-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 8, fleet 2 perturbation 1, ..., fleet 2 perturbation 8, final")
            }
          if(length(par) != 24)
            {
             stop("For a 2-fleet 8-perturbation 8-perturbation model par must be a vector of length 24")
            }
          results <- do.call(.CDMN8P8P, parlist);
          }
        else if(sum(p == c(9,9)) == length(p))
          {
          if(length(dates) != 20)
            {
             stop("For a 2-fleet 9-perturbation 9-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 9, fleet 2 perturbation 1, ..., fleet 2 perturbation 9, final")
            }
          if(length(par) != 26)
            {
             stop("For a 2-fleet 9-perturbation 9-perturbation model par must be a vector of length 26")
            }
          results <- do.call(.CDMN9P9P, parlist);
          }
        else if(sum(p == c(10,10)) == length(p))
          {
          if(length(dates) != 22)
            {
             stop("For a 2-fleet 10-perturbation 10-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 10, fleet 2 perturbation 1, ..., fleet 2 perturbation 10, final")
            }
          if(length(par) != 28)
            {
             stop("For a 2-fleet 10-perturbation 10-perturbation model par must be a vector of length 28")
            }
          results <- do.call(.CDMN10P10P, parlist);
          }
        else if(sum(p == c(11,11)) == length(p))
          {
          if(length(dates) != 24)
            {
             stop("For a 2-fleet 11-perturbation 11-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 11, fleet 2 perturbation 1, ..., fleet 2 perturbation 11, final")
            }
          if(length(par) != 30)
            {
             stop("For a 2-fleet 11-perturbation 11-perturbation model par must be a vector of length 30")
            }
          results <- do.call(.CDMN11P11P, parlist);
          }
        else if(sum(p == c(12,12)) == length(p))
          {
          if(length(dates) != 26)
            {
             stop("For a 2-fleet 12-perturbation 12-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 12, fleet 2 perturbation 1, ..., fleet 2 perturbation 12, final")
            }
          if(length(par) != 32)
            {
             stop("For a 2-fleet 12-perturbation 12-perturbation model par must be a vector of length 32")
            }
          results <- do.call(.CDMN12P12P, parlist);
          }
        else if(sum(p == c(13,13)) == length(p))
          {
          if(length(dates) != 28)
            {
             stop("For a 2-fleet 13-perturbation 13-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 13, fleet 2 perturbation 1, ..., fleet 2 perturbation 13, final")
            }
          if(length(par) != 34)
            {
             stop("For a 2-fleet 13-perturbation 13-perturbation model par must be a vector of length 34")
            }
          results <- do.call(.CDMN13P13P, parlist);
          }
        else if(sum(p == c(14,14)) == length(p))
          {
          if(length(dates) != 30)
            {
             stop("For a 2-fleet 14-perturbation 14-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 14, fleet 2 perturbation 1, ..., fleet 2 perturbation 14, final")
            }
          if(length(par) != 36)
            {
             stop("For a 2-fleet 14-perturbation 14-perturbation model par must be a vector of length 36")
            }
          results <- do.call(.CDMN14P14P, parlist);
          }
        else if(sum(p == c(15,15)) == length(p))
          {
          if(length(dates) != 32)
            {
             stop("For a 2-fleet 15-perturbation 15-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 15, fleet 2 perturbation 1, ..., fleet 2 perturbation 15, final")
            }
          if(length(par) != 38)
            {
             stop("For a 2-fleet 15-perturbation 15-perturbation model par must be a vector of length 38")
            }
          results <- do.call(.CDMN15P15P, parlist);
          }
        else if(sum(p == c(16,16)) == length(p))
          {
          if(length(dates) != 34)
            {
             stop("For a 2-fleet 16-perturbation 16-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 16, fleet 2 perturbation 1, ..., fleet 2 perturbation 16, final")
            }
          if(length(par) != 40)
            {
             stop("For a 2-fleet 16-perturbation 16-perturbation model par must be a vector of length 40")
            }
          results <- do.call(.CDMN16P16P, parlist);
          }
        else if(sum(p == c(17,17)) == length(p))
          {
          if(length(dates) != 36)
            {
             stop("For a 2-fleet 17-perturbation 17-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 17, fleet 2 perturbation 1, ..., fleet 2 perturbation 17, final")
            }
          if(length(par) != 42)
            {
             stop("For a 2-fleet 17-perturbation 17-perturbation model par must be a vector of length 42")
            }
          results <- do.call(.CDMN17P17P, parlist);
          }
        else if(sum(p == c(18,18)) == length(p))
          {
          if(length(dates) != 38)
            {
             stop("For a 2-fleet 18-perturbation 18-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 18, fleet 2 perturbation 1, ..., fleet 2 perturbation 18, final")
            }
          if(length(par) != 44)
            {
             stop("For a 2-fleet 18-perturbation 18-perturbation model par must be a vector of length 44")
            }
          results <- do.call(.CDMN18P18P, parlist);
          }
        else if(sum(p == c(19,19)) == length(p))
          {
          if(length(dates) != 40)
            {
             stop("For a 2-fleet 19-perturbation 19-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 19, fleet 2 perturbation 1, ..., fleet 2 perturbation 19, final")
            }
          if(length(par) != 46)
            {
             stop("For a 2-fleet 19-perturbation 19-perturbation model par must be a vector of length 46")
            }
          results <- do.call(.CDMN19P19P, parlist);
          }
        else if(sum(p == c(20,20)) == length(p))
          {
          if(length(dates) != 42)
            {
             stop("For a 2-fleet 20-perturbation 20-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 20, fleet 2 perturbation 1, ..., fleet 2 perturbation 20, final")
            }
          if(length(par) != 48)
            {
             stop("For a 2-fleet 20-perturbation 20-perturbation model par must be a vector of length 48")
            }
          results <- do.call(.CDMN20P20P, parlist);
          }
        else if(sum(p == c(21,21)) == length(p))
          {
          if(length(dates) != 44)
            {
             stop("For a 2-fleet 21-perturbation 21-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 21, fleet 2 perturbation 1, ..., fleet 2 perturbation 21, final")
            }
          if(length(par) != 50)
            {
             stop("For a 2-fleet 21-perturbation 21-perturbation model par must be a vector of length 50")
            }
          results <- do.call(.CDMN21P21P, parlist);
          }
        else if(sum(p == c(22,22)) == length(p))
          {
          if(length(dates) != 46)
            {
             stop("For a 2-fleet 22-perturbation 22-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 22, fleet 2 perturbation 1, ..., fleet 2 perturbation 22, final")
            }
          if(length(par) != 52)
            {
             stop("For a 2-fleet 22-perturbation 22-perturbation model par must be a vector of length 52")
            }
          results <- do.call(.CDMN22P22P, parlist);
          }
        else if(sum(p == c(23,23)) == length(p))
          {
          if(length(dates) != 48)
            {
             stop("For a 2-fleet 23-perturbation 23-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 23, fleet 2 perturbation 1, ..., fleet 2 perturbation 23, final")
            }
          if(length(par) != 54)
            {
             stop("For a 2-fleet 23-perturbation 23-perturbation model par must be a vector of length 54")
            }
          results <- do.call(.CDMN23P23P, parlist);
          }
        else if(sum(p == c(24,24)) == length(p))
          {
          if(length(dates) != 50)
            {
             stop("For a 2-fleet 24-perturbation 24-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 24, fleet 2 perturbation 1, ..., fleet 2 perturbation 24, final")
            }
          if(length(par) != 56)
            {
             stop("For a 2-fleet 24-perturbation 24-perturbation model par must be a vector of length 56")
            }
          results <- do.call(.CDMN24P24P, parlist);
          }
        else if(sum(p == c(25,25)) == length(p))
          {
          if(length(dates) != 52)
            {
             stop("For a 2-fleet 25-perturbation 25-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 25, fleet 2 perturbation 1, ..., fleet 2 perturbation 25, final")
            }
          if(length(par) != 58)
            {
             stop("For a 2-fleet 25-perturbation 25-perturbation model par must be a vector of length 58")
            }
          results <- do.call(.CDMN25P25P, parlist);
          }
        else if(sum(p == c(26,26)) == length(p))
          {
          if(length(dates) != 54)
            {
             stop("For a 2-fleet 26-perturbation 26-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 26, fleet 2 perturbation 1, ..., fleet 2 perturbation 26, final")
            }
          if(length(par) != 60)
            {
             stop("For a 2-fleet 26-perturbation 26-perturbation model par must be a vector of length 60")
            }
          results <- do.call(.CDMN26P26P, parlist);
          }
        else if(sum(p == c(27,27)) == length(p))
          {
          if(length(dates) != 56)
            {
             stop("For a 2-fleet 27-perturbation 27-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 27, fleet 2 perturbation 1, ..., fleet 2 perturbation 27, final")
            }
          if(length(par) != 62)
            {
             stop("For a 2-fleet 27-perturbation 27-perturbation model par must be a vector of length 62")
            }
          results <- do.call(.CDMN27P27P, parlist);
          }
        else if(sum(p == c(28,28)) == length(p))
          {
          if(length(dates) != 58)
            {
             stop("For a 2-fleet 28-perturbation 28-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 28, fleet 2 perturbation 1, ..., fleet 2 perturbation 28, final")
            }
          if(length(par) != 64)
            {
             stop("For a 2-fleet 28-perturbation 28-perturbation model par must be a vector of length 64")
            }
          results <- do.call(.CDMN28P28P, parlist);
          }
        else if(sum(p == c(29,29)) == length(p))
          {
          if(length(dates) != 60)
            {
             stop("For a 2-fleet 29-perturbation 29-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 29, fleet 2 perturbation 1, ..., fleet 2 perturbation 29, final")
            }
          if(length(par) != 66)
            {
             stop("For a 2-fleet 29-perturbation 29-perturbation model par must be a vector of length 66")
            }
          results <- do.call(.CDMN29P29P, parlist);
          }
        else if(sum(p == c(30,30)) == length(p))
          {
          if(length(dates) != 62)
            {
             stop("For a 2-fleet 30-perturbation 30-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 30, fleet 2 perturbation 1, ..., fleet 2 perturbation 30, final")
            }
          if(length(par) != 68)
            {
             stop("For a 2-fleet 30-perturbation 30-perturbation model par must be a vector of length 68")
            }
          results <- do.call(.CDMN30P30P, parlist);
          }
        else if(sum(p == c(31,31)) == length(p))
          {
          if(length(dates) != 64)
            {
             stop("For a 2-fleet 31-perturbation 31-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 31, fleet 2 perturbation 1, ..., fleet 2 perturbation 31, final")
            }
          if(length(par) != 70)
            {
             stop("For a 2-fleet 31-perturbation 31-perturbation model par must be a vector of length 70")
            }
          results <- do.call(.CDMN31P31P, parlist);
          }
        else if(sum(p == c(32,32)) == length(p))
          {
          if(length(dates) != 66)
            {
             stop("For a 2-fleet 32-perturbation 32-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 32, fleet 2 perturbation 1, ..., fleet 2 perturbation 32, final")
            }
          if(length(par) != 72)
            {
             stop("For a 2-fleet 32-perturbation 32-perturbation model par must be a vector of length 72")
            }
          results <- do.call(.CDMN32P32P, parlist);
          }
        else if(sum(p == c(33,33)) == length(p))
          {
          if(length(dates) != 68)
            {
             stop("For a 2-fleet 33-perturbation 33-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 33, fleet 2 perturbation 1, ..., fleet 2 perturbation 33, final")
            }
          if(length(par) != 74)
            {
             stop("For a 2-fleet 33-perturbation 33-perturbation model par must be a vector of length 74")
            }
          results <- do.call(.CDMN33P33P, parlist);
          }
        else if(sum(p == c(34,34)) == length(p))
          {
          if(length(dates) != 70)
            {
             stop("For a 2-fleet 34-perturbation 34-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 34, fleet 2 perturbation 1, ..., fleet 2 perturbation 34, final")
            }
          if(length(par) != 76)
            {
             stop("For a 2-fleet 34-perturbation 34-perturbation model par must be a vector of length 76")
            }
          results <- do.call(.CDMN34P34P, parlist);
          }
        else if(sum(p == c(35,35)) == length(p))
          {
          if(length(dates) != 72)
            {
             stop("For a 2-fleet 35-perturbation 35-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 35, fleet 2 perturbation 1, ..., fleet 2 perturbation 35, final")
            }
          if(length(par) != 78)
            {
             stop("For a 2-fleet 35-perturbation 35-perturbation model par must be a vector of length 78")
            }
          results <- do.call(.CDMN35P35P, parlist);
          }
        else if(sum(p == c(35,35)) == length(p))
          {
          if(length(dates) != 72)
            {
             stop("For a 2-fleet 35-perturbation 35-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 35, fleet 2 perturbation 1, ..., fleet 2 perturbation 35, final")
            }
          if(length(par) != 78)
            {
             stop("For a 2-fleet 35-perturbation 35-perturbation model par must be a vector of length 78")
            }
          results <- do.call(.CDMN35P35P, parlist);
          }
        else if(sum(p == c(36,36)) == length(p))
          {
          if(length(dates) != 74)
            {
             stop("For a 2-fleet 36-perturbation 36-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 36, fleet 2 perturbation 1, ..., fleet 2 perturbation 36, final")
            }
          if(length(par) != 80)
            {
             stop("For a 2-fleet 36-perturbation 36-perturbation model par must be a vector of length 80")
            }
          results <- do.call(.CDMN36P36P, parlist);
          }
        else if(sum(p == c(37,37)) == length(p))
          {
          if(length(dates) != 76)
            {
             stop("For a 2-fleet 37-perturbation 37-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 37, fleet 2 perturbation 1, ..., fleet 2 perturbation 37, final")
            }
          if(length(par) != 82)
            {
             stop("For a 2-fleet 37-perturbation 37-perturbation model par must be a vector of length 82")
            }
          results <- do.call(.CDMN37P37P, parlist);
          }
        else if(sum(p == c(38,38)) == length(p))
          {
          if(length(dates) != 78)
            {
             stop("For a 2-fleet 38-perturbation 38-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 38, fleet 2 perturbation 1, ..., fleet 2 perturbation 38, final")
            }
          if(length(par) != 84)
            {
             stop("For a 2-fleet 38-perturbation 38-perturbation model par must be a vector of length 84")
            }
          results <- do.call(.CDMN38P38P, parlist);
          }
        else if(sum(p == c(39,39)) == length(p))
          {
          if(length(dates) != 80)
            {
             stop("For a 2-fleet 39-perturbation 39-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 39, fleet 2 perturbation 1, ..., fleet 2 perturbation 39, final")
            }
          if(length(par) != 86)
            {
             stop("For a 2-fleet 39-perturbation 39-perturbation model par must be a vector of length 86")
            }
          results <- do.call(.CDMN39P39P, parlist);
          }
        else if(sum(p == c(40,40)) == length(p))
          {
          if(length(dates) != 82)
            {
             stop("For a 2-fleet 40-perturbation 40-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 40, fleet 2 perturbation 1, ..., fleet 2 perturbation 40, final")
            }
          if(length(par) != 88)
            {
             stop("For a 2-fleet 40-perturbation 40-perturbation model par must be a vector of length 88")
            }
          results <- do.call(.CDMN40P40P, parlist);
          }
        else if(sum(p == c(41,41)) == length(p))
          {
          if(length(dates) != 84)
            {
             stop("For a 2-fleet 41-perturbation 41-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 41, fleet 2 perturbation 1, ..., fleet 2 perturbation 41, final")
            }
          if(length(par) != 90)
            {
             stop("For a 2-fleet 41-perturbation 41-perturbation model par must be a vector of length 90")
            }
          results <- do.call(.CDMN41P41P, parlist);
          }
        else if(sum(p == c(42,42)) == length(p))
          {
          if(length(dates) != 86)
            {
             stop("For a 2-fleet 42-perturbation 42-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 42, fleet 2 perturbation 1, ..., fleet 2 perturbation 42, final")
            }
          if(length(par) != 92)
            {
             stop("For a 2-fleet 42-perturbation 42-perturbation model par must be a vector of length 92")
            }
          results <- do.call(.CDMN42P42P, parlist);
          }
        else if(sum(p == c(43,43)) == length(p))
          {
          if(length(dates) != 88)
            {
             stop("For a 2-fleet 43-perturbation 43-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 43, fleet 2 perturbation 1, ..., fleet 2 perturbation 43, final")
            }
          if(length(par) != 94)
            {
             stop("For a 2-fleet 43-perturbation 43-perturbation model par must be a vector of length 94")
            }
          results <- do.call(.CDMN43P43P, parlist);
          }
        else if(sum(p == c(44,44)) == length(p))
          {
          if(length(dates) != 90)
            {
             stop("For a 2-fleet 44-perturbation 44-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 44, fleet 2 perturbation 1, ..., fleet 2 perturbation 44, final")
            }
          if(length(par) != 96)
            {
             stop("For a 2-fleet 44-perturbation 44-perturbation model par must be a vector of length 96")
            }
          results <- do.call(.CDMN44P44P, parlist);
          }
        else if(sum(p == c(45,45)) == length(p))
          {
          if(length(dates) != 92)
            {
             stop("For a 2-fleet 45-perturbation 45-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 45, fleet 2 perturbation 1, ..., fleet 2 perturbation 45, final")
            }
          if(length(par) != 98)
            {
             stop("For a 2-fleet 45-perturbation 45-perturbation model par must be a vector of length 98")
            }
          results <- do.call(.CDMN45P45P, parlist);
          }
        else if(sum(p == c(46,46)) == length(p))
          {
          if(length(dates) != 94)
            {
             stop("For a 2-fleet 46-perturbation 46-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 46, fleet 2 perturbation 1, ..., fleet 2 perturbation 46, final")
            }
          if(length(par) != 100)
            {
             stop("For a 2-fleet 46-perturbation 46-perturbation model par must be a vector of length 100")
            }
          results <- do.call(.CDMN46P46P, parlist);
          }
        else if(sum(p == c(47,47)) == length(p))
          {
          if(length(dates) != 96)
            {
             stop("For a 2-fleet 47-perturbation 47-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 47, fleet 2 perturbation 1, ..., fleet 2 perturbation 47, final")
            }
          if(length(par) != 102)
            {
             stop("For a 2-fleet 47-perturbation 47-perturbation model par must be a vector of length 102")
            }
          results <- do.call(.CDMN47P47P, parlist);
          }
        else if(sum(p == c(48,48)) == length(p))
          {
          if(length(dates) != 98)
            {
             stop("For a 2-fleet 48-perturbation 48-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 48, fleet 2 perturbation 1, ..., fleet 2 perturbation 48, final")
            }
          if(length(par) != 104)
            {
             stop("For a 2-fleet 48-perturbation 48-perturbation model par must be a vector of length 104")
            }
          results <- do.call(.CDMN48P48P, parlist);
          }
        else if(sum(p == c(49,49)) == length(p))
          {
          if(length(dates) != 100)
            {
             stop("For a 2-fleet 49-perturbation 49-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 49, fleet 2 perturbation 1, ..., fleet 2 perturbation 49, final")
            }
          if(length(par) != 106)
            {
             stop("For a 2-fleet 49-perturbation 49-perturbation model par must be a vector of length 106")
            }
          results <- do.call(.CDMN49P49P, parlist);
          }
        else if(sum(p == c(50,50)) == length(p))
          {
          if(length(dates) != 102)
            {
             stop("For a 2-fleet 50-perturbation 50-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 50, fleet 2 perturbation 1, ..., fleet 2 perturbation 50, final")
            }
          if(length(par) != 108)
            {
             stop("For a 2-fleet 50-perturbation 50-perturbation model par must be a vector of length 108")
            }
          results <- do.call(.CDMN50P50P, parlist);
          }
        else if(sum(p == c(51,51)) == length(p))
          {
          if(length(dates) != 104)
            {
             stop("For a 2-fleet 51-perturbation 51-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 51, fleet 2 perturbation 1, ..., fleet 2 perturbation 51, final")
            }
          if(length(par) != 110)
            {
             stop("For a 2-fleet 51-perturbation 51-perturbation model par must be a vector of length 110")
            }
          results <- do.call(.CDMN51P51P, parlist);
          }
        else if(sum(p == c(52,52)) == length(p))
          {
          if(length(dates) != 106)
            {
             stop("For a 2-fleet 52-perturbation 52-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 52, fleet 2 perturbation 1, ..., fleet 2 perturbation 52, final")
            }
          if(length(par) != 112)
            {
             stop("For a 2-fleet 52-perturbation 52-perturbation model par must be a vector of length 112")
            }
          results <- do.call(.CDMN52P52P, parlist);
          }
        else if(sum(p == c(53,53)) == length(p))
          {
          if(length(dates) != 108)
            {
             stop("For a 2-fleet 53-perturbation 53-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 53, fleet 2 perturbation 1, ..., fleet 2 perturbation 53, final")
            }
          if(length(par) != 114)
            {
             stop("For a 2-fleet 53-perturbation 53-perturbation model par must be a vector of length 114")
            }
          results <- do.call(.CDMN53P53P, parlist);
          }
        else if(sum(p == c(54,54)) == length(p))
          {
          if(length(dates) != 110)
            {
             stop("For a 2-fleet 54-perturbation 54-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 54, fleet 2 perturbation 1, ..., fleet 2 perturbation 54, final")
            }
          if(length(par) != 116)
            {
             stop("For a 2-fleet 54-perturbation 54-perturbation model par must be a vector of length 116")
            }
          results <- do.call(.CDMN54P54P, parlist);
          }
        else if(sum(p == c(55,55)) == length(p))
          {
          if(length(dates) != 112)
            {
             stop("For a 2-fleet 55-perturbation 55-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 55, fleet 2 perturbation 1, ..., fleet 2 perturbation 55, final")
            }
          if(length(par) != 118)
            {
             stop("For a 2-fleet 55-perturbation 55-perturbation model par must be a vector of length 118")
            }
          results <- do.call(.CDMN55P55P, parlist);
          }
        else if(sum(p == c(56,56)) == length(p))
          {
          if(length(dates) != 114)
            {
             stop("For a 2-fleet 56-perturbation 56-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 56, fleet 2 perturbation 1, ..., fleet 2 perturbation 56, final")
            }
          if(length(par) != 120)
            {
             stop("For a 2-fleet 56-perturbation 56-perturbation model par must be a vector of length 120")
            }
          results <- do.call(.CDMN56P56P, parlist);
          }
        else if(sum(p == c(57,57)) == length(p))
          {
          if(length(dates) != 116)
            {
             stop("For a 2-fleet 57-perturbation 57-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 57, fleet 2 perturbation 1, ..., fleet 2 perturbation 57, final")
            }
          if(length(par) != 122)
            {
             stop("For a 2-fleet 57-perturbation 57-perturbation model par must be a vector of length 122")
            }
          results <- do.call(.CDMN57P57P, parlist);
          }
        else if(sum(p == c(58,58)) == length(p))
          {
          if(length(dates) != 118)
            {
             stop("For a 2-fleet 58-perturbation 58-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 58, fleet 2 perturbation 1, ..., fleet 2 perturbation 58, final")
            }
          if(length(par) != 124)
            {
             stop("For a 2-fleet 58-perturbation 58-perturbation model par must be a vector of length 124")
            }
          results <- do.call(.CDMN58P58P, parlist);
          }
        else if(sum(p == c(59,59)) == length(p))
          {
          if(length(dates) != 120)
            {
             stop("For a 2-fleet 59-perturbation 59-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 59, fleet 2 perturbation 1, ..., fleet 2 perturbation 59, final")
            }
          if(length(par) != 126)
            {
             stop("For a 2-fleet 59-perturbation 59-perturbation model par must be a vector of length 126")
            }
          results <- do.call(.CDMN59P59P, parlist);
          }
        else if(sum(p == c(60,60)) == length(p))
          {
          if(length(dates) != 122)
            {
             stop("For a 2-fleet 60-perturbation 60-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 60, fleet 2 perturbation 1, ..., fleet 2 perturbation 60, final")
            }
          if(length(par) != 128)
            {
             stop("For a 2-fleet 60-perturbation 60-perturbation model par must be a vector of length 128")
            }
          results <- do.call(.CDMN60P60P, parlist);
          }
        else if(sum(p == c(61,61)) == length(p))
          {
          if(length(dates) != 124)
            {
             stop("For a 2-fleet 61-perturbation 61-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 61, fleet 2 perturbation 1, ..., fleet 2 perturbation 61, final")
            }
          if(length(par) != 130)
            {
             stop("For a 2-fleet 61-perturbation 61-perturbation model par must be a vector of length 130")
            }
          results <- do.call(.CDMN61P61P, parlist);
          }
        else if(sum(p == c(62,62)) == length(p))
          {
          if(length(dates) != 126)
            {
             stop("For a 2-fleet 62-perturbation 62-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 62, fleet 2 perturbation 1, ..., fleet 2 perturbation 62, final")
            }
          if(length(par) != 132)
            {
             stop("For a 2-fleet 62-perturbation 62-perturbation model par must be a vector of length 132")
            }
          results <- do.call(.CDMN62P62P, parlist);
          }
        else if(sum(p == c(63,63)) == length(p))
          {
          if(length(dates) != 128)
            {
             stop("For a 2-fleet 63-perturbation 63-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 63, fleet 2 perturbation 1, ..., fleet 2 perturbation 63, final")
            }
          if(length(par) != 134)
            {
             stop("For a 2-fleet 63-perturbation 63-perturbation model par must be a vector of length 134")
            }
          results <- do.call(.CDMN63P63P, parlist);
          }
        else if(sum(p == c(64,64)) == length(p))
          {
          if(length(dates) != 130)
            {
             stop("For a 2-fleet 64-perturbation 64-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 64, fleet 2 perturbation 1, ..., fleet 2 perturbation 64, final")
            }
          if(length(par) != 136)
            {
             stop("For a 2-fleet 64-perturbation 64-perturbation model par must be a vector of length 136")
            }
          results <- do.call(.CDMN64P64P, parlist);
          }
        else if(sum(p == c(65,65)) == length(p))
          {
          if(length(dates) != 132)
            {
             stop("For a 2-fleet 65-perturbation 65-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 65, fleet 2 perturbation 1, ..., fleet 2 perturbation 65, final")
            }
          if(length(par) != 138)
            {
             stop("For a 2-fleet 65-perturbation 65-perturbation model par must be a vector of length 138")
            }
          results <- do.call(.CDMN65P65P, parlist);
          }
        else if(sum(p == c(66,66)) == length(p))
          {
          if(length(dates) != 134)
            {
             stop("For a 2-fleet 66-perturbation 66-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 66, fleet 2 perturbation 1, ..., fleet 2 perturbation 66, final")
            }
          if(length(par) != 140)
            {
             stop("For a 2-fleet 66-perturbation 66-perturbation model par must be a vector of length 140")
            }
          results <- do.call(.CDMN66P66P, parlist);
          }
        else if(sum(p == c(66,66)) == length(p))
          {
          if(length(dates) != 134)
            {
             stop("For a 2-fleet 66-perturbation 66-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 66, fleet 2 perturbation 1, ..., fleet 2 perturbation 66, final")
            }
          if(length(par) != 140)
            {
             stop("For a 2-fleet 66-perturbation 66-perturbation model par must be a vector of length 140")
            }
          results <- do.call(.CDMN66P66P, parlist);
          }
        else if(sum(p == c(67,67)) == length(p))
          {
          if(length(dates) != 136)
            {
             stop("For a 2-fleet 67-perturbation 67-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 67, fleet 2 perturbation 1, ..., fleet 2 perturbation 67, final")
            }
          if(length(par) != 142)
            {
             stop("For a 2-fleet 67-perturbation 67-perturbation model par must be a vector of length 142")
            }
          results <- do.call(.CDMN67P67P, parlist);
          }
        else if(sum(p == c(68,68)) == length(p))
          {
          if(length(dates) != 138)
            {
             stop("For a 2-fleet 68-perturbation 68-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 68, fleet 2 perturbation 1, ..., fleet 2 perturbation 68, final")
            }
          if(length(par) != 144)
            {
             stop("For a 2-fleet 68-perturbation 68-perturbation model par must be a vector of length 144")
            }
          results <- do.call(.CDMN68P68P, parlist);
          }
        else if(sum(p == c(69,69)) == length(p))
          {
          if(length(dates) != 140)
            {
             stop("For a 2-fleet 69-perturbation 69-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 69, fleet 2 perturbation 1, ..., fleet 2 perturbation 69, final")
            }
          if(length(par) != 146)
            {
             stop("For a 2-fleet 69-perturbation 69-perturbation model par must be a vector of length 146")
            }
          results <- do.call(.CDMN69P69P, parlist);
          }
        else if(sum(p == c(70,70)) == length(p))
          {
          if(length(dates) != 142)
            {
             stop("For a 2-fleet 70-perturbation 70-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 70, fleet 2 perturbation 1, ..., fleet 2 perturbation 70, final")
            }
          if(length(par) != 148)
            {
             stop("For a 2-fleet 70-perturbation 70-perturbation model par must be a vector of length 148")
            }
          results <- do.call(.CDMN70P70P, parlist);
          }
        else if(sum(p == c(71,71)) == length(p))
          {
          if(length(dates) != 144)
            {
             stop("For a 2-fleet 71-perturbation 71-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 71, fleet 2 perturbation 1, ..., fleet 2 perturbation 71, final")
            }
          if(length(par) != 150)
            {
             stop("For a 2-fleet 71-perturbation 71-perturbation model par must be a vector of length 150")
            }
          results <- do.call(.CDMN71P71P, parlist);
          }
        else if(sum(p == c(72,72)) == length(p))
          {
          if(length(dates) != 146)
            {
             stop("For a 2-fleet 72-perturbation 72-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 72, fleet 2 perturbation 1, ..., fleet 2 perturbation 72, final")
            }
          if(length(par) != 152)
            {
             stop("For a 2-fleet 72-perturbation 72-perturbation model par must be a vector of length 152")
            }
          results <- do.call(.CDMN72P72P, parlist);
          }
        else if(sum(p == c(73,73)) == length(p))
          {
          if(length(dates) != 148)
            {
             stop("For a 2-fleet 73-perturbation 73-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 73, fleet 2 perturbation 1, ..., fleet 2 perturbation 73, final")
            }
          if(length(par) != 154)
            {
             stop("For a 2-fleet 73-perturbation 73-perturbation model par must be a vector of length 154")
            }
          results <- do.call(.CDMN73P73P, parlist);
          }
        else if(sum(p == c(74,74)) == length(p))
          {
          if(length(dates) != 150)
            {
             stop("For a 2-fleet 74-perturbation 74-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 74, fleet 2 perturbation 1, ..., fleet 2 perturbation 74, final")
            }
          if(length(par) != 156)
            {
             stop("For a 2-fleet 74-perturbation 74-perturbation model par must be a vector of length 156")
            }
          results <- do.call(.CDMN74P74P, parlist);
          }
        else if(sum(p == c(75,75)) == length(p))
          {
          if(length(dates) != 152)
            {
             stop("For a 2-fleet 75-perturbation 75-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 75, fleet 2 perturbation 1, ..., fleet 2 perturbation 75, final")
            }
          if(length(par) != 158)
            {
             stop("For a 2-fleet 75-perturbation 75-perturbation model par must be a vector of length 158")
            }
          results <- do.call(.CDMN75P75P, parlist);
          }
        else if(sum(p == c(76,76)) == length(p))
          {
          if(length(dates) != 154)
            {
             stop("For a 2-fleet 76-perturbation 76-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 76, fleet 2 perturbation 1, ..., fleet 2 perturbation 76, final")
            }
          if(length(par) != 160)
            {
             stop("For a 2-fleet 76-perturbation 76-perturbation model par must be a vector of length 160")
            }
          results <- do.call(.CDMN76P76P, parlist);
          }
        else if(sum(p == c(77,77)) == length(p))
          {
          if(length(dates) != 156)
            {
             stop("For a 2-fleet 77-perturbation 77-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 77, fleet 2 perturbation 1, ..., fleet 2 perturbation 77, final")
            }
          if(length(par) != 162)
            {
             stop("For a 2-fleet 77-perturbation 77-perturbation model par must be a vector of length 162")
            }
          results <- do.call(.CDMN77P77P, parlist);
          }
        else if(sum(p == c(78,78)) == length(p))
          {
          if(length(dates) != 158)
            {
             stop("For a 2-fleet 78-perturbation 78-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 78, fleet 2 perturbation 1, ..., fleet 2 perturbation 78, final")
            }
          if(length(par) != 164)
            {
             stop("For a 2-fleet 78-perturbation 78-perturbation model par must be a vector of length 164")
            }
          results <- do.call(.CDMN78P78P, parlist);
          }
        else if(sum(p == c(79,79)) == length(p))
          {
          if(length(dates) != 160)
            {
             stop("For a 2-fleet 79-perturbation 79-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 79, fleet 2 perturbation 1, ..., fleet 2 perturbation 79, final")
            }
          if(length(par) != 166)
            {
             stop("For a 2-fleet 79-perturbation 79-perturbation model par must be a vector of length 166")
            }
          results <- do.call(.CDMN79P79P, parlist);
          }
        else if(sum(p == c(80,80)) == length(p))
          {
          if(length(dates) != 162)
            {
             stop("For a 2-fleet 80-perturbation 80-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 80, fleet 2 perturbation 1, ..., fleet 2 perturbation 80, final")
            }
          if(length(par) != 168)
            {
             stop("For a 2-fleet 80-perturbation 80-perturbation model par must be a vector of length 168")
            }
          results <- do.call(.CDMN80P80P, parlist);
          }
        else if(sum(p == c(81,81)) == length(p))
          {
          if(length(dates) != 164)
            {
             stop("For a 2-fleet 81-perturbation 81-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 81, fleet 2 perturbation 1, ..., fleet 2 perturbation 81, final")
            }
          if(length(par) != 170)
            {
             stop("For a 2-fleet 81-perturbation 81-perturbation model par must be a vector of length 170")
            }
          results <- do.call(.CDMN81P81P, parlist);
          }
        else if(sum(p == c(82,82)) == length(p))
          {
          if(length(dates) != 166)
            {
             stop("For a 2-fleet 82-perturbation 82-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 82, fleet 2 perturbation 1, ..., fleet 2 perturbation 82, final")
            }
          if(length(par) != 172)
            {
             stop("For a 2-fleet 82-perturbation 82-perturbation model par must be a vector of length 172")
            }
          results <- do.call(.CDMN82P82P, parlist);
          }
        else if(sum(p == c(83,83)) == length(p))
          {
          if(length(dates) != 168)
            {
             stop("For a 2-fleet 83-perturbation 83-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 83, fleet 2 perturbation 1, ..., fleet 2 perturbation 83, final")
            }
          if(length(par) != 174)
            {
             stop("For a 2-fleet 83-perturbation 83-perturbation model par must be a vector of length 174")
            }
          results <- do.call(.CDMN83P83P, parlist);
          }
        else if(sum(p == c(84,84)) == length(p))
          {
          if(length(dates) != 170)
            {
             stop("For a 2-fleet 84-perturbation 84-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 84, fleet 2 perturbation 1, ..., fleet 2 perturbation 84, final")
            }
          if(length(par) != 176)
            {
             stop("For a 2-fleet 84-perturbation 84-perturbation model par must be a vector of length 176")
            }
          results <- do.call(.CDMN84P84P, parlist);
          }
        else if(sum(p == c(85,85)) == length(p))
          {
          if(length(dates) != 172)
            {
             stop("For a 2-fleet 85-perturbation 85-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 85, fleet 2 perturbation 1, ..., fleet 2 perturbation 85, final")
            }
          if(length(par) != 178)
            {
             stop("For a 2-fleet 85-perturbation 85-perturbation model par must be a vector of length 178")
            }
          results <- do.call(.CDMN85P85P, parlist);
          }
        else if(sum(p == c(86,86)) == length(p))
          {
          if(length(dates) != 174)
            {
             stop("For a 2-fleet 86-perturbation 86-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 86, fleet 2 perturbation 1, ..., fleet 2 perturbation 86, final")
            }
          if(length(par) != 180)
            {
             stop("For a 2-fleet 86-perturbation 86-perturbation model par must be a vector of length 180")
            }
          results <- do.call(.CDMN86P86P, parlist);
          }
        else if(sum(p == c(87,87)) == length(p))
          {
          if(length(dates) != 176)
            {
             stop("For a 2-fleet 87-perturbation 87-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 87, fleet 2 perturbation 1, ..., fleet 2 perturbation 87, final")
            }
          if(length(par) != 182)
            {
             stop("For a 2-fleet 87-perturbation 87-perturbation model par must be a vector of length 182")
            }
          results <- do.call(.CDMN87P87P, parlist);
          }
        else if(sum(p == c(88,88)) == length(p))
          {
          if(length(dates) != 178)
            {
             stop("For a 2-fleet 88-perturbation 88-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 88, fleet 2 perturbation 1, ..., fleet 2 perturbation 88, final")
            }
          if(length(par) != 184)
            {
             stop("For a 2-fleet 88-perturbation 88-perturbation model par must be a vector of length 184")
            }
          results <- do.call(.CDMN88P88P, parlist);
          }
        else if(sum(p == c(89,89)) == length(p))
          {
          if(length(dates) != 180)
            {
             stop("For a 2-fleet 89-perturbation 89-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 89, fleet 2 perturbation 1, ..., fleet 2 perturbation 89, final")
            }
          if(length(par) != 186)
            {
             stop("For a 2-fleet 89-perturbation 89-perturbation model par must be a vector of length 186")
            }
          results <- do.call(.CDMN89P89P, parlist);
          }
        else if(sum(p == c(90,90)) == length(p))
          {
          if(length(dates) != 182)
            {
             stop("For a 2-fleet 90-perturbation 90-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 90, fleet 2 perturbation 1, ..., fleet 2 perturbation 90, final")
            }
          if(length(par) != 188)
            {
             stop("For a 2-fleet 90-perturbation 90-perturbation model par must be a vector of length 188")
            }
          results <- do.call(.CDMN90P90P, parlist);
          }
        else if(sum(p == c(91,91)) == length(p))
          {
          if(length(dates) != 184)
            {
             stop("For a 2-fleet 91-perturbation 91-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 91, fleet 2 perturbation 1, ..., fleet 2 perturbation 91, final")
            }
          if(length(par) != 190)
            {
             stop("For a 2-fleet 91-perturbation 91-perturbation model par must be a vector of length 190")
            }
          results <- do.call(.CDMN91P91P, parlist);
          }
        else if(sum(p == c(92,92)) == length(p))
          {
          if(length(dates) != 186)
            {
             stop("For a 2-fleet 92-perturbation 92-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 92, fleet 2 perturbation 1, ..., fleet 2 perturbation 92, final")
            }
          if(length(par) != 192)
            {
             stop("For a 2-fleet 92-perturbation 92-perturbation model par must be a vector of length 192")
            }
          results <- do.call(.CDMN92P92P, parlist);
          }
        else if(sum(p == c(93,93)) == length(p))
          {
          if(length(dates) != 188)
            {
             stop("For a 2-fleet 93-perturbation 93-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 93, fleet 2 perturbation 1, ..., fleet 2 perturbation 93, final")
            }
          if(length(par) != 194)
            {
             stop("For a 2-fleet 93-perturbation 93-perturbation model par must be a vector of length 194")
            }
          results <- do.call(.CDMN93P93P, parlist);
          }
        else if(sum(p == c(94,94)) == length(p))
          {
          if(length(dates) != 190)
            {
             stop("For a 2-fleet 94-perturbation 94-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 94, fleet 2 perturbation 1, ..., fleet 2 perturbation 94, final")
            }
          if(length(par) != 196)
            {
             stop("For a 2-fleet 94-perturbation 94-perturbation model par must be a vector of length 196")
            }
          results <- do.call(.CDMN94P94P, parlist);
          }
        else if(sum(p == c(95,95)) == length(p))
          {
          if(length(dates) != 192)
            {
             stop("For a 2-fleet 95-perturbation 95-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 95, fleet 2 perturbation 1, ..., fleet 2 perturbation 95, final")
            }
          if(length(par) != 198)
            {
             stop("For a 2-fleet 95-perturbation 95-perturbation model par must be a vector of length 198")
            }
          results <- do.call(.CDMN95P95P, parlist);
          }
        else if(sum(p == c(96,96)) == length(p))
          {
          if(length(dates) != 194)
            {
             stop("For a 2-fleet 96-perturbation 96-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 96, fleet 2 perturbation 1, ..., fleet 2 perturbation 96, final")
            }
          if(length(par) != 200)
            {
             stop("For a 2-fleet 96-perturbation 96-perturbation model par must be a vector of length 200")
            }
          results <- do.call(.CDMN96P96P, parlist);
          }
        else if(sum(p == c(97,97)) == length(p))
          {
          if(length(dates) != 196)
            {
             stop("For a 2-fleet 97-perturbation 97-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 97, fleet 2 perturbation 1, ..., fleet 2 perturbation 97, final")
            }
          if(length(par) != 202)
            {
             stop("For a 2-fleet 97-perturbation 97-perturbation model par must be a vector of length 202")
            }
          results <- do.call(.CDMN97P97P, parlist);
          }
        else if(sum(p == c(98,98)) == length(p))
          {
          if(length(dates) != 198)
            {
             stop("For a 2-fleet 98-perturbation 98-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 98, fleet 2 perturbation 1, ..., fleet 2 perturbation 98, final")
            }
          if(length(par) != 204)
            {
             stop("For a 2-fleet 98-perturbation 98-perturbation model par must be a vector of length 204")
            }
          results <- do.call(.CDMN98P98P, parlist);
          }
        else if(sum(p == c(99,99)) == length(p))
          {
          if(length(dates) != 200)
            {
             stop("For a 2-fleet 99-perturbation 99-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 99, fleet 2 perturbation 1, ..., fleet 2 perturbation 99, final")
            }
          if(length(par) != 206)
            {
             stop("For a 2-fleet 99-perturbation 99-perturbation model par must be a vector of length 206")
            }
          results <- do.call(.CDMN99P99P, parlist);
          }
        else if(sum(p == c(100,100)) == length(p))
          {
          if(length(dates) != 2002)
            {
             stop("For a 2-fleet 100-perturbation 100-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 100, fleet 2 perturbation 1, ..., fleet 2 perturbation 100, final")
            }
          if(length(par) != 208)
            {
             stop("For a 2-fleet 99-perturbation 99-perturbation model par must be a vector of length 208")
            }
          results <- do.call(.CDMN100P100P, parlist);
          }
       }
     if(any(results$Model$Results[,4] < 0))
       {stop("Predicted catch at some time steps is negative; consider decreasing natural mortality or increasing abundance. \n In transit fisheries you may also increase vulnerable abundance by delaying exit times")}
     return(results);
    }
