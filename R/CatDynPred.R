CatDynPred <-
function(x, method)
   {
     fleet.name <- x$Data$Properties$Fleets$Fleet;
     p          <- x$Model[[method]]$Type;
     options(warn=-1)
     if(class(x) != "catdyn")
       {stop("Pass an object of class 'catdyn' to CatDynPred")}
     if(class(method) != "character")
       {stop("method must be a character corresponding to one of the numerical methods passed to the CatDyn function")}  
     if(length(method) != 1)
       {stop("Provide the name of just one the numerical methods used to fit the model")}  
     if(sum(method == names(x$Model)) == 0)
       {stop("The method provided in 'method' was not used to numerically fit the model")}
     if(is.na(x$Model[[method]]$AIC))
       {stop("The selected method failed. Consider trying a different method")}
     if(length(fleet.name) == 1)
       {
        parlist <- list(par=log(as.numeric(x$Model[[method]]$bt.par)), 
                        dates=x$Model[[method]]$Dates, 
                        obscat1=x$Data$Data[[fleet.name]][,5], 
                        obseff1=x$Data$Data[[fleet.name]][,2], 
                        obsmbm1=x$Data$Data[[fleet.name]][,4], 
                        distr=x$Model[[method]]$Distr, 
                        properties=x$Data$Properties,
                        output="predict")
        if(p == 0)
          {
          results <- do.call(.CDMN0P, parlist);
          }
        else if(p == 1)
          {
          results <- do.call(.CDMN1P, parlist);
          }
        else if(p == -1)
          {
          results <- do.call(.CDMNT1P, parlist);
          }
        else if(p == 2)
          {
          results <- do.call(.CDMN2P, parlist);
          }
        else if(p == -2)
          {
          results <- do.call(.CDMNT2P, parlist);
          }
        else if(p == 3)
          {
          results <- do.call(.CDMN3P, parlist);
          }
        else if(p == -3)
          {
          results <- do.call(.CDMNT3P, parlist);
          }
        else if(p == 4)
          {
          results <- do.call(.CDMN4P, parlist);
          }
        else if(p == -4)
          {
          results <- do.call(.CDMNT4P, parlist);
          }
        else if(p == 5)
          {
          results <- do.call(.CDMN5P, parlist);
          }
        else if(p == -5)
          {
          results <- do.call(.CDMNT5P, parlist);
          }
        else if(p == 6)
          {
          results <- do.call(.CDMN6P, parlist);
          }
        else if(p == -6)
          {
          results <- do.call(.CDMNT6P, parlist);
          }
        else if(p == 7)
          {
          results <- do.call(.CDMN7P, parlist);
          }
        else if(p == -7)
          {
          results <- do.call(.CDMNT7P, parlist);
          }
        else if(p == 8)
          {
          results <- do.call(.CDMN8P, parlist);
          }
        else if(p == -8)
          {
          results <- do.call(.CDMNT8P, parlist);
          }
        else if(p == 9)
          {
          results <- do.call(.CDMN9P, parlist);
          }
        else if(p == -9)
          {
          results <- do.call(.CDMNT9P, parlist);
          }
        else if(p == 10)
          {
          results <- do.call(.CDMN10P, parlist);
          }
        else if(p == -10)
          {
          results <- do.call(.CDMNT10P, parlist);
          }
        else if(p == 11)
          {
          results <- do.call(.CDMN11P, parlist);
          }
        else if(p == -11)
          {
          results <- do.call(.CDMNT11P, parlist);
          }
        else if(p == 12)
          {
          results <- do.call(.CDMN12P, parlist);
          }
        else if(p == -12)
          {
          results <- do.call(.CDMNT12P, parlist);
          }
        else if(p == 13)
          {
          results <- do.call(.CDMN13P, parlist);
          }
        else if(p == -13)
          {
          results <- do.call(.CDMNT13P, parlist);
          }
        else if(p == 14)
          {
          results <- do.call(.CDMN14P, parlist);
          }
        else if(p == -14)
          {
          results <- do.call(.CDMNT14P, parlist);
          }
        else if(p == 15)
          {
          results <- do.call(.CDMN15P, parlist);
          }
        else if(p == -15)
          {
          results <- do.call(.CDMNT15P, parlist);
          }
        else if(p == 16)
          {
          results <- do.call(.CDMN16P, parlist);
          }
        else if(p == -16)
          {
          results <- do.call(.CDMNT16P, parlist);
          }
        else if(p == 17)
          {
          results <- do.call(.CDMN17P, parlist);
          }
        else if(p == -17)
          {
          results <- do.call(.CDMNT17P, parlist);
          }
        else if(p == 18)
          {
          results <- do.call(.CDMN18P, parlist);
          }
        else if(p == -18)
          {
          results <- do.call(.CDMNT18P, parlist);
          }
        else if(p == 19)
          {
          results <- do.call(.CDMN19P, parlist);
          }
        else if(p == -19)
          {
          results <- do.call(.CDMNT19P, parlist);
          }
        else if(p == 20)
          {
          results <- do.call(.CDMN20P, parlist);
          }
        else if(p == -20)
          {
          results <- do.call(.CDMNT20P, parlist);
          }
        else if(p == 21)
          {
          results <- do.call(.CDMN21P, parlist);
          }
        else if(p == 22)
          {
          results <- do.call(.CDMN22P, parlist);
          }
        else if(p == 23)
          {
          results <- do.call(.CDMN23P, parlist);
          }
        else if(p == 24)
          {
          results <- do.call(.CDMN24P, parlist);
          }
        else if(p == 25)
          {
          results <- do.call(.CDMN25P, parlist);
          }
        else if(p == 26)
          {
          results <- do.call(.CDMN26P, parlist);
          }
        else if(p == 27)
          {
          results <- do.call(.CDMN27P, parlist);
          }
        else if(p == 28)
          {
          results <- do.call(.CDMN28P, parlist);
          }
        else if(p == 29)
          {
          results <- do.call(.CDMN29P, parlist);
          }
        else if(p == 30)
          {
          results <- do.call(.CDMN30P, parlist);
          }
        else if(p == 31)
          {
          results <- do.call(.CDMN31P, parlist);
          }
        else if(p == 32)
          {
          results <- do.call(.CDMN32P, parlist);
          }
        else if(p == 33)
          {
          results <- do.call(.CDMN33P, parlist);
          }
        else if(p == 34)
          {
          results <- do.call(.CDMN34P, parlist);
          }
        else if(p == 35)
          {
          results <- do.call(.CDMN35P, parlist);
          }
        else if(p == 36)
          {
          results <- do.call(.CDMN36P, parlist);
          }
        else if(p == 37)
          {
          results <- do.call(.CDMN37P, parlist);
          }
        else if(p == 38)
          {
          results <- do.call(.CDMN38P, parlist);
          }
        else if(p == 39)
          {
          results <- do.call(.CDMN39P, parlist);
          }
        else if(p == 40)
          {
          results <- do.call(.CDMN40P, parlist);
          }
        else if(p == 41)
          {
          results <- do.call(.CDMN41P, parlist);
          }
        else if(p == 42)
          {
          results <- do.call(.CDMN42P, parlist);
          }
        else if(p == 43)
          {
          results <- do.call(.CDMN43P, parlist);
          }
        else if(p == 44)
          {
          results <- do.call(.CDMN44P, parlist);
          }
        else if(p == 45)
          {
          results <- do.call(.CDMN45P, parlist);
          }
        else if(p == 46)
          {
          results <- do.call(.CDMN46P, parlist);
          }
        else if(p == 47)
          {
          results <- do.call(.CDMN47P, parlist);
          }
        else if(p == 48)
          {
          results <- do.call(.CDMN48P, parlist);
          }
        else if(p == 49)
          {
          results <- do.call(.CDMN49P, parlist);
          }
        else if(p == 50)
          {
          results <- do.call(.CDMN50P, parlist);
          }
        else if(p == 51)
          {
          results <- do.call(.CDMN51P, parlist);
          }
        else if(p == 52)
          {
          results <- do.call(.CDMN52P, parlist);
          }
        else if(p == 53)
          {
          results <- do.call(.CDMN53P, parlist);
          }
        else if(p == 54)
          {
          results <- do.call(.CDMN54P, parlist);
          }
        else if(p == 55)
          {
          results <- do.call(.CDMN55P, parlist);
          }
        else if(p == 56)
          {
          results <- do.call(.CDMN56P, parlist);
          }
        else if(p == 57)
          {
          results <- do.call(.CDMN57P, parlist);
          }
        else if(p == 58)
          {
          results <- do.call(.CDMN58P, parlist);
          }
        else if(p == 59)
          {
          results <- do.call(.CDMN59P, parlist);
          }
        else if(p == 60)
          {
          results <- do.call(.CDMN60P, parlist);
          }
        else if(p == 61)
          {
          results <- do.call(.CDMN61P, parlist);
          }
        else if(p == 62)
          {
          results <- do.call(.CDMN62P, parlist);
          }
        else if(p == 63)
          {
          results <- do.call(.CDMN63P, parlist);
          }
        else if(p == 64)
          {
          results <- do.call(.CDMN64P, parlist);
          }
        else if(p == 65)
          {
          results <- do.call(.CDMN65P, parlist);
          }
        else if(p == 66)
          {
          results <- do.call(.CDMN66P, parlist);
          }
        else if(p == 67)
          {
          results <- do.call(.CDMN67P, parlist);
          }
        else if(p == 68)
          {
          results <- do.call(.CDMN68P, parlist);
          }
        else if(p == 69)
          {
          results <- do.call(.CDMN69P, parlist);
          }
        else if(p == 70)
          {
          results <- do.call(.CDMN70P, parlist);
          }
        else if(p == 71)
          {
          results <- do.call(.CDMN71P, parlist);
          }
        else if(p == 72)
          {
          results <- do.call(.CDMN72P, parlist);
          }
        else if(p == 73)
          {
          results <- do.call(.CDMN73P, parlist);
          }
        else if(p == 74)
          {
          results <- do.call(.CDMN74P, parlist);
          }
        else if(p == 75)
          {
          results <- do.call(.CDMN75P, parlist);
          }
        else if(p == 76)
          {
          results <- do.call(.CDMN76P, parlist);
          }
        else if(p == 77)
          {
          results <- do.call(.CDMN77P, parlist);
          }
        else if(p == 78)
          {
          results <- do.call(.CDMN78P, parlist);
          }
        else if(p == 79)
          {
          results <- do.call(.CDMN79P, parlist);
          }
        else if(p == 80)
          {
          results <- do.call(.CDMN80P, parlist);
          }
        else if(p == 81)
          {
          results <- do.call(.CDMN81P, parlist);
          }
        else if(p == 82)
          {
          results <- do.call(.CDMN82P, parlist);
          }
        else if(p == 83)
          {
          results <- do.call(.CDMN83P, parlist);
          }
        else if(p == 84)
          {
          results <- do.call(.CDMN84P, parlist);
          }
        else if(p == 85)
          {
          results <- do.call(.CDMN85P, parlist);
          }
        else if(p == 86)
          {
          results <- do.call(.CDMN86P, parlist);
          }
        else if(p == 87)
          {
          results <- do.call(.CDMN87P, parlist);
          }
        else if(p == 88)
          {
          results <- do.call(.CDMN88P, parlist);
          }
        else if(p == 89)
          {
          results <- do.call(.CDMN89P, parlist);
          }
        else if(p == 90)
          {
          results <- do.call(.CDMN90P, parlist);
          }
        else if(p == 91)
          {
          results <- do.call(.CDMN91P, parlist);
          }
        else if(p == 92)
          {
          results <- do.call(.CDMN92P, parlist);
          }
        else if(p == 93)
          {
          results <- do.call(.CDMN93P, parlist);
          }
        else if(p == 94)
          {
          results <- do.call(.CDMN94P, parlist);
          }
        else if(p == 95)
          {
          results <- do.call(.CDMN95P, parlist);
          }
        else if(p == 96)
          {
          results <- do.call(.CDMN96P, parlist);
          }
        else if(p == 97)
          {
          results <- do.call(.CDMN97P, parlist);
          }
        else if(p == 98)
          {
          results <- do.call(.CDMN98P, parlist);
          }
        else if(p == 99)
          {
          results <- do.call(.CDMN99P, parlist);
          }
        else if(p == 100)
          {
          results <- do.call(.CDMN100P, parlist);
          }
       }
     else if(length(fleet.name) == 2) 
       {
        parlist <- list(par=log(as.numeric(x$Model[[method]]$bt.par)), 
                        dates=x$Model[[method]]$Dates, 
                        obscat1=x$Data$Data[[fleet.name[1]]][,5], 
                        obseff1=x$Data$Data[[fleet.name[1]]][,2], 
                        obsmbm1=x$Data$Data[[fleet.name[1]]][,4],
                        obscat2=x$Data$Data[[fleet.name[2]]][,5], 
                        obseff2=x$Data$Data[[fleet.name[2]]][,2], 
                        obsmbm2=x$Data$Data[[fleet.name[2]]][,4], 
                        distr=x$Model[[method]]$Distr, 
                        properties=x$Data$Properties,
                        output="predict")
        if(sum(p==c(0,0)) == length(p))
          {
          results <- do.call(.CDMN0P0P, parlist);
          }
        else if(sum(p==c(0,1)) == length(p))
          {
          results <- do.call(.CDMN0P1P, parlist);
          }
        else if(sum(p==c(0,2)) == length(p))
          {
          results <- do.call(.CDMN0P2P, parlist);
          }
        else if(sum(p==c(0,3)) == length(p))
          {
          results <- do.call(.CDMN0P3P, parlist);
          }
        else if(sum(p==c(0,4)) == length(p))
          {
          results <- do.call(.CDMN0P4P, parlist);
          }
        else if(sum(p==c(0,5)) == length(p))
          {
          results <- do.call(.CDMN0P5P, parlist);
          }
        else if(sum(p==c(1,1)) == length(p))
          {
          results <- do.call(.CDMN1P1P, parlist);
          }                                                                                                                                            
        else if(sum(p==c(1,2)) == length(p))
          {
          results <- do.call(.CDMN1P2P, parlist);
          }
        else if(sum(p==c(1,3)) == length(p))
          {
          results <- do.call(.CDMN1P3P, parlist);
          }
        else if(sum(p==c(1,4)) == length(p))
          {
          results <- do.call(.CDMN1P4P, parlist);
          }
        else if(sum(p==c(1,5)) == length(p))
          {
          results <- do.call(.CDMN1P5P, parlist);
          }
        else if(sum(p==c(2,2)) == length(p))
          {
          results <- do.call(.CDMN2P2P, parlist);
          }
        else if(sum(p==c(2,3)) == length(p))
          {
          results <- do.call(.CDMN2P3P, parlist);
          }
        else if(sum(p==c(2,4)) == length(p))
          {
          results <- do.call(.CDMN2P4P, parlist);
          }
        else if(sum(p==c(2,5)) == length(p))
          {
          results <- do.call(.CDMN2P5P, parlist);
          }
        else if(sum(p==c(3,3)) == length(p))
          {
          results <- do.call(.CDMN3P3P, parlist);
          }
        else if(sum(p==c(3,4)) == length(p))
          {
          results <- do.call(.CDMN3P4P, parlist);
          }
        else if(sum(p==c(3,5)) == length(p))
          {
          results <- do.call(.CDMN3P5P, parlist);
          }
        else if(sum(p==c(4,4)) == length(p))
          {
          results <- do.call(.CDMN4P4P, parlist);
          }
        else if(sum(p==c(4,5)) == length(p))
          {
          results <- do.call(.CDMN4P5P, parlist);
          }
        else if(sum(p==c(5,5)) == length(p))
          {
          results <- do.call(.CDMN5P5P, parlist);
          }
        else if(sum(p==c(6,6)) == length(p))
          {
          results <- do.call(.CDMN6P6P, parlist);
          }
        else if(sum(p==c(7,7)) == length(p))
          {
          results <- do.call(.CDMN7P7P, parlist);
          }
        else if(sum(p==c(8,8)) == length(p))
          {
          results <- do.call(.CDMN8P8P, parlist);
          }
        else if(sum(p==c(9,9)) == length(p))
          {
          results <- do.call(.CDMN9P9P, parlist);
          }
        else if(sum(p==c(10,10)) == length(p))
          {
          results <- do.call(.CDMN10P10P, parlist);
          }
        else if(sum(p==c(11,11)) == length(p))
          {
          results <- do.call(.CDMN11P11P, parlist);
          }
        else if(sum(p==c(12,12)) == length(p))
          {
          results <- do.call(.CDMN12P12P, parlist);
          }
        else if(sum(p==c(13,13)) == length(p))
          {
          results <- do.call(.CDMN13P13P, parlist);
          }
        else if(sum(p==c(14,14)) == length(p))
          {
          results <- do.call(.CDMN14P14P, parlist);
          }
        else if(sum(p==c(15,15)) == length(p))
          {
          results <- do.call(.CDMN15P15P, parlist);
          }
        else if(sum(p==c(16,16)) == length(p))
          {
          results <- do.call(.CDMN16P16P, parlist);
          }
        else if(sum(p==c(17,17)) == length(p))
          {
          results <- do.call(.CDMN17P17P, parlist);
          }
        else if(sum(p==c(18,18)) == length(p))
          {
          results <- do.call(.CDMN18P18P, parlist);
          }
        else if(sum(p==c(19,19)) == length(p))
          {
          results <- do.call(.CDMN19P19P, parlist);
          }
        else if(sum(p==c(20,20)) == length(p))
          {
          results <- do.call(.CDMN20P20P, parlist);
          }
        else if(sum(p==c(21,21)) == length(p))
          {
          results <- do.call(.CDMN21P21P, parlist);
          }
        else if(sum(p==c(22,22)) == length(p))
          {
          results <- do.call(.CDMN22P22P, parlist);
          }
        else if(sum(p==c(23,23)) == length(p))
          {
          results <- do.call(.CDMN23P23P, parlist);
          }
        else if(sum(p==c(24,24)) == length(p))
          {
          results <- do.call(.CDMN24P24P, parlist);
          }
        else if(sum(p==c(25,25)) == length(p))
          {
          results <- do.call(.CDMN25P25P, parlist);
          }
        else if(sum(p==c(26,26)) == length(p))
          {
          results <- do.call(.CDMN26P26P, parlist);
          }
        else if(sum(p==c(27,27)) == length(p))
          {
          results <- do.call(.CDMN27P27P, parlist);
          }
        else if(sum(p==c(28,28)) == length(p))
          {
          results <- do.call(.CDMN28P28P, parlist);
          }
        else if(sum(p==c(29,29)) == length(p))
          {
          results <- do.call(.CDMN29P29P, parlist);
          }
        else if(sum(p==c(30,30)) == length(p))
          {
          results <- do.call(.CDMN30P30P, parlist);
          }
        else if(sum(p==c(31,31)) == length(p))
          {
          results <- do.call(.CDMN31P31P, parlist);
          }
        else if(sum(p==c(32,32)) == length(p))
          {
          results <- do.call(.CDMN32P32P, parlist);
          }
        else if(sum(p==c(33,33)) == length(p))
          {
          results <- do.call(.CDMN33P33P, parlist);
          }
        else if(sum(p==c(34,34)) == length(p))
          {
          results <- do.call(.CDMN34P34P, parlist);
          }
        else if(sum(p==c(35,35)) == length(p))
          {
          results <- do.call(.CDMN35P35P, parlist);
          }
        else if(sum(p==c(36,36)) == length(p))
          {
          results <- do.call(.CDMN36P36P, parlist);
          }
        else if(sum(p==c(37,37)) == length(p))
          {
          results <- do.call(.CDMN37P37P, parlist);
          }
        else if(sum(p==c(38,38)) == length(p))
          {
          results <- do.call(.CDMN38P38P, parlist);
          }
        else if(sum(p==c(39,39)) == length(p))
          {
          results <- do.call(.CDMN39P39P, parlist);
          }
        else if(sum(p==c(40,40)) == length(p))
          {
          results <- do.call(.CDMN40P40P, parlist);
          }
        else if(sum(p==c(41,41)) == length(p))
          {
          results <- do.call(.CDMN41P41P, parlist);
          }
        else if(sum(p==c(42,42)) == length(p))
          {
          results <- do.call(.CDMN42P42P, parlist);
          }
        else if(sum(p==c(43,43)) == length(p))
          {
          results <- do.call(.CDMN43P43P, parlist);
          }
        else if(sum(p==c(44,44)) == length(p))
          {
          results <- do.call(.CDMN44P44P, parlist);
          }
        else if(sum(p==c(45,45)) == length(p))
          {
          results <- do.call(.CDMN45P45P, parlist);
          }
        else if(sum(p==c(46,46)) == length(p))
          {
          results <- do.call(.CDMN46P46P, parlist);
          }
        else if(sum(p==c(47,47)) == length(p))
          {
          results <- do.call(.CDMN47P47P, parlist);
          }
        else if(sum(p==c(48,48)) == length(p))
          {
          results <- do.call(.CDMN48P48P, parlist);
          }
        else if(sum(p==c(49,49)) == length(p))
          {
          results <- do.call(.CDMN49P49P, parlist);
          }
        else if(sum(p==c(50,50)) == length(p))
          {
          results <- do.call(.CDMN50P50P, parlist);
          }
        else if(sum(p==c(51,51)) == length(p))
          {
          results <- do.call(.CDMN51P51P, parlist);
          }
        else if(sum(p==c(52,52)) == length(p))
          {
          results <- do.call(.CDMN52P52P, parlist);
          }
        else if(sum(p==c(53,53)) == length(p))
          {
          results <- do.call(.CDMN53P53P, parlist);
          }
        else if(sum(p==c(54,54)) == length(p))
          {
          results <- do.call(.CDMN54P54P, parlist);
          }
        else if(sum(p==c(55,55)) == length(p))
          {
          results <- do.call(.CDMN55P55P, parlist);
          }
        else if(sum(p==c(56,56)) == length(p))
          {
          results <- do.call(.CDMN56P56P, parlist);
          }
        else if(sum(p==c(57,57)) == length(p))
          {
          results <- do.call(.CDMN57P57P, parlist);
          }
        else if(sum(p==c(58,58)) == length(p))
          {
          results <- do.call(.CDMN58P58P, parlist);
          }
        else if(sum(p==c(59,59)) == length(p))
          {
          results <- do.call(.CDMN59P59P, parlist);
          }
        else if(sum(p==c(60,60)) == length(p))
          {
          results <- do.call(.CDMN60P60P, parlist);
          }
        else if(sum(p==c(61,61)) == length(p))
          {
          results <- do.call(.CDMN61P61P, parlist);
          }
        else if(sum(p==c(62,62)) == length(p))
          {
          results <- do.call(.CDMN62P62P, parlist);
          }
        else if(sum(p==c(63,63)) == length(p))
          {
          results <- do.call(.CDMN63P63P, parlist);
          }
        else if(sum(p==c(64,64)) == length(p))
          {
          results <- do.call(.CDMN64P64P, parlist);
          }
        else if(sum(p==c(65,65)) == length(p))
          {
          results <- do.call(.CDMN65P65P, parlist);
          }
        else if(sum(p==c(66,66)) == length(p))
          {
          results <- do.call(.CDMN66P66P, parlist);
          }
        else if(sum(p==c(67,67)) == length(p))
          {
          results <- do.call(.CDMN67P67P, parlist);
          }
        else if(sum(p==c(68,68)) == length(p))
          {
          results <- do.call(.CDMN68P68P, parlist);
          }
        else if(sum(p==c(69,69)) == length(p))
          {
          results <- do.call(.CDMN69P69P, parlist);
          }
        else if(sum(p==c(70,70)) == length(p))
          {
          results <- do.call(.CDMN70P70P, parlist);
          }
        else if(sum(p==c(71,71)) == length(p))
          {
          results <- do.call(.CDMN71P71P, parlist);
          }
        else if(sum(p==c(72,72)) == length(p))
          {
          results <- do.call(.CDMN72P72P, parlist);
          }
        else if(sum(p==c(73,73)) == length(p))
          {
          results <- do.call(.CDMN73P73P, parlist);
          }
        else if(sum(p==c(74,74)) == length(p))
          {
          results <- do.call(.CDMN74P74P, parlist);
          }
        else if(sum(p==c(75,75)) == length(p))
          {
          results <- do.call(.CDMN75P75P, parlist);
          }
        else if(sum(p==c(76,76)) == length(p))
          {
          results <- do.call(.CDMN76P76P, parlist);
          }
        else if(sum(p==c(77,77)) == length(p))
          {
          results <- do.call(.CDMN77P77P, parlist);
          }
        else if(sum(p==c(78,78)) == length(p))
          {
          results <- do.call(.CDMN78P78P, parlist);
          }
        else if(sum(p==c(79,79)) == length(p))
          {
          results <- do.call(.CDMN79P79P, parlist);
          }
        else if(sum(p==c(80,80)) == length(p))
          {
          results <- do.call(.CDMN80P80P, parlist);
          }
        else if(sum(p==c(81,81)) == length(p))
          {
          results <- do.call(.CDMN81P81P, parlist);
          }
        else if(sum(p==c(82,82)) == length(p))
          {
          results <- do.call(.CDMN82P82P, parlist);
          }
        else if(sum(p==c(83,83)) == length(p))
          {
          results <- do.call(.CDMN83P83P, parlist);
          }
        else if(sum(p==c(84,84)) == length(p))
          {
          results <- do.call(.CDMN84P84P, parlist);
          }
        else if(sum(p==c(85,85)) == length(p))
          {
          results <- do.call(.CDMN85P85P, parlist);
          }
        else if(sum(p==c(86,86)) == length(p))
          {
          results <- do.call(.CDMN86P86P, parlist);
          }
        else if(sum(p==c(87,87)) == length(p))
          {
          results <- do.call(.CDMN87P87P, parlist);
          }
        else if(sum(p==c(88,88)) == length(p))
          {
          results <- do.call(.CDMN88P88P, parlist);
          }
        else if(sum(p==c(89,89)) == length(p))
          {
          results <- do.call(.CDMN89P89P, parlist);
          }
        else if(sum(p==c(90,90)) == length(p))
          {
          results <- do.call(.CDMN90P90P, parlist);
          }
        else if(sum(p==c(91,91)) == length(p))
          {
          results <- do.call(.CDMN91P91P, parlist);
          }
        else if(sum(p==c(92,92)) == length(p))
          {
          results <- do.call(.CDMN92P92P, parlist);
          }
        else if(sum(p==c(93,93)) == length(p))
          {
          results <- do.call(.CDMN93P93P, parlist);
          }
        else if(sum(p==c(94,94)) == length(p))
          {
          results <- do.call(.CDMN94P94P, parlist);
          }
        else if(sum(p==c(95,95)) == length(p))
          {
          results <- do.call(.CDMN95P95P, parlist);
          }
        else if(sum(p==c(96,96)) == length(p))
          {
          results <- do.call(.CDMN96P96P, parlist);
          }
        else if(sum(p==c(97,97)) == length(p))
          {
          results <- do.call(.CDMN97P97P, parlist);
          }
        else if(sum(p==c(98,98)) == length(p))
          {
          results <- do.call(.CDMN98P98P, parlist);
          }
        else if(sum(p==c(99,99)) == length(p))
          {
          results <- do.call(.CDMN99P99P, parlist);
          }
        else if(sum(p==c(100,100)) == length(p))
          {
          results <- do.call(.CDMN100P100P, parlist);
          }
       }
     results$Model$Method     <- method;
     results$Model$AIC        <- x$Model[[method]]$AIC;
     class(results)           <- "CatDynMod";
     return(results);
    }
