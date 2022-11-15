function mdop,nus,gd,nu0

dop=(1D/gd)*sqrt(alog(2D)/3.1415926535D)*exp(-alog(2D)*((nus-nu0)/gd)^2)

return,dop

end