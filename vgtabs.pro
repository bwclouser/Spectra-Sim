function vgtabs,nus,gd,gc,conc,fp,nu0,pl,s

dop=(1D/gd)*sqrt(alog(2D)/3.1415926535D)*exp(-alog(2D)*((nus-nu0-fp)/gd)^2)
lor=(1D/(gc*3.1415926536D))*gc^2/((nus-nu0-fp)^2+gc^2)

vgt=convol(dop,lor,/edge_truncate,center=1)

nvgt=vgt/int_tabulated(nus,vgt)

shape=exp(-s*nvgt*conc*pl)

return,shape

end