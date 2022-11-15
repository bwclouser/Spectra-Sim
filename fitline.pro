function fitline,X,P,line

;P[0] is fine position of line
;P[1] is doppler broadening
;P[2] is lorentz broadening
;P[3] is concentration
;P[4] is HITRAN line center
;P[5] is Path Length
;P[6] is Temp. dependent line strength


dop=exp(-alog(2D)*((X-p[4]-P[0])/P[1])^2)
lor=1D/(((X-P[4]-P[0])/P[2])^2+1D)

vgt=convol(dop,lor,/edge_truncate,center=1)

nvgt=vgt/int_tabulated(X,vgt)

shape=exp(-P[6]*nvgt*P[3]*P[5])

return,shape

end



