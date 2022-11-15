function vgtsim,a,x,lim,res

;This profile is normalized to one. The profile is sometimes normalized to sqrt(pi) in the literature.

;a is gl/gd
;x is the domain over which the voigt function will be calculated
;lim is the domain over which the convolution integral will be allowed to run
;res is the resolution

t1=systime(1)

ts=dindgen(2*lim/res)*res-lim

x=dindgen(2*x/res+1)*res-x

dop=exp(-ts^2)
vgt=dblarr(3,n_elements(x))
vgt[0,*]=x

lor=1D/(a^2+(x[0]-res-ts)^2)
dt=2D*ts*res
dx2=res^2
dx=2D*(x-res)*res

for i=0,n_elements(x)-1 do begin
  
   lor=1D/(1/temporary(lor)-dt+dx[i]+dx2)
  ;lor=1D/(a^2+(x[i]-ts)^2)
   ;if i mod 5000 eq 0 then p=plot(ts,lor-1D/(a^2+(x[i]-ts)^2))
  
  ;p=plot(ts,lor)
  ;stop
  
  vgt[1,i]=a/sqrt(!pi)^3*int_tabulated(ts,dop*lor)
  
endfor
  
  ;dop=exp(-x^2)
  ;lor=1D/(a^2+x^2)
  ;vgt[2,*]=a/sqrt(!pi)^3*convol(dop,lor,/edge_truncate)*res

print,systime(1)-t1

return,vgt

end

  