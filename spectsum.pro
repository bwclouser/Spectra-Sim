function  spectsum,alpha,beta,pl,dl,steps,conc

t0=systime(1)

i0=10000.
nui=3789.5
nuf=3790
smin=1e-23
T=235.
P=200.

line=hitranreader(nui,nuf,smin,1,1)



sp=spectrum(nui,nuf,smin,t,p,pl,.001,h2o_161=conc,lines=line)

v=fltarr(n_elements(sp[1,*]))

sp[1,*]=sp[1,*]*exp(-alpha*pl)*i0

for i=0,steps-1 do begin
  
  z=spectrum(nui,nuf,smin,t,p,2D*i*dl,.001,h2o_161=conc,lines=line)
  ;plz=plot(z)
  int=beta*z[1,*]*exp(-2D*alpha*i*dl)*i0*dl
  v[*]+=int[*]
  ;plz=plot(z[0,*],v)
  ;gtb=get_kbrd()
  ;plz.close
  
endfor

sum=fltarr(3,n_elements(z[1,*]))

sum[0,*]=sp[0,*]
sum[1,*]=sp[1,*]
sum[2,*]=v[*]

print,systime(1)-t0

return, sum

end