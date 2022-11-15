pro icosscanapprox,outs,conc,wns,p,t,wres,range=range,plt=plt,hdata=hdata,iend=iend,norm=norm,out=out,passes=passes,methane=methane,outrange=outrange,refl=refl,length=length,n2o=n2o

if keyword_set(refl) then r=refl else r=0.99988d0
if keyword_set(length) then l=length else l=90d0
c=299792458d0

if not keyword_set(range) then ran=[3774d0,3780d0] else ran=range

if keyword_set(passes) then begin
  wres=passes*wns*l*1d-2/c
endif else begin
  passes=double(ulong((c*wres)/(wns*l*1d-2)))
endelse
print,passes

if not keyword_set(hdata) then begin
  ch2o=[conc/1e6*catm(p,t),1e-29]
  ch18=[conc/1d6*catm(p,t),1e-29]
  ch17=ch18
  chdo=ch18
  chd18o=ch18
  chd17o=ch18
  cd2o=ch18
  cco2=[370d0/1d6*catm(p,t),1e-28]
  co3=[.5d0/1d6*catm(p,t),1e-28]
  cn2o=[.31d0/1d6*catm(p,t),1e-26]
  ;cn2o=[600000d0/1d6*catm(p,t),1d-32]
  cco1=[.04d0/1d6*catm(p,t),1e-26]
  ;cco1=[10d0/1d6*catm(p,t),1e-26]
  cch4=[1.68d0/1d6*catm(p,t),1e-24]
  cno=[1d0/1d9*catm(p,t),1e-26]
  chf=[1d0/1d12*catm(p,t),1e-26]
  ;chf=[10000d0/1d12*catm(p,t),1e-26]
  chocl=[1d0/1d9*catm(p,t),1e-22]
  ch2s=[1d0/1d9*catm(p,t),1e-26]
 
  ;stop
  o21=spectrum(ran[0],ran[1],1e-28,t,p,600000d0,wres,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,/pexp,/nocon,/sbc)
endif else begin
    ch2o=[conc/1e6*catm(p,t),1e-29]
    if keyword_set(methane) then begin
      cch4=[methane/1e6*catm(p,t),1e-28]
      if keyword_set(n2o) then cn2o=[n2o/1d6*catm(p,t),1e-28] else cn2o=[0d0,1d-28]
      o21=spectrum(ran[0],ran[1],1e-28,t,p,600000d0,wres,lines=hdata,h2o_161=ch2o,h2o_181=ch2o,h2o_171=ch2o,h2o_162=ch2o,h2o_182=ch2o,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,ch4_312=cch4,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,/pexp,/nocon,/sbc)
      ;stop
    endif else begin
      o21=spectrum(ran[0],ran[1],1e-28,t,p,600000d0,wres,lines=hdata,h2o_161=ch2o,/pexp,/nocon,/sbc)
    endelse
    ;stop
endelse

out=dblarr(n_elements(o21[0,*]))

print,passes

wnsact=wres/passes

if not keyword_set(iend) then iend=80000d0/passes
norm=0d0
;stop
for i=0,iend do begin
  
  out+=shift(r^(2d0*passes*i)*exp(-l*2d0*passes*i*o21[2,*]),-i)
  norm+=r^(2d0*passes*i)
  ;stop
  
endfor

if keyword_set(outrange) then begin
  j=where(o21[0,*] ge outrange[0] and o21[0,*] le outrange[1])
  nj=n_elements(j)
  wnums=dindgen(nj)/double(nj-1)*(outrange[1]-outrange[0])+outrange[0]
endif else begin
  j=where(o21[0,*] ge 3777.3 and o21[0,*] le 3778.4)
  nj=n_elements(j)
  wnums=dindgen(nj)/double(nj-1)*(3778.4-3777.3)+3777.3
endelse


;absp=interpol(out[j],o21[0,j],wnums)/norm

;stop

;plt=plot(wnums,absp)

outs=dblarr(3,nj)
outs[0,*]=wnums
;outs[1,*]=absp
outs[1,*]=out[j]/norm
outs[2,*]=interpol(o21[1,j],o21[0,j],wnums)

end
