function spectrum,nui,nuf,smin,T,P,pl,res,h2o_161=h2o_161,h2o_181=h2o_181,h2o_171=h2o_171,h2o_162=h2o_162,h2o_182=h2o_182,h2o_172=h2o_172,h2o_262=h2o_262,co2_626=co2_626,co2_636=co2_636,co2_628=co2_628,co2_627=co2_627,o3_666=o3_666,o3_668=o3_668,o3_686=o3_686,o3_667=o3_667,o3_676=o3_676,n2o_446=n2o_446,n2o_456=n2o_456,n2o_546=n2o_546,n2o_448=n2o_448,n2o_447=n2o_447,co_26=co_26,co_36=co_36,co_28=co_28,co_27=co_27,ch4_211=ch4_211,ch4_311=ch4_311,ch4_212=ch4_212,ch4_312=ch4_312,no_46=no_46,no_56=no_56,no_48=no_48,hf_19=hf_19,hf_29=hf_29,hocl_165=hocl_165,hocl_167=hocl_167,h2s_121=h2s_121,h2s_141=h2s_141,h2s_131=h2s_131,nh3_411=nh3_411,sbc=sbc,lines=lines,nexp=nexp,pexp=pexp,nocon=nocon,voigt=voigt,writeln=writeln,setw=setw,nus=nus

time=systime(1)

c=299792458D                  ;speed of light in m/s
R=8.3144621D                  ;Universal Gas Constant
T0=296D                       ;reference temperature of HITRAN
P0=1013.25D                   ;1 Atm in hPa

;This program generates a customized spectrum for an arbitrary range of wavenumbers based on Hitran Data
;nui         -  Lowest wavenumber to simulate in the spectrum
;nuf         -  Highest wavenumber to simulate in the spectrum
;smin        -  When searching the Hitran database, only accept lines with linestrength greater than this value
;T           -  Temperature in Kelvin at which to simulate lines
;P           -  Pressure in mb at which to simulate lines
;pl          -  Path length in cm at which to simulate lines
;res         -  Resolution in cm^-1 at which to generate spectrum
;h2o_161     -  Concentration of H2 16O in molecules/cm^3
;h2o_181     -  Concentration of H2 18O in molecules/cm^3
;h2o_171     -  Concentration of H2 17O in molecules/cm^3
;h2o_162     -  Concentration of HD 18O in molecules/cm^3
;sbc         -  Use self broadening in calculation of lorentz width
;lines       -  Array of line structures (see hitranreader.pro for format). These are the lines that will be simulated if included.
;nexp        -  In third column of output array, include vgt lineshape times pl times strength
;pexp        -  In third column of output array, include vgt lineshape times conc times strength
;vgt         -  In third column of output array, inclued normalized voigt lineshape

;This program returns a 3 column array
;Column zero contains the wavenumbers in the spectrum
;Column one contains the transmittance of the spectrum
;Column two contains nothing if no flag is used, or one of the options above if a flag is called

;Setup
;sample call:
;h2o_161=hitranreader(3780,3800,1d-28,1,1)
;h2o_181=hitranreader(3780,3800,1d-28,1,2)
;h2o_171=hitranreader(3780,3800,1d-28,1,3)
;h2o_162=hitranreader(3780,3800,1d-28,1,4)
;lines=[h2o_161,h2o_181,h2o_171,h2o_162]
;oo=spectrum(3760,3800,1e-21,296,10,8600d0,.001,h2o_161=[60d0/1d6*catm(40d0,296d0),1d-23],h2o_181=[60d0/1d6*catm(40d0,296d0),1d-23],h2o_171=[60d0/1d6*catm(40d0,296d0),1d-23],h2o_162=[60d0/1d6*catm(40d0,296d0),1d-23],lines=lines,/pexp,/nocon)
;p=plot(oo[0,*],1d0-oo[1,*])



alines=[]

;If lines are supplied by the user, use those
;If lines are not supplied by the user, then use nui, nuf, and smin to find the appropriate lines
;The program only searches for lines that have a user supplied concentration

if keyword_set(lines) then begin
   alines=lines
endif else begin
 print,'***'
  if keyword_set(h2o_161) then begin
    lines161=hitranreader(nui,nuf,h2o_161[1],1,1)
    alines=[alines,lines161]
  endif
  if keyword_set(h2o_181) then begin
    lines181=hitranreader(nui,nuf,h2o_181[1],1,2)
    alines=[alines,lines181]
  endif
  if keyword_set(h2o_171) then begin
    lines171=hitranreader(nui,nuf,h2o_171[1],1,3)
    alines=[alines,lines171]
  endif
  if keyword_set(h2o_162) then begin
    lines162=hitranreader(nui,nuf,h2o_162[1],1,4) 
    alines=[alines,lines162]
  endif
  if keyword_set(h2o_182) then begin
    lines182=hitranreader(nui,nuf,h2o_182[1],1,5)
    alines=[alines,lines182]
  endif
  if keyword_set(h2o_172) then begin
    lines172=hitranreader(nui,nuf,h2o_172[1],1,6)
    alines=[alines,lines172]
  endif
  if keyword_set(h2o_262) then begin
    lines262=hitranreader(nui,nuf,h2o_262[1],1,7)
    alines=[alines,lines262]
  endif
  if keyword_set(co2_626) then begin
    lines626=hitranreader(nui,nuf,co2_626[1],2,1)
    alines=[alines,lines626]
  endif
  if keyword_set(co2_636) then begin
    lines636=hitranreader(nui,nuf,co2_636[1],2,2)
    alines=[alines,lines636]
  endif
  if keyword_set(co2_628) then begin
    lines628=hitranreader(nui,nuf,co2_628[1],2,3)
    alines=[alines,lines628]
  endif
  if keyword_set(co2_627) then begin
    lines627=hitranreader(nui,nuf,co2_627[1],2,4)
    alines=[alines,lines627]
  endif
  if keyword_set(o3_666) then begin
    lines666=hitranreader(nui,nuf,o3_666[1],3,1)
    alines=[alines,lines666]
  endif
  if keyword_set(o3_668) then begin
    lines668=hitranreader(nui,nuf,o3_668[1],3,2)
    alines=[alines,lines668]
  endif
  if keyword_set(o3_686) then begin
    lines686=hitranreader(nui,nuf,o3_686[1],3,3)
    alines=[alines,lines686]
  endif
  if keyword_set(o3_667) then begin
    lines667=hitranreader(nui,nuf,o3_667[1],3,4)
    alines=[alines,lines667]
  endif
  if keyword_set(o3_676) then begin
    lines676=hitranreader(nui,nuf,o3_676[1],3,5)
    alines=[alines,lines676]
  endif
  if keyword_set(n2o_446) then begin
    lines446=hitranreader(nui,nuf,n2o_446[1],4,1)
    alines=[alines,lines446]
  endif
  if keyword_set(n2o_456) then begin
    lines456=hitranreader(nui,nuf,n2o_456[1],4,2)
    alines=[alines,lines456]
  endif  
  if keyword_set(n2o_546) then begin
    lines546=hitranreader(nui,nuf,n2o_546[1],4,3)
    alines=[alines,lines546]
  endif
  if keyword_set(n2o_448) then begin
    lines448=hitranreader(nui,nuf,n2o_448[1],4,4)
    alines=[alines,lines448]
  endif
  if keyword_set(n2o_447) then begin
    lines447=hitranreader(nui,nuf,n2o_447[1],4,5)
    alines=[alines,lines447]
  endif
  if keyword_set(co_26) then begin
    lines26=hitranreader(nui,nuf,co_26[1],5,1)
    alines=[alines,lines26]
  endif
  if keyword_set(co_36) then begin
    lines36=hitranreader(nui,nuf,co_36[1],5,2)
    alines=[alines,lines36]
  endif
  if keyword_set(co_28) then begin
    lines28=hitranreader(nui,nuf,co_28[1],5,3)
    alines=[alines,lines28]
  endif  
  if keyword_set(co_27) then begin
    lines27=hitranreader(nui,nuf,co_27[1],5,4)
    alines=[alines,lines27]
  endif
  if keyword_set(ch4_211) then begin
    lines211=hitranreader(nui,nuf,ch4_211[1],6,1)
    alines=[alines,lines211]
  endif
  if keyword_set(ch4_311) then begin
    lines311=hitranreader(nui,nuf,ch4_311[1],6,2)
    alines=[alines,lines311]
  endif
  if keyword_set(ch4_212) then begin
    lines212=hitranreader(nui,nuf,ch4_212[1],6,3)
    alines=[alines,lines212]
  endif
  if keyword_set(ch4_312) then begin
    lines312=hitranreader(nui,nuf,ch4_312[1],6,4)
    alines=[alines,lines312]
  endif
  if keyword_set(no_46) then begin
    lines46=hitranreader(nui,nuf,no_46[1],8,1)
    alines=[alines,lines46]
  endif
  if keyword_set(no_56) then begin
    lines56=hitranreader(nui,nuf,no_56[1],8,2)
    alines=[alines,lines56]
  endif  
  if keyword_set(no_48) then begin
    lines48=hitranreader(nui,nuf,no_48[1],8,3)
    alines=[alines,lines48]
  endif
  if keyword_set(hf_19) then begin
    lines19=hitranreader(nui,nuf,hf_19[1],14,1)
    alines=[alines,lines19]
  endif
  if keyword_set(hf_29) then begin
    lines29=hitranreader(nui,nuf,hf_29[1],14,2)
    alines=[alines,lines29]
  endif
  if keyword_set(hocl_165) then begin
    lines165=hitranreader(nui,nuf,hocl_165[1],21,1)
    alines=[alines,lines165]
  endif
  if keyword_set(hocl_167) then begin
    lines167=hitranreader(nui,nuf,hocl_167[1],21,2)
    alines=[alines,lines167]
  endif
  if keyword_set(h2s_121) then begin
    lines121=hitranreader(nui,nuf,h2s_121[1],31,1)
    alines=[alines,lines121]
  endif
  if keyword_set(h2s_141) then begin
    lines141=hitranreader(nui,nuf,h2s_141[1],31,2)
    alines=[alines,lines141]
  endif
  if keyword_set(h2s_131) then begin
    lines131=hitranreader(nui,nuf,h2s_131[1],31,3)
    alines=[alines,lines131]
  endif
  if keyword_set(nh3_411) then begin
    lines411=hitranreader(nui,nuf,nh3_411[1],11,1)
    alines=[alines,lines411]
  endif
endelse

if alines eq !null then begin
  out=dblarr(3,2)
  out[0,*]=[nui,nuf]
  out[1,*]=[1d0,1d0]
endif else begin

if keyword_set(writeln) then begin
  
  srt=sort(alines.nu)
  
  openw,lun,writeln,/get_lun
  printf,lun,alines[srt],format='(I02,I01,F22,2G10.4,2F12,F22,F12,F12,A60,A6,A12,A1,F12,F12)'
  free_lun,lun
  
endif

;Housekeeping to make sure the lines are within the simulated spectrum

sep=max(alines.nu)-min(alines.nu)
if sep ge nuf-nui then nusep=sep else nusep=nuf-nui



nume=mean(alines.nu)


if not keyword_set(nus) then begin
  Nnu=3.*nusep/res
  Nnu=Nnu - Nnu mod 1
  nus=nume+(dindgen(Nnu)-Nnu/2)*res

endif else begin
  Nnu=n_elements(nus)
endelse

n=n_elements(alines.nu)

;This is the 'x-axis', ie the spectral range over which simulation will occur

i0=nus*0d0+1d0
sh0=nus*0d0

;simulate the lines in the order in which they appear in the array of line structures

for i=0,n-1 do begin
  
  ;select the correct concentration and mass based on the isotopologue
  
  mois=[alines[i].mol,alines[i].iso]
  
  if mois[0] eq 1 AND mois[1] eq 1 then begin
    mass=.018D
    conc=h2o_161[0]
  endif
  if mois[0] eq 1 AND mois[1] eq 2 then begin
    mass=.02D
    conc=h2o_181[0]
  endif
  if mois[0] eq 1 AND mois[1] eq 3 then begin
    mass=.019D
    conc=h2o_171[0]
  endif
  if mois[0] eq 1 AND mois[1] eq 4 then begin
    mass=.019D
    conc=h2o_162[0]
  endif
  if mois[0] eq 1 AND mois[1] eq 5 then begin
    mass=.021D
    conc=h2o_182[0]
  endif
  if mois[0] eq 1 AND mois[1] eq 6 then begin
    mass=.020D
    conc=h2o_172[0]
  endif
  if mois[0] eq 1 AND mois[1] eq 7 then begin
    mass=.020D
    conc=h2o_262[0]
  endif
  if mois[0] eq 2 and mois[1] eq 1 then begin
    mass=.044d
    conc=co2_626[0]
  endif
  if mois[0] eq 2 and mois[1] eq 2 then begin
    mass=.045d
    conc=co2_636[0]
  endif
  if mois[0] eq 2 and mois[1] eq 3 then begin
    mass=.046d
    conc=co2_628[0]
  endif
  if mois[0] eq 2 and mois[1] eq 4 then begin
    mass=.045d
    conc=co2_627[0]
  endif
  if mois[0] eq 3 and mois[1] eq 1 then begin
    mass=.048d
    conc=o3_666[0]
  endif
  if mois[0] eq 3 and mois[1] eq 2 then begin
    mass=.050d
    conc=o3_668[0]
  endif
  if mois[0] eq 3 and mois[1] eq 3 then begin
    mass=.050d
    conc=o3_686[0]
  endif
  if mois[0] eq 3 and mois[1] eq 4 then begin
    mass=.049d
    conc=o3_667[0]
  endif
  if mois[0] eq 3 and mois[1] eq 5 then begin
    mass=.049d
    conc=o3_676[0]
  endif
  if mois[0] eq 4 and mois[1] eq 1 then begin
    mass=.044d
    conc=n2o_446[0]
  endif
  if mois[0] eq 4 and mois[1] eq 2 then begin
    mass=.045d
    conc=n2o_456[0]
  endif
  if mois[0] eq 4 and mois[1] eq 3 then begin
    mass=.045d
    conc=n2o_546[0]
  endif
  if mois[0] eq 4 and mois[1] eq 4 then begin
    mass=.046d
    conc=n2o_448[0]
  endif
  if mois[0] eq 4 and mois[1] eq 5 then begin
    mass=.045d
    conc=n2o_447[0]
  endif
  if mois[0] eq 5 and mois[1] eq 1 then begin
    mass=.028d
    conc=co_26[0]
  endif
  if mois[0] eq 5 and mois[1] eq 2 then begin
    mass=.029d
    conc=co_36[0]
  endif
  if mois[0] eq 5 and mois[1] eq 3 then begin
    mass=.030d
    conc=co_28[0]
  endif
  if mois[0] eq 5 and mois[1] eq 4 then begin
    mass=.029d
    conc=co_27[0]
  endif
  if mois[0] eq 6 and mois[1] eq 1 then begin
    mass=.016d
    conc=ch4_211[0]
  endif
  if mois[0] eq 6 and mois[1] eq 2 then begin
    mass=.017d
    conc=ch4_311[0]
  endif
  if mois[0] eq 6 and mois[1] eq 3 then begin
    mass=.017d
    conc=ch4_212[0]
  endif
  if mois[0] eq 6 and mois[1] eq 4 then begin
    mass=.018d
    conc=ch4_312[0]
  endif
  if mois[0] eq 8 and mois[1] eq 1 then begin
    mass=.030d
    conc=no_46[0]
  endif
  if mois[0] eq 8 and mois[1] eq 2 then begin
    mass=.031d
    conc=no_56[0]
  endif
  if mois[0] eq 8 and mois[1] eq 3 then begin
    mass=.032d
    conc=no_48[0]
  endif
  if mois[0] eq 11 and mois[1] eq 1 then begin
    mass=.017d
    conc=nh3_411[0]
  endif
  if mois[0] eq 14 and mois[1] eq 1 then begin
    mass=.020d
    conc=hf_19[0]
  endif
  if mois[0] eq 14 and mois[1] eq 2 then begin
    mass=.021d
    conc=hf_29[0]
  endif
  if mois[0] eq 21 and mois[1] eq 1 then begin
    mass=.052d
    conc=hocl_165[0]
  endif
  if mois[0] eq 21 and mois[1] eq 2 then begin
    mass=.054d
    conc=hocl_167[0]
  endif
  if mois[0] eq 31 and mois[1] eq 1 then begin
    mass=.034d
    conc=h2s_121[0]
  endif
  if mois[0] eq 31 and mois[1] eq 2 then begin
    mass=.036d
    conc=h2s_141[0]
  endif
  if mois[0] eq 31 and mois[1] eq 3 then begin
    mass=.035d
    conc=h2s_131[0]
  endif
  
  shif=alines[i].nu+alines[i].shift*p/p0
  o=where(nus ge shif)
  m=min([o[0],Nnu-1-o[0]])
  s=o[0]-m
  f=o[0]+m
  nu=nus[s:f]
  partp=conc*100D0^2*R*T/6.0221415d23
  if not keyword_set(setw) then begin
    gd=(alines[i].nu/c)*((2d0*R*T/mass)^(.5d0))                                        ;20130819 - Removed sqrt(alog(2)) from this expression. I believe it is now correct for a gaussian?
    if keyword_set(sbc) then gc=(alines[i].fbc*(p-partp)/P0+alines[i].sbc*partp/P0)*(T0/T)^alines[i].bctdep else gc=(alines[i].fbc*(p-partp)/P0)*(T0/T)^alines[i].bctdep
  endif else begin
    gd=setw[0]
    gc=setw[1]
  endelse
  ;print,p0,p,partp
  ;stop
  print,gd,gc
  print,alines[i].nu
  ;print,shif
  nun=nu-shif
  
  dop=(1D0/(gd*sqrt(3.1415926536D0)))*exp(-((nun)/gd)^2)            ;Calculate the doppler lineshape
  lor=(gc/3.1415926536D0)/((nun)^2+gc^2)                        ;Calculate the lorentz lineshape
  
  a=gc/gd
  print,a
  ;print,gc,gd,a
  ;stop
  if conc ne 0 then begin
    if not keyword_set(nocon) then begin
      vgt=convol(dop,lor,/edge_truncate,center=1)                     ;Estimate the voigt lineshape by convoluting the doppler and lorentz. !!This is not the best way to do this, but it is quick!!
      nvgt=vgt*res  
    endif else begin
      ;a=gc/(4d0*3.14159d0*gd)                                        ;not sure why idl's voigt function recommends this
      a=gc/gd
      u=(nus-shif)/gd
      
      nvgt=voigt(a,u)/(gd*sqrt(3.14159d0))
    endelse
  endif
  ;stop
                                                      ;This may result in a few percent error at worst. Deviation seems to be mostly in the wings.
  
  ;ts=nun/gd
  ;dp=exp(-ts^2)
  ;xvgt=dblarr(n_elements(ts))
  ;x=ts
  ;for i=0, n_elements(nun)-1 do begin
  ;  lr=1D/(a^2+(x[i]-ts)^2)
  ;  xvgt[i]=(a/(3.1415926535D)^(1.5))*int_tabulated(ts,dp*lr)
    
  ;endfor
  
  ;p=plot(nun,nvgt)
  ;p=plot(nun,xvgt,'b')
  ;stop
    
  if keyword_set(nocon) then print,'*',int_tabulated(nus,nvgt) else print,'*',int_tabulated(nun,nvgt)
  ;stop
  
  st=stpar(alines[i].mol,alines[i].iso,T,alines[i].en,alines[i].nu,alines[i].s)
  
  ;print,st
  
  shape=exp(-st*nvgt*conc*pl)                                       ;Calculate the transmittance based on the Beers-Lambert law
  ;stop
    if keyword_set(nocon) then begin
    s=0
    f=nnu-1
  endif
  print,s,f
  if keyword_set(voigt) then sh0[s:f]=sh0[s:f]+nvgt
  if keyword_set(pexp) then sh0[s:f]=sh0[s:f]+st*nvgt*conc
  if keyword_set(nexp) then sh0[s:f]=sh0[s:f]+st*nvgt*pl
  
  i0[s:f]=i0[s:f]*shape
  ;stop
  ;p=plot(nu,i0[s:f])
  
  ;print,systime(1)-time

endfor

out=dblarr(3,nnu)

out[0,*]=nus
out[1,*]=i0
out[2,*]=sh0

endelse
;print,alines

print,systime(1)-time

return,out


end
