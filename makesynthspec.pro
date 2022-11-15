pro makesynthspec,wnbase,spec,bline,pres,temp,conc,nscans,nsamps,tc,mloss,len,range,tol,rate,approx=approx,polysmp=polysmp,lines=lines,icos=icos,direct=direct,writedir=writedir,linear=linear,breath=breath,sine=sine,setco2=setco2,basind=basind,noise=noise,hf=hf

;Generates synthetic ICOS spectra

;requires spectrum, writesbasefile, array2cpci, catm

;wnbase       -       output vector of wavenumbers for the given tuning curve
;spec         -       nsamps x nscans array containing spectra
;bline        -       nsamps x nscans array containing baselines
;pres         -       pressure in mb. Can be a scalar or vector of length nscans
;temp         -       temperature in K. Can be a scalar or vector of length nscans
;conc         -       H2O mixing ratio in ppm. Can be a scalar or vector of length nscans
;nscans       -       Number of scans that will be output. If pres, temp, or conc are entered as scalars, the program assumes a constant value over all scans
;nsamps       -       Number of samples per scan
;tc           -       Tuning curve parameters
;mloss        -       Fractional mirror loss
;len          -       Cavity length in meters
;range        -       [lowest, highest] in wavenumbers
;tol          -       ICOS tolerance, 1e-6 is a safe number
;rate         -       Sampling rate in Hz

;approx       -       Number of passes per calculation time interval. Defaults to 1, which means every mirror reflection is modelled... this is unnecesarily expensive. Values around 25 seem safe.
;polysmp      -       Polynomial coefficients as a function of wavenumber/1000
;lines        -       Structure of HITRAN parameters corresponding to a set of lines. Defaults to finding all lines within the range that are above a certain threshhold.
;icos         -       Models ICOS output
;direct       -       Models direct absorption output. Direct defaults to one, otherwise treated as number of passes through a Herriot cell of length len. 
;writedir     -       Writes CPCI structure, baseline file, and pte file to the specified directory
;linear       -       Calculates a linear tuning curve between the range[0] and range[1]
;noise        -       Adds fractional noise with value noise to the data


c=299792458d0
etfsr=-0.0172798d0  

if n_elements(pres) eq 1 then begin
  pconst=pres
  pres=dblarr(nscans)
  pres[*]=pconst
endif

if n_elements(temp) eq 1 then begin
  tconst=temp
  temp=dblarr(nscans)
  temp[*]=tconst
endif

if n_elements(conc) eq 1 then begin
  cconst=conc
  conc=dblarr(nscans)
  conc[*]=cconst
endif

if keyword_set(hf) then hfppt=hf else hfppt=1d0



;-----------------------------------------;
;-----------Generate Time/WN base---------;
;-----------------------------------------;

sampdt=1d0/rate                           ;sampling time
sampnum=dindgen(nsamps)                   ;vector of sample
tbase=sampnum*sampdt                      ;vector of sampling times

if keyword_set(linear) then begin          ;Assumes a linear tuning curve
  dwn=((range[0]-range[1])/etFSR)/(nsamps*1d-3)
  tc=dblarr(8)
  tc[2]=dwn
  tc[5]=1d0
  tc[7]=1d0
  ;stop
endif

                                          ;calculate tuning curve

snred=(sampnum-tc[0])/1000d0
frin=tc[1]+tc[2]*snred+tc[3]*snred^2+tc[4]*EXP(-snred/tc[5])+tc[6]*EXP(-snred/tc[7])
wnbase=frin*etFSR+range[1]

;-----------------------------------------;
;-----------Generate Input Ramp-----------;
;-----------------------------------------;


if keyword_set(polysmp) then begin
  
  ;calculates the shape of the input ramp given the polynomial coefficients
  ;in ICOS mode this ramp is further summed to get the output ramp shape
  ;in direct absorption mode this is taken to be the absorption-free output ramp as well.
  
  I0=dblarr(nsamps)
  npoly=n_elements(polysmp)
  for i=0,npoly-1 do I0+=polysmp[i]*(sampnum/1d3)^i
  
  if keyword_set(sine) then I0*=(1d0+sine[0]*sin(sine[1]*sampnum))
  
endif

if keyword_set(polywn) then begin
  ;not implemented yet
endif

;-----------------------------------------;
;-----------Generate Calc Vars------------;
;-----------------------------------------;


if keyword_set(icos) then begin
 

  refl=1d0-mloss                          ;mirror reflectivity calculated from mirror loss
  
  passdt=len/c                            ;time it takes for for light to traverse cell once
  
  passes=sampdt/passdt                    ;number of passes light makes through the cell in one sampling interval
  
  if not keyword_set(approx) then begin
    nn=1d0
    pscale=1d0
  endif else begin                        ;sets the number of passes per calculation time interval, and calculates the power scale. The latter is a factor that accounts for summing over the nn passes.
    nn=approx
    pscale=(1d0-refl^(2d0*nn))/(1d0-refl^2d0)
  endelse

  nend=long(round((alog(tol)/(2d0*alog(refl)))/nn)) ;the number of calculation intervals the program models

  csamps=ceil(nsamps*passes/nn)           ;samples in the calculation base
  ;print,csamps,nsamps*passes/nn 
  psamps=dindgen(csamps)                  ;vector of sample number in calculation base

  ctbase=dindgen(csamps)*passdt*nn        ;time of calculation intervals
  
  cI0=interpol(I0,tbase,ctbase)           ;interpolates baseline into calculation base
  
  cwnbase=interpol(wnbase,tbase,ctbase)   ;interpolates tuning curve into calculation base
  ;stop
  spec=fltarr(nsamps,nscans)
  bline=fltarr(nsamps,nscans)
  ;stop
endif


;-----------------------------------------;
;-----------Generate Spectra--------------;
;-----------------------------------------;


;;Simulates ICOS spectra;;

if keyword_set(icos) then begin
  
  cwnbase=reverse(cwnbase)                ;reverses tuning curve to account for ramp direction
  
  for i=0,nscans-1 do begin
    
    blinec=dblarr(csamps)
    outc=dblarr(csamps)
    
    ;for non-H2O species below I've used typical atmospheric values
    
    if keyword_set(breath) then begin
      ch2o=[conc[i]/1e6*catm(pres[i],temp[i]),1e-29]
      ch18=[conc[i]/1d6*catm(pres[i],temp[i]),1e-29]
      ch17=ch18
      chdo=ch18
      chd18o=ch18
      chd17o=ch18
      cd2o=ch18
      cco2=[5d4/1d6*catm(pres[i],temp[i]),1e-28]
      co3=[.5d0/1d6*catm(pres[i],temp[i]),1e-28]
      cn2o=[.31d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cn2o=[600000d0/1d6*catm(p,t),1d-32]
      cco1=[.04d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cco1=[10d0/1d6*catm(pres[i],temp[i]),1e-26]
      cch4=[1.68d0/1d6*catm(pres[i],temp[i]),1e-24]
      cno=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
      chf=[hfppt/1d12*catm(pres[i],temp[i]),1e-26]
      ;chf=[10000d0/1d12*catm(pres[i],temp[i]),1e-26]
      chocl=[1d0/1d9*catm(pres[i],temp[i]),1e-22]
      ch2s=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
    endif else begin
      ch2o=[conc[i]/1e6*catm(pres[i],temp[i]),1e-29]
      ch18=[conc[i]/1d6*catm(pres[i],temp[i]),1e-29]
      ch17=ch18
      chdo=ch18
      chd18o=ch18
      chd17o=ch18
      cd2o=ch18
      if keyword_set(setco2) then cco2=[setco2/1d6*catm(pres[i],temp[i]),1e-28] else cco2=[400d0/1d6*catm(pres[i],temp[i]),1e-28]
      co3=[.5d0/1d6*catm(pres[i],temp[i]),1e-28]
      cn2o=[.31d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cn2o=[600000d0/1d6*catm(p,t),1d-32]
      cco1=[.04d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cco1=[10d0/1d6*catm(pres[i],temp[i]),1e-26]
      cch4=[1.68d0/1d6*catm(pres[i],temp[i]),1e-24]
      cno=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
      chf=[hfppt/1d12*catm(pres[i],temp[i]),1e-26]
      ;chf=[10000d0/1d12*catm(pres[i],temp[i]),1e-26]
      chocl=[1d0/1d9*catm(pres[i],temp[i]),1e-22]
      ch2s=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
    endelse
    wres=0                        ;not used in this implementation of spectrum.pro since we're supplying our own wavenumber base
    ;stop
    
    if (i ge 1 and i le nscans-1) then begin
      if conc[i-1] eq conc[i] and pres[i-1] eq pres[i] and temp[i-1] eq temp[i] then begin
        spec[*,i]=spec[*,i-1]/spnoise
        bline[*,i]=bline[*,i-1]
      endif else begin
        ;o21 contains the single pass absorption of the cell
        if keyword_set(lines) then begin
          o21=spectrum(range[0],range[1],1e-28,temp[i],pres[i],len*1d2,wres,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,nus=cwnbase,lines=lines,/pexp,/nocon)
        endif else begin
          o21=spectrum(range[0],range[1],1e-28,temp[i],pres[i],len*1d2,wres,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,nus=cwnbase,/pexp,/nocon)
        endelse
        
        oneabs=reverse(reform(o21[1,*],n_elements(o21[0,*])))
    
        ;stop
        ;This for loop sums the modeled spectra of different path length to generate the observed spectra and theoretical absorption-free baseline
        for j=0l,nend-1 do begin
          outc+=pscale*shift(refl^(2d0*nn*j)*cI0*oneabs^(2d0*nn*j+1),2l*j)
          blinec+=pscale*shift(refl^(2d0*nn*j)*cI0,2l*j)
          ;if (j mod 10)  eq 0 then stop
        endfor
        
        spec[*,i]=interpol(outc,ctbase,tbase)
        bline[*,i]=interpol(blinec,ctbase,tbase)
      endelse
    endif else begin
      if keyword_set(lines) then begin
        o21=spectrum(range[0],range[1],1e-28,temp[i],pres[i],len*1d2,wres,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,nus=cwnbase,lines=lines,/pexp,/nocon)
      endif else begin
        o21=spectrum(range[0],range[1],1e-28,temp[i],pres[i],len*1d2,wres,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,nus=cwnbase,/pexp,/nocon)
      endelse
      ;stop
      oneabs=reverse(reform(o21[1,*],n_elements(o21[0,*])))

      ;stop
      ;This for loop sums the modeled spectra of different path length to generate the observed spectra and theoretical absorption-free baseline
      for j=0l,nend-1 do begin
        outc+=pscale*shift(refl^(2d0*nn*j)*cI0*oneabs^(2d0*nn*j+1),2l*j)
        blinec+=pscale*shift(refl^(2d0*nn*j)*cI0,2l*j)
        ;if (j mod 10)  eq 0 then stop
      endfor
      ;stop
      spec[*,i]=interpol(outc,ctbase,tbase)
      bline[*,i]=interpol(blinec,ctbase,tbase)
    endelse
    if keyword_set(noise) then begin
      spnoise=(1.0+noise*randomn(seed,nsamps))
    endif else begin
      spnoise=1.0
    endelse
    spec[*,i]*=spnoise
    
  endfor
  
  if keyword_set(writedir) then begin
    pte=dblarr(11,nscans)
    pte[0,*]=indgen(nscans)+1
    pte[1,*]=pres*.75006d0
    pte[2,*]=temp
    for i=0,nscans-1 do begin
      pte[3:10,i]=tc
    endfor
    spec[0:124,*]=0d0
    array2cpci,float(spec),writedir,0,nscans
    openw,wlun,writedir+'dir0/PTE.txt',/get_lun
    printf,wlun,pte,format='(I5,10F20.10)'
    free_lun,wlun
    writesbasefile,writedir+'dir0/sbase.'+string(npoly,format='(I02)')+'.ptb',0l,1l,1d3,range[1],0d0,fix(0),fix(1751),fix(npoly),fix(0),0,polysmp/1d3,0d0;,mirloss=mloss
  endif
  
endif


;;Simulates direct absorption spectra;;

if keyword_set(direct) then begin
  
  npasses=direct
  wnbase=reverse(wnbase)
  fI0=float(I0)
  
  spec=fltarr(nsamps,nscans)
  bline=fltarr(nsamps,nscans)

  
  for i=0,nscans-1 do begin
    
    if keyword_set(breath) then begin
      ch2o=[conc[i]/1e6*catm(pres[i],temp[i]),1e-29]
      ch18=[conc[i]/1d6*catm(pres[i],temp[i]),1e-29]
      ch17=ch18
      chdo=ch18
      chd18o=ch18
      chd17o=ch18
      cd2o=ch18
      cco2=[5d4/1d6*catm(pres[i],temp[i]),1e-28]
      co3=[.5d0/1d6*catm(pres[i],temp[i]),1e-28]
      cn2o=[.31d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cn2o=[600000d0/1d6*catm(p,t),1d-32]
      cco1=[.04d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cco1=[10d0/1d6*catm(pres[i],temp[i]),1e-26]
      cch4=[1.68d0/1d6*catm(pres[i],temp[i]),1e-24]
      cno=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
      chf=[hfppt/1d12*catm(pres[i],temp[i]),1e-26]
      ;chf=[10000d0/1d12*catm(pres[i],temp[i]),1e-26]
      chocl=[1d0/1d9*catm(pres[i],temp[i]),1e-22]
      ch2s=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
    endif else begin
      ch2o=[conc[i]/1e6*catm(pres[i],temp[i]),1e-29]
      ch18=[conc[i]/1d6*catm(pres[i],temp[i]),1e-29]
      ch17=ch18
      chdo=ch18
      chd18o=ch18
      chd17o=ch18
      cd2o=ch18
      cco2=[400d0/1d6*catm(pres[i],temp[i]),1e-28]
      co3=[.5d0/1d6*catm(pres[i],temp[i]),1e-28]
      cn2o=[.31d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cn2o=[600000d0/1d6*catm(p,t),1d-32]
      cco1=[.04d0/1d6*catm(pres[i],temp[i]),1e-26]
      ;cco1=[10d0/1d6*catm(pres[i],temp[i]),1e-26]
      cch4=[1.68d0/1d6*catm(pres[i],temp[i]),1e-24]
      cno=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
      chf=[hfppt/1d12*catm(pres[i],temp[i]),1e-26]
      ;chf=[10000d0/1d12*catm(pres[i],temp[i]),1e-26]
      chocl=[1d0/1d9*catm(pres[i],temp[i]),1e-22]
      ch2s=[1d0/1d9*catm(pres[i],temp[i]),1e-26]
    endelse
    if keyword_set(lines) then begin
      ;stop
      o21=spectrum(range[0],range[1],1e-28,temp[i],pres[i],len*1d2*npasses,wres,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,nus=wnbase,/sbc,lines=lines,/pexp,/nocon)
    endif else begin
      o21=spectrum(range[0],range[1],1e-28,temp[i],pres[i],len*1d2*npasses,wres,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,nus=wnbase,/sbc,/pexp,/nocon)
    endelse

     absorp=reverse(reform(o21[1,*],n_elements(o21[0,*])))
     
     absorp[0:124]=0d0
     
     spec[0,i]=float(absorp*I0)
     bline[0,i]=fI0

  endfor
  
  wnbase=reverse(wnbase)
  
  if keyword_set(noise) then spec=temporary(spec)*(1.0+float(noise))
  
  if keyword_set(writedir) then begin
    pte=dblarr(11,nscans)
    pte[0,*]=indgen(nscans)+1
    pte[1,*]=pres*.75006d0
    pte[2,*]=temp
    for i=0,nscans-1 do begin
      pte[3:10,i]=tc
    endfor
    spec[0:124,*]=0d0
    array2cpci,float(spec),writedir,0,nscans
    openw,wlun,writedir+'dir0/PTE.txt',/get_lun
    printf,wlun,pte,format='(I5,10F20.10)'
    free_lun,wlun
    writesbasefile,writedir+'dir0/sbase.'+string(npoly,format='(I02)')+'.ptb',0l,1l,1d3,range[1],0d0,fix(0),fix(1751),fix(npoly),fix(0),0,polysmp,0d0;,mirloss=mloss
  endif


endif  
  
end
