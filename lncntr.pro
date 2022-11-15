function lncntr,T,P,pl,hdata,conc=conc,mixr=mixr,acc=acc,noise=noise,serr=serr,berr=berr

;Uses approximation to voigt lineshape to estimate the line center absorption
;
;T    -   Gas Temperature
;P    -   Gas Pressure
;pl   -   Path Length
;hdata-   HiTran Structure
;conc -   Concentration in mol/cm^3
;mixr -   Mixing ratio in ppm
;acc  -   Use high accuracy formula; Set this equal to the limits of integration in the domain
;10 is a safe starting value for acc
;noise-   if this flag is set to a value, that value is treated as the noise on a scan and the
;output of the function is a concentration in molecules/cm^3

T0=296D0                           ;reference temperature of HITRAN
P0=1013.25D0                       ;1 Atm in hPa
c=299792458D0                      ;speed of light in m/s
R=8.3144621D0                      ;Universal Gas Constant



if keyword_set(conc) then con=conc
if keyword_set(mixr) then con=mixr/1000000D0*catm(P,T)

if hdata.mol eq 1 then begin
  if hdata.iso eq 1 then m=0.018d0
  if hdata.iso eq 2 then m=0.020d0
  if hdata.iso eq 3 then m=0.019d0
  if hdata.iso eq 4 then m=0.019d0
  if hdata.iso eq 5 then m=0.021d0
endif

if hdata.mol eq 2 then begin
  if hdata.iso eq 1 then m=0.044d0
endif

if hdata.mol eq 4 then begin
  if hdata.iso eq 1 then m=0.044d0
endif

if hdata.mol eq 6 then begin
  if hdata.iso eq 1 then m=0.016d0
  if hdata.iso eq 2 then m=0.017d0
  if hdata.iso eq 3 then m=0.017d0
  if hdata.iso eq 4 then m=0.018d0
endif

if keyword_set(berr) then gcor=(hdata.fbc*p/P0)*(T0/T)^(berr*hdata.bctdep) else gcor=(hdata.fbc*p/P0)*(T0/T)^hdata.bctdep      ;Estimates lorentz broadening ** assumes p_abs is small enough that self-broadening is negligeable
                                               ;and that there is no need to write (p-p_abs)/P0                              

gdop=((hdata.nu)/c)*sqrt(2D0*R*T/m)               ;m in kg/mol, T in kelvin, nu in cm-1

gv=0.5346*gcor+sqrt(0.2166*gcor^2+alog(2d0)*gdop^2)
;print,gv
dv=gcor/gv
;print,gdop,gcor,gv
if not keyword_set(acc) then begin
  lv=1./(2*gv*(1.065+0.477*dv+0.058*dv^2))
endif else begin
  a=gcor/gdop
  print,size(a)
  lv=exp(a^2)/(gdop*sqrt(3.1415926535d0))*(erf(acc+a)-erf(a))
endelse
;print,stpar(hdata.mol,hdata.iso,T,hdata.en,hdata.nu,hdata.s)*lv
if keyword_set(serr) then abs=1.-exp(-stpar(hdata.mol,hdata.iso,T,hdata.en,hdata.nu,hdata.s*serr)*lv*pl*con) else abs=1.-exp(-stpar(hdata.mol,hdata.iso,T,hdata.en,hdata.nu,hdata.s)*lv*pl*con)

if keyword_set(noise) then abs=alog(1.-noise)/(-stpar(hdata.mol,hdata.iso,T,hdata.en,hdata.nu,hdata.s)*lv*pl)

return,abs

end
