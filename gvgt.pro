function gvgt,p,T,fbc,bcdep,m,nu

T0=296D                           ;reference temperature of HITRAN
P0=1013.25D                       ;1 Atm in hPa
c=299792458D                      ;speed of light in m/s
R=8.3144621D                      ;Universal Gas Constant

gcor=(fbc*p/P0)*(T0/T)^bcdep      ;Estimates lorentz broadening ** assumes p_abs is small enough that self-broadening is negligeable
                                  ;and that there is no need to write (p-p_abs)/P0
                              
gdop=(nu/c)*sqrt(2D*R*T/m)        ;m in kg/mol, T in kelvin, nu in cm-1

gv=0.5346*gcor+sqrt(0.2166*gcor^2+alog(2D)*gdop^2)  ;this should have a factor of ln2 since gdop here is the 'sigma' of the gaussian.

return,gv

end