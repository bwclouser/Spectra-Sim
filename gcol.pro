function gcol,p,T,fbc,bcdep

;Estimates the lorentz broadening coefficient in cm-1 given
;p, pressure in hPa
;T, temperature in Kelvin
;fbc, foreign broadening coefficient in cm-1/atm
;bcdep, exponent governing temperature dependence

T0=296D                           ;reference temperature of HITRAN
P0=1013.25D                          ;1 Atm in hPa

gcor=(fbc*p/P0)*(T0/T)^bcdep      ;Estimates lorentz broadening ** assumes p_abs is small enough that self-broadening is negligeable
                                  ;and that there is no need to write (p-p_abs)/P0

return,gcor

end