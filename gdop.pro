function gdop,m,T,nu,sigma=sigma

;Calculates the doppler broadening coefficient in cm-1

c=299792458D0                  ;speed of light in m/s
R=8.3144621D0                  ;Universal Gas Constant

;m in kg/mol, T in kelvin, nu in cm-1

gdop=(nu/c)*sqrt(2D0*alog(2D0)*R*T/m)

if keyword_set(sigma) then gdop=gdop/sqrt(alog(2D0))

return,gdop

end