function g2gv,gcol,gdop

;takes doppler and lorentz coefficients and estimates voigt coefficient

gv=0.5346*gcol+sqrt(0.2166*gcol^2+gdop^2)

return,gv

end