function vgtaprx,nu,nu0,gv,gc

;nu equals domain
;nu0 is the linecenter
;gv is the voigt gamma
;gc is the collisional gamma


dv=gc/gv

lv=1./(2*gv*(1.065+0.477*dv+0.058*dv^2))

vgt=lv*((1-dv)*exp(-2.722*((nu-nu0)/(2*gv))^2)+gc*gv/(gv^2+(nu-nu0)^2)+0.016*dv*(1-dv)*(exp(-0.4*((abs(nu-nu0)/(2*gv))^2.25))-10./(10.+(abs(nu-nu0)/(2*gv))^2.25)))

return,vgt

end