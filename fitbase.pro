function fitbase,wn,data,lines,p,t,mult

T0=296D                           ;reference temperature of HITRAN
P0=1013.25D                       ;1 Atm in hPa

nlines=n_elements(lines.nu)

widths=dblarr(nlines)

widths=mult*(lines.fbc*p/P0)*(T0/T)^lines.bctdep

ind=findgen(n_elements(wn))

;stop

for i=0,nlines-1 do begin

  o=where(wn[ind] lt lines[i].nu-widths[i] OR wn[ind] gt lines[i].nu+widths[i])
  
  ind=ind[o]
  
  ;print,o
  
endfor



return, ind

end