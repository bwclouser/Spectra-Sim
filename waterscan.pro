function waterscan,smin,water,nui,nuf

hdo=[]
h18o=[]

str={w:ulong(0),h18:intarr(20),hd:intarr(20)}

i0=where(water.nu ge nui)
i0=i0[0]
ifin=where(water.nu ge nuf)
ifin=ifin[0]

n=n_elements(water)

out=replicate(str,ifin-i0+1)

for i=i0,ifin-1 do begin

  if water[i].s ge smin then begin
  
    ;look left
    
    j=i
    
    while water[j].nu ge water[i].nu-2. do begin
      
      print,j
      if water[j].iso eq 2 AND water[j].s ge smin then h18o=[h18o,[i,j]]
      if water[j].iso eq 4 AND water[j].s ge smin then hdo=[hdo,[[i],[j]]]
      j--
      
    endwhile
    
    ;look right
    
    j=i
    
    while water[j].nu le water[i].nu+2. do begin
      
      if water[j].iso eq 2 AND water[j].s ge smin then h18o=[h18o,[i,j]]
      if water[j].iso eq 4 AND water[j].s ge smin then hdo=[hdo,[[i],[j]]]
      j++
      
    endwhile
    
    ;print,n_elements(h18o),n_elements(hdo)
    
    ;out[i-i0].w=i
    ;if h18o ne !null then out[i-i0].h18=h18o
    ;if hdo ne !null then out[i-i0].hd=hdo
    
    endif
    
 endfor
 
 return,hdo
 
 end
    