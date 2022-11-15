pro edittiff,in,out,rv,gv,bv,ro,go,bo,rr,gr,br,text=text

inf=size(in)

cn=inf[2]
rn=inf[3]

out=in

for i=0,rn-1 do begin

  if keyword_set(text) then begin
    
    row=in[0,*,i]
    
    xx=where(abs(row-smooth(double(row),text[0])) ge text[1])
    
    out[0,xx,i]=255
    
    row=in[1,*,i]
    
    xx=where(abs(row-smooth(double(row),text[0])) ge text[1])
    
    out[1,xx,i]=255
    
    row=in[2,*,i]
    
    xx=where(abs(row-smooth(double(row),text[0])) ge text[1])
    
    out[2,xx,i]=255
  endif else begin

  redrow=in[0,*,i]
  greenrow=in[1,*,i]
  bluerow=in[2,*,i]
  
  k=where(abs(redrow)-rv le rr and abs(greenrow)-gv le gr and abs(bluerow)-bv le br)
  
  out[0,k,i]=ro
  out[1,k,i]=go
  out[2,k,i]=bo
  
  ;kk=where(redrow ge 255 or greenrow ge 255 or bluerow ge 255)

  ;out[0,kk,i]=254
  ;out[1,kk,i]=254
  ;out[2,kk,i]=254 
  endelse
  
endfor


end