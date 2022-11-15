pro edittiff2,in,out,sm,cut

inf=size(in)

cn=inf[2]
rn=inf[3]

out=in

ind=indgen(cn)

for i=0,rn-1 do begin
    ;print,i
    row=in[0,*,i]
    
    dv=deriv(smooth(double(row),sm))
    
    bit=1
    
    j=sm
    ;if i eq 1225 then stop
    while(j le cn-1-sm) do begin
      
      if dv[j] ge cut and bit eq 1 then bit=0;and dv[j-1] ge cut and dv[j-2] ge cut and dv[j-3] ge cut then bit=0
      if dv[j] le -cut and bit eq 0 then bit=1; and dv[j-1] le -cut and dv[j-2] le -cut and dv[j-3] le -cut then bit=1
      
      ind[j]=j*bit
      j++
      ;if i eq 1225 then begin
      ;  print,dv[j],cut,bit,j
      ;  stop
     ; endif
    endwhile
    ind[0:sm-1]=indgen(sm)
    ind[cn-sm:cn-1]=sm-indgen(sm)-1
    
    out[0,ind,i]=255
    out[1,ind,i]=255   
    out[2,ind,i]=255
 
endfor


end