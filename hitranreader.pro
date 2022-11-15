function hitranReader,nui,nuf,smin,molnum,isonum,ht2008=ht2008,getall=getall

;returns all lines of a given isotopologue in a target wavenumber range that are above a certain line strength
;sample call: h2o=hitranreader(3787.5,3790.5,1d-26,1,1) 

;nui    -     start of target range
;nuf    -     end of target range
;smin   -     linestrength threshold
;molnum -     HITRAN molecule number
;isonum -     HITRAN isotopologue number

z=systime(1)

htln2={mol:fix(0),iso:fix(0),nu:double(0),s:double(0),r:double(0),fbc:float(0),sbc:float(0),en:double(0),bctdep:float(0),shift:float(0),junk:string(0),errcodes:intarr(6),ref:string(0),linemix:string(0),upperstatwt:float(0),lowerstatwt:float(0)}
;Liz's IDL structure for HITRAN data

hlines=ulong64(2713968)

niso=n_elements(isonum)

num=molnum*10.+isonum

if keyword_set(HT2008) then begin
  openr,lun,'~/Data/Hitran/Hitran08.par',/get_lun 
  hlines=ulong64(2713968)
endif else begin  
  openr,lun,'~/Data/Hitran/5a4d65e0.par',/get_lun
  hlines=ulong64(5399562.)
endelse

;find start of range
;

;lines=[]
lines=LIST()

testnu1=float(0)
testnu2=float(0)

i=-1D

while testnu1 lt nui OR testnu1 ge nuf do begin

  i++
  point_lun,lun,i*162.+3.
  readf,lun,testnu1,format='(F12)'
  i++
  point_lun,lun,i*162.+3.
  readf,lun,testnu2,format='(F12)'
  
  deli=(nui-testnu1)/(testnu2+.1-testnu1)
  ;print,testnu1,testnu2,deli
  ;print,deli
  
  newi=i+long64(deli)
  ;print,'-',newi
  if newi ge hlines then newi=hlines-3
  if newi le 0 then newi=0
  ;print,'-',newi
  i=newi
  
  ;print,i
  
  ;print,i,testnu
  ;stop
endwhile
;print,'*'
while testnu1 ge nui do begin

  i--
  point_lun,lun,i*162.+3
  readf,lun,testnu1
  ;print,'*',testnu1
  ;print,i
  
endwhile

;print,i,testnu1

point_lun,lun,i*162.+3.
readf,lun,testnu1

;print,i,testnu1

while testnu1 le nuf do begin
  
  readf,lun,htln2,format='(I02,I01,F12,2G10.4,2F5,F10,F4,F8,A60,6I1,A12,A1,F7,F7)'
  for i=0,niso-1 do begin
    if 10*htln2.mol+htln2.iso eq num[i] then begin
      ;if htln2.s ge smin then lines=[lines,htln2]
      if htln2.s ge smin then lines.add,htln2
    endif
  endfor
  
  testnu1=htln2.nu
  
endwhile

free_lun,lun

print,systime(1)-z

return,lines.toarray()

end
 
  
  


  