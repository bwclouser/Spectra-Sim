function testspec,hdo,fname,R,top=top,bottom=bottom,ground=ground,species=species,range=range,logabs=logabs,p0=p0



if keyword_set(top) then begin
  
  t=296d0
  p=70d0
  
  ch2o=[3.5d0/1e6*catm(70d0,296d0),1e-24]
  ch18=[.3*3.5d0/1d6*catm(70d0,296d0),1e-24]
  ch17=ch18
  chdo=ch18
  chd18o=ch18
  chd17o=ch18
  cd2o=ch18
  cco2=[370d0/1d6*catm(70d0,296d0),1e-27]
  co3=[.5d0/1d6*catm(70d0,296d0),1e-24]
  cn2o=[.31d0/1d6*catm(70d0,296d0),1e-24]
  cco1=[.04d0/1d6*catm(70d0,296d0),1e-23]
  cch4=[1.68d0/1d6*catm(70d0,296d0),1e-24]
  chf=[1d0/1d12*catm(70d0,296d0),1e-20]
  chocl=[1d0/1d9*catm(70d0,296d0),1e-21]
  ch2s=[1d0/1d9*catm(70d0,296d0),1e-21]
  
endif

if keyword_set(bottom) then begin

  t=296d0
  p=70d0
  
  ch2o=[60d0/1e6*catm(70d0,296d0),1e-24]
  ch18=[.3*60d0/1d6*catm(70d0,296d0),1e-24]
  ch17=ch18
  chdo=ch18
  chd18o=ch18
  chd17o=ch18
  cd2o=ch18
  cco2=[380d0/1d6*catm(70d0,296d0),1e-27]
  co3=[.05d0/1d6*catm(70d0,296d0),1e-24]
  cn2o=[.32d0/1d6*catm(70d0,296d0),1e-24]
  cco1=[.09d0/1d6*catm(70d0,296d0),1e-23]
  cch4=[1.7d0/1d6*catm(70d0,296d0),1e-24]
  chf=[1d0/1d12*catm(70d0,296d0),1e-20]
  chocl=[1d0/1d9*catm(70d0,296d0),1e-21]
  ch2s=[1d0/1d9*catm(70d0,296d0),1e-21]

endif

if keyword_set(ground) then begin
  
  p=20d0
  t=295d0
  
  ch2o=[ground/1e6*catm(p,t),1e-28]
  ch18=[ground/1d6*catm(p,t),1e-28]
  ch17=ch18
  chdo=ch18
  chd18o=ch18
  chd17o=ch18
  cd2o=ch18
  cco2=[370d0/1d6*catm(p,t),1e-28]
  co3=[.5d0/1d6*catm(p,t),1e-28]
  cn2o=[.31d0/1d6*catm(p,t),1e-26]
  cco1=[.04d0/1d6*catm(p,t),1e-25]
  cch4=[1.68d0/1d6*catm(p,t),1e-24]
  cno=[1d0/1d9*catm(p,t),1e-26]
  chf=[1d0/1d12*catm(p,t),1e-23]
  chocl=[1d0/1d9*catm(p,t),1e-22]
  ch2s=[1d0/1d9*catm(p,t),1e-26]
  
endif


pz=[]
pp=[]

n=n_elements(hdo)

if not keyword_set(species) then begin

for i=0,n-1 do begin
  
fnamei=fname+string(i,format='(I03)')+'.dat'

G=R/(1d0-R)

o1=spectrum(hdo[i].nu-5d0,hdo[i].nu+5d0,1e-23,t,p,90d0,0.002d0,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,writeln=fnamei)
o2=spectrum(hdo[i].nu-5d0,hdo[i].nu+5d0,1e-23,t,p,90d0,0.002d0,h2o_162=ch18)

ab1=1d0-o1[1,*]
t1=1d0-G*ab1/(1d0+G*ab1)

ab2=1d0-o2[1,*]
t2=G*ab2/(1d0+G*ab2)
mm=max(t2)
t2=1d0-t2

p=plot(o1[0,*],t1)
p=plot(o2[0,*],t2,'r',/overplot)

p.xrange=[hdo[i].nu-2.5,hdo[i].nu+2.5]
p.yrange=[1.-4*mm,1.]

p.xtitle='Wavenumber $cm^{-1}$'
p.ytitle='Normalized ICOS Signal'
p.title=string(float(hdo[i].nu-2.5),format='(F12)')+' to '+string(float(hdo[i].nu+2.5),format='(F12)')

pp=[pp,p]

endfor

pz=pp

endif

if keyword_set(species) then begin
  
  pl=600000d0
  
  for i=0,n-1 do begin
  
  ;fnamei=fname+string(i,format='(I03)')+'.dat'

  G=R/(1d0-R)
  if keyword_set(range) then ran=range else ran=[hdo[i].nu-5d0,hdo[i].nu+5d0]
  o0=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_161=ch2o,h2o_181=ch18,h2o_171=ch17,h2o_162=chdo,h2o_182=chd18o,h2o_172=chd17o,h2o_262=cd2o,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,hf_19=chf,hf_29=chf,hocl_165=chocl,hocl_167=chocl,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,writeln=fnamei,/nocon)
  o1=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_161=ch2o,/nocon)
  o2=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_181=ch18,/nocon)  
  o3=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_171=ch17,/nocon)
  o4=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_162=chdo,/nocon)
  on=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_182=chd18o,/nocon)  
  on7=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_172=chd17o,/nocon)
  on2=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2o_262=cd2o,/nocon)
  o5=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,co2_626=cco2,co2_636=cco2,co2_628=cco2,co2_627=cco2,/nocon)
  o6=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,o3_666=co3,o3_668=co3,o3_686=co3,o3_667=co3,o3_676=co3,/nocon)
  o7=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,n2o_446=cn2o,n2o_456=cn2o,n2o_546=cn2o,n2o_448=cn2o,n2o_447=cn2o,/nocon)
  o8=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,co_26=cco1,co_36=cco1,co_28=cco1,co_27=cco1,/nocon)
  o9=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,ch4_211=cch4,ch4_311=cch4,ch4_212=cch4,/nocon)
  o10=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,hf_19=chf,hf_29=chf,/nocon)
  o11=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,hocl_165=chocl,hocl_167=chocl,/nocon)
  o12=spectrum(ran[0],ran[1],1e-23,t,p,pl,0.002d0,h2s_121=ch2s,h2s_141=ch2s,h2s_131=ch2s,/nocon)
  
  ;print,o0[0,0],o0[0,n_elements(o0[0,*])-1]
  ;print,o1[0,0],o1[0,n_elements(o1[0,*])-1]
  ;print,o4[0,0],o4[0,n_elements(o4[0,*])-1]
  
  
  
  ;ab0=1d0-o0[1,*]
  ;t0=1d0-G*ab0/(1d0+G*ab0)
  
  ;ab1=1d0-o1[1,*]
  ;t1=1d0-G*ab1/(1d0+G*ab1)

  ;ab2=1d0-o2[1,*]
  ;t2=1d0-G*ab2/(1d0+G*ab2)
  
  ;ab3=1d0-o3[1,*]
  ;t3=1d0-G*ab3/(1d0+G*ab3)
  
  ;ab4=1d0-o4[1,*]
  ;t4=G*ab4/(1d0+G*ab4)
  ;mm=max(t4)
  ;t4=1d0-t4
  
  ;abn=1d0-on[1,*]
  ;tn=G*abn/(1d0+G*abn)
  ;tn=1d0-tn
  
  ;ab5=1d0-o5[1,*]
  ;t5=1d0-G*ab5/(1d0+G*ab5)
  
  ;ab6=1d0-o6[1,*]
  ;t6=1d0-G*ab6/(1d0+G*ab6)
  
  ;ab7=1d0-o7[1,*]
  ;t7=1d0-G*ab7/(1d0+G*ab7) 
  
  ;ab8=1d0-o8[1,*]
  ;t8=1d0-G*ab8/(1d0+G*ab8)
  
  ;ab9=1d0-o9[1,*]
  ;t9=1d0-G*ab9/(1d0+G*ab9)     
  
  ;ab10=1d0-o10[1,*]
  ;t10=1d0-G*ab10/(1d0+G*ab10)
  
  ;ab11=1d0-o11[1,*]
  ;t11=1d0-G*ab11/(1d0+G*ab11)  
  
  ;ab12=1d0-o12[1,*]
  ;t12=1d0-G*ab12/(1d0+G*ab12)
 
  if keyword_set(logabs) then begin
         
    p0=plot(o0[0,*],1d0-o0[1,*],xthick=3,ythick=3,thick=1,name='All',ylog=1)
    p1=plot(o1[0,*],1d0-o1[1,*],xthick=3,ythick=3,thick=2,'b',name='$H_2^{16}O$',/overplot)
    p2=plot(o2[0,*],1d0-o2[1,*],xthick=3,ythick=3,thick=2,'c',name='$H_2^{18}O$',/overplot)
    p3=plot(o3[0,*],1d0-o3[1,*],xthick=3,ythick=3,thick=2,'gray',name='$H_2^{17}O$',/overplot)
    p4=plot(o4[0,*],1d0-o4[1,*],xthick=3,ythick=3,thick=2,'g',name='$HD^{16}O$',/overplot)
    pn=plot(on[0,*],1d0-on[1,*],xthick=3,ythick=3,thick=2,'indigo',name='$HD^{18}O$',/overplot)
    pn7=plot(on7[0,*],1d0-on7[1,*],xthick=3,ythick=3,thick=2,'violet',name='$HD^{17}O$',/overplot)
    pn2=plot(on2[0,*],1d0-on2[1,*],xthick=3,ythick=3,thick=2,'dark slate gray',name='$D_2O$',/overplot)
    p5=plot(o5[0,*],1d0-o5[1,*],xthick=3,ythick=3,thick=1,'r',name='$CO_2$',/overplot)
    p6=plot(o6[0,*],1d0-o6[1,*],xthick=3,ythick=3,thick=1,'orange',name='$O_3$',/overplot)
    p7=plot(o7[0,*],1d0-o7[1,*],xthick=3,ythick=3,thick=1,'brown',name='$N_2O$',/overplot)
    p8=plot(o8[0,*],1d0-o8[1,*],xthick=3,ythick=3,thick=1,'indigo',name='CO',/overplot)
    p9=plot(o9[0,*],1d0-o9[1,*],xthick=3,ythick=3,thick=1,'m',name='$CH_4$',/overplot)
    p10=plot(o10[0,*],1d0-o10[1,*],xthick=3,ythick=3,thick=1,'y',name='HF',/overplot)
    p11=plot(o11[0,*],1d0-o11[1,*],xthick=3,ythick=3,thick=1,'gray',name='HOCl',/overplot)
    p12=plot(o12[0,*],1d0-o12[1,*],xthick=3,ythick=3,thick=1,'olive',name='$H_2S$',/overplot)
  endif else begin
    
    p0=plot(o0[0,*],o0[1,*],xthick=3,ythick=3,thick=1,name='All')
    p1=plot(o1[0,*],o1[1,*],xthick=3,ythick=3,thick=2,'b',name='$H_2^{16}O$',/overplot)
    p2=plot(o2[0,*],o2[1,*],xthick=3,ythick=3,thick=2,'c',name='$H_2^{18}O$',/overplot)
    p3=plot(o3[0,*],o3[1,*],xthick=3,ythick=3,thick=2,'gray',name='$H_2^{17}O$',/overplot)
    p4=plot(o4[0,*],o4[1,*],xthick=3,ythick=3,thick=2,'g',name='$HD^{16}O$',/overplot)
    pn=plot(on[0,*],on[1,*],xthick=3,ythick=3,thick=2,'indigo',name='$HD^{18}O$',/overplot)
    pn7=plot(on7[0,*],on7[1,*],xthick=3,ythick=3,thick=2,'violet',name='$HD^{17}O$',/overplot)
    pn2=plot(on2[0,*],on2[1,*],xthick=3,ythick=3,thick=2,'dark slate gray',name='$D_2O$',/overplot)
    p5=plot(o5[0,*],o5[1,*],xthick=3,ythick=3,thick=1,'r',name='$CO_2$',/overplot)
    p6=plot(o6[0,*],o6[1,*],xthick=3,ythick=3,thick=1,'orange',name='$O_3$',/overplot)
    p7=plot(o7[0,*],o7[1,*],xthick=3,ythick=3,thick=1,'brown',name='$N_2O$',/overplot)
    p8=plot(o8[0,*],o8[1,*],xthick=3,ythick=3,thick=1,'indigo',name='CO',/overplot)
    p9=plot(o9[0,*],o9[1,*],xthick=3,ythick=3,thick=1,'m',name='$CH_4$',/overplot)
    p10=plot(o10[0,*],o10[1,*],xthick=3,ythick=3,thick=1,'y',name='HF',/overplot)
    p11=plot(o11[0,*],o11[1,*],xthick=3,ythick=3,thick=1,'gray',name='HOCl',/overplot)
    p12=plot(o12[0,*],o12[1,*],xthick=3,ythick=3,thick=1,'olive',name='$H_2S$',/overplot)

  endelse
  ;stop

  p1.xrange=[hdo[i].nu-2.5,hdo[i].nu+2.5]
  p1.yrange=[.995,1.]

  p1.xtitle='Wavenumber ($cm^{-1}$)'
  p1.ytitle='Normalized ICOS Signal'
  p1.title=string(float(hdo[i].nu-2.5),format='(F12)')+' to '+string(float(hdo[i].nu+2.5),format='(F12)')

  !null=legend(target=[p0,p1,p2,p3,p4,pn,pn7,pn2,p5,p6,p7,p8,p9,p10,p11,p12],shadow=0,linestyle=6)

  pp=[p0,p1,p2,p3,p4,pn,pn7,pn2,p5,p6,p7,p8,p9,p10,p11,p12]

  pz=[[pz],[pp]]

endfor
  
endif

return,pp

end
