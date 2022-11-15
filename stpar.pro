function stpar,mol,iso,T,E,nu,s0

;Calculates the HITRAN linestrength as a function of temperature
;Uses 'parsum.dat', which is a compilation of partition functions Q(T) for molecules in HITRAN

k=1.3806503d-23                 ;Boltzmann's Contsant
hc=1.98644568d-23               ;h*c in units of J*cm

names=''

q=fltarr(109,2931)

openr,lun1,'/media/benjamin/Data/Hitran/parsum.dat',/get_lun
readf,lun1,names
readf,lun1,q
free_lun,lun1

;print,s0

if mol eq 1 then co=iso
if mol eq 2 then co=6+iso
if mol eq 3 then co=14+iso
if mol eq 4 then co=32+iso
if mol eq 5 then co=37+iso
if mol eq 6 then co=42+iso
if mol eq 8 then co=48+iso
if mol eq 11 then co=55+iso
if mol eq 14 then co=62 
if mol eq 21 then co=77+iso
if mol eq 31 then co=92+iso

q0=q[co,226]

qt=interpol(q[co,*],q[0,*],T)

kbT=k*T
kbT0=k*296D0
hcE=hc*E

S=s0*(exp(-hcE/kbT)*q0*(1D0-exp(-hc*nu/kbT)))/(exp(-hcE/kbT0)*qt*(1D - exp(-hc*nu/kbT0)))
;print,qt
return,S

end
