function dsdtpar,mol,iso,T,E,nu,s0

;Calculates the HITRAN linestrength as a function of temperature
;Uses 'parsum.dat', which is a compilation of partition functions Q(T) for molecules in HITRAN

k=1.3806503d-23                 ;Boltzmann's Contsant
hc=1.98644568d-23               ;h*c in units of J*cm

names=''

q=fltarr(109,2931)

openr,lun1,'~/Documents/HITRAN08/parsum.dat',/get_lun
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
if mol eq 14 then co=62 
if mol eq 21 then co=77+iso
if mol eq 31 then co=92+iso

q0=q[co,226]

tstar=dindgen(11)+T-5

i1=where(q[0,*] ge T-5)
i2=where(q[0,*] ge T+5)

qt=interpol(q[co,i1[0]:i2[0]],q[0,i1[0]:i2[0]],tstar)
;p=plot(q[0,*],q[co,*])
ft=poly_fit(tstar,qt,2)
;print,ft
dqtdt=2d0*ft[2]*tstar[5]+ft[1]

kbT=k*T
kbT0=k*296D0
hcE=hc*E
print,dqtdt,qt[5]
;dSdt=s0*q0*hc/(k*T^2)*exp(-hcE/kbT)*(E-(E+nu)*exp(-hc*nu/kbT))/(exp(-hcE/kbT0)*qt*(1D - exp(-hc*nu/kbT0)))

S=s0*(exp(-hcE/kbT)*q0*(1D - exp(-hc*nu/kbT)))/(exp(-hcE/kbT0)*qt[5]*(1D - exp(-hc*nu/kbT0)))

dSdT=s0*q0*exp(-hcE/kbT)/(qt[5]*exp(-hcE/kbT0)*(1D0-exp(-hc*nu/kbT0)))*((hc/(k*T^2))*(E-(E+nu)*exp(-hc*nu/kbT))-dqtdt/qt[5]*(1d0-exp(-hc*nu/kbT)))

frac=dsdT/S
print,dSdt,S,frac

return,frac
end
