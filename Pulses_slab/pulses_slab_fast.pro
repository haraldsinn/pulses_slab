; -------------------------------------------------------------
; --- Calculates time evolution of a single XFEL-pulse --------
; --- and adds subsequent pulses ------------------------------
; --- materials: natural silicon, silicon28, natural diamond,
; --- and diamond12 -------------------------------------------
; --- to run program in idl/gdl: '.run pulse_slab_gdl' --------
; --- needs subroutines: read_dat.pro, colors_pro, 1plot.pro, 
; --- additional for gdl from idl: spline.pro, interpol_idl.pro
; --- needs inputs files: format ASCII: -----------------------
; --- col1: temperature, col2: cp [J/g*K] or lambda [W/cm*K]---
; -------------------------------------------------------------
; --- Harald Sinn, 27.7.2007, DESY ----------------------------
; --- fast version: Harald Sinn, 10.4.2012 E.XFEL--------------
; -------------------------------------------------------------
;
; crystal geometry:    beam direction
;                _______    |    ______
;               |       \   |   /      |
;               |        \_____/       |
;               |         _____        |
;               |        /  |  \       |
;               |_______/   |   \______|
;                          \|/  
;                           V
;
; --- adjustable parameters ------------- ---------------------
pulse_energy=2e-3 	; energy per photon pulse in J
abso=703.0	    	; absorption length in microns 
pulse_sep=0.22       	; separation between pulses in usec.
hot_pulses=10     	; number of subsequent pulses
cold_pulses=2       	; number of pulse-length at end to cool down
FWHM_beam=1000.      	; FWHM of round gaussian beam in micron
material='C_nat'    	; materials: Si_nat, Si28, C_nat, C12
T_base= 100.         	; base crystal temperature in K
xtal_rad=2000.      	; radius of crystal disk in micron
xtal_thick_c=200    	; thickness of crystal in center
xtal_rad_c=500     	; radius of thin region in the center 
xtal_thick_e=200    	; thickness of crystal at edge
;
;
; --- program parameters -------------------------------------
maxit=100000L   	; number of iterations per pulse
maxstep=1.0 		; maximum time step in ns
minstep=pulse_sep*1e3/maxit	; minimum stepsize in ns
itdivide=5e4  		; the higher the smaller the min.time step
x_points=100            ; number of points along crystal
T_min=1.                ; minimum Temp. for spline interpol.
T_max=1800              ; max. Temp. for spline interpol.
T_points=18000          ; points for spline interpolations
pulse_length=1e-4     	; pulse length in usec in output file, no impact on calc. 
pulse_number=hot_pulses+cold_pulses+1
;
;
; --- definition of input file names and variables -----------
case material of 
  'Si_nat': begin
  rho=2.33               ; density in g/cm^3
  cp_file='cp_silicon_debye.dat'
  lambda_file='lambda_si_nat.dat'
  end
  'Si28': begin 
  rho=2.33
  cp_file='cp_silicon_debye.dat'
  lambda_file='lambda_si_28.dat'
  end 
  'C_nat': begin
  rho=3.52
  cp_file='cp_diamond_debye.dat'
  lambda_file='lambda_C_nat.dat'
  end
  'C12': begin
  rho=3.52 
  cp_file='cp_diamond_debye.dat'
  lambda_file='lambda_C_12.dat'
  end
  else: print, 'material not recognized'
endcase
;
;
; --- read in data files --------------------------------------
read_dat, cp_file, data
cp_T=data(*,0)
cp_raw=data(*,1)
read_dat, lambda_file, data
lambda_T=data(*,0)
lambda_raw=data(*,1)
;
;
; --- interpolate data -------------------- -------------------
T_axis=findgen(T_points)/(T_points-1)*(T_max-T_min)+T_min 
cp=spline(cp_T,cp_raw,T_axis,0)
lambda=spline(lambda_T,lambda_raw,T_axis,0)
;
;
; --- calculate initial temperature profile -----------------
; -- integrate cp over T
cp_int=T_axis
sum=0
for i=1, T_points-1 do begin
  sum=sum+cp[i]
  cp_int[i]=sum
endfor
cp_int=cp_int*(T_axis[1]-T_axis[0])
; -- find index of T_min
i_t_min=0
while T_base-T_axis[i_t_min] gt 1e-2 do i_t_min=i_t_min+1
cp0=cp_int[i_T_min]
; -- generate x-axis,thickness profile, mass of each crystal ring
x_axis=findgen(x_points)/x_points*xtal_rad
dx=x_axis[1]-x_axis[0]
x_thick=x_axis
for i=0,x_points-1 do begin 
  if x_axis[i]- xtal_rad_c lt 1e-3 then begin 
    x_thick[i]=xtal_thick_c
  endif else begin 
    x_thick[i]=(x_axis[i]-xtal_rad_c)*(xtal_thick_e-xtal_thick_c) $
     /(xtal_rad-xtal_rad_c)+xtal_thick_c
; --quadratic variation of thickness profile gives smoother gradient 
;    x_thick[i]=(x_axis[i]-xtal_rad_c)^2*(xtal_thick_e-xtal_thick_c) $
;     /(xtal_rad-xtal_rad_c)^2+xtal_thick_c
  end 
endfor
x_area=x_axis  ; area facing the X-rays
x_area[0]=!pi*dx^2
for i=1, x_points-1 do begin 
  x_area[i]=2*!pi*(x_axis[i]+dx/2)*dx
end
x_volume=x_axis
x_volume[0]=x_area[0]*x_thick[0]
for i=1, x_points-2 do begin 
  x_volume[i]=x_area[i]*(x_thick[i]+x_thick[i+1])/2
end
i=x_points-1
x_volume[i]=x_area[i]*(3*x_thick[i]-x_thick[i-1])/2 
x_mass=x_volume*rho*1e-12 ; mass of each ring in g
; -- gaussian energy distribution of beam
sigma=FWHM_beam/2.35482
gauss_beam=1/(2*!pi*sigma^2)$
        *exp(-(x_axis+dx/2)^2/(2*sigma^2)>(-50))*pulse_energy $
        *x_area
; -- absorbed energy in each ring
gauss_abs=gauss_beam*(1-exp(-x_thick/abso))
q_per_pulse=total(gauss_abs)
q_pulse=gauss_abs/x_mass
T_pulse=interpol_idl(T_axis,cp_int-cp0,q_pulse)
T_final_pulse=$
 interpol_idl(T_axis,cp_int-cp0,q_per_pulse/total(x_mass)*hot_pulses)
;
;
j_left=x_axis
j_right=x_axis
dT_dt=x_axis
train_time=findgen(pulse_number)*pulse_sep
T_max=train_time
T_min=train_time
T_edge=train_time
T_max[0]=max(T_pulse)
T_min[0]=T_base
T_edge[0]=T_base
p_hwhm=fix(FWHM_beam/xtal_rad*(x_points-1)/2.0)
;
;
; --- start pulse train ------------------------------------------
;
;
@1plot
aa=systime(1)
for i=1,pulse_number-1 do begin 
  ; --- calculate thermal diffusion equation for T-distr.----------
  plot, x_axis, T_pulse
  oplot, x_axis, T_pulse*0+T_final_pulse, line=2, color=green
  print, 'pulse number:', i 
  print, ' step,    time step [10-8],  total time [10-6s],'$
                ,'    T_max,     Q[pulses],    Q_peak/Q_total'
  time=0L
  j=1L
  kappa=lambda/cp
  ; --- start individual pulse -----------------------------------
  while (j lt maxit) and time lt pulse_sep do begin 
    kappa_pulse=interpol(kappa,T_axis,T_pulse)  ; 
    ;
    ; --- calculate heat flows left and right of each element
    T_pulse_m1=shift(T_pulse,1)
    j_left=-(T_pulse-T_pulse_m1)/dx
    j_left=j_left*2*!pi*x_axis*x_thick
    j_right=shift(j_left,-1)
    j_left[0]=0
    j_right[x_points-1]=0
    dT_dt=(j_left-j_right)*kappa_pulse/(rho*x_volume) 
    dt=1./max(dT_dt/(T_pulse-T_base+1.))/itdivide<maxstep/10>minstep/10
    time=time+dt/100.  ; time in microseconds
    T_pulse=T_pulse+dT_dt*dt
    if j mod 500 eq 0 or j eq 1 then begin 
      if j ne 1 then  oplot,x_axis, T_pulse, color=blue
      print,j,dt,time,max(T_pulse)
    endif 
    j=j+1
  endwhile 
  ; calculate total energy in crystal 
  cp_int_pulse=interpol(cp_int,T_axis,T_pulse)-cp0
  q_total_j=total(cp_int_pulse*x_mass)    
  ; calculate energy in FHHM area
  q_fwhm_j=total((cp_int_pulse*x_mass)[0:p_hwhm])
  print,j,dt,time,max(T_pulse),q_total_j/q_per_pulse,q_fwhm_j/q_total_j
  oplot,x_axis, T_pulse, color=red
  print, 'sim. T_final:', T_pulse[0]
  T_min[i]=max(T_pulse)
  T_edge[i]=T_pulse[x_points-1]
  ; --- renormalize pulse energy ------------------------------
  if i lt hot_pulses then begin 
    cp_int_pulse=interpol(cp_int,T_axis,T_pulse)-cp0
    q_pulse_i=cp_int_pulse*i/q_total_j*q_per_pulse ; normalize energy
    q_pulse_i=q_pulse_i+q_pulse  ; add new pulse
    T_pulse=interpol_idl(T_axis,cp_int-cp0,q_pulse_i)
  endif else begin 
    cp_int_pulse=interpol(cp_int,T_axis,T_pulse)-cp0
    q_pulse_i=cp_int_pulse*hot_pulses/q_total_j*q_per_pulse ; normalize energy
    T_pulse=interpol_idl(T_axis,cp_int-cp0,q_pulse_i)
  endelse
  T_max[i]=max(T_pulse)
  wait, 0.1 ; -- needed for graphics buffer to follow up 
endfor 
print, 'total time:', systime(1)-aa, 'sec'
;
;
; --- plot results for complete pulse train ----------------------
plot, train_time, T_edge, yrange=[T_base-10, max(T_max*1.1)]$
      , title='blue=T_max, black=T_edge, green=T_inf'$
     , xtitle='microsec', ytitle='Kelvin'
train_time_long=[train_time,train_time+pulse_length]
T_curve=[T_min,T_max]
T_curve=T_curve(sort(train_time_long))
train_time_long=train_time_long(sort(train_time_long))
oplot, train_time_long, T_curve, color=blue
oplot, train_time, T_min*0+T_final_pulse, color=green
print, 'integrated beam power  :', total(gauss_beam)*1000, ' mJ'
print, 'absorbed in one shot   :', q_per_pulse*1000, ' mJ = ',$
                q_per_pulse/total(gauss_beam)*100,'%'
print, 'absorbed in pulse train:', q_per_pulse*hot_pulses, ' J'
print, 'average temperature after pulse train(calc):', T_final_pulse
openw,1, 'pulses_slab.dat'
  for i=0, pulse_number-1 do begin 
    printf, 1, train_time(i), T_min(i)
    printf, 1, train_time(i)+pulse_length, T_max(i)  
  endfor 
close, 1
end
