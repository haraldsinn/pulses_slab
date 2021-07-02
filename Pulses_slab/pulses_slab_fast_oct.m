% -------------------------------------------------------------
% --- Calculates time evolution of a single XFEL-pulse --------
% --- and adds subsequent pulses ------------------------------
% --- materials: natural silicon, silicon28, natural diamond,
% --- and diamond12 -------------------------------------------
% --- needs inputs files: format ASCII: -----------------------
% --- col1: temperature, col2: cp [J/g*K] or lambda [W/cm*K]---
% -------------------------------------------------------------
% --- Harald Sinn, 27.7.2007, DESY ----------------------------
% --- Harald Sinn, 10.4.2012, EXFEL: converted to Octave -----
% -------------------------------------------------------------
%
% crystal geometry:    beam direction
%                _______    |    ______
%               |       \   |   /      |
%               |        \_____/       |
%               |         _____        |
%               |        /  |  \       |
%               |_______/   |   \______|
%                          \|/  
%                           V
%
% --- adjustable parameters ------------- ---------------------
% function pulse_slab()
pulse_energy=2e-3;   	% Energy per photon pulse in J
abso=703.0;             % absorption length in microns
pulse_sep=0.22;       	% separation between pulses in usec.
hot_pulses=10;     	% number of subsequent pulses+1
cold_pulses=2;       	% number of pulse-length at end to cool down
FWHM_beam=1000.0;      	% FWHM of round gaussian beam in micron
material='C_nat';    	% materials: Si_nat, Si28, C_nat, C12
T_base= 100.0;         	% base crystal temperature in K
xtal_rad=2000.0;      	% radius of crystal disk in micron
xtal_thick_c=200;    	% thickness of crystal in center
xtal_rad_c=500;     	% radius of thin region in the center 
xtal_thick_e=200;    	% thickness of crystal at edge
%
%
% --- program parameters -------------------------------------
maxit=100000;           % number of iterations per pulse
maxstep=1.0;            % maximum time step in ns
minstep=pulse_sep*1e3/maxit;  	% minimum time step in ns
itdivide=5e4;   	% the higher the smaller the dynamical time step
x_points=100;     	% number of points along crystal
T_min=1;     		% minimum Temp. for spline interpol.
T_max=1800;            	% max. Temp. for spline interpol.
T_points=18000;         % points for spline interpotations 50000
pulse_length=1e-4;	% pulse length (us) in output file, no impact on calc.
pulse_number=hot_pulses+cold_pulses+1;
%
%
% --- definition of input file names and variables -----------
switch material  
  case {'Si_nat'}
  rho=2.33;               % density in g/cm^3
  cp_file='cp_silicon_debye.dat';
  lambda_file='lambda_si_nat.dat';
  case {'Si28'}
  rho=2.33;
  cp_file='cp_silicon_debye.dat';
  lambda_file='lambda_si_28.dat';
  case{'C_nat'} 
  rho=3.52;
  load 'cp_diamond_debye.dat';
  cp_dat=cp_diamond_debye;
  load 'lambda_C_nat.dat';
  lambda_dat=lambda_C_nat;
  case {'C12'} 
  rho=3.52; 
  cp_file='cp_diamond_debye.dat';
  lambda_file='lambda_C_12.dat';
  otherwise 
    print('material not recognized');
endswitch
%
%
% --- read in data files --------------------------------------
points=size(cp_dat)(1);
cp_T=cp_dat(1:points,1);
cp_raw=cp_dat(1:points,2);
points=size(lambda_dat)(1);
lambda_T=lambda_dat(1:points,1);
lambda_raw=lambda_dat(1:points,2);
%
%
% --- interpolate data -------------------- -------------------
T_axis=linspace(T_min, T_max, T_points);
cp=interp1(cp_T,cp_raw,T_axis,'spline');
lambda=interp1(lambda_T,lambda_raw,T_axis,'spline');
%
%
% --- calculate initial temperature profile -----------------
% -- integrate cp over T
cp_int=T_axis;
sum_cp=0;
for i=1:T_points
  sum_cp=sum_cp+cp(i);
  cp_int(i)=sum_cp;
endfor
cp_int=cp_int*(T_axis(2)-T_axis(1));
% -- find index of T_min
i_t_min=1;
while (T_base-T_axis(i_t_min)>1e-2)
  i_t_min=i_t_min+1;
endwhile
cp0=cp_int(i_t_min);
% -- generate x-axis,thickness profile, mass of each crystal ring
x_axis=linspace(0.0,x_points-1.0,x_points)/x_points*xtal_rad;
dx=x_axis(2)-x_axis(1);
x_thick=x_axis;
for i=1:x_points 
  if (x_axis(i)-xtal_rad_c<1e-3) 
    x_thick(i)=xtal_thick_c;
  else 
    x_thick(i)=(x_axis(i)-xtal_rad_c)*(xtal_thick_e-xtal_thick_c) ...
       /(xtal_rad-xtal_rad_c)+xtal_thick_c;
  endif 
endfor
x_area=x_axis;  		% area facing the X-ray
x_area(1)=pi*dx^2;
for i=2:x_points 
  x_area(i)=2*pi*(x_axis(i)+dx/2)*dx;
endfor
x_volume=x_axis;
x_volume(1)=x_area(1)*x_thick(1);
for i=1:x_points-1
  x_volume(i)=x_area(i)*(x_thick(i)+x_thick(i+1))/2;
endfor
i=x_points;
x_volume(i)=x_area(i)*(3*x_thick(i)-x_thick(i-1))/2 ;
x_mass=x_volume*rho*1e-12; % mass of each ring in g
% -- gaussian energy distribution of beam
sigma=FWHM_beam/2.35482;
gauss_beam=1./(2*pi*sigma^2)*exp(max(-50,-(x_axis.+dx/2).^2./ ... 
   (2*sigma^2)))*pulse_energy.*x_area;
% -- absorbed energy in each ring
gauss_abs=gauss_beam.*(1.-exp(-x_thick/abso));
q_per_pulse=sum(gauss_abs);
q_pulse=gauss_abs./x_mass;
T_pulse=interp1(cp_int-cp0,T_axis, q_pulse);
T_final_pulse=interp1(cp_int-cp0,T_axis, q_per_pulse/sum(x_mass)*hot_pulses);
%
%
j_left=x_axis;
j_right=x_axis;
dT_dt=x_axis;
train_time=linspace(0, pulse_number-1, pulse_number)*pulse_sep;
T_max=train_time;
T_min=train_time;
T_edge=train_time;
T_max(1)=max(T_pulse);
T_min(1)=T_base;
T_edge(1)=T_base;
p_hwhm=fix(FWHM_beam/xtal_rad*(x_points-1)/2.0);
%
%
% --- start pulse train ------------------------------------------
%
%
tic()
for i=1:pulse_number-1 
  % --- calculate thermal diffusion equation for T-distr.----------
  plot(x_axis, T_pulse)
  line(x_axis, T_pulse*0+T_final_pulse,'color','g')
  pause(0.01) % -- needed for graphics buffer to follow up 
  printf('pulse number: %g \n', i) 
  printf('  step, time step [10-8], total time [10-6s], T_max,  Q[pulses], Q_peak/Q_total \n')
  time=0.0;
  j=1; %  j=int32(1);
  kappa=lambda./cp;
  kappa_pulse=interp1(T_axis,kappa,T_pulse,'*nearest')./(rho.*x_volume);  % linear 2.5 x slower
  jfactor=2*pi.*x_axis.*x_thick;
  % --- start individual pulse -----------------------------------
  while ((j< maxit)&(time< pulse_sep)) 
    kappa_pulse=interp1(T_axis,kappa,T_pulse,'*nearest')./(rho.*x_volume);  % linear 2.5 x slower
    %
    % --- calculate heat flows left and right of each element
    T_pulse_m1=shift(T_pulse,1); 
    j_left=-(T_pulse-T_pulse_m1)/dx.*jfactor;
    j_right=shift(j_left,-1);  
    j_left(1)=0;
    j_right(x_points)=0;
    dT_dt=(j_left.-j_right).*kappa_pulse ;
    dt=1.0/max(dT_dt./(T_pulse-T_base+1.))/itdivide;
    dt=min(dt,maxstep/10.0);  
    dt=max(dt,minstep/10.0);
    time=time+dt/100. ; % time in microseconds
    T_pulse=T_pulse.+dT_dt*dt;
    if (mod(j,500)==0)|(j==1) 
      if j != 1  
        line(x_axis, T_pulse,'color','b')
        pause(0.01) % -- needed for graphics buffer to follow up 
       endif
      fprintf(1,'%8d   %5f       %10f       %5f\n', j,dt,time,max(T_pulse))
      fflush(1);
    endif 
    j=j+1;
  endwhile 
  cp_int_pulse=interp1(T_axis,cp_int,T_pulse,'*linear')-cp0;  %linear gives much better quality   
  q_total_j=sum(cp_int_pulse.*x_mass);
  % calculate energy in FHHM area
  q_fwhm_j=sum((cp_int_pulse.*x_mass)(1:p_hwhm+1));
  fprintf(1,'%8d   %5f       %10f       %5f    %5f    %5f\n', j,dt,time,max(T_pulse),q_total_j/q_per_pulse,q_fwhm_j/q_total_j)
  line(x_axis,T_pulse,'color','r')
  fprintf(1,'sim. T_final: %12.3f \n', T_pulse(1))
  pause(0.01) % -- needed for graphics buffer to follow up 
  T_min(i+1)=max(T_pulse);
  T_edge(i+1)=T_pulse(x_points);
  % --- renormalize pulse energy ------------------------------
  if i < hot_pulses 
    cp_int_pulse=interp1(T_axis,cp_int,T_pulse,'*spline')-cp0;
    q_pulse_i=cp_int_pulse*i./q_total_j*q_per_pulse; % normalize energy
    q_pulse_i=q_pulse_i+q_pulse;  % add new pulse
    T_pulse=interp1(cp_int-cp0,T_axis, q_pulse_i,'spline');
  else 
    cp_int_pulse=interp1(T_axis,cp_int,T_pulse,'*spline')-cp0;
    q_pulse_i=cp_int_pulse*hot_pulses./q_total_j*q_per_pulse; % normalize energy
    T_pulse=interp1(cp_int-cp0,T_axis, q_pulse_i,'spline');
  endif
  T_max(i+1)=max(T_pulse);
  pause(0.2) % -- needed for graphics buffer to follow up 
endfor
toc()
%
%
% --- plot results for complete pulse train --------------------
train_time_long=[train_time,train_time+pulse_length];
T_curve=[T_min,T_max];
[train_time_long,index]=sort(train_time_long);
T_curve=T_curve(index);
plot(train_time_long, T_curve, 'color', 'b')
line(train_time,T_edge)
line(train_time,T_min*0+T_final_pulse, 'color', 'g')
title('blue=T_{max}, black=T_{edge}, green=T_{inf}')
xlabel('microseconds')
ylabel('Kelvin')
axis([0, max(train_time_long),T_base-10.0, max(T_max*1.1)])
fprintf(1,'integrated beam power  : %8f mJ\n', sum(gauss_beam)*1000)
fprintf(1,'absorbed in one shot   : %8f mJ =  %8f percent\n', q_per_pulse*1000, q_per_pulse/sum(gauss_beam)*100)
fprintf(1,'absorbed in pulse train: %8f J\n', q_per_pulse*hot_pulses)
fprintf(1,'final temperature after pulse train(calc): %8f \n', T_final_pulse)
filename='pulses_slab_mtl.dat';
fid=fopen(filename,'w');
for i=1:pulse_number 
  fprintf(fid, '%8.5f  %8.5f \n', train_time(i), T_min(i));
  fprintf(fid, '%8.5f  %8.5f \n', train_time(i)+pulse_length, T_max(i));
endfor
fclose(fid);
% endfunction
%

