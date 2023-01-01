function [delV_total] = compute_deltaV(V,h,h_f,RP,A_solarpanel,A_bus,m_dry,Cd,rho_800,MU,a,n_orbit,lifetime,delV_margin)
%This function is used to calculate delta V for the given subsystem. It
%takes velocity, current altitude, deorbit altitude, solarpanel area, dry
%mass to compute delta V
delV_deorbit = V*(h - h_f)/(4*(RP+h_f)); %Delta V required for deorbiting in km/s
A_drag = A_solarpanel + A_bus; %To compute drag force area is required in m2
perturb_drag = -(2*pi*Cd*A_drag*(a^2)*rho_800/m_dry)*1000; %Drag force
V1a = sqrt(MU/(a+perturb_drag)); %Hohmann transfer calculations
a2 = a+perturb_drag;
V1b = sqrt(MU*((2/a2)-(2/(a+a2))));
V2b = sqrt(MU/a);
V2a = sqrt(MU*((2/a)-(2/(a+a2))));
delV_maintain = ((V2b-V2a)+(V1b-V1a))*n_orbit*365*lifetime/1000; %km/s Drag causes change in semimajor axis which needs correction
delV_total = (delV_maintain + delV_deorbit)*(1+delV_margin); %Total delta V