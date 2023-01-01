%% Instructions
%This code is written by Sweekar Bengaluru. Please reach out at
%ssweekar@stanford.edu for clarifications.
%This code is intended to help rapidly size the satellite with minimum
%inputs. Mission requirements, design variables, initial masses needs to be
%specified.
%This code implements following logic, initially dry mass is computed based
%on initial masses specified, with the initial masses and design values,
%delta V is computed. Delta V is also dependent on solarpanel area and dry
%mass.With this Delta V, mass of propellant can be calculated which serves
%as an input to Propulsion subsystem sizing. To compute ADCS system,
%disturance torques needs to be known. To compute disturbance toruqes,
%inertias needs to be known for which masses of elements needs to be known. 
%But the size of the solar panels are
%not known and neither its mass. Hence, with initial assumption these
%values are computed, initial power requirements is obtained and used as an
%input for the power subsystem computation, masses and areas of solar panel
%from power subsystem is used to iterate until there is a convergence
%This thus forms the inner loop for convergence. However, with this change
%and with new masses computed for ADCS, communication, Power, propulsion,
%dry mass would be different from initial value which would affect delta V,
%mass of propellant and hence all the other subsystems. Hence, the outer
%loop is set to iterate on dry mass. During these iterations, computed
%values can go outside the specifications of vendor components especially
%in ADCS and Propulsion subsystem category. Hence, quality checks are
%implemented to ascertain performance values are within the specifications.

%NOTE: Several inputs are set in each individual functions, please check to
%confirm if acceptable, without change inputs runs as default. Used to
%prevent user from inputting large values at the start since lot of these
%values are unknown.
%Default sensor and actuator data are used directly in the ADCS
%Default Transmitter, antenna and receiver used directly in Comm
%Default Power profile imported directly in Compute power secction

%% Cleanup
clear all
close all
clc

%% Global constants
MU = 3.986e5; %km3/s
RP = 6370; %km
rho_800 = 2.27814498 * 10^(-9); %Compute a model
g = 9.81; %m/s2 acceleration due to gravity
R = 8.314; %J/mol-K Universal gas constant
MW_air = 28*10^(-3); %kg/mol Molecular weight of air
Mag_Earth = 7.8*10^(15); %Tm3 Magnetic flux of earth
PHI = 1368.4; %W/m2 solar flux at 1AU

%% Mission requirements
h = 800; %Altitude of operation in Km
h_f = 120; %Final altitude to which satellite needs to be deorbited in Km
lifetime = 5; %Mission duration in years
RANN = 20.688; %Deg for 1030LTAN
slew_rate = 0.059; %Slewrate deg/s
slew_rate_deg = 30; % Slew rate in degs
slew_rate_time = 8.5 ; %Slew duration in mins

%% Design variables This specifies basic geometry and propellant choice
Isp = 1390; %in s. Change based on propellant choice, value entered for Xenon
MW_Xe = 0.132;%in kg/mol Molecular weight of Xenon
bus_l = 2.4; %in m Satellite bus length
bus_w = 2.4; %in m Satellite bus width
bus_h = 2.7; %in m Satellite bus height
Cd = 2.2; %Coefficient of discharge of satellite body
form_factor_sp = 9/7; %Form factor of solar panel dimension ratio of length to height

%% Margin
delV_margin = 0.25; %Delta V margin for design
comm_contingency = 200/1000; % in Kg Communication subsystem contingency
comm_margin = 0.25; %Communication subsystem margin
Power_contingency = 0.2035; %Power subsystem contingency on power mode calc

%% Initial mass and geometry estimation
mpay = 243.5; %kg
m_prop_init = 28.445; %kg
m_ADCS_init = 55.7679; %kg
m_TTC_init = 21.422; %kg
m_OBC_init = 50.06;%kg
m_power_init = 219.36;%kg
m_struct_init = 262.749;%kg
m_therm_init = 28.73;%kg
m_misc_init = 98.12;%kg
m_dry_init = mpay + m_prop_init + m_ADCS_init + m_TTC_init + m_OBC_init...
              + m_power_init + m_struct_init + m_therm_init + m_misc_init;%kg
sp_h = 1.8; %in m Initial solar panel height
sp_l = 1.4; %in m Initial solar panel length
A_solarpanel = sp_h*sp_l*2; %in m2 Initial solar panel area
A_bus = bus_l*bus_h;   %in m2 Satellite bus area 
m_sp_init = 126.0498; %in kg Initial solar panel mass
A_solarpanel_old = 100; %m random initialisation for starting inner while loop

%% Orbit elements and orbit events
a = RP + h; %km Semimajor axis
e = 0; %Eccentricity =0 means a circular orbit
i = acosd(-360/(365.25*2.06247*10^(14)*a^(-7/2))); %Inclination angle to obtain sun synch in degrees
Omega = 0; %Arg of latitude in degs
Gamma = 0; %True anomaly
T = 2*pi*sqrt(a^3/MU); %Period in s
n = 2*pi/T; %Mean motion in 1/s
lambda = asind(RP/a); %To calculate eclipse half angle in deg
T_eclipse = (2*lambda/360)*T; %Eclipse duration in s
V = sqrt(MU/a); %Velocity in km/s
n_orbit = floor(24*3600/T); %Number of orbits in a day
T_day = T - T_eclipse; %ins
%% Convergence criterion
tol_mass = 10; %in kg tolerance for convergence of dry mass
tol_area = 0.5; % in m2 tolerance for convergence of solar array area
count = 0; %Counter for debugging while loop
i = 0; %Counter for deugging inner while loop

%% DeltaV Budget
delV_total = compute_deltaV(V,h,h_f,RP,A_solarpanel,A_bus,m_dry_init,Cd,rho_800,MU,a,n_orbit,lifetime,delV_margin); %in km/s To compute initial delta V
m_propellant = (exp((delV_total*1000)/(Isp*g))*m_dry_init)-m_dry_init; %To compute initial propellant mass in kg
m_prop = compute_propulsion(m_propellant,R,MW_Xe,MW_air);     %To compute initial propulsion subsystem mass
m_dry_final = mpay + m_prop + m_ADCS_init + m_TTC_init + m_OBC_init...
              + m_power_init + m_struct_init + m_therm_init + m_misc_init; %To compute initial dry mass, only for intialising while loop
T_init = m_dry_final*delV_total*1000/(.217*3600*24*365*5); %Initial thrust level for deugging and verification
m_sp = m_sp_init; 
while abs(m_dry_final - m_dry_init)>tol_mass  %Outer while loop to set convergence on Dry mass
    count = count+1;
    while abs(A_solarpanel - A_solarpanel_old)>tol_area || abs(m_sp_old - m_sp)>tol_mass %Inner while loop to set convergence on Solar panel designs
        i = i+1;
        delV_total = compute_deltaV(V,h,h_f,RP,A_solarpanel,A_bus,m_dry_final,Cd,rho_800,MU,a,n_orbit,lifetime,delV_margin); %in km/s computes delta V
        m_propellant = (exp((delV_total*1000)/(Isp*g))*m_dry_final)-m_dry_final; %in Kg computes mass of propellant
        [m_prop,P_prop] = compute_propulsion(m_propellant,R,MW_Xe,MW_air);     %in Kg returns mass and power needs of propulsion subsystem
        sp_h = sqrt(form_factor_sp*A_solarpanel); %Computes new solar panel dimensions
        sp_l = sqrt((1/form_factor_sp)*A_solarpanel); %Computes new solar panel dimensions
        [I,Ix,Iy,Iz] = computeInertia(m_sp,m_dry_final,bus_l,bus_w,bus_h, sp_h, sp_l); %Computes inertia, returns maximum inertia and inertia on all 3 axis
        Td = ((3*MU)/(2*a^3))*(Iz - Iy)*sind(2);%Gravity gradient Disturbance torque in Nm
        [m_adcs,P_adcs] = computeADCS(Td,I,Mag_Earth,a,T,slew_rate_deg,slew_rate_time); %Returns ADCS mass and power requirements
        A_solarpanel_old = A_solarpanel; %Resets solarpanel area for next iteration
        m_sp_old = m_sp; %Resets solarpanel mass for next iteration
        [A_solarpanel,m_power,m_solararray] = computePower(lifetime,P_prop,m_dry_final,PHI,T_day,Power_contingency); %Computes new solarpanel area, mass of solar panel and mass of power subsystem
        m_sp = m_solararray; %Assigns mass of solarpanel computed
    end
    delV_total = compute_deltaV(V,h,h_f,RP,A_solarpanel,A_bus,m_dry_final,Cd,rho_800,MU,a,n_orbit,lifetime,delV_margin); %Computes Delta V
    m_propellant = (exp((delV_total*1000)/(Isp*g))*m_dry_final)-m_dry_final; %Computes mass of propellant
    [m_prop,P_prop] = compute_propulsion(m_propellant,R,MW_Xe,MW_air);  %Returns propulsion subsystem mass and power   
    m_dry_init = m_dry_final; %Prepares dry mass for next iteration
    sp_h = sqrt((form_factor_sp)*A_solarpanel); %New dimension of solar panel
    sp_l = sqrt((1/form_factor_sp)*A_solarpanel);
    [I,Ix,Iy,Iz] = computeInertia(m_sp,m_dry_final,bus_l,bus_w,bus_h, sp_h, sp_l); %Returns inertias and max inertia
    Td = ((3*MU)/(2*a^3))*(Iz - Iy)*sind(2); %Computes disturbance torque
    [m_adcs,P_adcs] = computeADCS(Td,I,Mag_Earth,a,T,slew_rate_deg,slew_rate_time); %Returns mass and power requirement of ADCS subsystem
    [m_comm,P_comm] = computeComm(comm_contingency,comm_margin); %Returns mass and power draw of communication subsystem. 
    [A_solarpanel,m_power] = computePower(lifetime,P_prop,m_dry_final,PHI,T_day,Power_contingency); %Returns mass and power draw of Power subsystem
    m_struct = (exp((1.0128*log(m_dry_final))-1.4332))*1.2; %Empirical fit for structural mass computation
    m_dry_final = (mpay) + m_prop + m_adcs + m_comm + m_OBC_init*1.2...
              + m_power + m_struct + m_therm_init*1.2 + m_misc_init; %Dry mass computation in kg
    P_adcs = P_adcs*1.2; %adding margin on ADCS power draw
    T_final = (m_dry_final*delV_total*1000)/(.217*3600*24*365*5); %Consistency check on chosen thruster
end

if T_final>20*10^(-3) || T_final < (m_dry_final*delV_total*1000)/(365*5*24*3600*0.25) %Thrust has to be greater than minimum thrust to meet perf req and must be less than max thrust allowed bu busek thruster
    disp('Thruster needs to be changed')
end
