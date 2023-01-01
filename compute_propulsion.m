function [m_prop,P_prop] = compute_propulsion(m_propellant,R,MW_Xe,MW_air)
%This function computes and returns mass and power of propulsion subsystem
%% Inputs
tank_strength = 950;%MPa Ultimate yield strength
burst_pressure = 30;%Mpa Burst pressure
density_tank = 2710; %kg/m3
Xe_pressure = 15;%MPa xenon storage pressure
init_pressurant_press = 27.6; %MPa Initial pressurant pressure
final_pressurant_press = Xe_pressure;
init_pressurant_temp = 323;%K
final_pressurant_temp = 293;%K
Thruster_mass = 0.98; %Kg from spec sheet
tank_contingency = 0.2; %Contingency for change of tank material or storage condition
m_power_manag_unit = 1.1; %kg mass of power management unit
m_prop_manag_unit = 5; %kg mass of propulsion management unit

%% Computation
V_tank = (m_propellant*R*293*1.05)/(15*10^6 * MW_Xe); %Computes volume required to store Xenon in m3 with 5perc expansion
r_tank = ((3*V_tank)/(4*pi))^(1/3); %Computes radius of tank assuming sphere
t_tank = (burst_pressure*r_tank)/(2*tank_strength); %Computes thickness by equating stresses and pressure
m_tank_1 = density_tank*4*pi*(r_tank^2)*t_tank*(1.2); %Mass of tank considering hollowsphere
m_tank_2 = (2.7086*10^(-8)*V_tank^3)-(6.1703*10^(-5)*V_tank^2) + (6.629*10^(-2)*V_tank) + 1.3192; %From empirical fit
m_tank = max(m_tank_1, m_tank_2);
V_pressur = (final_pressurant_press*V_tank)/((final_pressurant_temp*init_pressurant_press/init_pressurant_temp)-15); %Volume of pressurant requried
m_pressur = (init_pressurant_press*10^6)*V_pressur*MW_air/(R*init_pressurant_temp); %Mass of pressurant required
PV_calc = init_pressurant_press*V_pressur;
m_press_tank = 0.7266*(PV_calc)^2+2.5119*PV_calc +2.9826 ; %Empirical fit
m_prop = Thruster_mass + (m_tank*(1+tank_contingency)) + (m_press_tank*(1+tank_contingency)) + m_pressur +m_power_manag_unit  + m_prop_manag_unit;
P_prop = 0.2*10^3; %W from spec sheet
disp('Busek BHT200 Hall thruster is chosen')

end