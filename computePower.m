function [A_solarpanel,m_power,m_solararray] = computePower(lifetime,P_prop,m_dry_final,PHI,Td,Power_contingency)
%% Inputs and assumption
XE = 0.65; %Efficiency of power source in eclipse mode
XD = 0.85; %Efficiency of power source in Sun lit mode 
n_batteries = 2;
Energy_density_battery = 60; %Wh/kg From Battery spec sheet
rho_solararray = 8.3366; %kg/m2 From solar array spec sheet
DoD = 0.2; %Depth of charge
eta_conv = 0.247; %Efficiency of conversion
eta_inher = 0.72; %Internal efficiency
eta_point = cosd(20); %Pointing efficiency
degrad_factor = 3.75/100; %Degrad factor
battery_eff = 0.9;


%% Calculation
eta_degrad = (1-degrad_factor)^lifetime; %Efficiency at EOL
P_saa = PHI*eta_conv*eta_inher*eta_point*eta_degrad; %Power per unit area of a solar array
Pe = ((150*2103.878)+(2.17*2103.878)+(1.2*600)+(42*1210*2)+(6*841.2))*(1+Power_contingency); %Power profile imported and calculated directly
Pd = ((P_prop*1311.5)+(1.2*600)+(42*1210*2)+(5.5*9.5*60)+(6*1575.2)+(2.17*2626.85)+(150*3938.82))*(1+Power_contingency); %Power profile imported and calculated directly
P_sa = ((Pe/XE)+(Pd/XD))/Td; %Power required by solar array in W
A_solarpanel = P_sa/P_saa; %in m2 Area of solar array
P_eclip = Pe*XE; %Power need in eclipse mode in W
C = P_eclip/(DoD*battery_eff*n_batteries*3600); %Battery capacity
m_battery = C*n_batteries/Energy_density_battery; %Mass of battery in kg
m_solararray = rho_solararray * A_solarpanel; %Mass of solar array in kg
PMAD_mass = 0.045*P_sa; %PMAD mass by scaling law
Wiring_mass = 0.02*m_dry_final; %Wiring mass by empirical fit
m_power = m_battery + m_solararray + Wiring_mass+PMAD_mass; %Total mass
