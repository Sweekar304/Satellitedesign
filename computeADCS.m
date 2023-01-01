function [m_adcs,P_adcs] = computeADCS(Td,I,Mag_Earth,a,T,slew_rate_deg,slew_rate_time)
%This function takes in inertias and disturance toruqe and validtes if the
%actuators still hold good. If they do it returns value.
%% Sensors
magnetometer_mass = 0.1; %Clydespace MAG3 in kg
sunsensor_mass = 0.3; %Bradford sunsensor in kg
startracker_mass = 0.158; %RocketlabUSA in kg
Memsgyros = 0.015; %Vectornav in kg
magnetometer_power = 0.45; %Clydespace MAG3 in W
sunsensor_power = 0.0; %Bradford sunsensor in W
startracker_power = 1; %RocketlabUSA in W
Memsgyros_power = 0.22;%Vectornav in W

%% Actuators
RW_mass = 0.6; %Kg RocketLaUSA reaction wheel mass
RW_power = 10; %W RocketlabUSA reaction wheel power
MT_mass = 60*10^(-3); %kg Newspace magtorq
MT_Power = 1.2; %W Newspace magtorq

%% Computation
Trw_1 = Td*1.2; %Required torque
Trw_2 = 4*pi*slew_rate_deg*I/(180*(60*slew_rate_time)^2); %Torque to oovercome slew
Trw = max(Trw_1,Trw_2);
Hrw = (Td*T)/(4*sqrt(2)); %Angular momentum capcity required in Nms
if (Trw<100*10^(-3)) && (Hrw < 200*10^(-3)) %Validates use of actuator
    m_rw = RW_mass*4; 
    P_rw = RW_power*4;
else
    disp('Check reaction wheel selected')
end
B = (Mag_Earth*2)/((a*1000)^3); %Computes magnetic field in T
D_mt = Td/B; %Computes dipole moment in Am2
if D_mt<100
    m_mt = MT_mass;
    P_mt = MT_Power;
else
    disp('Check magentic torquer used')
end    
   
%% 
m_adcs = magnetometer_mass + sunsensor_mass + startracker_mass + Memsgyros + m_mt + m_rw;
P_adcs = P_mt + P_rw*1.05 + magnetometer_power + startracker_power*1.5 + Memsgyros_power + sunsensor_power;

disp('Clydespace Mag3 magnetometer is selected')
disp('Bradford Sunsensor is selected')
disp('RocketlabUSA startracker is selected')
disp('VectorNav mems is selected')
disp('RocketlabUSA reaction wheel selected')
disp('Newspace NMTRX mag torquers selected')

end