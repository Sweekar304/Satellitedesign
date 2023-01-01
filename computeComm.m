function [m_comm,P_comm] = computeComm(comm_contingency, comm_margin)
%This function computes communication subsystem mass and power draw
%Based on data link transmitter and antennas are already chosen

%% High gain antenna and transmitter
m_data_trans = (59*1.05)/1000; %Innfolight SDR in kg
m_HG_antenna = 290/1000; %Syrlinks Span XT3 Antenna in kg
P_data_trans = 2.5;%W
P_HG_antenna = 3; %W

%% Low gain antenna and transreceivers
m_transreceivers = 132*1.03/1000; % in Kg ISISpace antenna and transreceivers
m_LG_antennas = 100/1000; %in kg
P_transreceivers = 132*1.03/1000; %in W
P_LG_antennas = 100/1000;% in W

m_comm = m_data_trans + m_HG_antenna + m_transreceivers + m_LG_antennas + comm_contingency; %in Kg
m_comm = m_comm*(1+comm_margin);


P_contingency = 10.1;
P_comm = P_data_trans + P_HG_antenna + P_transreceivers + P_LG_antennas + P_contingency; %in W
P_comm = P_comm*(1+comm_margin);