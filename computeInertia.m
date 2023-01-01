function [I,Ix,Iy,Iz] = computeInertia(m_sp,m_dry_final,bus_l,bus_w,bus_h, sp_h, sp_l)
%This function computes inertia for the given mass distribution
m_sb = (m_dry_final - m_sp); %Bus Mass
m_sp = m_sp/6; %Solar panel divided into 6 cells
l = bus_l;%Reassignment for easier computation
w = bus_w;
h = bus_h;
I_sbx = (m_sb/12)*(l^2 + h^2); %Compouting inertia assuming bus as a cube
I_sby = (m_sb/12)*(w^2 + h^2);
I_sbz = (m_sb/12)*(l^2 + w^2);
I_spy = m_sp*(sp_h^2)/12; %Computing inertia assuming solar panel as a plate
I_spz = m_sp*(sp_l^2)/12;
I_spx = I_spy + I_spz; %Perpendicular axis theorem
Izspt = (I_spz+(m_sp*(2.4^2)))+(I_spz+(m_sp * (3.8^2)))+(I_spz+(m_sp *(5.2^2))); %Parallel axis theorem
Izspt = 2*Izspt;
Ixspt = (I_spx+(m_sp*(2.4^2)))+(I_spx+(m_sp * (3.8^2)))+(I_spx+(m_sp *(5.2^2))); %Parallel axis theorem
Ixspt = 2*Ixspt;
Iyspt = 6*I_spy;
Ix = Ixspt + I_sbx; %Total inertia along x kg-m2
Iy = Iyspt + I_sby; %Total inertia along y kg-m2
Iz = Izspt + I_sbz; %Total inertia along z kg-m2
I = max(Ix,Iy);
I = max(I,Iz);

