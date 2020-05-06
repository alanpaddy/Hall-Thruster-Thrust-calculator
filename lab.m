 group_1=[11	0	0	10.52	1.769	95.2	0.82
11	0	0	11.22	1.769	98.3	1.012
11	0	0	11.81	1.769	96.9	1.365
11	0	0	12.45	1.77	104.3	1.562
11	0	0	13.21	1.769	110.8	1.758
11	0	0	13.87	1.77	114.2	1.797
11	0	0	14.67	1.77	114.5	1.836
];

group_2=[14	0	0	13.77	1.505	80.8	1.014
14	0	0	13.82	1.505	79	1.208
14	0	0	13.88	1.505	80	1.403
14	0	0	13.93	1.505	78.4	1.6
14	0	0	13.97	1.505	80.4	1.798
14	0	0	14	1.505	82.3	1.913
14	0	0	14.02	1.505	85.4	2.11
14	0	0	14.05	1.505	91.8	2.306
14	0	0	14.1	1.505	99.2	2.501

];

T_laser_1 =[0.001315628	0.001782536	0.002657029	0.00347071	0.003701098	0.003928469	0.003888689];

T_laser_2 =[0.001088082	0.001532437	0.00215997	0.001957768	0.002177051	0.002268518	0.002429657	0.002782394	0.003415658];




theta=30; %Plume divergence half angle
nu_cos=1/2*(1+cos(deg2rad(theta))^2); 
Ipp_Ip=0.2; %Proportion of doubly charged ions
alpha=(1+1/sqrt(2)*Ipp_Ip)/(1+Ipp_Ip);
nu_m=0.88; % mass utilization efficiency
g_0=9.80665; %Gravity [m/s^2]
q=1.6*10^-19; %Fundamental charge (C)
M_a =83.798; %Atomic mass of Kripton Ion [amu]
V_plasma=10; %10 [V]
V_cathode_ground=30; %Potential of cathode below ground[V]
M_ion = 1.46*10^-25; %Mass of Kripton ion [kg]
nu_d =0.77; %Beam of current fraction in discharge current
m_cathode = 7.436*10^-4*M_a*10^-6*8; %8 is cathode flow rate sscm VARIABLE!

I_cathode_1 = 1.0028; % 1.001 [A]
V_cathode_1 = 50.8; % 44.7 [V]
P_cathode_1 = V_cathode_1*I_cathode_1;
I_cathode_2 = 1.028; % 1[A]
V_cathode_2 = 50.7; % 46.9 [A]
P_cathode_2 = V_cathode_2*I_cathode_2;

b=size(group_1);
c=size(group_2);

%Group 1 (anode flow =7)
for i=1:b(1)
    V_anode_1(i) = group_1(i,6);
    I_anode_1(i) = group_1(i,7);
    V_b_1(i) = V_anode_1(i)-(V_plasma+V_cathode_ground);
    I_b_1(i) = I_anode_1(i)*nu_d;
    
    m_anode_1(i)=7.436*10^-4*M_a*10^-6*group_1(i,1);
    P_anode_1(i) =V_anode_1(i)* I_anode_1(i);
    P_electromag_1(i) = group_1(i,2)*group_1(i,3)+group_1(i,4)*group_1(i,5); % P=V*I+V*I
      
    T_2_1(i) = alpha*nu_cos*sqrt((2*M_ion)/q)*I_b_1(i)*sqrt(V_b_1(i));
    
    Isp_1(i) = alpha*nu_cos*nu_m/g_0*sqrt((2*q*V_b_1(i))/M_ion);
    Isp_T_2_1(i) = T_2_1(i)/(m_anode_1(i)*g_0);
    Isp_T_Laser_1(i) = T_laser_1(i)/(m_anode_1(i)*g_0);
    
    T_1_1(i) = m_anode_1(i)*Isp_1(i)*g_0;
    
    nu_anode_1_T_1_1(i) = T_1_1(i).^2/(2*m_anode_1(i)*P_anode_1(i));
    nu_anode_1_T_2_1(i) = T_2_1(i).^2/(2*m_anode_1(i)*P_anode_1(i));
    nu_anode_1_T_laser_1(i) = T_laser_1(i).^2/(2*m_anode_1(i)*P_anode_1(i));
    
    nu_total_1(i) = T_2_1(i).^2/(2*(m_anode_1(i)+m_cathode)*(P_anode_1(i)+P_cathode_1+P_electromag_1(i)));
end

%Group 2 (anode flow = 12)
for j=1:c(1)
    V_anode_2(j) = group_2(j,6);
    I_anode_2(j) = group_2(j,7);
    V_b_2(j) = V_anode_2(j)-(V_plasma+V_cathode_ground);
    I_b_2(j) = I_anode_2(j)*nu_d;
    
    m_anode_2(j)=7.436*10^-4*M_a*10^-6*group_2(j,1);
    P_anode_2(j) = V_anode_2(j)*I_anode_2(j);
    P_electromag_2(j) = group_2(j,2)*group_2(j,3)+group_2(j,4)*group_2(j,5); % P=V*I
    
    T_2_2(j) = alpha*nu_cos*sqrt((2*M_ion)/q)*I_b_2(j)*sqrt(V_b_2(j));
    
    Isp_2(j) = alpha*nu_cos*nu_m/g_0*sqrt((2*q*V_b_2(j))/M_ion);
    Isp_T_Laser_2(j) = T_laser_2(j)/(m_anode_2(j)*g_0);
    Isp_T_2_2(j) = T_2_2(j)/(m_anode_2(j)*g_0);
    
    T_1_2(j) = m_anode_2(j)*Isp_2(j)*g_0;
    
    nu_anode_2_T_1_2(j) = T_1_2(j).^2/(2*m_anode_2(j)*P_anode_2(j));
    nu_anode_2_T_2_2(j) = T_2_2(j).^2/(2*m_anode_2(j)*P_anode_2(j));
    nu_anode_2_T_laser_2(j) = T_laser_2(j).^2/(2*m_anode_2(j)*P_anode_2(j));
    nu_total_2(j) =T_2_2(j).^2/(2*(m_anode_2(j)+m_cathode)*(P_anode_2(j)+P_cathode_2+P_electromag_2(j)));
end




%Plot thrusts v Power
plot(P_anode_1, T_1_1,':k*'); 
hold on; 
plot(P_anode_2, T_1_2,'--ko'); 
hold on;
plot( P_anode_1,T_2_1,':b*');
hold on;
plot( P_anode_2,T_2_2,'--bo');
hold on;
plot( P_anode_1,T_laser_1,':m*');
hold on;
plot(P_anode_2, T_laser_2,'--mo');
title('Thrust vs Anode Power')
xlabel('Anode Power [W]') % x-axis label
ylabel('Thrust [N]') % y-axis label
legend('Thrust 1, 11 sscm', 'Thrust 1, 14 sscm','Thrust 2, 11 sscm ', 'Thrust 2, 14 sscm', 'Thrust Laser, 11 sscm', 'Thrust Laser, 14 sscm')
grid on






%Plot anode efficiency v Power
figure
plot(P_anode_1, nu_anode_1_T_1_1,':k*');
hold on
plot(P_anode_2, nu_anode_2_T_1_2,'--ko');
hold on;
plot(P_anode_1, nu_anode_1_T_2_1,':b*');
hold on
plot(P_anode_2, nu_anode_2_T_2_2,'--bo');
hold on;
plot(P_anode_1, nu_anode_1_T_laser_1,':m*');
hold on
plot(P_anode_2, nu_anode_2_T_laser_2,'--mo');
title('Anode Efficiency vs Anode Power')
xlabel('Anode Power [W]') % x-axis label
ylabel('Anode Efficiency') % y-axis label
legend('Thrust 1, 11 sscm', 'Thrust 1, 14 sscm','Thrust 2, 11 sscm', 'Thrust 2, 14 sscm','Thrust Laser, 11 sscm', 'Thrust Laser, 14 sscm')
grid on



%Plot Isp v Power
figure
plot(P_anode_1, Isp_1,':k*');
hold on
plot(P_anode_2, Isp_2,'--ko');
hold on;
plot(P_anode_1, Isp_T_2_1,':b*');
hold on
plot(P_anode_2, Isp_T_2_2,'--bo');
hold on;
plot(P_anode_1, Isp_T_Laser_1,':m*');
hold on
plot(P_anode_2, Isp_T_Laser_2,'--mo');
title('Isp vs Anode Power')
xlabel('Anode Power [W]') % x-axis label
ylabel('Isp (s)') % y-axis label
legend('Direct, 11 sscm', 'Direct, 14 sscm','Thrust 2, 11 sscm', 'Thrust 2, 14 sscm','Thrust Laser, 11 sscm', 'Thrust Laser, 14 sscm')
grid on

%Plot I discharge vs V discharge
figure
plot(I_anode_1, V_anode_1,':k*');
hold on
plot(I_anode_2, V_anode_2,'--ko');
hold on;
title('Discharge current vs discharge voltage')
xlabel('Anode current [A]') % x-axis label
ylabel('Anode Voltage [V]') % y-axis label
legend(' 11 sscm', '14 sscm')
grid on



