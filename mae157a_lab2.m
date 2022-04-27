%{
 MAE 157A Lab 2, Trajectory Analysis
 Written by David Tran
%}

%% Clear Cache
clc; close all; clear all;

%% Load flight.mat variables
load('Flight1.mat');
load('Flight2.mat');
load('Flight3.mat');
load('Flight4.mat');

%% Constants
rho_air = 1.225;                %density of air in [kg/m^3]
g0 = 9.81;                       %acceleration of gravity [m/s^2]
launch_angle = 5 * pi / 180;    %initial launch angle [rad]
para_Cd = 1;                    %parachute coefficient of drag
rocket_Cd = 0.55;               %rocket coefficient of drag
para_D = 12 / 39.37;            %parachute diameter [m]
rocket_D1 = 0.98 / 39.37;         %rocket 1 diameter [m]
rocket_D2 = 1.38 / 39.37;         %rocket 2 diameter [m]
A_rocket1 = pi * rocket_D1 ^2 /4; %cross-sectional area of rocket 1 in [m2]
A_rocket2 = pi * rocket_D2 ^2 /4; %cross-sectional area of rocket 2 in [m2]
A_para = pi * para_D^2 / 4;
t_delay = 4;                    %ejection charge delay in [s]

%rocket launch masses
M_f1 = 0.0892;                        %final mass L1 (dry mass) [kg]
M_f2 = 0.0892;                        %final mass L2 (dry mass) [kg]
M_f3 = 0.0892;                        %final mass L3 (dry mass) [kg]
M_f4 = 0.0504;                        %final mass L4 (dry mass) [kg]

M_B4 = 0.019;                      %mass of estes b4 in [kg]
M_A8 = 0.016;                      %mass of estes b6 in [kg]

M_o1 = M_f1 + M_B4;                %initial mass L1 (wet mass) [kg]
M_o2 = M_f2 + M_A8;                %initial mass L2 (wet mass) [kg]
M_o3 = M_f3 + M_A8;                %initial mass L3 (wet mass) [kg]
M_o4 = M_f4 + M_B4;                %initial mass L4 (wet mass) [kg]
M_o5 = M_f5 + M_B4;                %initial mass L5 (wet mass) [kg]
M_o6 = M_f6 + M_A8;                %initial mass L6 (wet mass) [kg]

%Load data from thrustcurve.org (or manually type it bc idek how to get it
%in csv format from that website)
T_B4 = [0.418 1.673 4.076 6.69 9.304 11.496 12.75 11.916 10.666 9.304 ...
    7.214 5.645 4.809 4.182 3.763 3.554 3.345 3.345 3.345 3.345 3.449 ...
    3.449 2.404 1.254 0 ]; %thrust array for b4-4
t_B4 = [0.02 0.04 0.065 0.085 0.105 0.119 0.136 0.153 0.173 0.187 0.198 ... 
    0.207 0.226 0.258 0.326 0.422 0.549 0.665 0.776 0.863 0.94 0.991 ...
    1.002 1.01 1.03 ]; %time array for b4-4

T_A8 = [0.01635 0.512 2.115 4.358 6.794 8.588 9.294 9.730 8.845 7.179 ...
    5.063 3.717 3.205 2.884 2.499 2.371 2.307 2.371 2.371 2.243 1.794 ...
    1.153 0.448 0.0]; %thrust array for A8-3
t_A8 = [0.0033 0.041 0.084 0.127 0.166 0.192 0.206 0.226 0.236 0.247 ...
    0.261 0.277 0.306 0.351 0.405 0.467 0.532 0.579 0.632 0.652 0.668 ...
    0.684 0.703 0.73]; %time array for b6-4

m_dot_b4 = 0.006 / t_B4(end);            %propellant burn rate of B4-4 in [kg/s]
m_dot_a8 = 0.003 / t_A8(end);           %propellant burn rate of B6-4 in [kg/s]

%% Loop and Eulerian step simulation
t = 0;              %Initial time [sec]
dstep = 0.01;          %time step [sec]
t_final = 7;        %final time [sec]
alt = [];            %Altitude [m]
vel = [];            %Velocity [m/s]
acc = [];            %Acceleration [m/s^2]

%Simulation Loop

%L1
ii = 2;
M_o1t = [];
z41 = [];
z41(1) = 0;
v41 = [];
v41(1) = 0;
a41 = [];
T41 = [];
T41(1) = 0;
a41(1) = 0;
t = [];
t(1) = 0;
delta_t = 0.02;
while t(ii-1) < t_final
    t(ii) = t(ii -1) + delta_t;
    if t(ii) < t_B4(end)
        %burn phase
        M_o1t(ii) = M_o1 - m_dot_b4 * t(ii);
        W41(ii) = M_o1t(ii) * g0;
        if t(ii) > 0
            T41(ii) = interp1(t_B4,T_B4,t(ii),'linear', 'extrap'); %thrust
%             disp(T41(ii));
        end
%         I_sp(ii) = T41(ii) / (m_dot_b6 * g);
        v41(ii) = v41(ii-1) + a41(ii-1)*delta_t;
        D41(ii) = 0.5 .* rho_air .* v41(ii)^2 .* rocket_Cd .* A_rocket1;
        a41(ii) = (T41(ii) - W41(ii) - D41(ii)) ./ M_o1t(ii);
        z41(ii) = z41(ii-1) + v41(ii-1)*delta_t; %in [m]
        
    elseif t(ii) >= t_B4(end) && t(ii) < (t_B4(end) + t_delay)
        %no more propellant but no delay charge yet
        T41(ii) = 0; %no more thrust :(
        M_o1t(ii) = M_f1; % mass does not change
        W41(ii) = M_o1t(ii) * g0;
        v41(ii) = v41(ii-1) + a41(ii-1)*delta_t;
        D41(ii) = 0.5 .* rho_air .* v41(ii)^2 .* rocket_Cd .* A_rocket1;
        a41(ii) = (T41(ii) - W41(ii) - D41(ii)) ./ M_o1t(ii);
        z41(ii) = z41(ii-1) + v41(ii-1)*delta_t; %in [m]
    elseif t(ii) >= (t_B4(end) + t_delay)
        T41(ii) = 0;  %no more thrust :(
        M_o1t(ii) = M_f1; % mass does not change
        W41(ii) = M_o1t(ii) * g0;
        v41(ii) = v41(ii-1) + a41(ii-1)*delta_t;
        D41_para = 0.5 * para_Cd * v41(ii)^2 * A_para * rho_air;
        a41(ii) = (T41(ii) - W41(ii) + D41_para) ./ M_o1t(ii);
        z41(ii) = z41(ii-1) + v41(ii-1)*delta_t; %in [m]
    end

ii = ii + 1;
end

[dummy_var, ind_max] = max(z41);
[dummy_var, ind_exp_max] = max(table2array(Flight1(:, 3)));

%L2
ii = 2;
M_o2t = [];
z42 = [];
z42(1) = 0;
v42 = [];
v42(1) = 0;
a42 = [];
a42(1) = 0;
T42 = [];
T42(1) = 0;
t42 = [];
t42(1) = 0;
while t42(ii-1) < t_final
    t42(ii) = t42(ii -1) + delta_t;
    if t42(ii) < t_B4(end)
        %burn phase
        M_o2t(ii) = M_o2 - m_dot_a8 * t42(ii);
        W42(ii) = M_o2t(ii) * g0;
        if t42(ii) > 0
            T42(ii) = interp1(t_B4,T_B4,t42(ii),'linear', 'extrap'); %thrust
        end
        v42(ii) = v42(ii-1) + a42(ii-1)*delta_t;
        D42(ii) = 0.5 .* rho_air .* v42(ii)^2 .* rocket_Cd .* A_rocket2;
        a42(ii) = (T42(ii) - W42(ii) - D42(ii)) ./ M_o2t(ii);
        z42(ii) = z42(ii-1) + v42(ii-1)*delta_t; %in [m]
        
    elseif t42(ii) >= t_B4(end) && t42(ii) < (t_B4(end) + t_delay)
        %no more propellant but no delay charge yet
        T42(ii) = 0; %no more thrust :(
        M_o2t(ii) = M_f2; % mass does not change
        W42(ii) = M_o2t(ii) * g0;
        v42(ii) = v42(ii-1) + a42(ii-1)*delta_t;
        D42(ii) = 0.5 .* rho_air .* v42(ii)^2 .* rocket_Cd .* A_rocket2;
        a42(ii) = (T42(ii) - W42(ii) - D42(ii)) ./ M_o2t(ii);
        z42(ii) = z42(ii-1) + v42(ii-1)*delta_t; %in [m]
    elseif t42(ii) >= (t_B4(end) + t_delay)
        T42(ii) = 0;  %no more thrust :(
        M_o2t(ii) = M_f2; % mass does not change
        W42(ii) = M_o2t(ii) * g0;
        v42(ii) = v42(ii-1) + a42(ii-1)*delta_t;
        D42_para = 0.5 * para_Cd * v42(ii)^2 * A_para * rho_air;
        a42(ii) = (T42(ii) - W42(ii) + D42_para) ./ M_o2t(ii);
        z42(ii) = z42(ii-1) + v42(ii-1)*delta_t; %in [m]
    end

ii = ii + 1;
end

[dummy_var, ind_max_2] = max(z42);
[dummy_var, ind_exp_max_2] = max(table2array(Flight2(:, 3)));

%L3
ii = 2;
M_o3t = [];
z8 = [];
z8(1) = 0;
v8 = [];
v8(1) = 0;
a8 = [];
a8(1) = 0;
T8 = [];
T8(1) = 0;
t8 = [];
t8(1) = 0;
while t8(ii-1) < t_final
    t8(ii) = t8(ii -1) + delta_t;
    if t8(ii) < t_A8(end)
        %burn phase
        M_o3t(ii) = M_o3 - m_dot_a8 * t42(ii);
        W8(ii) = M_o3t(ii) * g0;
        if t8(ii) > 0
            T8(ii) = interp1(t_A8,T_A8,t8(ii),'linear', 'extrap'); %thrust
        end
        v8(ii) = v8(ii-1) + a8(ii-1)*delta_t;
        D8(ii) = 0.5 .* rho_air .* v8(ii)^2 .* rocket_Cd .* A_rocket1;
        a8(ii) = (T8(ii) - W8(ii) - D8(ii)) ./ M_o3t(ii);
        z8(ii) = z8(ii-1) + v8(ii-1)*delta_t; %in [m]
        
    elseif t8(ii) >= t_A8(end) && t8(ii) < (t_A8(end) + t_delay)
        %no more propellant but no delay charge yet
        T8(ii) = 0; %no more thrust :(
        M_o3t(ii) = M_f3; % mass does not change
        W8(ii) = M_o3t(ii) * g0;
        v8(ii) = v8(ii-1) + a8(ii-1)*delta_t;
        D8(ii) = 0.5 .* rho_air .* v8(ii)^2 .* rocket_Cd .* A_rocket1;
        a8(ii) = (T8(ii) - W8(ii) - D8(ii)) ./ M_o3t(ii);
        z8(ii) = z8(ii-1) + v8(ii-1)*delta_t; %in [m]
    elseif t8(ii) >= (t_A8(end) + t_delay)
        T8(ii) = 0;  %no more thrust :(
        M_o3t(ii) = M_f3; % mass does not change
        W8(ii) = M_o3t(ii) * g0;
        v8(ii) = v8(ii-1) + a42(ii-1)*delta_t;
        D8_para = 0.5 * para_Cd * v8(ii)^2 * A_para * rho_air;
        a8(ii) = (T8(ii) - W8(ii) + D8_para) ./ M_o3t(ii);
        z8(ii) = z8(ii-1) + v8(ii-1)*delta_t; %in [m]
    end

ii = ii + 1;
end

[dummy_var, ind_max_3] = max(z8);
[dummy_var, ind_exp_max_3] = max(table2array(Flight3(:, 3)));

%L4
ii = 2;
M_o4t = [];
z43 = [];
z43(1) = 0;
v43 = [];
v43(1) = 0;
a43 = [];
a43(1) = 0;
T43 = [];
T43(1) = 0;
t43 = [];
t43(1) = 0;
while t43(ii-1) < t_final
    t43(ii) = t43(ii -1) + delta_t;
    if t43(ii) < t_B4(end)
        %burn phase
        M_o4t(ii) = M_o4 - m_dot_b4 * t43(ii);
        W43(ii) = M_o4t(ii) * g0;
        if t43(ii) > 0
            T43(ii) = interp1(t_B4,T_B4,t43(ii),'linear', 'extrap'); %thrust
        end
        v43(ii) = v43(ii-1) + a43(ii-1)*delta_t;
        D43(ii) = 0.5 .* rho_air .* v43(ii)^2 .* rocket_Cd .* A_rocket1;
        a43(ii) = (T43(ii) - W43(ii) - D43(ii)) ./ M_o4t(ii);
        z43(ii) = z43(ii-1) + v43(ii-1)*delta_t; %in [m]
        
    elseif t43(ii) >= t_B4(end) && t43(ii) < (t_B4(end) + t_delay)
        %no more propellant but no delay charge yet
        T43(ii) = 0; %no more thrust :(
        M_o4t(ii) = M_f4; % mass does not change
        W43(ii) = M_o4t(ii) * g0;
        v43(ii) = v43(ii-1) + a43(ii-1)*delta_t;
        D43(ii) = 0.5 .* rho_air .* v43(ii)^2 .* rocket_Cd .* A_rocket1;
        a43(ii) = (T43(ii) - W43(ii) - D43(ii)) ./ M_o4t(ii);
        z43(ii) = z43(ii-1) + v43(ii-1)*delta_t; %in [m]
    elseif t43(ii) >= (t_B4(end) + t_delay)
        T43(ii) = 0;  %no more thrust :(
        M_o4t(ii) = M_f4; % mass does not change
        W43(ii) = M_o4t(ii) * g0;
        v43(ii) = v43(ii-1) + a43(ii-1)*delta_t;
        D43_para = 0.5 * para_Cd * v43(ii)^2 * A_para * rho_air;
        a43(ii) = (T43(ii) - W43(ii) + D43_para) ./ M_o4t(ii);
        z43(ii) = z43(ii-1) + v43(ii-1)*delta_t; %in [m]
    end

ii = ii + 1;
end

[dummy_var, ind_max_4] = max(z43);
[dummy_var, ind_exp_max_4] = max(table2array(Flight4(:, 3)));

%plot only up to apogee
figure;
plot(t(1:ind_max), z41(1:ind_max).*3.28084,'r', 'LineWidth', 1.5);
hold on;
plot(t42(1:ind_max_2), z42(1:ind_max_2).*3.28084, 'm' ,'LineWidth', 1.5);
plot(t8(1:ind_max_3), z8(1:ind_max_3).*3.28084, 'b' ,'LineWidth', 1.5);
plot(t43(1:ind_max_4), z43(1:ind_max_4).*3.28084, 'k' ,'LineWidth', 1.5);
scatter(table2array(Flight1(1:ind_exp_max, 1)), table2array(Flight1(1:ind_exp_max, 3)), 12, 'filled' ,'r');
scatter(table2array(Flight2(1:ind_exp_max_2, 1)), table2array(Flight2(1:ind_exp_max_2, 3)), 12, 'filled', 'm');
scatter(table2array(Flight3(1:ind_exp_max_3, 1)), table2array(Flight3(1:ind_exp_max_3, 3)), 12, 'filled', 'b');
scatter(table2array(Flight4(1:ind_exp_max_4, 1)), table2array(Flight4(1:ind_exp_max_4, 3)), 12, 'filled', 'k');
xlabel('Time $t$ (s)','Interpreter','latex');
ylabel('Altitude $h$ (ft)', 'Interpreter','latex');
legend('S1', 'S2', 'S3', 'S4', 'L1', 'L2', 'L3', 'L4', 'location','best');
set(gca,'FontSize',16)
xlim([0 6]);
ylim([0 500]);


%% Experimental Data
% I clicked on the xlsx file of interest in the "Current Folder" tab in 
% MATLAB and imported only the numbers

t_burnout = 1.03;    %burnout time for estes b4-4 in [s]

% obtain apogee, burnout, and landing
apogee = [];
ind_ap = [];
burnout = [];
ind_burn = [];

[apogee(1), ind_ap(1)] = max(table2array(Flight1(:,3)));
for ii = 1:length(table2array(Flight1(:,7)))
    if table2array(Flight1(ii,7)) <= 0.25
        burnout(1) = table2array(Flight1(ii,7));
        ind_burn(1) = ii;
        break;
    end
end

%Flight 1
%alt. vs. time
figure;
scatter(table2array(Flight1(:,1)), table2array(Flight1(:,3)), 12, 'filled');
hold on;
xline(table2array(Flight1(64,1)), '-.', 'Burnout');
xline(table2array(Flight1(ind_ap(1),1)),'-.', 'Apogee');
xline(table2array(Flight1((end - 5),1)),'-.', 'Landing');
xlabel('Time $t$ (s)', 'Interpreter','latex');
ylabel('Altitude $h$ (ft)', 'Interpreter','latex');
set(gca,'FontSize',16);
grid on;
legend('Measured', 'location', 'best');
xlim([0 10]);

%vel. vs. time
figure;
dt = table2array(Flight1(2,1)) - table2array(Flight1(1,1));
vel = movmean(gradient(table2array(Flight1(:,3))) ./ gradient(table2array(Flight1(:,1))),4);
scatter(table2array(Flight1(:,1)), vel, 12, 'filled');
hold on;
xline(table2array(Flight1(64,1)), '-.', 'Burnout');
xline(table2array(Flight1(ind_ap(1),1)),'-.', 'Apogee');
xline(table2array(Flight1((end - 5),1)),'-.', 'Landing');
xlabel('Time $t$ (s)', 'Interpreter','latex');
ylabel('Velocity $v_z$ (ft/s)', 'Interpreter','latex');
set(gca,'FontSize',16);
grid on;
legend('Measured', 'location', 'best');
hold off;
xlim([0 10]);

%accel. vs. time
a = movmean(gradient(vel) ./ gradient(table2array(Flight1(:,1))),8);
figure;
hold on;
scatter(table2array(Flight1(:,1)), a, 12, 'filled');
xline(table2array(Flight1(64,1)), '-.', 'Burnout');
xline(table2array(Flight1(ind_ap(1),1)),'-.', 'Apogee');
xline(table2array(Flight1((end - 5),1)),'-.', 'Landing');
xlabel('Time $t$ (s)', 'Interpreter','latex');
ylabel('Acceleration $a_z$ (ft/s$^2$)','Interpreter','latex');
set(gca,'FontSize',16);
grid on;
legend('Measured', 'location', 'best');
xlim([0 10]);


%% Print Quantities
vel2 = movmean(gradient(table2array(Flight2(:,3))) ./ gradient(table2array(Flight2(:,1))),4);
vel3 = movmean(gradient(table2array(Flight3(:,3))) ./ gradient(table2array(Flight3(:,1))),4);
vel4 = movmean(gradient(table2array(Flight4(:,3))) ./ gradient(table2array(Flight4(:,1))),4);

disp('%======================SIMULATION RESULTS=========================');
disp('---MEASURED APOGEE---');
fprintf('L1: %d [ft]\n', table2array(Flight1(ind_exp_max, 3)));
fprintf('L2: %d [ft]\n', table2array(Flight2(ind_exp_max_2, 3)));
fprintf('L3: %d [ft]\n', table2array(Flight3(ind_exp_max_3, 3)));
fprintf('L4: %d [ft]\n\n', table2array(Flight4(ind_exp_max_4, 3)));

disp('---SIMULATED APOGEE---');
fprintf('S1: %d [ft]\n', max(z41)*3.28084);
fprintf('S2: %d [ft]\n', max(z42)*3.28084);
fprintf('S3: %d [ft]\n', max(z8)*3.28084);
fprintf('S4: %d [ft]\n\n', max(z43)*3.28084);

disp('---BURNOUT VEL EXP---');
% L3 is flying in the range 206-310, calculate all values here
% A8-3 burns for 0.7s which is index 219
fprintf('L1: %d [ft/s]\n', vel(64)); %burnout at around 1 s
fprintf('L2: %d [ft/s]\n', vel2(68)); %burnout at around 1 s
fprintf('L3: %d [ft/s]\n', vel3(219)); %burnout at around 0.55s
fprintf('L4: %d [ft/s]\n\n', vel4(64)); %burnout at around 1s

disp('---BURNOUT VEL SIM---');
fprintf('S1: %d [ft/s]\n', max(v41)*3.28084);
fprintf('S2: %d [ft/s]\n', max(v42)*3.28084);
fprintf('S3: %d [ft/s]\n', max(v8)*3.28084);
fprintf('S4: %d [ft/s]\n\n', max(v43)*3.28084);

disp('---TIME TO APOGEE EXP---');
fprintf('L1: %d [s]\n', table2array(Flight1(ind_exp_max, 1)));
fprintf('L2: %d [s]\n', table2array(Flight2(ind_exp_max_2, 1)));
fprintf('L3: %d [s]\n', table2array(Flight3(ind_exp_max_3, 1)));
fprintf('L4: %d [s]\n\n', table2array(Flight4(ind_exp_max_4, 1)));

disp('---TIME TO APOGEE SIM---');
fprintf('S1: %d [s]\n', t(ind_max));
fprintf('S2: %d [s]\n', t42(ind_max_2));
fprintf('S3: %d [s]\n', t8(ind_max_3));
fprintf('S4: %d [s]\n\n', t43(ind_max_4));


%% Altimeter Analysis
p_unc = 0.2; %plus-minus uncertainty of pressure reading in [psi]
p_unc_si = 6895 * p_unc; %convert pressure uncertainty to [Pa]

load('atmos.mat'); %load file

dPdh = gradient(P_Pa(:)) ./ gradient(h_m(:)); %pressure gradient wrt height [Pa/m]

h_unc = p_unc_si ./ dPdh; %height uncertainty in [m]
h_unc_ft = h_unc .* 3.281; %height uncertainty in [ft]

disp('%============================ALTIMETER ANALYSIS=====================');
for ii = 1:length(h_unc_ft)
    if h_unc_ft(ii) <= -500
        fprintf('Height where uncertainty is +/- 500 ft: %d [m] = %d [ft]\n\n', h_m(ii), h_m(ii)*3.281);
        break;
    end
end

%% Space Shot

%Constants
A_exit = 0.025;                     %exit area in [m^2]
P_exit = 101325;                    %exit pressure in [Pa]
T_jet = 4500 * 4.448;               %jet thrust in [N]
I_sp_SL = 220;                      %jet sea-level specific impulse [s]
m_prop = 850 / 2.205;               %propellant mass in [kg]
m_payload = 400 / 2.205;            %payload mass in [kg]
D = 1.25 / 3.281;                   %tank diameter in [m]
A_tank = pi * D^2 / 4;              %tank area in [m^2]
C_D_tank = 0.75;                    %tank drag coefficient
p_chamber = 200 * 6895;             %chamber pressure in [Pa]
rho_prop = 48.3 * 16.018;           %propellant density in [kg/m^3]
m_l_tank = 0.1905 * 17.858;         %mass per unit length of the tank in [kg/m]
m_dot_ss = T_jet / (I_sp_SL * g0);  %mass flow rate in [kg/s]
t_burn = m_prop / m_dot_ss;         %burn time in [s]
R_e = 6.38e6;                       %radius of earth in [m]
L_tank = m_prop / (rho_prop * A_tank);  %length of tank in [m]
m_tank = L_tank * m_l_tank; %tank mass in [kg]
m_tot_i = m_prop + m_payload + m_tank;  %total initial mass in [kg]
m_empty = m_tank + m_payload; %empty mass in [kg]


%account for variation of atm properties with altitude
P_amb = P_Pa;
T_press = (P_exit - P_amb) .* A_exit;
T_tot = T_jet + T_press;

%loop through time


%calc traj
ii = 2;
M_ss = [];
zss = [];
zss(1) = 0;o90p21
vss = [];
vss(1) = 0;
a_ss = []; %lol
a_ss(1) = 0;
Tss = [];
Tss(1) = 0;
tss = [];
tss(1) = 0;
g = [];
g(1) = g0;
g(2) = g0;
t_end_ss = 300;
delta_t_ss = 0.5;
M_ss = [];
M_ss(1) = m_tot_i;
Wss = [];
Wss(1) = m_tot_i * g0;
while tss(ii-1) < t_end_ss
    tss(ii) = tss(ii -1) + delta_t_ss;
    if tss(ii) < t_burn
        %burn phase
        M_ss(ii) = m_tot_i - m_dot_ss * tss(ii);
        zss(ii) = zss(ii-1) + vss(ii-1)*delta_t_ss; %in [m]
        g(ii) = g0 * (R_e / (R_e + zss(ii)))^2;
        vss(ii) = vss(ii-1) + a_ss(ii-1)*delta_t_ss;
        Wss(ii) = M_ss(ii) * g(ii);
        if tss(ii) > 0
            Tss(ii) = interp1(h_m,T_tot,zss(ii),'linear', 'extrap'); %thrust
            rho_h(ii) = interp1(h_m, rho_kg_m3,zss(ii),'linear','extrap');
        end
        Dss(ii) = 0.5 .* rho_h(ii) .* vss(ii)^2 .* C_D_tank .* A_tank; %change properties
        a_ss(ii) = (Tss(ii) - Wss(ii) - Dss(ii)) ./ M_ss(ii);
        
    else
        Tss(ii) = 0;  %no more thrust :(
        M_ss(ii) = m_empty; % mass does not change
        zss(ii) = zss(ii-1) + vss(ii-1)*delta_t_ss; %in [m]
        g(ii) = g0 * (R_e / (R_e + zss(ii)))^2;
        Wss(ii) = M_ss(ii) * g(ii);
        if tss(ii) > 0
            rho_h(ii) = interp1(h_m, rho_kg_m3,zss(ii),'linear','extrap');
        end
        vss(ii) = vss(ii-1) + a_ss(ii-1)*delta_t_ss;
        Dss(ii) = 0.5 * C_D_tank * vss(ii)^2 * A_tank * rho_h(ii);
        a_ss(ii) = (Tss(ii) - Wss(ii) - Dss(ii)) ./ M_ss(ii);

    end
ii = ii + 1;
end

%Print apogee, burnout vel and alt
[max_ss, max_ss_ind] = max(vss);
[max_zss, max_zss_ind] = max(zss);
disp('%==============================SPACE SHOT==========================');
fprintf('Apogee: %d [ft]\n', max_zss*3.28084);
fprintf('Burnout velocity: %d [ft/s]\n', max_ss*3.28084);
fprintf('Burnout altitude: %d [ft]\n', zss(max_ss_ind)*3.28084);
fprintf('Burnout thrust: %d [N] = %d [lbf]\n\n', T_tot(max_ss_ind), T_tot(max_ss_ind) / 4.448);

%Plot altitude vs time
figure;
plot(tss, zss.*3.28084, 'LineWidth', 2);
xlabel('Time $t$ (s)', 'Interpreter','latex');
ylabel('Altitude $h$ (ft)','Interpreter','latex');
xline(tss(max_ss_ind), '-.', 'Burnout');
xline(tss(max_zss_ind), '-.', 'Apogee');
set(gca,'FontSize',16);