%% ME 315 Project Simulation
% Written by John McCaw
% 3/4/2019
% This program computes transient heat transfer for a flat plate (with
% bi-axial symmetry); for the 315 project, this represents one-fourth of a
% 3D printer bed. 
%% Data: 
% Material Properties:
k_al = 237; % W/mK
Power = 12*5; % W
vol = (1.8/1000)*(120/1000)*(120/1000); % m^3
% qdot = Power/vol; % W/m^3
qdot = 9.8561*10^6; % Experimental Value
h_top =319.822; % W/m^2K Experimental Value solved from system of equations from experimental data
h_bottom = h_top; % W/m^2K
h = h_top; % W/m^2K
Cp_al = 900; % J/kgK 
rho_al = 2700; % kg/m^3
alpha = k_al/(rho_al*Cp_al); % m^2/s

% Mesh Properties:
dx = 0.02667; % m 
dy = 0.02167; % m
dz = 1.8/1000; % m
dt = 0.001; % s
Fo_x = alpha*dt/(dx^2); % [1]
Fo_y = alpha*dt/(dy^2); % [1]
end_time = 100; % sec
time_steps = ceil(end_time/dt);
data = zeros(4,4,time_steps);

% Environment:
T_amb = 23.6; % degC

% Node Areas and initial Temperatures
A11 = dx*dy/4; A12 = dx*dy/2; A13 = dx*dy/2; A14 = dx*dy/4; % m
T11 = T_amb; T12 = T_amb; T13 = T_amb; T14 = T_amb; % degC

A21 = dx*dy/2; A22 = dx*dy; A23 = dx*dy; A24 = dx*dy/2; % m
T21 = T_amb; T22 = T_amb; T23 = T_amb; T24 = T_amb; % degC

A31 = dx*dy/2; A32 = dx*dy; A33 = dx*dy; A34 = dx*dy/2; % m
T31 = T_amb; T32 = T_amb; T33 = T_amb; T34 = T_amb; % degC

A41 = dx*dy/4; A42 = dx*dy/2; A43 = dx*dy/2; A44 = dx*dy/4; % m
T41 = T_amb; T42 = T_amb; T43 = T_amb; T44 = T_amb; % degC

T = [T11, T12, T13, T14; T21, T22, T23, T24; T31, T32, T33, T34; T41, T42, T43, T44];
%% Transient Calculation:
for p = 1:time_steps
    % use forward definition to calc new time_step. 
    % store temp at timestep p in data array
    data(:,:,p) = T;
%     if(max(T) > 50)
%         qdot = 7500000;
%     elseif(max(T) < 47)
%         qdot = 9.8561*10^6; % Experimental Value
%     end
    T(1,1) = 2*Fo_y*T(2,1) + 2*Fo_x*T(1,2) + ((2*h*dt)/(rho_al*Cp_al*dz))*T_amb + (qdot*dt)/(rho_al*Cp_al) + (1-2*Fo_y - 2*Fo_x - (2*h*dt)/(dz*Cp_al*rho_al))*T(1,1);
    T(1,2) = Fo_y*(T(1,1) + T(1,3)) + 2*Fo_x*T(2,2) + ((2*h*dt)/(Cp_al*rho_al*dz))*T_amb + qdot*dt/(Cp_al*rho_al) + (1-2*Fo_y - 2*Fo_x - (2*h*dt)/(rho_al*Cp_al*dz))*T(1,2);
    T(1,3) = Fo_y*(T(1,2) + T(1,4)) + 2*Fo_x*T(2,3) + ((2*h*dt)/(Cp_al*rho_al*dz))*T_amb + qdot*dt/(Cp_al*rho_al) + (1-2*Fo_y - 2*Fo_x - (2*h*dt)/(rho_al*Cp_al*dz))*T(1,3);
    T(1,4) = 2*Fo_y*T(2,4) + 2*Fo_x*T(1,3) + ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz)*T_amb + qdot*dt/(rho_al*Cp_al) + (1-2*Fo_x - 2*Fo_y - ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz))*T(1,4);
    T(2,1) = Fo_y*(T(1,1) + T(3,1)) + 2*Fo_x*T(2,2) + ((2*h*dt)/(Cp_al*rho_al*dz))*T_amb + qdot*dt/(Cp_al*rho_al) + (1-2*Fo_y - 2*Fo_x - (2*h*dt)/(rho_al*Cp_al*dz))*T(2,1);
    T(2,2) = Fo_y*(T(3,2) + T(1,2)) + Fo_x*(T(2,3) + T(2,1)) + ((2*h*dt)/(Cp_al*rho_al*dz))*T_amb + (dt*qdot/(Cp_al*rho_al)) + (1-2*Fo_x -2*Fo_y - ((2*h*dt)/(Cp_al*rho_al*dz)))*T(2,2);
    T(2,3) = Fo_y*(T(3,3) + T(1,3)) + Fo_x*(T(2,4) + T(2,2)) + ((2*h*dt)/(Cp_al*rho_al*dz))*T_amb + (dt*qdot/(Cp_al*rho_al)) + (1-2*Fo_x -2*Fo_y - ((2*h*dt)/(Cp_al*rho_al*dz)))*T(2,3);
    T(2,4) = 2*Fo_x*T(2,3) + 2*Fo_y*(T(1,4) + T(3,4)) + ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz)*T_amb + (qdot*dt)/(rho_al*Cp_al) + (1 - 2*Fo_x -4*Fo_y - ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz))*T(2,4);
    T(3,1) = Fo_y*(T(2,1) + T(4,1)) + 2*Fo_x*T(3,2) + ((2*h*dt)/(Cp_al*rho_al*dz))*T_amb +qdot*dt/(Cp_al*rho_al) + (1-2*Fo_y - 2*Fo_x - (2*h*dt)/(rho_al*Cp_al*dz))*T(3,1);
    T(3,2) = Fo_y*(T(4,2) + T(2,2)) + Fo_x*(T(3,3) + T(3,1)) + ((2*h*dt)/(Cp_al*rho_al*dz))*T_amb + (dt*qdot/(Cp_al*rho_al)) + (1-2*Fo_x -2*Fo_y - ((2*h*dt)/(Cp_al*rho_al*dz)))*T(3,2);
    if (p < 5000) % This is for corner conditioning; finite difference can produce a peak at the node one element in from a corner for certain simulation conditions. This may be commented out for certain situations. 
        conditioner = 0.999;
    else
        conditioner = 0.99;
    end
    T(3,3) = Fo_y*(T(4,3) + T(2,3)) + Fo_x*(T(3,4) + T(3,2)) + ((2*h*dt)/(Cp_al*rho_al*dz))*(T_amb) + (dt*qdot/(Cp_al*rho_al)) + conditioner*(1-2*Fo_x -2*Fo_y - ((2*h*dt)/(Cp_al*rho_al*dz)))*T(2,2);
    T(3,4) = 2*Fo_x*T(3,3) + 2*Fo_y*(T(4,4) + T(2,4)) + ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz)*T_amb + (qdot*dt)/(rho_al*Cp_al) + (1 - 2*Fo_x -4*Fo_y - ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz))*T(3,4);
    T(4,1) = 2*Fo_y*T(4,2) + 2*Fo_x*T(3,1) + ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz)*T_amb + qdot*dt/(rho_al*Cp_al) + (1-2*Fo_x - 2*Fo_y - ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz))*T(4,1);
    T(4,2) = 2*Fo_x*T(3,2) + 2*Fo_y*(T(4,1) + T(4,3)) + ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz)*T_amb + (qdot*dt)/(rho_al*Cp_al) + (1 - 2*Fo_x -4*Fo_y - ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz))*T(4,2);
    T(4,3) = 2*Fo_x*T(3,3) + 2*Fo_y*(T(4,2) + T(4,4)) + ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz)*T_amb + (qdot*dt)/(rho_al*Cp_al) + (1 - 2*Fo_x -4*Fo_y - ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dz))*T(4,3);
    T(4,4) = 2*Fo_x*T(4,3) + 2*Fo_y*T(3,4) + ((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dy + 1/dz)*T_amb + (qdot*dt)/(rho_al*Cp_al) + (1-2*Fo_x - 2*Fo_y -((2*h*dt)/(rho_al*Cp_al))*(1/dx + 1/dy + 1/dz))*T(4,4);
    % loop through matrix T calculating each new temp at timestep
    % p+1
    % Assign new temp array to matrix T
end
%% Plotting: 
surf(data(:,:, 1));
set(gca,'nextplot','replacechildren'); 
vid = VideoWriter('test.avi');
open(vid);
[xx, yy] = meshgrid(1:0.1:10); 
for i =1:400:time_steps
    surf(interp2(data(:,:, i), xx, yy), 'EdgeColor', 'interp');
    view(-37.5+180,30)
    zlabel('^\circC')
    xlabel('X, m')
    ylabel('Y, m')
    tstring = strcat('t = ', num2str(round(i*dt, 2)));
    title({'Transient Conduction for Quarter Plate, Raised-Edge Insulation',tstring})
    tops = max(data(:,:,i));
    bottoms = min(data(:,:,i));
    zlim([bottoms(end)-2 (tops(1)+2)])
    frame = getframe(gcf);
    writeVideo(vid, frame);
end
close(vid);