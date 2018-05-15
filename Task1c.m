clc
clear all
close all

%% Initialization

s = tf('s');
dampz=1;
switch dampz
    case 1  %for underdamped
        cp=0.4;
    case 2  %for critically damped
        cp=2.0112
    case 3  %for overdamped
        cp=3
end
% cp =0.4;                   %damping coefficient
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cc = 2*sqrt(mp*kp);         %critical damping coefficient
zeta = cp/cc;               %damping ratio
omega = sqrt(kp/mp);        %natural frequency
omega_d = (1-zeta^2)*omega; %damped natural frequency
% Transfer function
G1 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));
G2 = tf([cp kp],[mp cp kp])

fr=0:0.001:10e2;
[mag,phase,wout]=bode(G1,fr);
omega_max=wout(find(mag==max(mag))); %find resonant frequency


%% Task 1d
f=0.5;
dt=0.1;
Tfinal=10;
t = 0:dt:Tfinal
u=0.05*sin(2*pi*f*t);
figure (3)
lsim(G2,u,t);        %sine input

figure (4)
impulse(0.05*G2,1);  %Impulse input

%% Plots
% 
% figure(1);
% bode(G,fr)
% hold on
% figure(2);
% step(G)
% xlim([0 5])
% hold on