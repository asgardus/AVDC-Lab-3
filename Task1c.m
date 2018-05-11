clc

%% Initialization

s = tf('s');
cp = 0.4;                   %damping coefficient
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cc = 2*sqrt(mp*kp);         %critical damping coefficient
zeta = cp/cc;               %damping ratio
omega = sqrt(kp/mp);        %natural frequency
omega_d = (1-zeta^2)*omega; %damped natural frequency
fr=0:0.001:10e2;
[mag,phase,wout]=bode(T,fr);
omega_max=wout(find(mag==max(mag))); %find resonant frequency

% Transfer function
T = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));

%% Plots

figure(1);
bode(T,fr)
hold on
figure(2);
step(T)
xlim([0 5])
hold on