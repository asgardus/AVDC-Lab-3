clear;
clc

%% Initialization 2DOF

s = tf('s');
cp = 0.8;                   %damping coefficient
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cs = 0.05;                  %damping coefficient
ks = 0.0632;                %spring constant
ms = 0.16;                  %mass

% Transfer functions
Num = (cp*cs*s^2+(kp*cs+ks*cp)*s+kp*ks);
Den = (mp*ms)*s^4+(mp*cs+ms*cp+ms*cs)*s^3+(mp*ks+kp*ms+cp*cs+ks*ms)*s^2+(cs*kp+cp*ks)*s+(kp*ks);
T1 = Num/Den;

fr=0:0.001:10e2;
[mag,phase,wout]=bode(T1,fr);
omega_max1=wout(find(mag==max(mag))); %find resonant frequency

%% Initialization 1DOF

cp1 = 0.4;                   %damping coefficient
kp1 = 6.32;                  %spring constant
mp1 = 0.16;                  %mass
cc = 2*sqrt(mp1*kp1);         %critical damping coefficient
zeta = cp1/cc;               %damping ratio
omega = sqrt(kp1/mp1);        %natural frequency
omega_d = (1-zeta^2)*omega; %damped natural frequency

% Transfer function
T2 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));

fr=0:0.001:10e2;
[mag,phase,wout]=bode(T2,fr);
omega_max2=wout(find(mag==max(mag))); %find resonant frequency

%% Plots

figure(1);
bode(T1,fr)
hold on
figure(2);
step(T1)
% xlim([0 15])
hold on
figure(1);
bode(T2,fr)
hold on
legend('show')
figure(2);
step(T2)
legend('show')
hold on