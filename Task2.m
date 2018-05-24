clc
close all
clear all

%% Initialization

s = tf('s');
cp = 0.4;                   %damping coefficient
dd_mat = 1.2;       %derivative gain
dp_mat = 0;                 %proportional gain
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cc = 2*sqrt(mp*kp);         %critical damping coefficient
zeta = cp/cc;               %damping ratio
omega = sqrt(kp/mp);        %natural frequency
omega_d = (1-zeta^2)*omega; %damped natural frequency

% Transfer function(damped passive)
T1 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));

for i=1:length(dd_mat)
    dd = dd_mat(i);
    for j=1:length(dp_mat)
            dp = dp_mat(j);

            % Transfer function
            T2 = kp/(mp*s^2+dd*s+kp+dp); 

            fr=0:0.001:10e2;
            [mag,phase,wout]=bode(T1,fr);
            omega_max=wout(find(mag==max(mag)));    %find resonant frequency

            [mag,phase,wout]=bode(T2,fr);
            omega_max2=wout(find(mag==max(mag)));   %find resonant frequency          

        figure(1);
        bode(T2,fr)
        hold on
        figure(2);
        step(T2)        
        hold on
    end
end

figure(1);
bode(T1,fr)
legend('show')
figure(2);
step(T1)
xlim([0 5])
legend('show')