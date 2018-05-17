clc
close all
clear all
%% Initialization and Model
src=1;
s = tf('s');
cp = 0.4;                   %damping coefficient
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cc = 2*sqrt(mp*kp);         %critical damping coefficient
zeta = cp/cc;               %damping ratio
omega = sqrt(kp/mp);        %natural frequency
omega_d = (1-zeta^2)*omega; %damped natural frequency
z_pin = 0;                  %starting position of the mass
hi = 1;                   %integrative gain
hp = 1;                    %proportional gain
hd = 1;                     %derivative gain
hd_mat = hd;                %derivative gain matrix
hp_mat = hp;                %proportional gain matrix
hi_mat = hi;                %integrative gain matrix
src = 3;                    % 1 - sine, 2 - impulse, 3 - step
simout = sim('extratask2'); %simulate model

%% Transfer function

T1 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s)); %simple mass spring damper system

for i=1:length(hd_mat)
    hd = hd_mat(i);
    for j=1:length(hp_mat)
        hp = hp_mat(j);
        for k=1:length(hi_mat)
            hi = hi_mat(k);

%           Transfer function
            T2 = kp*s/(mp*s^3+hd*s^2+(kp+hp)*s+hi);            %actuated system 

            fr=0:0.001:10e2;
            [mag,phase,wout]=bode(T1,fr);
            omega_maxT1=wout(find(mag==max(mag)));      %find resonant frequency

            [mag,phase,wout]=bode(T2,fr);
            omega_maxT2=wout(find(mag==max(mag)));      %find resonant frequency          

        figure(1);
        bode(T2,fr)
        hold on
        figure(2);
        step(T2)        
        hold on
        end
    end
end

% figure(1);
% bode(T1,fr)
% legend('show')
% figure(2);
% step(T1)
xlim([0 5])
% legend('show')

clc