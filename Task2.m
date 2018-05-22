clc
close all
clear all
%% Initialization (PD)
src=2;                      %define source (1=Sine,2=Impulse, 3=PSD)
s = tf('s');
cp = 0.4;                   %damping coefficient
hd_mat = 1.2;               %derivative gain
hp_mat = 0;                 %proportional gain
hi_mat= 0;                  %integral gain
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cc = 2*sqrt(mp*kp);         %critical damping coefficient
zeta = cp/cc;               %damping ratio
omega = sqrt(kp/mp);        %natural frequency
omega_d = (1-zeta^2)*omega; %damped natural frequency

%% PD
% % Transfer function
% T1 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));
% hi=hi_mat;
% for i=1:length(hd_mat)
%     hd = hd_mat(i);
%     for j=1:length(hp_mat)
%             hp = hp_mat(j);
% 
%             % Transfer function
%             T2 = kp/(mp*s^2+hd*s+kp+hp); 
% 
%             fr=0:0.001:10e2;
%             [mag,phase,wout]=bode(T1,fr);
%             omega_max=wout(find(mag==max(mag)));    %find resonant frequency
% 
%             [mag,phase,wout]=bode(T2,fr);
%             omega_max2=wout(find(mag==max(mag)));   %find resonant frequency          
% 
%         figure(1);
%         bode(T2,fr)
%         hold on
%         figure(2);
%         step(T2)        
%         hold on
%         sim('Main')
%     end
% end
% 
% % figure(1);
% % bode(T1,fr)
% % legend('show')
% % figure(2);
% % step(T1)
% % xlim([0 5])
% % legend('show')
%% PID
a=[20];
% z=[0.8];
% hd=mp*(2*w*z+a);
% hi=a*w^2;
% hp=w^2+2waz-kp;
% Transfer function
T1 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));
w=omega_d;
z=0.8;

for i=1:length(a)
   
    for j=1:length(z)
           hd=mp*(2*w*z+a);
           hi=a*w^2;
           hp=w^2+2*w*a*z-kp;

            % Transfer function
            T2 = kp/(mp*s^2+hd*s+kp+hp); 

            fr=0:0.001:10e2;
            [mag,phase,wout]=bode(T1,fr);
            omega_max=wout(find(mag==max(mag)));    %find resonant frequency

            [mag,phase,wout]=bode(T2,fr);
            omega_max2=wout(find(mag==max(mag)));   %find resonant frequency          

%         figure(1);
%         bode(T2,fr)
%         hold on
%         figure(2);
%         step(T2)        
%         hold on
        sim('Main')
    end
end
