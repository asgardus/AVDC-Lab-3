clc
clear all


%% Initialization

s = tf('s');
damping=1;
switch damping
    case 1  %for underdamped
        cp=0.4;
    case 2  %for critically damped
        cp=2.0112
    case 3  %for overdamped
        cp=8;
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


%% Sine Excitation
f=0.5;
dt=0.02;
Tfinal=10;
t = 0:dt:Tfinal
u=0.05*sin(2*pi*f*t);

figure (1)
lsim(G1,u,t);        %sine input
hold on;
plot(t,u)
hold on
ylabel('Amplitude (m)')
set(findall(gcf,'type','line'), 'LineWidth', 1.5);
grid on
legend('Response','Excitation')
%% Impulse Excitation


 stepinput = 0.05*ones([1 length(t)]);
for i = 51:length(t)
   stepinput(i) = 0; 
end

figure(4)
lsim(G1,stepinput,t) % impulse excitation
hold on
plot(t,stepinput)
hold on;
set(findall(gcf,'type','line'), 'LineWidth', 1.5);
grid on
ylabel('Amplitude(m)')
% legend('Response','Excitation')
xlim([0 6])
%% PSD
w=0:25;
G3=freqresp(G1,w)
s=(4.028*10^(-7))./(2.88*10^(-4)+0.68*w.^2+w.^4)
G4=abs(G3(:))';
PSD=s.*(G4.^2);

figure (3)
semilogy(w,PSD)
hold on;
set(findall(gcf,'type','line'), 'LineWidth', 1.5);
grid on
xlabel('Frequency (rad/s)')
ylabel('PSD (W/Hz)')
%% Plots
% 
% figure(1);
% bode(G,fr)
% hold on
% figure(2);
% step(G)
% xlim([0 5])
% hold on