

%% Initialization
s = tf('s');

cp =0.4;                   %damping coefficient
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cc = 2*sqrt(mp*kp);         %critical damping coefficient
zeta = cp/cc;               %damping ratio
omega = sqrt(kp/mp);        %natural frequency
omega_d = (1-zeta^2)*omega; %damped natural frequency


G1 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));
G2 = tf([cp kp],[mp cp kp])
%% PD 
d_d=1.2;
d_p=0;
H1=(kp)/(mp*s^2+d_d*s+kp+d_p);
%% PID
w=omega;
z=0.5;
a=15;
hd=mp*(2*w*z+a);
hi=a*w^2;
hp=w^2+2*w*a*z-kp;

%Transfer function PID
K1 = (kp*s)/(mp*s^3+hd*s^2+s*(kp+hp)+hi); 
%% Skyhook
T=1.5;
%Transfer function skyhook
L1=(kp)/(mp*s^2+T*s+kp);

%% Bode (Task 5.1)
%% Sine Excitation
f=0.5;
dt=0.02;
Tfinal=10;
t = 0:dt:Tfinal
u=0.05*sin(2*pi*f*t);

figure (1)
plot(t,u)
hold on
lsim(G1,u,t);        %sine input
hold on;

lsim(H1,u,t);
hold on
lsim(K1,u,t)
hold on
lsim(L1,u,t)
hold on
ylabel('Amplitude (m)')
set(findall(gcf,'type','line'), 'LineWidth', 1.5);
grid on
legend('Excitation','Passive','PD','PID','Skyhook')
%% Impulse Excitation


 stepinput = 0.05*ones([1 length(t)]);
for i = 51:length(t)
   stepinput(i) = 0; 
end

figure(4)
plot(t,stepinput)
hold on
lsim(G1,stepinput,t) % impulse excitation
hold on
lsim(H1,stepinput,t)
hold on
lsim(K1,stepinput,t)
hold on;
lsim(L1,stepinput,t)
hold on
set(findall(gcf,'type','line'), 'LineWidth', 1.5);
grid on
ylabel('Amplitude(m)')
% legend('Response','Excitation')
xlim([0 6])
legend('Excitation','Passive','PD','PID','Skyhook')
%% PSD
w=0:0.1:25;
G3=freqresp(G1,w)
s=(4.028*10^(-7))./(2.88*10^(-4)+0.68*w.^2+w.^4)
G4=abs(G3(:))';
PSD=s.*(G4.^2);

H2=freqresp(H1,w)
H3=abs(H2(:))';
PSDpd=s.*(H3.^2);

K2=freqresp(K1,w)
K3=abs(K2(:))';
PSDpid=s.*(K3.^2);


L2=freqresp(L1,w)
L3=abs(L2(:))';
PSDsky=s.*(L3.^2);

figure (3)
semilogy(w,PSD)
hold on;
semilogy(w,PSDpd)
hold on;
semilogy(w,PSDpid)
hold on;
semilogy(w,PSDsky)
hold on;
set(findall(gcf,'type','line'), 'LineWidth', 1.5);
grid on
xlabel('Frequency (rad/s)')
ylabel('PSD (W/Hz)')
legend('Passive','PD','PID','Skyhook')