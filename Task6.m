clear;
clc

%% Initialization 2DOF
f=2;                      %frequency (rad/s)
dt=0.02;
Tfinal=10;
t = 0:dt:Tfinal;
s = tf('s');
cp = 0.8;                   %damping coefficient
kp = 6.32;                  %spring constant
mp = 0.16;                  %mass
cs = 0.05;                  %damping coefficient
ks = 0.0632;                %spring constant
ms = 0.16;                  %mass
T_mat = 0.15;
for i=1:length(T_mat)
    T = T_mat(i);                      %derivative gain
    src = 3;                    % 1 - impulse, 2 - step, 3 - sine

    % State-Space matrices
    A = [0,1,0,0;-ks/ms,-T/ms,ks/ms,0;0,0,0,1;ks/mp,T/mp,-(ks+kp)/mp,-cp/mp];
    B = [0,0;0,0;0,0;kp/mp,cp/mp];
    o = [1 1 1 1];
    C = diag(o);
    D = zeros(4,2);
    simout = sim('task6ss');
    track(i) = T;
end

% State space to TF Derivative controlled
[b,a]=ss2tf(A,B,C,D,1);
skyhook = tf(b(1,:),a);

% Transfer functions
Num = (cp*cs*s^2+(kp*cs+ks*cp)*s+kp*ks);
Den = (mp*ms)*s^4+(mp*cs+ms*cp+ms*cs)*s^3+(mp*ks+kp*ms+cp*cs+ks*ms)*s^2+(cs*kp+cp*ks)*s+(kp*ks);
T2 = Num/Den;

fr=0:0.001:10e2;
[mag,phase,wout]=bode(T2,fr);
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
T1 = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));

fr=0:0.001:10e2;
[mag,phase,wout]=bode(T1,fr);
omega_max2=wout(find(mag==max(mag))); %find resonant frequency

%% Sine excitation
u=0.05*sin(f*t);
figure(3);
lsim(T2,u,t);
hold on
plot(z_s,'r');
legend('show')


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
bode(skyhook,fr)
legend('show')
figure(2);
step(T2)
hold on
step(skyhook)
legend('show')
hold on