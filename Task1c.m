clc
s = tf('s');
cp = 0.4; 
kp = 6.32;
mp = 0.16;
cc = 2*sqrt(mp*kp);
zeta = cp/cc;
omega = sqrt(kp/mp);
T = (omega^2+(2*zeta*omega*s))/(s^2+omega^2+(2*zeta*omega*s));
omega_d = (1-zeta^2)*omega;
[mag,phase,wout]=bode(T,fr);
omega_max=wout(find(mag==max(mag)));
figure(1);
fr=0:0.001:10e4;
bode(T,fr);
hold on
figure(2);
step(T);
xlim([0 5])
hold on