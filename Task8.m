clear all
close all
clc

m=22000;
J=700000;
c=40000;
k=600000;
L=6;
k1=k;
k2=k;
c1=c;
c2=c;
l1=L;
l2=L;
% state space of the system
A=[0 1 0 0;(-2*k/m) (-2*c/m) 0 0; 0 0 0 1;0 0 (-2*k*L^2)/J (-2*c*L^2)/J]
B=[0 0 0 0;k/m c/m k/m c/m;0 0 0 0;-k*L/J -c*L/J k*L/J c*L/J]
C=[1 0 0 0;0 0 1 0];
D=[0 0 0 0;0 0 0 0];
sys=ss(A,B,C,D)
% for natural freq: c=0
An=[0 1 0 0;(-2*k/m) 0 0 0; 0 0 0 1;0 0 (-2*k*L^2)/J 0];
Bn=[0 0 0 0;k/m 0 k/m 0;0 0 0 0;-k*L/J 0 k*L/J 0];
Cn=[1 0 0 0;0 0 0 0;0 0 1 0;0 0 0 0];
Dn=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];

src=1;%1 for impulse,2 for sine
f=8;%set frequency in Hz (1 or 8)

sys1=ss(An,Bn,Cn,Dn)
G=ss2tf(An,Bn,Cn,Dn,1)
[num,den]=ss2tf(An,Bn,Cn,Dn,3);
bounce=tf(num(1,:),den);
pitch=tf(num(3,:),den);
% bode(bounce,pitch);
simout8=sim('task8ss')
%% Task 9.1
A1=[0 1 0 0;-(k1+k2)/m 0 (l1*k1-l2*k2)/m 0; 0 0 0 1;-(l1*k1-l2*k2)/J 0 -(k1*l1^2+k2*l2^2)/J 0]
B1=[0 0 0 0;k1/m k2/m 1/m 1/m;0 0 0 0;-k1*l1/J k2*l2/J -l1/J l2/J]
C1=eye(4);
D1=zeros(4);

time=1:10
c_z=190000;
% cx=[1e6 2e6 3e6 4e6 5e6 6e6];
cx=5e6;
figure(1)
plot(pitch_passive.time,pitch_passive.data,'LineWidth',1.5)
hold on
legend('pitch_passive')
figure(2)
plot(input.time, input.data,'LineWidth',1.5)
hold on
figure(2)
plot(z_passive.time, z_passive.data,'LineWidth',1.5)
hold on
for i=1:length(cx)
    c_x=cx(i);
    simout9=sim('Task9ss')
    figure(1)
plot(pitch.time, pitch.data,'LineWidth',1.5)
hold on
figure(2)
plot(bounce.time, bounce.data,'LineWidth',1.5)
hold on
disp(Fa1)
disp(Fa2)
% figure(2)
f1(i)=max(Fa1.data);
f2(i)=max(Fa2.data);
figure(3)
plot(Fa1.time, Fa1.data,'LineWidth',1.5)
hold on
plot(Fa2.time, Fa2.data,'LineWidth',1.5)
hold on
end
figure(1)%pitch graph
xlim([0 5])
grid on
legend('Passive','Active')
xlabel('Time(s)')
ylabel('Pitch angle(rad)')
figure(2)
xlim([0 5])
xlabel('Time(s)')
ylabel('Bounce amplitude(m)')
grid on
legend('Excitation','Passive','Active')
figure(3)
xlim([0 5])
legend('Fa1','Fa2')
ylabel('Actuator Force(N)')
xlabel('Time(s)')
grid on
