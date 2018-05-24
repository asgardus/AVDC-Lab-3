%% This is a Matlab file for designing H_infinity controller (assignment 3
%% of SD2231)
clear all
clc
s=tf('s');

% systme parameters
m=22000;   %kg
j=700e3;   %kgm^2
c=40e3;    %Ns/m
k=600e3;%600e3; %N/m
L=6;       %m

%% State space model for skyhook contorl
Ask=[0 1 0 0
    -2*k/m 0 0 0
    0 0 0 1
    0 0 -2*k*L^2/j 0];
Bsk=[0 0 0 0
    k/m k/m -1/m -1/m
    0 0 0 0
    -L*k/j L*k/j L/j -L/j];
Csk=[0 1 0 0
    0 0 0 1];
Dsk=zeros(2,4);

%% H_inf using linmod syntax

%state space: The same as skyhook

      
%Weighting functions

%For penalizing actuator force
Wa1=(0.00175*s+1)/(0.00025*s+1);
Wa2=Wa1;

%For penalizing bounce and pitch motions
eps=1;
wnb=7.39;            %Find the right equation or value for wnb
wnchi=7.86;          %Find the right equation or value for wnchi
s1b=-eps+1i*sqrt(wnb^2-eps^2);
s2b=-eps-1i*sqrt(wnb^2-eps^2);
s1chi=-eps+1i*sqrt(wnchi^2-eps^2);
s2chi=-eps-1i*sqrt(wnchi^2-eps^2);
kb=5000;%input('Enter the gain for Wb = '); 
kchi=40000;%input('Enter the gain for Wchi = ');
Wb=(kb*s1b*s2b)/((s-s1b)*(s-s2b));
Wchi=(kchi*s1chi*s2chi)/((s-s1chi)*(s-s2chi));

%Extracting the extended model
[A_Pe,B_Pe,C_Pe,D_Pe] = linmod('Extended_model');% state space parameters of the extended system: Pe
Pe=ss(A_Pe,B_Pe,C_Pe,D_Pe);

%Calculating the controller
ncont = 2;%Number of control inputs
nmeas = 2;%Number of measured outputs provided to the controller
Pe=minreal(Pe);%This syntax cancels pole-zero pairs in transfer
%functions. The output system has minimal order and the same response
%characteristics as the original model.
[K,Pec,gamma,info]=hinfsyn(Pe,nmeas,ncont,'method','lmi'); % for working with the error
[Ainf, Binf, Cinf, Dinf]=ssdata(K);

%Now use the controller K in your simulation
%% initialization
src=1;
J=j;
f=8; %set frequency as 1Hz or 8Hz

k1=k;
k2=k;
c1=c;
c2=c;
l1=L;
l2=L;
%% Passive damped
A=[0 1 0 0;(-2*k/m) (-2*c/m) 0 0; 0 0 0 1;0 0 (-2*k*L^2)/J (-2*c*L^2)/J];
B=[0 0 0 0;k/m c/m k/m c/m;0 0 0 0;-k*L/J -c*L/J k*L/J c*L/J];
C=[1 0 0 0;0 0 1 0];
D=[0 0 0 0;0 0 0 0];
simout8=sim('task8ss')
% figure(1)%pitch
% plot(pitch_passive.time,pitch_passive.data,'LineWidth',1.5)
% hold on
% 
% figure(2)%excitation
% plot(input.time, input.data,'LineWidth',1.5)
% hold on
% figure(2)%bounce
% plot(z_passive.time, z_passive.data,'LineWidth',1.5)
% hold on
%% Skyhook (comment in when comparing)
% A1=[0 1 0 0;-(k1+k2)/m 0 (l1*k1-l2*k2)/m 0; 0 0 0 1;-(l1*k1-l2*k2)/J 0 -(k1*l1^2+k2*l2^2)/J 0]
% B1=[0 0 0 0;k1/m k2/m 1/m 1/m;0 0 0 0;-k1*l1/J k2*l2/J -l1/J l2/J]
% C1=eye(4);
% D1=zeros(4);
% 
% 
% c_z=190000;
% % cx=[1e6 2e6 3e6 4e6 5e6 6e6];
% c_x=5e6;
% simout9=sim('Task9ss')
% figure(1)%pitch skyhook
% plot(pitch.time, pitch.data,'LineWidth',1.5)
% hold on
% figure(2)%bounce skyhook
% plot(bounce.time, bounce.data,'LineWidth',1.5)
% hold on

% figure(2)
% f1(i)=max(Fa1.data);
% f2(i)=max(Fa2.data);
% figure(3)
% plot(Fa1.time, Fa1.data,'LineWidth',1.5)
% hold on
% plot(Fa2.time, Fa2.data,'LineWidth',1.5)
% hold on

%% H infinity

% A1=[0 1 0 0;-(k1+k2)/m 0 (l1*k1-l2*k2)/m 0; 0 0 0 1;-(l1*k1-l2*k2)/J 0 -(k1*l1^2+k2*l2^2)/J 0];
% B1=[0 0 0 0;k1/m k2/m 1/m 1/m;0 0 0 0;-k1*l1/J k2*l2/J -l1/J l2/J];
% C1=eye(4);
% D1=zeros(4);
%parameter variation (change in parameters)
mx=22000;   %kg
jx=700e3;   %kgm^2 
cx=40e3;    %Ns/m
kx=600e3; %N/m
L=6;       %m
% A1=Ask;
% B1=Bsk;
% C1=Csk;
% D1=Dsk;

A1=[0 1 0 0
    -2*kx/mx 0 0 0
    0 0 0 1
    0 0 -2*kx*L^2/jx 0];
B1=[0 0 0 0
    kx/mx kx/mx -1/mx -1/mx
    0 0 0 0
    -L*kx/jx L*kx/jx L/jx -L/jx];
C1=[0 1 0 0
    0 0 0 1];
D1=zeros(2,4);

%% plot()
simout10=sim('task10')
figure(1)
plot(pitch_inf.time,pitch_inf.data,'LineWidth',1.5)
hold on
grid on
% legend('passive','Hinf')%add 'skyhook b/w passive and Hinf when comparing with skyhook'
legend('original','Inc K','dec K')
xlim([0 5])
xlabel('Time(s)')
ylabel('Pitch angle(rad)')
hold on
figure(2)
plot(bounce_inf.time,bounce_inf.data,'LineWidth',1.5)
xlim([0 5])
xlabel('Time(s)')
ylabel('Bounce amplitude(m)')
% legend('Excitation','passive','Hinf')%add 'skyhook b/w passive and Hinf when comparing with skyhook'
legend('original','Inc J','dec J')% comment in for mass variation
hold on
grid on
figure(3)
plot(Fa1_inf.time, Fa1_inf.data,'LineWidth',1.5)
hold on
plot(Fa2_inf.time, Fa2_inf.data,'LineWidth',1.5)
hold on
legend('Fa1','Fa2','Inc J Fa1','Inc J Fa2','Dec J Fa1','Dec J Fa2')
% legend('Fa1','Fa2')
xlabel('Time(s)')
ylabel('Actuator Force (N)')
grid on
hold on
f1=max(Fa1_inf.data)
f2=max(Fa2_inf.data)
%% bode for weighting functions
fr=0:0.1:1e4
figure(7)
bode(Wa1,fr)
hold on

bode(Wb,fr)
hold on
bode(Wchi,fr)
hold on
legend('Wa1','Wb','Wchi')
