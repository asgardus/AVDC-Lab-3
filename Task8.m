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
%% Task 9.1
A1=[0 1 0 0;-(k1+k2)/m 0 (l1*k1-l2*k2)/m 0; 0 0 0 1;-(l1*k1-l2*k2)/J 0 -(k1*l1^2+k2*l2^2)/J 0]
B1=[0 0 0 0;k1/m k2/m 1/m 1/m;0 0 0 0;-k1*l1/J k2*l2/J -l1/J l2/J]
C1=eye(4);
D1=zeros(4);
src2=1;

c_z=230000;
c_x=12e6;
