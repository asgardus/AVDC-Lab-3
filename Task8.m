clear all
close all
clc

m=22000;
J=700000;
c=40000;
k=600000;
L=6;
% state space of the system
A=[0 1 0 0;(-2*k/m) (-2*c/m) 0 0; 0 0 0 1;0 0 (-2*k*L^2)/J (-2*c*L^2)/J]
B=[0 0 0 0;k/m c/m k/m c/m;0 0 0 0;-k*L/J -c*L/J k*L/J c*L/J]
C=[1 0 0 0;0 0 1 0];
D=[0 0 0 0;0 0 0 0];
sys=ss(A,B,C,D)
% for natural freq: c=0
An=[0 0 1 0;(-2*k/m) 0 0 0; 0 0 0 1;0 0 (-2*k*L^2)/J 0];
Bn=[0 0 0 0;k/m 0 k/m 0;0 0 0 0;-k*L/J 0 k*L/J 0];
Cn=[1 0 0 0;0 0 1 0];
Dn=0;
src=1;%1 for impulse,2 for sine
f=8; %set frequency in Hz (1 or 8)
