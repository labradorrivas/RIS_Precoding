clc
clear all

Nm=4;      % Number of Antennas per UT
B=8;       % Number of Antennas at BS
Nr=8;      % Number of RIS elements
UT=2;       % Number of Users

syms IB [B B] 
syms HB [B Nr] 
syms Phi [Nr Nr] 
syms HA1 [Nr Nm] 
syms Q1 [Nm Nm] 

SE=log(det(IB+(HB*Phi*HA1)*Q1*(HB*Phi*HA1)'));
syms sigma