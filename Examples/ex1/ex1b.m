% ex1b.m
%
% Script for producing PostScript file ex1a.ps
% Second figure to example 1 in the VertEgg manual
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 1 October 1996

% Set up grid
H = 100;
dz0 = 2.0;
ve_init;
ve_grid(H,dz0);

% Compute "srcsact" solution
K = 0.01;
W = 0.001;
P = 1/(24*3600);
alpha = - log(0.01) / 3600;
A = srcsact(K,W,P,alpha);

clf;
plot(A,ZE)
title('Stationary solution with source term')
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')

print -deps ex1b.ps   % Make PostScript file
