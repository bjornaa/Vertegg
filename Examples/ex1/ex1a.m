% ex1a.m
%
% Script for producing PostScript file ex1a.ps
% Figure to example 1 in the VertEgg manual
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

% Space discretization
ve_init;
ve_grid(100, 0.5);

% Compute "eggsact" solution
A = eggsact(1000, 0.01, 0.001);

% Make plot
clf;
plot(A,ZE)
title('Stationary solution')
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')

print -deps ex1a.ps   % Make PostScript file

