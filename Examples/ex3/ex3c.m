% ex3c
%
% Script for making PostScript file ex3c.ps
%
% Figure 3 to example 3 in the VertEgg Manual
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

% Space discretization
ve_init
ve_grid(100,2)

% Plot only part of the depth range
II = 41:Ncell;

% Compute the eggsact solution
B = eggsact(1000, 0.0005, -0.001);

% make the plot
clf
plot(B(II), ZE(II))
title('Eggsact solution for sinking eggs')
xlabel('Egg concentration   [eggs/m^3]');
ylabel('Depth   [m]');


print -deps ex3c.ps