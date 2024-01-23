% ex3d
%
% Script for making PostScript file ex3d.ps
%
% Figure 4 to example 3 in the VertEgg Manual
% Note: runsens must have been run first. 
% Note: To reproduce the figure, the parameters in 
%       runsens must have the following values
%       K = 0.0005, W = -0.001, dz = 2 and dt = 120.
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

II = 41:Ncell;
clf
plot(X(II,1),ZE(II),'-.',  X(II,2), ZE(II), ':', ...
     X(II,3),ZE(II),'--',  X(II,5), ZE(II), '-')
legend('ftcs', 'lwendrof', 'upstream', 'minlim')
title('Computational stationary solution')
xlabel('Egg concentration   [eggs/m^3]');
ylabel('Depth   [m]');

print -deps ex3d.ps