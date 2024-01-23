% ex3a
%
% Script for making PostScript file ex3a.ps
%
% Figure 1 to example 3 in the VertEgg Manual
% Note: runsens must have been run first. 
% Note: To reproduce the figure, the parameters in 
%       runsens must have the following values
%       K = 0.01, W = 0.001, dz = 2 and dt = 120.
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

cla
plot(X(:,1),ZE,'-.',  X(:,2), ZE, ':', ...
     X(:,3),ZE,'--',  X(:,5), ZE, '-')
legend('ftcs', 'lwendrof', 'upstream', 'minlim')
title('Computational stationary solution')
xlabel('Egg concentration   [eggs/m^3]');
ylabel('Depth   [m]');

print -deps ex3a.ps