% ex4b
%
% Script for making PostScript file ex3b.ps
%
% Figure 3 to example 3 in the VertEgg Manual
% Note: runsens must have been run first. 
% Note: To reproduce the figure, the parameters in 
%       runsens must have the following values
%       K = 0.01, W = 0.001, dz = 2 and dt = 120.
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

clf
plot(Y(:,1),ZE,'-.',  Y(:,2), ZE, ':', ...
     Y(:,3),ZE,'--',  Y(:,5), ZE, '-')
legend('ftcs', 'lwendrof', 'upstream', 'minlim')
title('Computational stationary solution')
xlabel('Egg concentration   [eggs/m^3]');
ylabel('Depth   [m]');


print -deps ex3b.ps