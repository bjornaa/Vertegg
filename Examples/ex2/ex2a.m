% ex2a
%
% Script for making PostScript file ex2a.ps
%
% Makes first figure in example 2 in the VertEgg manual
% Note: One of the run-scripts must have been run first
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

% Compute "eggsact" stationary solution
B = eggsact(M0,K0,W0);
% Compute Root Mean Square Deviation
R = ve_rmsd(X,B);

% Time labels
T = [1:nout];
TL = T*outstep;

clf
semilogy(TL,R)
title('Root mean square deviation from exact solution')
xlabel('Time [hours]')
ylabel('RMS deviation [eggs/m^3]')

print -deps ex2a.ps   % Make PostScript file




