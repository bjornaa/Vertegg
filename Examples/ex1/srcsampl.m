% srcsampl
%
% Example script for testing the "srcsact" solution
%
% srcsampl.m is part of the ex1 example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 01 October 1996

clc;

disp(' ')
disp(' Example 1.1,  Stationary solution with constant coefficients')
disp('    including source terms')
disp(' ')

% Diffusion koefficient K = 0.01 m^2/s correspond to windy conditions
%   in the mixed layer.
K = 0.01;
% W = 0.001 m/s is a typical ascending velocity for cod eggs.
W = 0.001;
% Constant spawning rate, 1 egg/day/m
P = 1/(24*3600);
% Mortality 1%/hour
%   exp(-alpha*3600) = 0.01
alpha = - log(0.01) / 3600;

fprintf(1,' K     = %f m^2/s\n', K)
fprintf(1,' W     = %f m/s\n', W)
fprintf(1,' P     = %f eggs/s/m\n', P)
fprintf(1,' alpha = %f \n', alpha)

H = 100;
dz0 = 2.0;
ve_init;
ve_grid(H,dz0);

disp(' ');
disp('  Press any key to continue'); pause
clc;
disp(' ')
disp('Compute the cell averages')
disp('   A = srcsact(K,W,P,alpha)')

A = srcsact(K,W,P,alpha);

disp(' ')
disp('Properties of A:')
fprintf(1,'Vertical integral,  ve_int(A)  = %f\n', ve_int(A))
fprintf(1,'Mean,               ve_mean(A) = %f\n', ve_mean(A))
fprintf(1,'Standard deviation, ve_std(A)  = %f\n', ve_std(A))
disp(' ')
disp(' Equation 1.x gives the vertical integral as Ptot/alpha')
Ptot = P*H;
fprintf(1,'Ptot = P*H = %f\n', Ptot)
fprintf(1,'Thuss the vertical integral = %f\n', Ptot/alpha)

disp(' ')
disp('Make a plot of the solution')
disp(' ')
disp('Command: plot(A,ZE)')

clf;
plot(A,ZE)
title('Stationary solution')
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')

disp(' ')
disp('  srcsampl is finished')
