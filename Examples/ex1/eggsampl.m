% eggsampl
%
% Example script for testing the "eggsact" solution
%
% eggsampl.m is part of the example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 13 August 1995

clc;

disp(' ')
disp(' Example 1,  Stationary solution with constant coefficients')

M = 1000;
K = 0.01;
W = 0.001;

fprintf(1,' M = %f eggs/m^3\n', M)
fprintf(1,' K = %f m^2/s\n', K)
fprintf(1,' W = %f m/s\n', W)

H = 100;
dz0 = 0.5;
ve_init;
ve_grid(H,dz0);

disp(' ');
disp('  Press any key to continue'); pause
clc;

disp(' ')
disp('Compute concentration values at specific depths:')
disp(' ')
fprintf(1,' At surface, eggsact(M,K,W,0)    = %f\n', eggsact(M,K,W,0))
fprintf(1,' At 10 m,    eggsact(M,K,W,10)   = %f\n', eggsact(M,K,W,10))
fprintf(1,' At bottom,  eggsact(M,K,W,Hcol) = %f\n', eggsact(M,K,W,-Hcol))
disp(' ')
disp('Both positive and negative depth values can be used:')
fprintf(1,' At 10 m,    eggsact(M,K,W,-10)  = %f\n', eggsact(M,K,W,-10))

disp(' ')
disp('  Press any key to continue'); pause
clc;

disp(' ')
disp('Compute the cell averages')
disp('   A = eggsact(M,K,W)')

A = eggsact(M,K,W);

disp(' ')
disp('Properties of A:')
fprintf(1,'Vertical integral,  ve_int(A)  = %f\n', ve_int(A))
fprintf(1,'Mean,               ve_mean(A) = %f\n', ve_mean(A))
fprintf(1,'Standard deviation, ve_std(A)  = %f\n', ve_std(A))

% Compute the exact values, formulas ...
m = W/K;
mean = - (1 / m) + H / (exp(m*Hcol) - 1);
var  = (2 - exp(-m*Hcol)*(m^2*Hcol^2 + 2*m*Hcol + 2)) / ...
                (m^2*(1-exp(-m*Hcol))) - mean^2;
std = sqrt(var);

disp(' ')
disp('The exact values can be computed')
fprintf(1,'Integral           = %f\n', M);
fprintf(1,'Mean               = %f\n', mean);
fprintf(1,'Standard deviation = %f\n', std);

disp(' ')
disp('The errors are due to the space discretisation')
disp(' ')
disp('Press any key to continue'); pause
clc

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
disp('  Eggsampl is finished')
disp('  Run "srcsampl" to see the effect of source terms')
disp('  Press any key to quit'); pause

clc
