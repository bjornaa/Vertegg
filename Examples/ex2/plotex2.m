% plotex2
%
% Script for plotting results from runlw.
%
% runlw must have been run first.
%
% plotex2 is part of the example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 13 August 1995


clc

disp(' ');
disp(' plotex2 -- Plotting results from runlw')
disp(' ')
disp(' The last result from the simulation')
disp(' can easily be plotted')

clf
plot(A,ZE)

disp(' ')
disp('Press a key, to compare with exact solution'); pause

B = eggsact(M0,K0,W0);
hold on
plot(B,ZE,'r')

disp(' ')
disp('Press a key, to make the figure more presentable'); pause

title('Egg distributions')
xlabel('Egg concentration  [eggs/m^3]')
ylabel('Depth [m]')
legend('Numerical solution', 'Exact solution')

disp(' '),
disp('Press a key to continue'); pause

clc
clf

disp(' ');
disp(' The time evolution of the mean');

T = [1:nout];     
TL = T*outstep;       % Time labels
mu = ve_mean(X);      % The means

plot(TL,mu)
title('Time evolution of center of gravity')
xlabel('Time [hours]')
ylabel('Mean [m]')

disp(' ')
disp('Press any key to continue'); pause
disp(' ')
disp(' A better way to look at the convergence of the solution')
disp(' is to plot the time evolution of the root mean square')
disp(' deviation from the exact steady state solution')
disp(' ')
disp(' A logarithmic Y-axis is used')

R = ve_rmsd(X,B);

clf
semilogy(TL,R)
title('Root mean square deviation from exact solution')
xlabel('Time [hours]')
ylabel('RMS deviation [eggs/m^3]')


disp(' ')
disp('Press any key to continue'); pause

clc;
disp(' ')
disp(' An isoplet diagram show the time evolution of the solution')

clf;
cl = [10:10:110];
c = contour(TL,ZE,X,cl);
clabel(c)
title('Time evolution of vertical egg concentratuion')
xlabel('Time [hours]');
ylabel('Depth [m]');
zlabel('Egg concentration [eggs/m^3]');

disp(' ')
disp('Press a key to continue'); pause
disp(' ')
disp(' The same surface viewed in 3D')

clf;

surfl(TL,ZE,X)
%shading interp
title('Time evolution of vertical egg concentration')
xlabel('Time [hours]');
ylabel('Depth [m]');
zlabel('Egg concentration [eggs/m^3]');

disp(' ')
disp('The end')




