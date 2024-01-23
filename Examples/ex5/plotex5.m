% plotex5.m
%
% Script for plotting the results in Ex. 5
%
% runex5 must have been run first
%
% plotex5 is part of the example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 14 August 1995

disp(' ');
disp(' plotex5 -- Plotting results in Example 5')

disp(' ')
disp('  Reproduce figure 1 in Westgård (1989)')

% Results after 12 hours

t12 = 12 / outstep;
clf;
hold on;
plot(A0, ZE, 'y');       % The start distribution   
plot(X(:,t12), ZE, 'g');   % After 12 hours
title('Figure 1 in Westgård (1989)')
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')
legend('Initial distribution', 'Distribution after 12 hours')

disp(' ')
disp('Press any key to continue'); pause
disp(' ')
disp(' Add in blue, the 2 egg groups separately')

plot([X1(:,t12) X2(:,t12)], ZE, 'b')
hold off;

disp(' ')
disp('Press any key to continue'); pause
clc

disp('Plot results from all time steps')

clf;
hold on;
plot(A0 , ZE, 'y');
plot(X, ZE, 'g');
title('Modell solutions at al times)')
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')

disp(' ')
disp('Press any key to continue'); pause
disp(' ')
disp(' Add the stationary solution in red')

% Compute stationary solution
Y1 = sstate(M1,K,W1);
Y2 = sstate(M2,K,W2);
Y = Y1 + Y2;

plot(Y, ZE, 'r');
hold off

disp(' ')
disp('Press any key to continue'); pause

disp(' ')
disp(' A nicer view of the time evolution')

clf;
T = [0:outstep:simtime];    % Time array in hours
surf(T,ZE,[A0 X]);
title('Time evolution of solution')
xlabel('Time [hours]')
ylabel('Depth [m]')
zlabel('Egg concentration [eggs/m^3]')

disp(' ')
disp('Press any key to continue');pause

disp(' ')
disp(' Look at the convergence to the stationary solution by plotting')
disp(' the time evolution of the root mean square deviation from the')
disp(' steady state solution')

RMS = ve_rmsd([A0 X], Y);

clf;
semilogy(T,RMS);
title('Root mean square deviation from steady state')
xlabel('Time')
ylabel('RMS deviation [eggs/m^3]')
pause;  % Press any key to contonue

disp(' ')
disp('The end')


