%
% plotres.m
% Plot results from modell run
%

res = fopen('result.dat', 'rt');

% Read info on first line
tmp = fscanf(res, '%*s %*s %f,',5);
depth   = tmp(1)
dz0     = tmp(2)
outstep = tmp(3)
simtime = tmp(4)
Ngrp    = tmp(5)
clear tmp
%
% Deduced variables
nout = simtime / outstep; % 

% Define grid
ve_init
ve_grid(depth,dz0)


A   = zeros(Ncell, Ngrp);
X   = zeros(Ncell, nout+1);
Mu  = zeros(nout+1, Ngrp);
Std = zeros(nout+1, Ngrp);
Max = zeros(nout+1, Ngrp);

tic;
for t = 1:nout+1   % +1 because initial distr.
  % Read the data for time t-1
  for k = 1:Ncell
    A(k,:) = fscanf(res,'%f',Ngrp)';
  end
  
  % X(:,t) is the sum of all egg groups at time t-1
  if (Ngrp == 1)
    X(:,t) = A(:,1);
  else
    X(:,t) = sum(A')';   % Row-sums
  end %if
  
  % Mu(t:grp) = mean depth, time t-1
  Mu(t,:) = ve_mean(A);
  Std(t,:) = ve_std(A);
  Max(t,:) = max(A);
  
end %for t
toc;

clf

disp(' ')
disp(' For control, plot the distributions at all times')

plot(X,ZE)
xlabel('Egg concentration [eggs/m^3]')
ylabel('Depth [m]')
title('Egg distribution at all times')

disp(' ')
disp('   Press any key to continue'); pause;
disp(' ')
disp('Plot final distributions for each group')

plot(A,ZE)
xlabel('Egg concentration [eggs/m^3]')
ylabel('Depth [m]')
title('Final distribution for each group')
legend('Grp 1', 'Grp 2', 'Grp 3', 'Grp 4', 'Grp 5')

print -deps ex6f.ps

disp(' ')
disp('   Press any key to continue'); pause;
disp(' ')
disp(' Plot the time evolution of the')
disp(' mean depth of each group')

TL = [0:nout]*outstep/24;    % Time labels


plot(TL, Mu)
axis([0 10 -Hcol 0]);
xlabel('Time [days]');
ylabel('Depth [m]');
title('Mean depth of egg group');
legend('Grp 1', 'Grp 2', 'Grp 3', 'Grp 4', 'Grp 5')

print -deps ex6g.ps

disp(' ')
disp('   Press any key to continue'); pause;
disp(' ')
disp('   Plot the time evolution of the standard deviation')
disp('   of depth within the groups')

plot(TL, Std)
xlabel('Time [days]');
ylabel('Standard deviation [m]');
title('Standard depth deviation in the groups');
legend('Grp 1', 'Grp 2', 'Grp 3', 'Grp 4', 'Grp 5')

disp(' ')
disp('   Press any key to continue'); pause;
disp(' ')
disp('   Plot the time evolution of the max concentration')

plot(TL, Max)
xlabel('Time [days]');
ylabel('Egg concentration [eggs/m^3]');
title('Maximum egg concentration for each group');
legend('Grp 1', 'Grp 2', 'Grp 3', 'Grp 4', 'Grp 5')

disp(' ')
disp('   The end')