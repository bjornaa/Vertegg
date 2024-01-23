%
% plotgrp3.m
% Plot results from group 3 (of 5) from modell run
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
X3   = zeros(Ncell, nout/2+1);

tic;
save = 1;  % Switch for saving
for t = 1:nout+1   % +1 because initial distr.
  % Read all data for time t-1
  for k = 1:Ncell
    A(k,:) = fscanf(res,'%f',Ngrp)';
  end
  
  % save the third group in X3, for even t
  if (save == 1)
    X3(:,t/2) = A(:,3); 
  end % if
  save = 1 - save;
  
end %for t
toc;

disp(' ')
disp(' Time evolution ')

plot(X3,ZE)
xlabel('Egg concentration [eggs/m^3]')
ylabel('Depth [m]')
title('Daily egg distributions')

print -deps ex6d.ps

disp(' ')
disp('   Press any key to continue'); pause;
disp(' ')
disp('Plot final distribution against sstate')

M = 25.0;  % Vertical integrated consentration [eggs/m^2]
SS3 = sstate(M, K, W(:,3));
plot(X3(:,10),ZE)
xlabel('Egg concentration [eggs/m^3]')
ylabel('Depth [m]')
title('Final distribution against steady state solution')
axis([0 5 -150 -100])
hold on
plot(SS3,ZE,'r:')
%plot(SS3,ZE,'r')
legend('posmet','sstate')
hold off

print -deps ex6e.ps

disp(' ')
disp('   The end')