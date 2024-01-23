% ex5b.m
%
% Script to produce PostScript film ex5b.ps
% Second figure to example 5 in the VertEgg manual
% Reproduces figure 1 in (Westgård 1989) + steady state
%
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

t12 = 12 / outstep;

% Stationary solution
Y1 = sstate(M1,K,W1);
Y2 = sstate(M2,K,W2);
Y = Y1 + Y2;


clf
hold on;
plot(A0, ZE, ':');         % The start distribution   
plot(X(:,t12), ZE,'--');   % After 12 hours
plot(Y,ZE);                % Steady state
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')
legend('Initial distribution', 'Distribution after 12 hours', ...
       'Stationary distribution')
hold off;


print -deps ex5b.ps



