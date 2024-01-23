% ex5a.m
%
% Script to produce PostScript file ex5b.ps
% First figure to example 5 in the VertEgg manual
% Reproduces figure 1 in (Westgård 1989).
%
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 10 October 1995

t12 = 12 / outstep;

clf;
hold on;
plot(A0, ZE, ':');         % The start distribution   
plot(X(:,t12), ZE);   % After 12 hours
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')
legend('Initial distribution', 'Distribution after 12 hours')
hold off;


print -deps ex5a.ps



