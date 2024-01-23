% ex8a.m
%
% Script for producing postscript file ex8a.ps
% Figure ... in the VertEgg report
%
% Requires that runex8 has already been runned
% in mode = 1.

B = X(:,20);  % Distribution after 20 days

clf

h = plot(B,ZE);
axis([0 12 -600 0])
box on
set(h, 'LineWidth', 1.2)
ylabel('Depth [m]')
xlabel('Egg concentration  [Eggs/m^3]')

print -deps ex8a.ps
