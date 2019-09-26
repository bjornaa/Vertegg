% mdepth.m
%
% Script to produce mdepth.ps
% First figure in theory chapter in the VertEgg manual
% Shows dependence of mean depth mu on m = w/K.
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 22 September 1995

clf;

H = 100;
m = [-1:0.01:-0.01 0.01:0.01:1];
mu = -(1 ./ m) + H ./ (exp(m*H) - 1);

plot(m,mu)

xlabel('m = w/K [1/m]');
ylabel('Depth [m]');

print mdepth.ps