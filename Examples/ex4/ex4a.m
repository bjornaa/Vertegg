% ex4a.m
%
% Script to produce ex4a.ps
% First figure to example 4 in the VertEgg manual
% Reproduction of Fig. 1 in Sundby (1983)
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 22 September 1995

clf;

mu = 1.6e-3;     % Dynamic molecular viscosity
rho = 1027;      % Density of sea water

d = [0:0.1:5];   % Range of diameters (in mm)

% Plot the velocity curves
hold on;
axis([0 5 0 5]);
title('Terminal velocity at different density differences')
xlabel('Diameter [mm]')
ylabel('Velocity [mm/s]')
for drho = [0.25 0.5 1:6]   % The values of drho used by Sundby
  % Calling eggvel, drho must be an array of size = size(d)
  % and d must be given in m.
  W = eggvel(drho * ones(size(d)), d/1000);
  plot(d, W*1000);
end

% Add Reynolds number curves
d(1) = [];   % Remove d(1) to prevent divison by zero;
W = 0.5 * mu ./ (rho * d/1000);
W = 1000 * W;       % Convert to mm/s
plot(d, W, ':');    % Re = 0.5
plot(d, 10*W, ':'); % Re = 5.0

print -deps ex4a.ps


