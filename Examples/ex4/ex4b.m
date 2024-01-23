% ex4b.m
%
% Script to produce PostScript file ex4b.ps
% Second figure to example 4 in the VertEgg manual
% Contours of egg velocities
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 22 September 1995

clf

d = [0:0.1:4];      % mm
drho = [0:0.2:6];   % kg/m^3
dtab = ones(size(drho))' * d / 1000;    % d constant in colums
drhotab = drho' * ones(size(d)); % drho constant in rows 
[W Re] = eggvel(drhotab, dtab);

% Velocity contours
cl = [0.1:0.1:0.5 1.0:0.5:3.0 4.0:9.0];
c = contour(d, drho, 1000*W, cl, '-');
clabel(c);
title('Terminal egg velocity [mm/s]');
xlabel('Diameter [mm]')
ylabel('Density difference [kg/m^3]')
zlabel('Velocity [mm/s]')

% Reynolds number contours
hold on
cl = [0.5 5.0];
contour(d, drho, Re, cl, ':');
hold off

print -deps ex4b.ps

