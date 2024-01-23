% velex
%
% Example script for testing the eggvel
%
% velex.m is part of the example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 13 August 1995

clc;
clf;

disp(' ')
disp(' -- Example 4,  Velocity computations')

% Necessary constants
mu = 1.6e-3;          % Dynamic molecular viscosity
g = 9.81;             % Acceleration due to gravity
rho = 1027;           % Standard density of sea water

disp(' ')
disp('   First we will plot the maximum diameter Dmax')
disp('   for which Stoke''s formula applies')
disp('   The formula is')
disp('   Dmax = ( 9 * mu^2 / (rho * g * drho) )^(1/3)')


drho = [0.1:0.1:5];   % Range of density differences
Dmax = ( 9 * mu^2 ./ (rho * g * drho) ).^(1/3);

plot(drho, Dmax * 1000)  % Use mm as unit
title('Maximum egg diameter for Stokes'' formula')
xlabel('Density difference')
ylabel('Diameter [mm]')

disp(' ')
disp('Press any key to continue')
pause
clc
clf

disp('   Next we will plot velocity as function of diameter')
disp('   for different density differences.')

d = [0:0.1:5];   % Range of diameters (in mm)

hold on;
axis([0 5 0 5]);
title('Terminal velocity at different density differences')
xlabel('Diameter [mm]')
ylabel('Velocity [mm/s]')
for drho = [0.25 0.5 1:6]   % The values of drho used by Sundby
  % Calling eggvel, drho must be an array of size = size(d)
  % d must be given in m.
  W = eggvel(drho * ones(size(d)), d/1000, mu);
  plot(d, W*1000, 'g');
  fprintf(1,'  drho = %5.2f kg/m^3\n', drho)
end

disp(' ')
disp('Press any key to complete the figure')
pause    % Press any key

disp(' ')
disp('  Add the curves Re = 0.5 and Re = 5.0 in red');
disp('  This gives figure 1 in Sundby (1983)')

d(1) = [];   % Remove d(1) to prevent divison by zero;
W = 0.5 * mu ./ (rho * d/1000);
W = 1000 * W;       % Convert to mm/s
plot(d, W, 'r');    % Re = 0.5
plot(d, 10*W, 'r'); % Re = 5.0

disp('Press any key to continue')
pause;  % Press any key
clc
clf

hold off

disp('   Make a velocity table')

d = [0:0.1:4];      % mm
drho = [0:0.2:6];   % kg/m^3
dtab = ones(size(drho))' * d / 1000;    % d constant in colums
drhotab = drho' * ones(size(d)); % drho constant in rows 
[W Re] = eggvel(drhotab, dtab);

disp('   and make a contour plot of it');

cl = [0.1:0.1:0.5 1.0:0.5:3.0 4.0:9.0];
c = contour(d, drho, 1000*W, cl, 'g');
clabel(c);
title('Terminal egg velocity [mm/s]');
xlabel('Diameter [mm]')
ylabel('Density difference [kg/m^3]')
zlabel('Velocity [mm/s]')

disp(' ');
disp('Press a key complete the figure');pause
disp(' ');
disp('   Add the red lines Re = 0.5, Re = 5');

hold on
cl = [0.5 5.0];
contour(d, drho, Re, cl, 'r');
hold off

disp(' ')
disp('Press any key to continue'); pause
disp(' ')
disp('   Make a surface plot of the velocity table\')

surf(d, drho, 1000*W)
title('Terminal velocity');
xlabel('Diameter [mm]')
ylabel('Density difference [kg/m^3]')
zlabel('Velocity [mm/s]')


disp(' ')
disp('The end')


