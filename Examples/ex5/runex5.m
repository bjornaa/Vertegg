% runex5.m
%
% Script for reproducing results from Westgård (1989)
%
% Non-constant coeffients, more egg-groups
%
% Some of the variables after execution;
%  S, T         : Hydrography
%  M1, M2       : Integrated concentration in the two groups
%  K            : Eddy diffusivity 
%  W1, W2       : Velocities in the two groups
%  A10, A20, A0 : Initial egg distribution, grp1, grp2, total.
%  X1, X2, X    : Time evolution, grp1, grp2, total.
%  A1, A2, A    : Last distribution from simulation
%  outstep      : Hours between columns in X    
%  nout         : Number of columns in X
%
%  runex5 is part of the example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 14 August 1995
%

disp(' ')
disp(' Example 5 -- Non constant coefficients, more egg groups')
disp(' ')

% Initialize VertEgg
% Depth = 100 m, dz = 2 m.
ve_init
ve_grid(100,2)

% Declare arrays
S = zeros(Ncell+1,1);    % Salinity
T = zeros(Ncell+1,1);    % Temperature
K = zeros(Ncell+1,1);    % Eddy diffusivity
A1 = zeros(Ncell,1);     % Egg concentration, group 1
A2 = zeros(Ncell,1);     % Egg concentration, group 2

disp(' ')
disp(' Pycnocline at 50 m')
disp(' Upper layer: S = 34, T = 7')
disp(' Lower layer: S = 35, T = 5')

% Set up salinity and temperature
I25 = ones(25,1);        
S(1:25) = 34*I25; S(26) = 34.5; S(27:51) = 35*I25;
T(1:25) =  7*I25; T(26) =  6  ; T(27:51) =  5*I25;

% Compute eddy diffusivity by the formula of Sundby (1983)
Wind = 5;
K0 = (76.1 + 2.26*Wind^2) * 1e-4;
K(1:25) = K0*I25; K(26) = 0.55*K0; K(27:51) = 0.1*K0*I25;

disp(' ')
disp(' Eddy diffusivity K corresponds to wind speed 5 m/s')
fprintf(1,' Upper layer K = %f, lower kayer K = %f\n', K0, 0.1*K0)
%
%

% The number of eggs, the diameter and neutral salinity in egg-goup 1.
% Computation of egg velocity array.
M1 = 55;        
d1 = 0.0015;
Se1 = 33;
W1 = eggvelst(S,T,d1,Se1); 
%
% Same for egg-group 2.
M2 = 55;
d2 = 0.0013;
Se2 = 36;
W2 = eggvelst(S,T,d2,Se2);

disp(' ')
disp(' Egg group I')
fprintf(1,'   Integrated concentration : %f eggs/m^2\n', M1)
fprintf(1,'   Diameter : %f mm\n', 1000*d1)
fprintf(1,'   Neutral salinity : %f psu\n', Se1)
disp(' ')
disp(' Egg group II')
fprintf(1,'   Integrated concentration : %f eggs/m^2\n', M2)
fprintf(1,'   Diameter : %f mm\n', 1000*d2)
fprintf(1,'   Neutral salinity : %f psu\n', Se2)


% Set up initial egg distribution in group 1
I10 = ones(10,1);
A1( 1:10) = 0*I10;
A1(11:20) = 0.75*I10;
A1(21:30) = 1.25*I10;
A1(31:40) = 0.75*I10;
A1(41:50) = 0*I10;
% Use the same initital distribution for group 2
A2 = A1;
% Save the intial distributions for later use.
A10 = A1; A20 = A2;
A0 = A10 + A20;

% Time control for numerical solution
dt = 120;      % Time step [s]
outstep = 2;   % Time between saving of model results [hours]
simtime = 96;   % Total simulation time [hours]
% Deduced time variables
nstep  = outstep*3600/dt;  % Number of steps between outputs
nout   = simtime/outstep;  % Number of output steps


disp(' ')
fprintf(1, 'Time step: %d s\n', dt);
fprintf(1, 'Time between output: %d hours\n', outstep);
fprintf(1, 'Simulation time: %d hours\n', simtime);

% Starting the numerical integration
disp(' ')
disp('Press any key to start time integration'); pause
disp(' ')
for t = 1:nout
  A1 = lwendrof(A1,K,W1,nstep,dt);
  X1(:,t) = A1;
  A2 = lwendrof(A2,K,W2,nstep,dt);
  X2(:,t) = A2;
  fprintf(1,'Simulated %d hours\n', t*outstep)
end

% Add the results from the two egg groups
X = X1 + X2;

% Example 5 is finished
disp(' ')
disp(' runex5 is finsihed')



