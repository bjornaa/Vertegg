% runex2.m
% 
% Script for runing on of the transient methods 
% with constant coefficients
%
% Some of the variables after execution;
%  M0      : Vertical integrated concentration
%  K0,  K  : Eddy diffusivity used 
%  W0,  W  : Terminal velocity used
%  A0      : Initial egg distribution
%  X       : Time evolution          
%  A       : Last distribution = X(:,nout)
%  outstep : Hours between columns in X    
%  nout    : Number of columns in X
%
%  runex2 is part of the example directory ex2 in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 13 August 1995

%
% The defining variables are collected here,
% for easy modification.

% -- Space --
H = 100;       % Depth [m]
dz0 = 2;        % Grid size [m]

% -- Time --
dt = 120;      % Time step [s]
outstep = 1;   % Time between saving of model results [hours]
simtime = 96;  % Total simulation time [hours]

% -- Coefficients --
K0 = 0.01;     % Eddy diffusivity [m^2/s]
W0 = 0.001;    % Egg velocity     [m/s]

% -- Initial conditions --
M0 = 1000;     % Vertical integral of concentration [eggs/m^2]

clc;  % Clear the command window

disp(' Running the  Lax-Wendroff scheme')
disp(' with constant coefficients')

% Space discretisation
ve_init;
ve_grid(H,dz0);

% Echo the rest of the variables to screen
fprintf(1,' Time step           = %d s\n', dt)
fprintf(1,' Simulation time     = %d hours\n', simtime)
fprintf(1,' Time between output = %d hours\n', outstep)
disp(' ')
fprintf(1,' Eddy diffusivity: K0 = %f m^2/s\n', K0)
fprintf(1,' Egg velocity:     W0 = %f m/s\n', W0)
disp(' ')
fprintf(1,' Total mass   = %5.0f eggs/m^2\n', M0)

% Make ready for new page,
disp(' ')
disp('    Press any key to check stability'); pause
clc

% Compute the characteristic values
C     = W0 * dt/dz;
S     = K0 * dt/(dz*dz);
Pcell = C/S;

disp(' ')
fprintf(1,' Courant number,       C     = %f\n', C)
fprintf(1,' Diffusion parameter,  S     = %f\n', S)
fprintf(1,' Cell Peclet number,   Pcell = %f\n', Pcell)
disp(' ')

% Test for positivity
% This test is for the Lax-Wendroff method
% and must be modified for the other methods
stab = C^2 + 2*S;
pos1 = (abs(C) <= stab);
pos2 = (stab <= 1);
pos = pos1 & pos2;
if (pos)
  fprintf(1,' Lax-Wendroff scheme is positive\n')
else
  fprintf(1,'***Warning: Lax-Wendroff is not positive\n')
  fprintf(1,'  abs(C)    = %f\n', abs(C))
  fprintf(1,'  C^2 + 2*S = %f\n', stab)
end
disp(' ')


% Make initital distribution
%A0 = ve_rand(M0);
A0 = spawn(M0,50);

% Make constant vectors at the flux points from K0 and W0
K = K0*ones(Ncell+1,1);   % Make a vector from K0
W = W0*ones(Ncell+1,1);   %   and W0

% Deduced time variables
nstep  = outstep*3600/dt;  % Number of steps between outputs
nout   = simtime/outstep;  % Number of output steps

disp(' ')
disp('Press any key to start time integration');pause

clc;
disp(' ')
disp(' Starting time integration');
disp(' ')
A = A0; 
X = []; % Clean up X

% Perform the time loop
% The call to lwendrof can be replaced with
% calls to ftcs, upstream, posmet, or minlim.
tic
for t = 1:nout
  A = lwendrof(A,K,W,nstep,dt);
  X(:,t) = A;
  fprintf(1,'Simulated %d hours\n', t*outstep)
end
toc

disp(' ')
disp(' Time integration is finished')











