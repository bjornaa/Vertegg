%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runex8
%
% Example script for time evolution with source term
%
% Author: Bjørn Ådlandsvik,     (bjorn@imr.no)
% Institute of Marine Research, 11 January 1999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% ---- Test for mode, and set if not set already
%        mode = 1 gives continuous egg production
%        mode = 2 gives egg production for first 24 hours
%        Defaults to mode = 1 if unset.
%  
if (exist('mode') ~= 1)
  mode = 1;
end

%
% --- Print start message
%
disp('Running example 8')
fprintf(1,'Mode = %d\n', mode)

%
% --- Define grid
%
depth   = 600;
dz0     = 5;

%
% --- Initialize VertEgg and set grid
%
ve_init;
ve_grid(depth, dz0);

%
% --- Read the profile-file 
%
profile  = '../ex6/testprof.dat';
fprintf(1, 'Profile file : %s\n', profile)
fid = fopen(profile, 'rt');
tmp = fscanf(fid, '%f');
fclose(fid);

% - Extract the data
K0      = 0.012;         % Eddy diffusivity scaling
inr = length(tmp);
iZ = tmp(1:4:inr, 1);
iZ = -abs(iZ);           % Make sure that depths are negative
iS = tmp(2:4:inr, 1);
iT = tmp(3:4:inr, 1);
iK = tmp(4:4:inr, 1);
iK = iK * K0;
%
% --- Interpolate the profile data to the flux points
%
S = interp1(iZ,iS,ZF);
T = interp1(iZ,iT,ZF);
K = interp1(iZ,iK,ZF);

%
% --- Egg diameters and salinities
%
d = 3.2 / 1000.;
Se = 34.5;

%
% --- Compute egg velocities
%
W = eggvelst(S,T,d,Se);

%
% --- Spawning term
%     Constant spawning rate: 1 egg/day/m in the lower 100 m
%     (500-600 m)
%
P = zeros(Ncell,1);
P(Ncell-19:Ncell) = 1/(24*3600) * ones(20,1);

%
% --- The loss term
%     Given in percent/day
%
%     Increased "mortality" near neutral level to account
%     for loss by hatching
%     
%     "Bakground" mortality = m1
%     Increased mortality at neutral level = m2
%     Linear transition is used between the two values  using 
%     5 grid cells above and below the neutral level
%
m1 = 10;
m2 = 50;
IN = max(find(W < 0));     % Find neutral level
mort = ones(Ncell,1)*m1;   % Background mortality
mort(IN) = m2;
for i = [IN-1:IN-5]
   mort(i) = m2*(i-IN+6)/6 + m1*(IN-i)/6;
end
for i = [IN+1:IN+5]
   mort(i) = m2*(IN+6-i)/6 + m1*(i-IN)/6;
end
% Compute exponential rate, exp(-alpha*24*3600) = 1-mort/100
alpha = -log(1-mort/100)/(24*3600);

%
% --- Settings for time steping
% 
outstep = 24;    % Output  each day
simtime = 30*24; % 30 days
dt = 600;
nstep  = outstep*3600/dt  % Number of steps between outputs
nout   = simtime/outstep  % Number of output steps

%
% --- Initialization of time loop
% 
A0 = zeros(Ncell,1);
A = A0; 
%X = []; % Clean up X
X = zeros(Ncell,nout);

%
% --- Time loop
%
tic
for t = 1:nout
  if (t == 2) & (mode == 2)
    P = zeros(Ncell,1);
  end
  A = minlim(A,K,W,nstep,dt,P,alpha);
%  A = upstream(A,K,W,nstep,dt,P,alpha);
%  A = lwendrof(A,K,W,nstep,dt,P,alpha);
%  A = fluxlim(A,K,W,nstep,dt,P,alpha);
  X(:,t) = A;
end
toc


disp(' runex8 is finished')



