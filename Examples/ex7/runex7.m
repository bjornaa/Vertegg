% runex7
%
% Example script for description of steady state
% with sstate and the Monte Carlo method
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 14 October 1996
%



%
% --- Define grid
%
Depth   = 600;
dz0     = 5;

% Number of sstate profiles
Nprof = 1000

% Eggdiameter and standard deviation 
eggdia = [3.2 0.2];  %  [mm]
eggdia = 0.001 * eggdia; % Convert to m
% Egg salinity and standart deviation
eggsal = [34.5 0.2];

disp(' ')
fprintf(1,'Number of profiles : %d\n', Nprof)
fprintf(1,'Eggdiameter with std deviation: %f  %f\n', eggdia(1), eggdia(2))
fprintf(1,'Eggsalinity with std deviation: %f  %f\n', eggsal(1), eggsal(2))
disp(' ')

%
% --- Initialize VertEgg and set grid
%
ve_init;
ve_grid(Depth, dz0);

%
% --- Read the profile-file 
%
profile  = '../ex6/testprof.dat';
fprintf(1, 'Profile file : %s\n', profile)
fid = fopen(profile, 'rt');
tmp = fscanf(fid, '%f');
fclose(fid);
% - Eddy diffusivity scaling
K0      = 0.012;
% - Extract the data
inr = length(tmp);
iZ = tmp(1:4:inr, 1);
iZ = -abs(iZ);           % Make sure that depths are negative
iS = tmp(2:4:inr, 1);
iT = tmp(3:4:inr, 1);
iK = tmp(4:4:inr, 1);
iK = iK * K0;
%
% --- Interpolate the profile data to the Flux points
%
S = interp1(iZ,iS,ZF);
T = interp1(iZ,iT,ZF);
K = interp1(iZ,iK,ZF);

%
% --- Compute egg diameters and salinities
%
d  = eggdia(1) + eggdia(2)*randn(1,Nprof);
Se = eggsal(1) + eggsal(2)*randn(1,Nprof); 


%
% --- Initiate Monte Carlo loop
%
M = 100;             % Vertical integrated consentration [eggs/m^2]
W = [];              % Velocity profile
A = zeros(Ncell,1);  % Vector for egg distribution
Anorm = A;           % Normalised egg distr.

tic;
for i = 1:Nprof
  if (rem(i,25) == 0)
    fprintf(1,' --- %d ----\n', i);
  end;

  W = eggvelst(S, T, d(i), Se(i));
  X = sstate(1,K,W);

  A = A + X;

% Normalize,
  Aold = Anorm;
  Anorm = (M / ve_int(A)) * A;
  
% - Compute rms increment from last distribution  
  R(i) = ve_rmsd(Aold,Anorm);    
  
end;

disp(' runex7 is finished, plot results with plotex7')



