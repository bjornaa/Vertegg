% modell.m  (Bedre navn ??)
%
% Egg simulation modell with File I/O
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 24 August 1995
%

%
% --- Read the setup file, modell.sup
%
sup = fopen('modell.sup', 'rt');

profile  = fscanf(sup, '%*s %s', 1);
initfile = fscanf(sup, '%*s %s', 1);
outfile  = fscanf(sup, '%*s %s', 1);

depth   = fscanf(sup, '%*s %f', 1);
dz0     = fscanf(sup, '%*s %f', 1);
dt      = fscanf(sup, '%*s %f', 1);
outstep = fscanf(sup, '%*s %f', 1);
simtime = fscanf(sup, '%*s %f', 1);

K0      = fscanf(sup, '%*s %f', 1);

eggdia = [];
eggsal = [];
Ngrp = fscanf(sup, '%*s %d', 1);
for i = 1:Ngrp
  eggdia(i) = fscanf(sup, '%f', 1);
  eggsal(i) = fscanf(sup, '%f', 1);
end
eggdia = eggdia / 1000;   % Konvert to m

fclose(sup);

%
% --- Print the contents of the setup file
%

fprintf(1,'Physical profile file : %s\n', profile)
fprintf(1,'Initial distrib. file : %s\n', initfile)
fprintf(1,'Output file           : %s\n', outfile)
disp(' ')
fprintf(1,'Depth            : %5.1f m\n', depth)
fprintf(1,'Space step       : %5.2f m\n', dz0)
fprintf(1,'Time step        : %5.0f s\n', dt)
fprintf(1,'Output time step : %5.1f hours\n', outstep)
fprintf(1,'Simulation time  : %5.1f hours\n', simtime)
disp(' ')
fprintf(1,'Vertical mixing in mixed layer %6.4f m^2/s', K0)
disp(' ')
fprintf(1,'Number of egg groups : %d\n', Ngrp)
for i = 1:Ngrp
  fprintf(1,'Group %2d, diameter %5.2f mm, salinity %6.3f\n',...
              i, eggdia(i)*1000, eggsal(i) )
end


%
% --- Initialize VertEgg and set grid
%
ve_init;
ve_grid(depth, dz0);
%
% --- Read the profile-file
%
fid = fopen(profile, 'rt');
tmp = fscanf(fid, '%f');
fclose(fid);
tmp = reshape(tmp, 4, length(tmp)/4)';  % Remake matrix from file
iZ = -abs(tmp(:,1));   % Depths must be negative
iS = tmp(:,2);
iT = tmp(:,3);
iK = K0*tmp(:,4);      % Scale the turbulence
%clear tmp;
%
% --- Interpolate the profile data to the Flux points
%
S = interp1(iZ,iS,ZF);
T = interp1(iZ,iT,ZF);
K = interp1(iZ,iK,ZF);
%
% --- Read the file with initial distribution
%
fid = fopen(initfile, 'rt');
tmp = fscanf(fid, '%f');
fclose(fid);
tmp = reshape(tmp, Ngrp+1, length(tmp)/(Ngrp+1))';  % Shape in file
iZ = -abs(tmp(:,1));
iA = tmp(:,2:Ngrp+1);
clear tmp
%
% --- Interpolate to egg points
%
A0 = interp1(iZ, iA, ZE);
%
% --- Compute the egg velocities
%
W = [];
for i = 1:Ngrp
  W(:,i) = eggvelst(S, T, eggdia(i), eggsal(i));
end

%
% --- Check stability
%
s = K * dt/(dz*dz) * ones(1,Ngrp);  % Same for all groups
C = W * dt/dz;
P = C ./ s;

lwstab = C .* C + 2 * s;
if any(lwstab > 1)
  disp('*** Warning: Lax-Wendroff is instable***')
end

%
% krav til upstream, Cp(i) - Cm(i+1) + S(i) + S(i+1) <= 1
% Garanterer positivitet => stabilitet
Cp = 0.5*dt*(W + abs(W))/dz;
Cm = 0.5*dt*(W - abs(W))/dz;
% Justere rendene
Cp(1,:) = zeros(1,Ngrp);
Cm(1,:) = zeros(1,Ngrp);
Cp(Ncell+1,:) = zeros(1,Ngrp);
Cm(Ncell+1,:) = zeros(1,Ngrp);
usstab(1:Ncell) = Cp(1:Ncell) - Cm(2:Ncell+1) + s(1:Ncell) + s(2:Ncell+1);
if any(usstab > 1)
  disp('*** Warning: upstream is instable***')
  disp('*** Warning: posmet  is instable***')
end

wigg = 2 ./ (1 - C);
if any(max(abs(P) - wigg) > 0)
  disp('*** Warning: Lax-Wendroff will develop wiggles***')
end

%
%  Initiate time loop
%
fid = fopen(outfile, 'wt');    % Open output file
nstep = outstep * 3600 / dt;   % # steps between output
nout  = simtime / outstep;     % # output times
A = A0;

% Write grid-info at start of file
fprintf(fid, 'depth = %6.1f, dz = %5.2f,', Hcol, dz);
% Write time-info at start of file
fprintf(fid, ' outstep = %d, simtime = %d,', outstep, simtime);
% Write nr of egg-groups
fprintf(fid, ' Ngrp = %d\n', Ngrp);

% Save initial distribution
for k = 1:Ncell
  for g = 1:Ngrp
    fprintf(fid, '  %6.2f', A(k,g));
  end
  fprintf(fid, '\n');
end
%
%
%
X = [];

for t = 1:nout
  fprintf(1,' --- t = %d --- \n', t*outstep);
  for i = 1:Ngrp
   A(:,i) = posmet(A(:,i), K, W(:,i), nstep, dt);
  end

% Save A as column t in X
  if (Ngrp == 1)
    X(:,t) = A;
  else
    X(:,t) = sum(A')';  % Row-sums (virker ikke hvis A = vektor)
  end

% Save A in output file
  for k = 1:Ncell
    for i = 1:Ngrp
      fprintf(fid, '  %6.2f', A(k,i));
    end
    fprintf(fid, '\n');
  end
  
end
%
% --- Clean up and finish
%
fclose(fid);
  


