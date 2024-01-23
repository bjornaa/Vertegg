% runsens.m
% 
% Script for testing the different numerical schemes with
% constant coefficients
%
% Some of the variables after execution;
%  M0      : Vertical integrated concentration
%  K0,  K  : Eddy diffusivity used 
%  W0,  W  : Terminal velocity used
%  A0      : Initial egg distribution
%  X       : The computed stationary distributions
%  Y       : The errors in the computed distributions
%
%  runsens is part of VertEgg example 3.
%  
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 13 September 1995

%
% The testparameters are put first
% for easy modifications.
%
% Diffusion coefficient [m^2/s]
%K0 = 0.01;
K0 = 0.0005;
% Vertical velocity     [m/s]
%W0 = 0.001;
W0 = -0.001;
% Space step            [m]
%dz0 = 2;
dz0 = 2;
% Time step             [s]
%dt = 120; 
dt = 120; 

%
% Other defining variables
%

% -- Space --
H = 100;       % Depth [m]

% -- Time --
simtime = 72;   % Total simulation time [hours]
 
% -- Initial conditions --
M0 = 1000;     % Vertical integral of concentration [eggs/m^2]

% Initialize VertEgg.
ve_init;
ve_grid(H,dz0);


% Compute the characteristic values
C     = W0 * dt/dz;
S     = K0 * dt/(dz*dz);
Pcell = abs(C)/S;

disp(' ')
fprintf(1, '  K        W       dz    dt  C     S     Pcell\n')
fprintf(1, ' %7.5f %7.4f %4.1f %5.0f %5.3f %5.3f %5.1f\n', ...
            K0, W0, dz, dt, C, S, Pcell)
disp(' ')


% Make initital distribution
A0 = spawn(M0,50);

% Make constant vectors at the flux points from K0 and W0
K = K0*ones(Ncell+1,1);   % Make a vector from K0
W = W0*ones(Ncell+1,1);   %   and W0
X = [];   % Remove old stuff in X

% Deduced time variables
nstep  = simtime*3600/dt;  % Number of steps 

% Perform the time loop
tic
A = ftcs(A0,K,W,nstep,dt);
tid(1) = toc;
X(:,1) = A;

tic
A = lwendrof(A0,K,W,nstep,dt);
tid(2) = toc;
X(:,2) = A;

tic
A = upstream(A0,K,W,nstep,dt);
tid(3) = toc;
X(:,3) = A;

tic
A = posmet(A0,K,W,nstep,dt);
tid(4) = toc;
X(:,4) = A;

tic
A = minlim(A0,K,W,nstep,dt);
tid(5) = toc;
X(:,5) = A;

% Compute eggsact solution
B = eggsact(M0,K0,W0);
Bmean = ve_mean(B);
Bstd  = ve_std(B);

% Root mean square error
R = ve_rmsd(X,B);

% Error in mean
Emean = ve_mean(X) - Bmean;
Estd  = ve_std(X) - Bstd;

%% Make table
disp(' ')
fprintf(1, '          RMSE    Emean   Estd     CPU-time\n');
fprintf(1, 'ftcs     %6.3f  %6.3f  %6.3f   %5.1f\n',...
              R(1), Emean(1), Estd(1), tid(1))
fprintf(1, 'lwendrof %6.3f  %6.3f  %6.3f   %5.1f\n',...
              R(2), Emean(2), Estd(2), tid(2))
fprintf(1, 'usptream %6.3f  %6.3f  %6.3f   %5.1f\n',...
              R(3), Emean(3), Estd(3), tid(3))
fprintf(1, 'posmet   %6.3f  %6.3f  %6.3f   %5.1f\n',...
              R(4), Emean(4), Estd(4), tid(4))
fprintf(1, 'minlim   %6.3f  %6.3f  %6.3f   %5.1f\n',...
              R(5), Emean(5), Estd(5), tid(5))

Y = X - B *ones(1,5);











