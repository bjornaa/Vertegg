function A = upstream(A0, K, W, nstep, dt, P, alpha)

% upstream --  Numerical integration of transport equation
% --------------------------------------------------------
% USAGE: A = upstream(A0, K, W, nstep, dt, P, alpha)
%
% INPUT:
%    A0          : Start concentration   [eggs/m^3]
%    K           : Eddy diffusivity      [m^2/s]
%    W           : Terminal velocity     [m/s]
%    nstep       : Number of integration steps
%    dt          : Time step             [s]
%    P     (opt) : Egg production        [eggs/m^3/s]
%    alpha (opt) : Egg loss rate         [1/s]
%
%    A0 lives at the egg-points, size(A0) = (Ncell x 1).
%    K and W live at the flux-points, size = (Ncell+1 x1).
%    If P and alpha are present they shall live at the egg-points.
%
% OUTPUT:
%    A       : Result concentration  [eggs/m^3]
%    
%    A lives at egg-points in the same way as A0
%
% DESCRIPTION
%   Integrates the convection-diffusion equation by the 
%   upstream method. Starting with the concentration A0,
%   nstep integration steps are performed. The result is 
%   saved in A.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 21 August 1995.
%   Last modified: 7 October 1995
% ----------------------------------------------------

global Ncell
global dz

if (nargin == 5)
  source = 0;    % No source
elseif (nargin == 7)
  source = 1;    % Source
else
  error('upstream must have 5 or 7 arguments');
end


Cp = 0.5*dt*(W + abs(W))/dz;
Cm = 0.5*dt*(W - abs(W))/dz;
S  = K * dt / (dz*dz);

II  = [2 : Ncell];      % Index array
Im = [1 : Ncell-1];     % Index array - 1

% Initiate the A and F
A = A0;
F = zeros(Ncell+1,1);

Bm = Cm(II)-S(II);   % Coeff. to A(Im)
B  = Cp(II)+S(II);   % Coeff. to A(II)

% Source terms
if (source == 1)
  QP = P .* (1 - exp(-alpha*dt)) ./ alpha;
  QA = exp(-alpha*dt) - 1;
else
  QP = zeros(Ncell,1);
  QA = zeros(Ncell,1);
end 


% Time integration loop
for t = 1:nstep
  % Compute the flux
  F(II) = Bm .* A(Im) + B .* A(II);
  % Compute the source term
  Q = QP + QA .* A;
  % Update the solution
  A = A + F(2:Ncell+1) - F(1:Ncell) + Q;
end
  
