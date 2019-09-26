function A = fluxlim(A0, K, W, nstep, dt, P, alpha)

% fluxlim     Numerical integration of transport equation
% ---------------------------------------------------------
% USAGE: A = fluxlim(A0, K, W, nstep, dt, P, alpha)
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
%   flux-limited method. Starting with the concentration
%   in A0 nstep integration steps are performed. The 
%   result is saved in A.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 3 September 1995.
%   Last modified: 7 October 1995
% ----------------------------------------------------

global Ncell
global dz

if (nargin == 5)
  source = 0;    % No source
elseif (nargin == 7)
  source = 1;    % Source
else
  error('fluxlim must have 5 or 7 arguments');
end


C  = W * dt / dz;       % Courant number
S  = K * dt / (dz*dz);  % Diffusion parameter
Cp = 0.5*dt*(W + abs(W))/dz;
Cm = 0.5*dt*(W - abs(W))/dz;

II  = [2 : Ncell];      % Index array

% Initiate  A and F
A = A0;
F = zeros(Ncell+1,1);
A1 = zeros(Ncell,1);

% Lax-Wendroff flux coefficients (+ diffusion)
BLWm = 0.5*(C(II)-C(II).*C(II)) - S(II);  % Coeff. to A(II-1)
BLW  = 0.5*(C(II)+C(II).*C(II)) + S(II);  % Coeff. to A(II)

% Low-order coefficients
% Upstream flux component (+ diffusion)
BUSm = Cm(II)-S(II);
BUS  = Cp(II)+S(II);

% Epsilon 
epsilon = 1.0e-30;

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

  % Compute the Lax-Wendroff flux
  F(II) = BLWm .* A(II-1) + BLW .* A(II);
  % Compute the source term
  Q = QP + QA .* A;
  
  % Flux check and limitation
  A1 = A + F(2:Ncell+1) - F(1:Ncell) + Q;
  if (any(A1 < epsilon))
    IJ = 1 + find(A1(2:Ncell-1) < epsilon); 
      F(IJ) = BUSm(IJ-1) .* A(IJ-1) + BUS(IJ-1) .* A(IJ);
      F(IJ+1) = BUSm(IJ) .* A(IJ) + BUS(IJ) .* A(IJ+1);
    if (A1(1) < eps)
      F(2) = BUSm(1) * A(1) + BUS(1) * A(2);
    end
    if (A1(Ncell) < epsilon)
      F(Ncell) = BUSm(Ncell-1) * A(Ncell-1) + BUS(Ncell-1) * A(Ncell);
    end
    % Use the limited fluxes
    A = A + F(2:Ncell+1) - F(1:Ncell) + Q;
  else
     % Use the unmodified Lax-Wendroff fluxes
     A = A1;
  end 
  
end
