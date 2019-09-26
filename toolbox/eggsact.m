function A = eggsact(M,K,W,Z)

% eggsact -- Exact stationary solution, const. coeff.
% ---------------------------------------------------
% USAGE: A = eggsact(M,K,W,Z)
%
% INPUT:
%    M       : Vertical integrated concentration  [eggs/m^2]
%    K       : Eddy diffusivity      [m^2/s]
%    W       : Terminal velocity     [m/s]
%    Z (opt) : Vertical coordinate   [m]
%
%    M, K, W are scalars. Z can be arbitrary array.
%    IF Z is ommitted, ZE is used as vertical coordinates.
%
%  OUTPUT:
%    A       : Concentration at depth Z. [eggs/m^3]
%
%    If Z is present, size(A) = size(Z),
%    otherwise, size(A) = (Ncell x 1).
%
%  DESCRIPTION
%    Computes the exact stationary solution of the
%    convection diffusion equation with constant
%    eddy diffusivity K and velocity W.
%
%    If Z is present, returns array of pointwise values.
%    If Z is not present, returns exact cell averages.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 30 August 1995.
% ----------------------------------------------------

global Ncell    % Number of grid cells
global dz       % Vertical step size  [m]
global Hcol     % Depth of water column [m]
global ZF       % N+1-vector of flux-point depths [m]

% Check number of arguments
if ~((nargin == 3) | (nargin == 4))
  error(' eggsact must have 3 or 4 arguments')
end

% Check that M,K,W are scalars
if (any(size(M) ~= [1 1]) | any(size(K) ~= [1 1]) | ...
    any(size(W) ~= [1 1]))
  error('M, K, W must be scalars')
end

% 

% Check positivity of K
if (K <= 0)
  error('Must have positive diffusion coefficient K')
end

m = W/K;

if (nargin == 3)
  % Compute 
  if (m == 0) 
    A = ones(Ncell,1) * M / Hcol;
  elseif (m > 0)
    faktor = M / ((1-exp(-m*Hcol)) * dz);
    A = faktor * (exp(m*ZF([1:Ncell])) - exp(m*ZF([2:Ncell+1])));
  else  % m < 0
    faktor = - M / ((1-exp(m*Hcol)) *dz);
    A = faktor * (exp(m*(Hcol+ZF(1:Ncell))) - exp(m*(Hcol+ZF(2:Ncell+1))));
  end

else  %  (nargin == 4)

  Z = -abs(Z);  % Tolerate positive depths

  if (m == 0) 
    faktor = M / Hcol;
    A = faktor * ones(size(Z));
  elseif (m > 0)
    faktor = M * m / (1 - exp(-m*Hcol));
    A = faktor * exp(m*Z);
  else % (m < 0)
    faktor = - M * m / (1 - exp(m*Hcol));
    A = faktor * exp(m*(Hcol+Z));
  end

end %if


