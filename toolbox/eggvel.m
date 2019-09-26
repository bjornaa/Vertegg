function [W, Re] = eggvel(drho,d,mu)

% eggvel    Terminal velocity of egg in sea water.
% -------------------------------------------------
% USAGE: [W, Re] = eggvel(drho, d, mu)
%
% INPUT:
%    drho     : Buoyancy of egg   [kg/m^3]
%               Density of water - density of egg.
%    d        : Diameter of egg  [m]
%    mu (opt) : Dynamic molecular viscosity [kgm^-1s^-1]
%
%    drho, d and mu can be matrices (of the same shape).
%    mu can be also be a scalar or omitted.
%    The sign of drho is positive if the egg is ascending.
%
% OUTPUT:
%   W        : terminal velocity   [m/s]
%   Re (opt) : Reynolds number
%
% DESCRIPTION:  
%   Computes the terminal velocity of a small sphere 
%   in sea water. The routine should be useful for
%   Reynolds number Re less than 5.
%   If Re < 0.5 Stokes' formula is used,
%   otherwise a formula of DallaValle is used
%   With only two arguments a default value 1.6e-3 
%   is used for mu.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 9 August 1995.
%  
% REFERENCES:
%   Sundby, S. 1983,
%   A one-dimensional model for the vertical ...
%   Deep-Sea Research, 30, 645-661.
% ---------------------------------------------------

do_fortran_indexing = 'true';

% Correct number of arguments ??
if (nargin ~= 2) & (nargin ~= 3)
  error('Usage: [W, Re] = eggvel(drho, d, mu)')
end

% Shape testing;
[r s] = size(drho);
if any(size(d) ~= [r s])
  error('drho and d must have the same shape');
end

% Test for availability of mu and sort out arrayness
if (nargin == 2)   % Use default value for mu
  mu = 1.6e-3 * ones(r,s);
else 
  if any(size(mu) ~= [r s])    % if shape of mu is wrong
    if all(size(mu) == [1 1])  %   if mu is scalar
      mu = mu * ones(r,s);     %     make it a constant array
     else
      error('mu must be a scalar or shaped like drho and d');
    end
  end
end

% Initiate arrays, with default values
W = zeros(r, s); 
Dmax = ones(r,s);

g = 9.81;      % Acceleration due to gravity
rho = 1025.0;  % Standard density of sea water

sgn = sign(drho);    % sgn = +1 for positiv oppdrift
drho = abs(drho);


% Maximum diameter for Stokes' formula
II = find(drho);    % Dont do anything if drho = 0 
Dmax(II) = ( 9 * mu(II).^2 ./ (rho * g * drho(II)) ).^(1/3);

% Stokes' formula
IS = find(d(II) <= Dmax(II));
W(II(IS)) =  (1/18) * g * d(II(IS)).^2 .* drho(II(IS)) ./ mu(II(IS));

% Dallavalle's formula
ID = find(d(II) > Dmax(II));
W(II(ID)) = 0.08825 * (d(II(ID))-0.4*Dmax(II(ID))) .* ...
            drho(II(ID)).^(2/3) .* mu(II(ID)).^(-1/3);

% correct the sign
W = sgn .* W;

% Compute the Reynolds number
Re = rho * W .* d ./mu;  
