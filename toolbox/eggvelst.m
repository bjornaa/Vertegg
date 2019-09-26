function [W, Re] = eggvelst(S,T,d,Se)

% eggvelst -- Egg velocity from salinity and temperature
% ------------------------------------------------------
% USAGE: [W, Re] = eggvelst(S,T,d,Se)
%
% INPUT:
%   S       : Salinity of the environment  [psu]
%   T       : Temperature  --- " ----      [deg C]
%   d       : Egg diameter                 [m]
%   Se      : Egg salinity                 [psu]
%
%  All arguments can be arrays of the same shape.
%  Alternatively d and/or Se may be scalars.
%
% OUTPUT:
%  W        : Terminal egg velocity   [m/s]
%  Re (opt) : Reynolds number 
%
% DESCRIPTION
%   Computes the terminal egg velocity given the hydrography
%   of the environment and the salinity Se where the egg is
%   neutral buoyant.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 11 August 1995.
% ----------------------------------------------------

% Check the number of arguments
if (nargin ~= 4)
  error('Usage: [W, Re] = eggvelst(S,T,d,Se)')
end

% Check array shapes
[mS nS] = size(S);
[mT nT] = size(T);
[md nd] = size(d);
[mSe nSe] = size(Se);

if ((mT ~= mS) | (nT ~= nS))
  error('S and T must have the same shape')
end

if ((md ~= mS) | (nd ~= nS))
  if ((md == 1) & (nd == 1))   % Scalar
    d = d * ones(mS, nS);
  else
    error('d must have same shape as S or be a scalar')
  end 
end %if

if ((mSe ~= mS) | (nSe ~= nS))
  if ((mSe == 1) & (nSe == 1))   % Scalar
    Se = Se * ones(mS, nS);
  else
    error('Se must have same shape as S or be a scalar')
  end 
end %if

% Compute density difference
drho = dens0(S,T) - dens0(Se,T);

% Temperature dependent molecular viscosity (dynamic)
%mu  = 1.854e-03 * exp(-0.02783*T);  %% Sjekk opp
% Er inlining raskere ???
mu = molvisk(S,T);

[W, Re] = eggvel(drho,d,mu);