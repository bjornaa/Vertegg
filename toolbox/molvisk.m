function mu = molvisk(S,T)
% 
% molvisk -- dynamic molecular viscosity of sea water
% ---------------------------------------------------
% USAGE: mu = molvisk(S,T)
%
% INPUT:
%   S      : Salinity        [psu]
%   T      : Temperature     [deg C]
%
%   S and T may be arrays of the same shape
%
% OUTPUT:
%   mu  : dynamic moldecular viscosity  [kg/(ms)]
%
%   mu is an array of the same shape as S and T
%
% DESCRIPTION:
%   Computes an approximation to the dynamic molecular
%   viscosity of sea water. The formula is 
%   mu = 0.001 * (1.7915 - 0.0538*T + 0.0007*T.^2 + 0.0023*S)
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 20 september 1995
% -------------------------------------------------

% Check number of arguments
if (nargin ~= 2)
  error('Usage: mu = molvisk(S,T)')
end

% Test for array shape mangler

mu = 0.001 * (1.7915 - 0.0538*T + 0.0007*T.^2 + 0.0023*S);
