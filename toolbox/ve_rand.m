function Y = ve_rand(M)

% ve_rand -- make random egg distribution
% ---------------------------------------
% USAGE: Y = ve_rand(M)
%
% INPUT:
%   M      : Vertical integrated concentration [eggs/m^2]
%
%   M is a scalar
%
% OUTPUT:
%   Y      : Random egg distribution  [eggs/m^3]
%
%   Y is vector of size [Ncell 1]
% 
% DESCRIPTION
%   A uniform distribution is used to generate random 
%   values between 0 and 1. The values are thereafter
%   scaled to make ve_int = M.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 16 August 1995.
% -----------------------------------------------

if (nargin ~= 1) 
  error('USAGE: Y = ve_rand(M)')
end

global Ncell;

rand('uniform');  % Needed in octave, legal in MATLAB

% Make a random vector,
Y = rand(Ncell,1); 

% and scale it.
m = ve_int(Y);
Y = (M / m) * Y;
