function  A = spawn(M,Z)

% spawn --  concentrated egg distribution
% ---------------------------------------
% USAGE: A = spawn(M,Z)
%
% INPUT:
%   M     : Vertical integrated concentration  [eggs/m^2]
%   Z     : Spawning depth                     [m]]
%
%   M and Z are scalars.
%
% OUTPUT:
%   A     : Egg distribution             [eggs/m^3]
%
%   A is a column vector, living at the egg points.
%
% DESCRIPTION
%   Returns a vertical egg distribution A with vertical
%   integral M, concentrated as much as possible around 
%   depth = Z.
%   If (ZE(Ncell) < Z < ZE(1), then Z is the mean of A.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 11 August 1995.
% ----------------------------------------------------

%% Should be extended to allow for vectors M and Z
%% Tolkning 
%% M,Z vektorer samme størrelse, skal ha summere bidragene
%% A = sum(spawn(M(i),Z(i)))
%%      i
   
% Check number of arguments
if (nargin ~= 2)
  error('Usage: A = spawn(M,Z)')
end

global Ncell
global ZE
global Hcol
global dz

A = zeros(Ncell,1);   % Initiate A

depth = abs(Z);      % Make sure depth is positive

if (depth > Hcol)
  error('Must have abs(Z) <= Hcol')
end

% Find i s.th. -ZE(i) <= depth < -ZE(i+1)
i = floor(0.5 + depth/dz);

if (i == 0)
  A(1) = M / dz;
elseif (i == Ncell)
  A(Ncell) = M / dz;
else
  l = (ZE(i) + depth)/dz;
  A(i) = (1-l) * M/dz;
  A(i+1) = l * M/dz;
end


