function s = ve_std(A)

% ve_std -- Standard deviation of egg distribution.
% -------------------------------------------------
% USAGE: s = ve_std(A)
%
% INPUT:
%    A     : Egg distribution    [eggs/m^3];
%
%    A must be a matrix where the columns live on egg-points,
%    but a row-vector at egg-points is also accepted,
%    size(A) = (Ncell x n) or (1 x Ncell)
%
%  OUTPUT:
%    s     : The standard deviation   [m];
%
%    If A is a matrix of size (Ncell x n), s is a
%    row-vector of length n.
%
%  DESCRIPTION
%    Computes the standard deviation of an egg 
%    distribution A,
%    ve_mean is a vector function.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 10 August 1995.
% -----------------------------------------------------

global Ncell;

% Correct number of arguments ??
if (nargin ~= 1)
  error('USAGE: s = ve_std(A)')
end

% Check shape
[m n] = size(A);
if ~((m == Ncell) | ((m == 1) & (n == Ncell)))
  error('Wrong shape of A')
end

M0 = eggmom(A,0);
M1 = eggmom(A,1);
M2 = eggmom(A,2);

mu = M1 ./ M0;
s  = sqrt(M2 ./ M0 - mu .* mu);


