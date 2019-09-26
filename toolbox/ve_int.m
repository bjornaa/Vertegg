function M = ve_int(A)

% ve_int -- Vertical integral of egg distribution.
% ------------------------------------------------
% USAGE: M = ve_int(A)
%
% INPUT:
%    A      : Egg distribution    [eggs/m^3];
%
%    A must be a matrix where the columns live on egg-points,
%    but a row-vector at egg-points is also accepted,
%    size(A) = (Ncell x n) or (1 x Ncell)
%
%  OUTPUT:
%    M      : The vertical integral   [eggs/m^2];
%
%    If A is a matrix of size (Ncell x n), M is a
%    row-vector of length n.
%
%  DESCRIPTION
%    Computes the vertical integral of an egg distribution.
%    ve_int is a vector function.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 10 August 1995.
% -----------------------------------------------------

global Ncell;

% Correct number of arguments ??
if (nargin ~= 1)
  error('USAGE: M = ve_int(A)')
end

% Check shape
[m n] = size(A);
if ~((m == Ncell) | ((m == 1) & (n == Ncell)))
  error('Wrong shape of A')
end

M = eggmom(A,0);


