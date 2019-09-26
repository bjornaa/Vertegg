function M = eggmom(A,p)

% eggmom     moment of egg distribution
% ---------------------------------------------------
% USAGE: M = eggmom(A, p)
%
% INPUT:
%    A      : Egg distribution   [eggs/m^3]
%    p      : Order of moment 
%
%    A must be a matrix where the columns live on egg-points,
%    but a row-vector at egg-points is also accepted,
%    size(A) = (Ncell x n) or (1 x Ncell)
%
%  OUTPUT:
%    M      : The p-th moment of A
%
%    If A is a matrix of size (Ncell x n), M becomes a 
%    row-vector of length m containing the moments of
%    the columns of A. p must not be negative.
%
%  DESCRIPTION
%    Computes the p-th moment of the egg-distribution A,
%
%         0
%    M = int(z^p a(z) dz) 
%        -H
%
%    where a(z) is the piecewice constant function 
%    a(z) = A(i) for ZF(i+1) < z < ZF(i).
%    
%    If A is a matrix, the moments of the columns are
%    calculated.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 10 August 1995.
% -----------------------------------------------------

global ZF
global Ncell

% Correct number of arguments ??
if (nargin ~= 2)
  error('USAGE: M = eggmom(A, p)')
end

% Check shape
[m n] = size(A);
if ~((m == Ncell) | ((m == 1) & (n == Ncell)))
  error('Wrong shape of A')
end
% Transpose a row vector
if (m == 1)  
  A = A'
end

% Check non-negativity of p
if (p < 0)
  error('p must not be negative')
end

% Index-arrays
II = [1:Ncell];
Ip = [2:Ncell+1];

M = 1 / (p+1) * sum(diag(ZF(II).^(p+1) - ZF(Ip).^(p+1)) * A);



