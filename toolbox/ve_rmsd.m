function R = ve_rmsd(X,Y)

% ve_rmsd -- Root mean square deviation 
% -------------------------------------
% USAGE: R = ve_rmsd(X,Y)
%
% INPUT:
%   X    : Egg concentration    [eggs/m^3]
%   Y    : Egg concentration    [eggs/m^3]
%
%   X and Y may be matrixes of the same size
%   with Ncell rows. X and/or Y may also be
%   vectors of length Ncell.
%
% OUTPUT:
%   R    : Root mean square deviation [eggs/m^3]
%
%   R is a row vector with length =
%   max(colums(X), columns(Y)).
%
% DESCRIPTION:
%   Computes the root mean square deviation R 
%   between the columns of X and Y. If one
%   argument is a vector, it is compared to
%   all columns of the other argument.
%   
%   If X and Y are matrices
%     R(j) = sqrt(sum((X(i,j)-Y(i,j))^2)/Ncell).
%   If Y is a vector
%     Y(j) = sqrt(sum((X(i,j)-Y(i))^2)/Ncell).
%                i
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 12 August 1995.
% ------------------------------------------------


%global dz
global Ncell

% Check number of arguments
if (nargin ~= 2)
  error('USAGE: R = ve_rmsd(X,Y)')
end

% Shape checking
[mx nx] = size(X);
if (mx ~= Ncell)
  if ((mx == 1) & (nx == Ncell)) %  Row-vector
    X = X'; mx = Ncell; nx = 1;
  else
    error('X must have Ncell rows')
  end
end

[my ny] = size(Y);
if (my ~= Ncell)
  if ((my == 1) & (ny == Ncell)) %  Row-vector
    Y = Y'; my = Ncell; ny = 1;
  else
    error('Y must have Ncell rows')
  end
end

n = max(nx, ny);
if (((nx ~= n) & (nx ~= 1)) | ((ny ~= n) & (ny ~= 1)))
  error('If X and Y are matrices, the shape must be equal');
end

if ((nx == 1) & (nx < n))
  X = X * ones(1,n);  % Repeat the columns
end

if ((ny == 1) & (ny < n))
  Y = Y * ones(1,n);  % Repeat the columns
end

% Root mean square deviation
R = sqrt(sum((X-Y) .* (X-Y))/Ncell);


