function int = ve_drint(A, z1, z2)

% ve_drint    integration over a depth range
% ------------------------------------------
% USAGE: int = ve_drint(A, z1, z2)
%
% INPUT:
%    A      : Egg distribution   [eggs/m^3]
%    z1     : First integration limit   [m]
%    z2     : Second integration limit  [m]
%
%    A must be a matrix where the columns are 
%    egg distributions. But a row-vector can also 
%    be accepted, size(A) = (Ncell x n) or (1 x Ncell)
%
%  OUTPUT:
%    int    : The integral of A over the depth range
%
%    If A is a matrix of size (Ncell x n), int becomes
%    a row-vector of length m containing the integrals
%    of the columns of A.
%
%  DESCRIPTION
%    Computes the vertical integral
%
%           z2
%    int = int(a(z) dz) 
%          z1
%
%    where a(z) is the piecewice constant function 
%    a(z) = A(i) for ZF(i+1) < z < ZF(i).
%    
%    If A is a matrix, the integral is computed for
%    each column
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 23 September 1995.
% -----------------------------------------------------


global Hcol
global dz
global ZF

% Make -Hcol <= z2 <= z1 <= 0
z1 = -abs(z1);
z2 = -abs(z2);
if (z1 < z2)
  tmp = z2;
  z2 = z1;
  z1 = tmp;
end
if (z2 < -Hcol) 
  disp('***Warning from ve_drint:')
  disp('   Depth range extends below bottom')
  z2 = -Hcol;
end

% Find the grid cells containing the points
i1 = ceil(-z1/dz);
if (i1 == 0)
  i1 = 1;
end
i2 = ceil(-z2/dz);
if (i2 == 0)
  i2 = 1;
end

% Perform the integration
if (i1 < i2)
  int = (z1-ZF(i1+1))*A(i1,:);        % Cell i1
  int = int + sum(A(i1+1:i2-1,:))*dz; % Cells between
  int = int + (ZF(i2)-z2)*A(i2,:);    % Cell i2
else
  int = (z1-z2)*A(i1,:)      % Cell i1 = Cell i2
end %if
  

