function A = sstate(M,K,W)

% sstate -- Steady state solution
% -------------------------------
% USAGE: A = sstate(M,K,W)
%
% INPUT:
%   M      : Vertical integrated concentration [eggs/m^2]
%   K      : Eddy diffusivity [m^2/s]
%   W      : Egg velocity     [m/s]
%
%   M is scalar, K and W are column vectors of length Ncell+1,
%   K and W live at the flux-points.
%
% OUTPUT:
%   A      : Stationary solution  [eggs/m^3]
%
%   A lives at the egg points, size = (Ncell x 1).
%
% DESCRIPTION
%   Computes the steady state solution of the convection-
%   diffusion equation, with eddy diffusivity K and egg
%   velocity W variable in the water column.
%
%   The solution is computed by the solution where K and W
%   are viewed as exact on the grid cells.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 29 August 1995.
% ----------------------------------------------------

global Ncell;
global dz;
global ZF;
global ZE;

% Cell averaged m
II = [1 : Ncell];
m = 0.5*(W(II) ./ K(II) + W(II+1) ./ K(II+1));

% Find prelinary constants in C by continuity
lnC = zeros(Ncell,1);
for i = 2:Ncell;
  lnC(i) = lnC(i-1) + (m(i-1)-m(i))*ZF(i);
end

% Adjust to prevent exp from overflowing
lnC = lnC - max(m.*ZE+lnC);

% prelinimary value of A;
A = zeros(Ncell,1);
II = find(m ~= 0);
  A(II) = (exp(m(II).*ZF(II)+lnC(II)) - exp(m(II).*ZF(II+1)+lnC(II))) ./ m(II);
I0 = find(m == 0);
  A(I0) = exp(lnC(I0)) .* dz;

%  Adjust for integral condition
A =  (M / (sum(A) * dz)) * A;

