function ve_grid(H, dz0)

% ve_grid    Set up vertical grid for VertEgg
% -----------------------------------------------------
% USAGE:  ve_grid(H, dz0)
%
% INPUT: 
%    H          : Depth of water column [m]
%    dz0        : Grid size             [m]
%
% OUTPUT:
%    None
%
% DESCRIPTION:
%  Sets up the vertical grid used in VertEgg, given the depth
%  H of the water column and the grid size dz0.
%  (Re)defines all global variables.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 9 August 1995.
% ----------------------------------------------------


global Ncell    % Number of grid cells
global dz       % Vertical step size  [m]
global Hcol     % Depth of water column [m]
global ZE       % N-vector of egg-point depths [m]
global ZF       % N+1-vector of flux-point depths [m]

% Test for correct number of arguments
if (nargin ~= 2)
  error('usage:  ve_grid(H, dz0)')
end

if (dz0 == 0) 
  error('dz0 must be nonzero')
end

Hcol = abs(H);  % Make sure the values are positive
dz = abs(dz0);

Ncell = ceil(Hcol/dz);

if (Hcol ~= Ncell*dz)
  error('H must be divisible with dz0')
end 

% Define the egg- and flux-points.
ZE  = [-0.5*dz : -dz : -(Ncell-0.5)*dz]';
ZF = [0 : -dz : -Hcol]';

% Write a message
fprintf(1,' --- Grid for VertEgg --- \n')
fprintf(1,'Hcol   = %f\n',Hcol)
fprintf(1,'dz     = %f\n',dz)
fprintf(1,'Ncell  = %d\n',Ncell)
fprintf(1,'\n')

