% ve_init.m    Initialize VertEgg
% ---------------------------------------------------
% USAGE: ve_init
% 
% DESCRIPTION
%  Script for initializing of VertEgg.
%  Declares the global variables, so they become
%  available in the workspace
%  The actual values are set by ve_grid.
%  
%  AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 9 August 1995.
% ----------------------------------------------------


global Ncell    % Number of grid cells
global dz       % Vertical step size  [m]
global Hcol     % Depth of water column [m]
global ZE       % N-vector of egg-point depths [m]
global ZF       % N+1-vector of flux-point depths [m]

fprintf(1,'\n --- VertEgg is initialized --- \n\n')

