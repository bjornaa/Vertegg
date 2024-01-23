% visklsq.m  95/09/22
% 
% Example script for least squares regression
% of dynamic molekular viscosity
%
% This is part of example 4 of VertEgg
%

clc;

disp('The table of dynamic molecular viscosity')
disp('from Sverdrup et al (1952).')

mutab = [
1.79   1.52  1.31  1.14  1.01  0.89  0.80
1.82   1.55  1.34  1.17  1.03  0.91  0.82
1.85   1.58  1.36  1.19  1.05  0.93  0.84
1.88   1.60  1.38  1.21  1.07  0.95  0.86
1.89   1.61  1.39  1.22  1.09  0.96  0.87
]

mutab = 0.001 * mutab;  % Convert to SI units

[M N] = size(mutab);

MN = M*N;

disp('The corresponding temperatures (horizontally)')
T = [0:5:30]
disp('and salinities (vertically)')
S = [0:10:30 35]

disp(' ')
disp('Press any key to perform the regression');pause

% Make T and S into MxN tables
T1 = ones(M,1)*T;
S1 = S' * ones(1,N);

% Reshape everything to column vectors
B = reshape(mutab,MN,1);
TI = reshape(T1,MN,1);
SI = reshape(S1,MN,1);

% The normal matrix
A = [ones(MN,1) TI TI.^2 SI];

% Equal weighting
V = eye(MN);

% Find the coefficients
X = lscov(A,B,V);
% and make them more visible
X * 1000

% Compute errors
E = (A*X -B);
disp('Maximum (absolute) error [kg/(ms)]')
max(abs(E))

% Relativ erros
Er = E ./ B;
disp('Makximum relative error [percent]')
max(abs(Er))

