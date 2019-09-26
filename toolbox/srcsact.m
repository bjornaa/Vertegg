function Y = srcsact(K,W,P,alpha,Z)

% srcsact -- Stationary solution, const. coeff., source term
% ----------------------------------------------------------
% USAGE: Y = srcsact(K,W,P,alphaZ)
%
% INPUT:
%    K       : Eddy diffusivity      [m^2/s]
%    W       : Terminal velocity     [m/s]
%    P       : Egg production        [eggs/m^3/s]
%    alpha   : Egg loss rate         [1/s]
%    Z (opt) : Vertical coordinate   [m]
%
%    K, W, P and alpha are scalars. Z can be arbitrary array.
%    IF Z is ommitted, ZE is used as vertical coordinates.
%
%  OUTPUT:
%    Y       : Concentration at depth Z.     [eggs/m^3]
%
%    If Z is present, size(Y) = size(Z),
%    otherwise, size(Y) = (Ncell x 1).
%
%  DESCRIPTION
%    Computes the exact stationary solution of the
%    convection diffusion equation with constant
%    eddy diffusivity K, velocity W, egg production P
%    and loss rate alpha.
%
%    If Z is present, returns array of pointwise values.
%    If Z is not present, returns exact cell averages.
%
% AUTHOR: Bjørn Ådlandsvik  (bjorn@imr.no)
%   Institute of Marine Research, 6 October 1995.
% ----------------------------------------------------

global Ncell    % Number of grid cells
global dz       % Vertical step size  [m]
global Hcol     % Depth of water column [m]
global ZF       % N+1-vector of flux-point depths [m]

% Check number of arguments
if ~((nargin == 4) | (nargin == 5))
  error(' srcsact must have 4 or 5 arguments')
end

% Check scalarity of arguments
if (any(size(K) ~= [1 1]) | any(size(W) ~= [1 1]) | ...
    any(size(P) ~= [1 1]) | any(size(alpha) ~= [1 1]))
  error('M, K, P and alpha must be scalars')
end

% Check positivity of arguments
if ((K <= 0) | (P < 0) | (alpha <= 0))
  error('K, P and alpha must be positive')
end


discr = W*W + 4*K*alpha;

a = (sqrt(discr) + abs(W)) / (2*K);
b = (sqrt(discr) - abs(W)) / (2*K);

if (W ~= 0)
  lnA = log(abs(W)*P/(alpha*b*K)) + ...
          log(1-exp(-b*Hcol)) - log(1-exp(-(a+b)*Hcol));
  lnB = log(abs(W)*P/(alpha*a*K)) - b*Hcol + ...
          log(1-exp(-a*Hcol)) - log(1-exp(-(a+b)*Hcol));
end

if (nargin == 4) % Z not present
  if (W == 0) 
    Y = P / alpha * ones(Ncell,1);
  else
    Z = ZF;
    if (W < 0)  
      Z = ZF(Ncell+1:-1:1);
    end  
    Y = sign(W) * (exp(lnA) * (1/a) * ...
                 (exp(a*Z(1:Ncell))-exp(a*Z(2:Ncell+1)))) / dz ...
        -  exp(lnB + log((1-exp(-b*dz))/b) - b*Z(2:Ncell+1)) / dz ...
        + P / alpha;
  end  % if (W==0)

else   % nargin == 5
  if (W == 0)
    Y = P / alpha * ones(size(Z));
  else
    Z = -abs(Z);   % Tolerate positive depth values
    if (W < 0)
      Z = -Hcol - Z;
    end
    Y = exp(lnA + a*Z) -  exp(lnB - b*Z) + P/alpha;
  end  % if (W == 0)
  
end  % if (nargin   

