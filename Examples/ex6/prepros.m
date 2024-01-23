% Prepros.m ---
% 
% Preprosessing script in example 6 in VertEgg
%
% Bjørn Ådlandsvik, bjorn@imr.no
% Havforskningsinstituttet, 12 October 1996
%

%
% --- Initialize VertEgg and set grid
%
disp('Initialize VertEgg');
depth = 600 % m
dz0   = 5   % m
ve_init;
ve_grid(depth, dz0);
%
% --- Read the profile-file
%
disp('Reading profile data');
load  testprof.dat

% --- Extract the input arrays
iZ = -abs(testprof(:,1));   % Depths shall be negative
iS = testprof(:,2);         % Second column is salinity
iT = testprof(:,3);         % Temperatue
K0 = 0.0120;                % Surface vertical mixing
iK = K0*testprof(:,4);      % Scaled mixing coefficient
clear testprof;             
%
% --- Interpolate the profile data to the flux points
%
S = interp1(iZ,iS,ZF);
T = interp1(iZ,iT,ZF);
K = interp1(iZ,iK,ZF);
%
disp('Plotting the physical profiles')
clf
subplot(1,3,1)
plot(S,ZF)
title('Salinity')
ylabel('Depth [m]')
xlabel('[psu]')
subplot(1,3,2)
plot(T,ZF)
title('Temperature')
xlabel('[deg C]')

subplot(1,3,3)
semilogx(K,ZF)
title('Vertical mixing')
xlabel('[m^2/s]')

%print 'ex6a.ps'

disp('Push a key to continue');pause

disp('Preliminary egg calculations')

eggsal = 34.5 
eggdiam = 3.2 
W = eggvelst(S,T,eggdiam/1000.0,eggsal);
W(100)

cla
clf

% Compute Steady state solution
M = 25; % eggs/m^2
ASS = sstate(M,K,W);

disp('Plotting egg velocity and steady state solution')

dt = 600 % 10 min

subplot(1,2,1)
plot(1000*W,ZF)
title('Egg velocity [mm/s]')
ylabel('Depth [m]')
hold on
plot([0 0],[0 -Hcol],'r:') % Add line W = 0
hold off

subplot(1,2,2)
plot(ASS,ZE)
title('Steady state solution')
xlabel('[Eggs/m^3]')

print ex6b.ps

disp('Push a key to continue');pause

disp('Plot numerical characteristics of the transient problem');
c = W * dt / dz;
s = K * dt / (dz*dz);
Pcell = abs(c) ./ s;

cla
clf
subplot(1,3,1)
plot(c,ZF)
axis([-1.0 1.0 -Hcol 0])
title('Courant number')
ylabel('Depth [m]')

subplot(1,3,2)
plot(s,ZF)
axis([0 0.5 -Hcol 0])
title('Diffusive parameter')

subplot(1,3,3)
plot(Pcell,ZF)
title('Cell Peclet number')

%
% --- Check stability ----
%
disp('Starting stability check')

if (any(2*s > 1)) 
  disp('*** FTCS, Lax-Wendroff and Upstream are all instable')
end

if (any(c.*c > 2*s))
  disp('*** FTCS is instable')
else
  if (any(abs(c) > 2*s))
    disp('*** FTCS is stable but not positive')
  else
    disp('*** FTCS is positive')
  end
end
  
if (any(c.*c + 2*s > 1))
  disp('*** Lax-Wendroff is instable')
else
  if (any(abs(c) > c.*c + 2*s))
    disp('*** Lax-Wendroff is stable but not positive')
  else
    disp('*** Lax-Wendroff is positive (and stable)')
  end
end

if (any(abs(c) + 2*s > 1))
  disp('*** Upstream is not stable')
else
  disp('*** Upstream is positive (and stable)')
end


wigg = 2 ./ (1 - C);
if any(max(abs(P) - wigg) > 0)
  disp('*** Warning: Lax-Wendroff will develop wiggles***')
end



print ex6c.ps
