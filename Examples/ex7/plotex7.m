% plotex7
%
% Example script for plotting results from runex7
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
% Institute of Marine Research, 14 October 1996
%
disp('Plotting the normalised solution')
clf;
plot(Anorm,ZE)
title('Stationary solution')
xlabel('Concentration [eggs/m^3]')
ylabel('Depth [m]')

disp(' ')
disp('Compute some statistics:')
disp(' ')
fprintf(1, 'Minimum egg diameter %4.2f mm\n', min(d)*1000)
fprintf(1, 'Maximum egg diameter %4.2f mm\n', max(d)*1000)
fprintf(1, 'Mean egg diameter    %4.2f mm\n', mean(d)*1000)
fprintf(1, 'Standard deviation   %4.2f mm\n', std(d)*1000)
disp(' ')
fprintf(1, 'Minimum egg salinity  %5.2f \n', min(Se))
fprintf(1, 'Maximum egg salinity  %5.2f \n', max(Se))
fprintf(1, 'Mean egg salinity     %5.2f \n', mean(Se))
fprintf(1, 'Standard deviation    %5.2f \n', std(Se))
disp(' ')
fprintf(1, 'Mean depth of distribution %5.1f m\n', -ve_mean(Anorm))
fprintf(1, 'Standard deviation         %5.1f m\n', ve_std(Anorm))

print -deps ex7a.ps

disp(' ')
disp('Press any key to continue'); pause
disp(' ')

disp('Use ve_drint to integrate in 100 m layers')
disp('and plot histogram')

for i = [1:6]
  B(i) = ve_drint(Anorm,(i-1)*100,i*100);
end
X = [-50 -150 -250 -350 -450 -550];
[XX,BB] = bar(X,B);

plot(BB,XX)
title('Histogram of the egg distribution')
xlabel('[Eggs/m^2]')
ylabel('Depth [m]')

disp(' ')
disp('Press any key to continue'); pause
disp(' ')

disp('The convergence of the Monte Carlo method')
disp('can be studied by plotting the "increment"')
disp('R(i) = rms difference between the normalised')
disp('Monte Carlo solutions with i and i-1 samples')

clf
semilogy(R)
title('Root mean square increment')
xlabel('Number of samples')

print -deps ex7b.ps

