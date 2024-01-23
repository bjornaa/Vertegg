Y = ve_int(X);

%plot(Y, 'LineWidth', 1.5)
semilogy(Y, 'LineWidth', 1.5)
xlabel('Egg age [days]')
ylabel('Number of eggs  [eggs/m^2]')
box on


print -deps ex8c.ps
