B1 = X(:,1);  % 0-1 days
B4 = X(:,4);  % 3-4
B7 = X(:,7); % 6-7 days
B10 = X(:,10);

%hold on
plot([B1 B4 B7 B10], ZE)
box on
text(1.1, -500.0, '0-1 days')
text(0.8, -250.0, '3-4 days')
text(2.8, -170.0, '6-7 days')
text(0.1, -90.0, '9-10 days')

print -deps ex8b.ps
