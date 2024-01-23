% animex5
%
% Script for animating results from runex5.
%
% runex5 must have been run first.
%
% animex5 is part of the example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 14 August 1995

clc;

disp(' ')
disp(' Animation of results from example 5')

clg

disp(' ')
disp(' Move or resize this window so that it does not')
disp(' cover the graphic window')
disp(' ')
disp(' Press any key to continue'); pause


T = [1:nout];


plot([A10 A20 A0], ZE)
axis([0 15 -Hcol 0])
title('Time = 0')
drawnow

for t = T
  title_string = sprintf('Time = %d', t*outstep);
  plot([X1(:,t) X2(:,t) X(:,t)], ZE)
  axis([0 15 -Hcol 0])
  title(title_string);
  drawnow;
end

disp(' ')
disp(' The end')
