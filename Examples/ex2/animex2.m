% animex2
%
% Script for animating results from runlw.
%
% runlw must have been run first.
%
% animex2 is part of the example directory in VertEgg
%
% Author: Bjørn Ådlandsvik,  (bjorn@imr.no)
%  Institute of Marine Research, 13 August 1995

clc;

disp(' ')
disp(' Animation of results from runlw')

clg

disp(' ')
disp(' Move or resize this window so that it does not')
disp(' cover the graphic window')
disp(' ')
disp(' Press any key to continue'); pause


T = [1:nout];
  
% Time = 0
  plot(A0,ZE)
  title('Time = 0')
  axis([0 100 -Hcol 0])
  drawnow

for t = T
  title_string = sprintf('Time = %d', t*outstep);

  plot(X(:,t), ZE)
  axis([0 100 -Hcol 0])
  title(title_string)
  drawnow
end

disp(' ')
disp(' The end')
