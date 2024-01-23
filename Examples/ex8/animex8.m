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

%
% Gjør det bedre med å putte data-ene på plass
%  

clc;

disp(' ')
disp(' Animation of results from example 8')

clf

disp(' ')
disp(' Move or resize this window so that it does not')
disp(' cover the graphic window')
disp(' ')
disp(' Press any key to continue'); pause

if mode == 1
  maxegg = 12;
else
  maxegg = 3;
end


%T = [1:nout];
T = [1:10];
  
% Time = 0
  h = plot(A0,ZE);
  title('Time = 0')
  axis([0 maxegg -Hcol 0])
  drawnow

for t = T
  title_string = sprintf('Time = %d', t*outstep);
  set(h, 'XData', X(:,t))
  title(title_string)
  pause(1)
end

disp(' ')
disp(' The end')
