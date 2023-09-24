close all;     clear all;     clc;

% The outputs of the actual plant and the FNN outputs
load datumtel.dat
plot(datumtel(:,1),datumtel(:,2));
hold on
plot(datumtel(:,1),datumtel(:,3),'r');

% The initial and final membership functions
load init_memb.dat;
load final_memb.dat;
x=size(init_memb);
no_rules=x(2);

figure(2)
title('Initial membership functions')
hold on
for i=1:no_rules
    plot(datumtel(:,1),init_memb(:,i));
    axis([0 360 0 1]);
end
figure(3)
title('Final membership functions')
hold on
for i=1:no_rules
    plot(datumtel(:,1),final_memb(:,i));
    axis([0 360 0 1]);
end
%The initial and final membership functions rule by rule
for i=1:no_rules
    figure(i+3)
    ti=sprintf('RULE %d: Blue line: Initial MF,   Red line: Final MF',i)
    title(ti)
    hold on
    plot(datumtel(:,1),init_memb(:,i),datumtel(:,1),final_memb(:,i),'r');
    axis([0 360 0 1]);
end

figure(no_rules+4)
load derror.dat
plot(derror(:,1),derror(:,2));
ti=sprintf('Evolution of the error curve')
title(ti)

figure(no_rules+5)
load derror.dat
plot(derror(:,1),derror(:,3));
ti=sprintf('Evolution of dpp')
title(ti)