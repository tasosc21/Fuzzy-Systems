function my_display_2(plot_no,x,y,tlt)

subplot(2,1,plot_no);
plot(x,y);
axis([min(x) max(x) 0 1.05]);
grid on;
title(tlt);


