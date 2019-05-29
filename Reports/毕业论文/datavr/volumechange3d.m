x=linspace(5,300);
y=exp(-(x-5)/50);
plot(x,y,'Color',[0.6,0.1,0.9]);
set(gca,'ticklength',[0,0]);
title('\fontsize{52}3D Scene');
xlabel('\fontsize{52}Distance From the Viewer(m)');
%set(gca,'fontsize',32);
ylabel('\fontsize{52}Volume');
legend('\fontsize{31}min distance=5m, max distance=300m');