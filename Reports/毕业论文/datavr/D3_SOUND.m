x=0:50:400;
y2=[6,6,4,4,4,2,5,2,6];   % buzhun you bug
y3=[8,8,8,7,5,1,7,5,3];
y4=[8,7,4,5,4,5,4,5,4];
%ax1=plot(x,y2);
%ax2=plot(x,y3);
%ax3=plot(x,y4);
plot(x,y2,x,y3,x,y4);
line([50,50],[6,8],'Marker','.');
title('3D-Scene');
xlabel('Time delay(ms)');
ylabel('Rating');
legend('No Sound','Only Collision Sound','Collision Sound+Rolling Sound');
legend('location','best');
axis([0,400,0,10]);