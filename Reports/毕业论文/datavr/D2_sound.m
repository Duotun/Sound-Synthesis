
x=0:50:400;
y2=[7,6,5,5,7,4,4,3,1];
y3=[8,6,7,7,6,8,4,8,5]; % bodong qiguai
y4=[7,7,8,6,6,5,9,7,6]  ;           %bi jinyou gun dong faner chale
plot(x,y2,x,y3,x,y4);
title('2D-Scene');
xlabel('Time delay(ms)');
ylabel('Rating');
legend('No Sound','Only Collision Sound','Collision Sound+Rolling Sound');
legend('location','best');
axis([0,400,0,10]);