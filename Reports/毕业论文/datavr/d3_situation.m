xdata=0:50:400;
ydata11=[10,10,8.6,6.6,7.2,4.1,5.7,3.7,2.4;
    10,8.7,5.6,5.9,5.7,4.7,5,2,6.5;
    0,0,0,1.4,0.27,0.35,0,0,0;
    0,2.3,0,0,0,0,0,0,0
    8.4,6,9.8,7.5,6.6,7,6.6,6.3,5
    8.7,9,7.7,7.9,6,7.5,5,7.2,7
     ]; %no sound
 
ydata1=mean(ydata11);

ydata21=[10,5.8,8.3,6.7,6.9,4.7,3.8,6.2,3.3;
    10,10,9.7,7.7,3.3,1.2,3.5,1.6,6
    9.5,8.3,7.6,4.4,4.1,7.1,6.6,2.6,2.7
    9.8,9.5,7.9,6,4.9,3.9,4.3,3.7,3.9
    9,10,9.6,8.4,7,7,7.1,7.1,6.9
    9,9,9.4,7.1,7.3,7.9,5.4,6.9,7]; %collision sound    
ydata2=mean(ydata21);
ydata31=[10,10,10,8.6,8.2,2.6,4.3,0,2.9;
    10,10,7.1,8.3,5.3,4.9,5.5,3.6,4.9
    10,7.7,8.6,8,6.4,3.9,3.3,5.6,1.6
    10,9.4,6.3,8.1,7.1,8.7,4.4,3.5,5.7
    10,8.4,8.1,7.9,6.6,7.3,7.2,5.3,7.1
    8.5,9.4,9.8,7.8,7.1,6.9,8,7.7,7.8];%both sound
ydata3=mean(ydata31); 
ydata4=[8.1,7.8,7.2,6.1,5.3,5.8,4.5,4.4,4];      %+sound
ydata5=[8,6.9,5.1,3.8,3.7,2.4,2,2.4,2.2]        %control
y1=[std(ydata1)/10
    std(ydata1)/10
    std(ydata1)/10
    std(ydata1)/10
    std(ydata1)/10
    std(ydata1)/10
    std(ydata1)/10
    std(ydata1)/10
    std(ydata1)/10];
y2=[std(ydata2)/10
    std(ydata2)/10
    std(ydata2)/10
    std(ydata2)/10
    std(ydata2)/10
    std(ydata2)/10
    std(ydata2)/10
    std(ydata2)/10
    std(ydata2)/10];
y3=[std(ydata3)/10
    std(ydata3)/10
    std(ydata3)/10
    std(ydata3)/10
    std(ydata3)/10
    std(ydata3)/10
    std(ydata3)/10
    std(ydata3)/10
    std(ydata3)/10];
y4=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
y5=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
a=errorbar(xdata,ydata1,y1,y1);hold on;
b=errorbar(xdata,ydata2,y2,y2); hold on;
c=errorbar(xdata,ydata3,y3,y3); hold on;
d=errorbar(xdata,ydata4,y4,y4);hold on;
e=errorbar(xdata,ydata5,y5,y5);
pl=legend('No Sound-Experiment','Only Collision Sound-Expeirment','Both Sound-Experiment','Sound situation','Control Situation');
pl.FontSize=14;
%a.Marker='o';
a.MarkerSize=3;
a.CapSize=3;
%b.Marker='o';
b.MarkerSize=3;
b.CapSize=3;
%c.Marker='o';
c.MarkerSize=3;
c.CapSize=3;
%d.Marker='o';
d.MarkerSize=3;
d.CapSize=3;
%e.Marker='o';
e.MarkerSize=3;
e.CapSize=3;
title('\fontsize{16}3D-Scene');
xlabel('\fontsize{16}Time delay(ms)');
ylabel('\fontsize{16}Rating£¨0-10£©');
axis([-50,450,0,15]);