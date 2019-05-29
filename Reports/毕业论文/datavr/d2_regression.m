x=0:50:400;
ydata11=[10,10,8.6,6.6,7.2,4.1,5.7,3.7,2.4;
    10,8.7,5.6,5.9,5.7,4.7,5,2,6.5;
    0,0,0,1.4,0.27,0.35,0,0,0;
    0,2.3,0,0,0,0,0,0,0
    8.4,6,9.8,7.5,6.6,7,6.6,6.3,5
    8.7,9,7.7,7.9,6,7.5,5,7.2,7]; %no sound

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
plot(x,ydata3,'.','color','red')
axis([-100,500,0,10]);
title({'\fontsize{12}2D-Scene';'\fontsize{12}Full Sound'})
xlabel('\fontsize{12}Time Delay(ms)')
ylabel('\fontsize{12}Rating')
x1=x';
y1=ydata3';
x2=[ones(9,1),x1];
[b,bint,r,rint,stats]=regress(y1,x2,1);
y=b(2)*x+b(1);
hold on,p1=plot(x,y,'r','color','green');
legend(p1,'\fontsize{10}slope= -0.0139 , error=0.4463');
%str={'slope= -0.0131 , error=0.6'};
dim=[.2 .2 .0 .0];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');