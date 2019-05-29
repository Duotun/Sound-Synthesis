y=[0 0 1;0 1 1; 1 1 1;1 1 1;1 1 1;1 1 1;1 1 1;1 1 1] %2d - distance and delay
[p,tble,stats] = anova2(y,4);
x=[1 0 1;1 1 1;1 1 1;0 0 0;1 1 1;1 1 1;0 0 0;1 0 1; 1 1 1;1 0 1;
   0 1 1;1 1 1 ;0 0 1;1 0 1; 1 1 1;0 1 1;1 0 1; 1 1 1]%3d mass ratiao a
                                                      %nd time delay
[p,tble,stats] = anova2(x,3);
z=[0 0 1;0 0 1;1 1 1;1 0 0;0 1 1;1 1 1;0 0 1;1 0 0;1 1 1;0 0 1 ;0 0 1;
    1 1 1;0 0 1;0 0 1;1 1 1 ]     %3d distance and time delay
[p,tble,stats] = anova2(z,3);

q=[1 0 1;0 0 0;0 0 0;0 0 0;1 0 1;0 0 1;0 1 1;0 1 0]
[p,tble,stats] = anova2(q,2);