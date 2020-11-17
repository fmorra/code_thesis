% Try to draw the iteration scheme 
clear all; close all; clc;
iter = linspace(0,3,4);
iter_y = 4*ones(1, length(iter));
step = linspace(0,3,12);
step_y = 2*ones(1, length(step));
cycle = linspace(0,3,4);
cycle_y = 0*ones(1, length(cycle));
op = linspace(0,3,6);
op_y = -2*ones(1, length(op));
figure()
hold on;

annotation('textbox',[.1 .82 .1 .2],'string','Job submission','EdgeColor','none');
annotation('textbox',[.28 .82 .1 .2],'string','Sph\_tools','EdgeColor','none');
annotation('textbox',[.45 .82 .1 .2],'string','New loads','EdgeColor','none');
xa = [.13 .13];
ya = [.83 .93];
annotation('doublearrow',xa,ya)
xa = [.17 .2];
ya = [.83 .93];
annotation('doublearrow',xa,ya)
xa = [.3 .26];
ya = [.83 .93];
annotation('doublearrow',xa,ya)
xa = [.38 .55];
ya = [.83 .93];
annotation('doublearrow',xa,ya)
plot(iter,iter_y,'k.-','markersize',12);
for i=1:(length(iter)-1)
    text(iter(i)+(iter(i+1)-iter(i))/2,4.5,num2str(i))
end
annotation('textbox',[.0001 .7 .1 .2],'string','Iteration','EdgeColor','none');

plot(step,step_y,'k.-','markersize',12);
xa = [.13 .17];
ya = [.65 .83];
annotation('doublearrow',xa,ya)
xa = [.9 .3];
ya = [.65 .83];
annotation('doublearrow',xa,ya)
for i=1:(length(step)-1)
    text(step(i)+(step(i+1)-step(i))/2,2.5,num2str(i))
end
annotation('textbox',[.05 .5 .1 .2],'string','Step','EdgeColor','none');

plot(cycle,cycle_y,'k.-','markersize',12);
xa = [.13 .13];
ya = [.45 .63];
annotation('doublearrow',xa,ya)
xa = [.9 .2];
ya = [.45 .63];
annotation('doublearrow',xa,ya)
for i=1:(length(cycle)-1)
    text(cycle(i)+(cycle(i+1)-cycle(i))/2,0.5,num2str(i))
end
annotation('textbox',[.05 .3 .1 .2],'string','Cycle','EdgeColor','none');

plot(op,op_y,'k.-','markersize',12);
xa = [.13 .13];
ya = [.25 .43];
annotation('doublearrow',xa,ya)
xa = [.9 .4];
ya = [.25 .43];
annotation('doublearrow',xa,ya)
for i=1:(length(op)-1)
    %text(op(i)+(op(i+1)-op(i))/2,-1.5,num2str(i))
    pos_vec(i) = op(i)+(op(i+1)-op(i))/2;
end
annotation('textbox',[.0001 .1 .1 .2],'string','Operations','EdgeColor','none');
text(pos_vec(1)-0.3,-2.5,'Reports generated')
text(pos_vec(2)-0.3,-1.5,'Diff calculated')
text(pos_vec(3)-0.2,-2.5,'csvs saved')
text(pos_vec(4)-0.3,-1.5,'stress coupling')
text(pos_vec(5)-0.3,-2.5,'job submission')

yline(-3.5,'color','w');
axis([0 3 -3.5 5])
%axes('Color','none','XColor','none');
set(gca,'ycolor','w')
set(gca,'xcolor','w')
axis off
hold off;
saveas(gcf, 'C:\Users\fabri\Desktop\TU Delft\iter_scheme.png');

