spring=readtable('figure3a.csv');
spring=sortrows(spring,'x');

autumn=readtable('figure3b.csv');
autumn=sortrows(autumn,'x');

dx=round(mean(diff(autumn.x)),2);
dx=round(mean(diff(spring.x)),2);
dx=.285;

x=6-dx/2+(1:numel(spring.x))'*dx;

y_sp=spring.y;
y_au=autumn.y;

figure; hold on
bar(x,y_sp)
plot(spring.x,spring.y,'.')

figure; hold on
bar(x,y_au)
plot(autumn.x,autumn.y,'.')


sum(x.*y_sp)/sum(y_sp)
sum(x.*y_au)/sum(y_au)
