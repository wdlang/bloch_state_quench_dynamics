% to plot the function beta
clear all; close all; clc; tic; myfont = 22;

a1 = 5; a2 = 7;
b1 = 1.25; b2 =1.75;

h1 = figure ;
hold on
x1 = [-4, -2 ,0, 2];
y1 = [1, 1, 1, 1];
plot(x1, y1, 'o')

x2 = [ -2 ,0, 2, 4];
y2 = -[1, 1, 1, 1];
plot(x2, y2, 'o')

x4 = [-4, -2 ,0, 2, 4];
y4 = 0*[1, 1, 1, 1, 1];
plot(x4, y4, '.', 'markersize',15)

hold on 
for s1 = 1: 4
    x3 = x1(s1)+0.02 : 0.01: x2(s1)-0.02;
    y3 = y1(s1)-0.02 : -0.01: y2(s1)+0.02;
    plot(x3, y3, 'linewidth', 2)
end

arrow([0 -b1],[0 b2],'width',2,'TipAngle',15,'BaseAngle',20,'FaceColor','b','EdgeColor','b')
arrow([-a1 0],[a2 0],'width',2,'TipAngle',15,'BaseAngle',20,'FaceColor','b','EdgeColor','b')

text(2-0.2 ,-0.2 ,'$ 2 $','fontsize',22,'Interpreter','latex')
text(4-0.2 ,-0.2 ,'$ 4 $','fontsize',22,'Interpreter','latex')
text(0.1 ,-0.2 ,'$ 0 $','fontsize',22,'Interpreter','latex')
text(-2-0.4 ,-0.2 ,'$ -2 $','fontsize',22,'Interpreter','latex')
text(-4-0.4 ,-0.2 ,'$ -4 $','fontsize',22,'Interpreter','latex')

text(0.3 ,1 ,'$ +1 $','fontsize',22,'Interpreter','latex')
text(0.3 ,-1 ,'$ -1 $','fontsize',22,'Interpreter','latex')
text(0.3 ,b2-0.2  ,'$ \beta /\pi  $','fontsize',22,'Interpreter','latex')
text(5.6 ,-0.2 ,'$ z/ \pi $','fontsize',22,'Interpreter','latex')

box off
axis off 

str = 'beta.eps';
print(h1, '-depsc', str)

h2 = figure;
hold on
a1 = -0.3; a2 = 1.3;
b1 = -2.7*pi, b2 = 2.7*pi;

inc = 0.1*pi;
a = 0.33;
x1 = 0: 0.01: a;
y1 = zeros(1, length(x1)) - inc;
plot(x1, y1, 'linewidth', 2)

x1 = a: 0.01: 1;
y1 = -2*pi*ones(1, length(x1)) - inc;
plot(x1, y1, 'linewidth', 2)

x1 = a * ones(1, 100);
y1 = -2*pi- inc  + 2*pi/100* (1:100);
plot(x1, y1, ':', 'linewidth', 2)

b = 1 - a;
x1 = 0: 0.01: b;
y1 = zeros(1, length(x1)) + inc;
plot(x1, y1, 'r', 'linewidth', 2)

x1 = b: 0.01: 1;
y1 = -2*pi*ones(1, length(x1)) + inc;
plot(x1, y1, 'r','linewidth', 2)

x1 = b * ones(1, 100);
y1 = -2*pi +  inc  + 2*pi/100* (1:100);
plot(x1, y1, ':r', 'linewidth', 2)

x1 = -0.02:0.01:0.02;
y1 = -2*pi*ones(1, length(x1));
plot(x1, y1, 'linewidth', 2)

y1 = (-0.07:0.01:0.07)*pi;
x1 = 1*ones(1, length(y1));
plot(x1, y1, 'linewidth', 2)

text(-0.08 , -0.3*pi  ,'$ 0 $','fontsize',22,'Interpreter','latex')
text(1.25 , -0.3*pi  ,'$ s $','fontsize',22,'Interpreter','latex')
text(0.05 , 2.65*pi  ,'$ \gamma  $','fontsize',22,'Interpreter','latex')
text(-0.3 , -2*pi  ,'$ -2 \pi   $','fontsize',22,'Interpreter','latex')
text(0.8 , -2.4*pi  ,'$ \gamma_+    $','fontsize',22,'Interpreter','latex')
text(0.8 , -1.6*pi  ,'$ \gamma_-    $','fontsize',22,'Interpreter','latex')
text(a -0.05 , 0.4*pi  ,'$ s_c    $','fontsize',22,'Interpreter','latex')
text(b-  0.18 , 0.4*pi  ,'$ T - s_c    $','fontsize',22,'Interpreter','latex')
text(0.97 , 0.4*pi  ,'$ T    $','fontsize',22,'Interpreter','latex')

arrow([0 b1],[0 b2],'width',2,'TipAngle',15,'BaseAngle',20,'FaceColor','b','EdgeColor','b')
arrow([a1 0],[a2 0],'width',2,'TipAngle',15,'BaseAngle',20,'FaceColor','b','EdgeColor','b')

axis off 
box off

print(h2, '-depsc', 'gamma.eps')