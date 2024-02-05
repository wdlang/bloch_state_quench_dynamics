% to plot the spectrum of the tight binding model 
clear all; close all; clc; 

x1 = (-1.5:0.01:1.5)*pi;
y1 = zeros(1, length(x1));

y2 = -6.5:0.01:4.5;
x2 = zeros(1, length(y2));

x3 = pi*(-1:0.01:1);
y3 = -2 * cos(x3);

x4 = pi*(-1:0.1:1);
y4 = -2*cos(x4);

y5 = 0:0.01:2;
x5 = -pi*ones(1, length(y5));

y6 = 0:0.01:2;
x6 = pi*ones(1, length(y5));

x7 = pi*(-0.3:0.1:1.1);
y7 = 2*(x7-pi/2);

x8 = -pi*(-0.3:0.1:1.1);
y8 = -2*(x8+pi/2);

x9 = pi*(-0.6:0.01: -0.35);
y9 = 2*(x9 - pi/2);

x10 = pi*(1.15: 0.01: 1.4);
y10 = 2*(x10 - pi/2);

x11 = -pi*(-0.6:0.01: -0.35);
y11 = -2*(x11 + pi/2);

x12 = -pi*(1.15: 0.01: 1.4);
y12 = -2*(x12 + pi/2);

h1 = figure;
plot (x1, y1, 'linewidth', 2)
hold on 
plot (x2, y2, 'linewidth', 2)
plot (x3, y3, 'linewidth', 2, 'color', 'r')
plot(x4, y4, '.','markersize',15, 'color', 'r')

plot (x5, y5, ':', 'linewidth', 2, 'color', 'r')
plot (x6, y6,':', 'linewidth', 2, 'color', 'r')

plot (x7, y7,'o', 'linewidth', 2, 'color', 'b')
plot (x8, y8,'o', 'linewidth', 2, 'color', 'b')

plot(x9, y9, ':', 'linewidth', 2,  'color', 'b' )
plot(x10, y10, ':', 'linewidth', 2,  'color', 'b' )
plot(x11, y11, ':', 'linewidth', 2,  'color', 'b' )
plot(x12, y12, ':', 'linewidth', 2,  'color', 'b' )
% 
% r= 0.9;
% plot(pi/2+r*cos(pi*(0:0.01:2)) , r*sin(pi*(0:0.01:2)), ':', 'linewidth', 2)
% plot(-pi/2+r*cos(pi*(0:0.01:2)) , r*sin(pi*(0:0.01:2)), ':', 'linewidth', 2)

arrow([-1.5*pi 0],[pi*1.5 0],'width',1.5,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')
arrow([0,-2.5],[ 0,4.5 ],'width',1.5,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')

text(-1.15*pi ,-0.55 ,'$- \pi $','fontsize',22,'Interpreter','latex')
text(0.85 *pi ,-0.55 ,'$ + \pi $','fontsize',22,'Interpreter','latex')
text(1.4*pi ,-0.7 ,'$ q $','fontsize',22,'Interpreter','latex')
text(0.1*pi ,4.3, '$ \varepsilon $','fontsize',24,'Interpreter','latex')
text(0.2*pi ,0.55 ,'$ + q_i $','fontsize',22,'Interpreter','latex')
text(-0.50*pi ,0.55 ,'$ -q_i $','fontsize',22,'Interpreter','latex')

text(-0.85*pi ,-2.55*pi ,'$ |R_n\rangle  $','fontsize',22,'Interpreter','latex')
text(0.55*pi ,-2.55*pi ,'$ |L_n\rangle  $','fontsize',22,'Interpreter','latex')

% text(0.1*pi ,2 ,'$ 2 $','fontsize',22,'Interpreter','latex')

axis off 
% axis equal 

str = strcat('luttinger.eps');
print(h1,'-depsc',str)