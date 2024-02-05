% to have a look of the function \sum_{n=1}^\infty sin (n \theta )/ n 
clear all; close all; clc; myfont = 22;

philist = pi* (0:0.001:4);

N = 10000;

flist = zeros(1, length(philist));

for s1 = 1: length(philist)
    phi = philist(s1);
    flist(s1) =  sum( sin(phi*(1:N))./(1:N));
end

h1 = figure;
plot(philist/pi  , flist/pi )
xlabel('$\theta/\pi$', 'fontsize', myfont,'Interpreter','latex')
ylabel('$summation $', 'fontsize', myfont,'Interpreter','latex')
set(gca,'fontsize', myfont)