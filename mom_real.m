% quench dynamics of the bloch state both in the mom space and in the real
% space; 2016.may.09 
clear all; close all; clc; myfont = 22;

L = 400;   N = 2*L+1;
ki = 200;
U = 1;
location = 4;
qi = 2*pi*ki/N;
deltaq = 2*pi/N;
deltaE = 2*sin(qi)*deltaq;
g = U/N;
T = 2*pi/deltaE;

steps = 1000;
loop = 15.5;
dt = T/steps;
tlist = dt*(0:steps*loop);
denlist = zeros(1, length(tlist));
plist = zeros(2, length(tlist));

xlist = -L:L;
xlist = xlist';
psii = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);
psif = (1/sqrt(N))*exp(i*(-2*pi*ki/N)*xlist);
psiif = [psii, psif]';

H = zeros(N, N);
for s= 1:(N-1)
    H(s,s+1) = -1;     H(s+1,s) = -1;
end
H(1,N) = -1;  H(N,1) = -1;
H(L+1, L+1) = U;

[VV,DD] = eig(H);
dd = diag(DD);

psi1 = VV'*psii;
for s = 1: length(tlist)
    time = tlist(s);
    psi = VV*(exp(-i*time*dd).*psi1);   
    
    plist(:,s) = abs(psiif*psi).^2;
    denlist(s) = abs(psi(L+1 + location ))^2;
end
denlist = N*denlist;

h1 = figure;
plot(tlist, denlist)
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2)
% set(gca, 'position', [0.15  0.15  0.8  0.8] )
set(gca, 'fontsize', myfont)
xlim([0 max(tlist)])
% ylim([0 1])
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2 ,':')
% plot(dt*(0:Tmax), plist2, dt*(0:Tmax), plist_ana,'--')
xlabel('$t$','fontsize',myfont,'Interpreter','latex');
ylabel('$ |\langle n | \Psi(t ) \rangle |^2  $','fontsize',myfont,'Interpreter','Latex')
str = 'real.eps';
print(h1,'-depsc',str)

h2 = figure;
plot(tlist, plist(1,:), tlist, plist(2,:),'--')
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2)
% set(gca, 'position', [0.15  0.15  0.8  0.8] )
set(gca, 'fontsize', myfont)
xlim([0 max(tlist)])
ylim([0 1])
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2 ,':')
% plot(dt*(0:Tmax), plist2, dt*(0:Tmax), plist_ana,'--')
xlabel('$t$','fontsize',myfont,'Interpreter','latex');
ylabel('$ P_i \;\& \; P_r $','fontsize',myfont,'Interpreter','Latex')
str = 'mom.eps';
print(h1,'-depsc',str)