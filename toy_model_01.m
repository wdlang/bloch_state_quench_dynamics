% reflection of a bloch state by a barrier; real space; looking for a
% better quantity suitable for experimental observation
% 2016.02.17
clear all; close all; clc; myfont = 22;

L = 100;   N = 2*L+1 ;
delta = 4/L;
ki = 50;
U = 0.4;
location = 6;
% cutoff = 50; 
dt = 0.1;     Tmax = 9000;
plist = zeros(1, 1+Tmax);
plist2 = zeros(1, 1+Tmax);
plist_ana = zeros(1, 1+Tmax);

xlist = -L:L;
xlist = xlist';
basis = exp(i*2*pi*xlist*xlist'/N)/sqrt(N);
% blochgroup1 = zeros(N,2*cutoff+1);
% blochgroup2 = zeros(N,2*cutoff+1);
% for s1 = -cutoff : cutoff
%     blochgroup1(:,s1 + cutoff + 1) = exp(i*2*pi*(ki+s1)/N*xlist)/sqrt(N);
%     blochgroup2(:,s1 + cutoff + 1) = exp(i*2*pi*(-ki+s1)/N*xlist)/sqrt(N);
% end
% blochgroup = [blochgroup1, blochgroup2 ];
% blochgroup2 = blochgroup';

H = zeros(N, N);
for s= -L : L
    H(s+L+1,s+L+1) = max ( delta* (s-L/2), -delta*(s+L/2));
end
H = H + (U/N)* ones(N ,N);

[VV,DD] = eig(H);
dd = diag(DD);

psii = zeros(N, 1);
psii(ki+L+1) = 1 ;
psi1 = VV'*psii;
% g = U/N;
% Delta = 4*pi*sin(2*pi*ki/N)/N;
% T = 2*pi/Delta;
% rotation = (1-i*g*T)/(1+i*g*T);

for s = 1:Tmax
    psi = basis*(VV*(exp(-i*dt*s*dd).*psi1));   
    plist(s+1) = abs(psi(L+1 + location ))^2;
    
%     proj = blochgroup2*psi;
%     1 - norm(proj)^2;
%     amp0 = blochgroup(L + 1 + location,:)*proj;
%     plist2(s+1) = N*abs(amp0)^2;
%     
%     p = floor(s*dt/T);
%     tdiff = s*dt - p*T;
%     amp0 = (1 -i*2*g*tdiff/(1+ i*g*T))*(rotation^p);
%     
%     amp = i*sin(2*pi*ki/N*location) + cos(2*pi*ki/N*location) * amp0 ;
%     for ss = 1:cutoff
%         amp = amp + cos(2*pi*(ki+ss)/N*location)* (2*g/ss/Delta )/(1+i *g*T)*(1 - exp(-i*ss*Delta*tdiff))*(rotation^p);
%         amp = amp + cos(2*pi*(ki-ss)/N*location)* (-2*g/ss/Delta )/(1+i *g*T)*(1 - exp(i*ss*Delta*tdiff))*(rotation^p);
%     end
%     
%     plist_ana(s+1) = abs(amp)^2 ;
end
plist = N*plist;

h1 = figure;
plot(dt*(0:Tmax), plist)
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2)
% set(gca, 'position', [0.15  0.15  0.8  0.8] )
set(gca, 'fontsize', myfont)
% ylim([0 3])
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2 ,':')
% plot(dt*(0:Tmax), plist2, dt*(0:Tmax), plist_ana,'--')
xlabel('$t$','fontsize',myfont,'Interpreter','latex');
ylabel('$|\psi_n|^2$','fontsize',myfont,'Interpreter','latex');
str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki),', n=',num2str(location));
title(str,'fontsize',myfont)
str = strcat('toy_U=', num2str(U),'_N=',num2str(N),'_ki=',num2str(ki),'_n=',num2str(location),'.jpg');
% print(h1,'-djpeg',str)