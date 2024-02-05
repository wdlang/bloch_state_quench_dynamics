% snapshots of the wave function
clear all; close all; clc; myfont = 22;

L = 200;   N = 2*L+1;
ki = 100;
U = 0.5;
levitation = 2 ;
location = 5;
deltaq = 2*pi/N;
deltaE = 2*sin(2*pi*ki/N)*deltaq;
T = 2* pi/ deltaE; 
steps = 8 ;
tlist = T*(0:steps*2)/steps;
denlist = zeros(N, length(tlist));

xlist = -L:L;
xlist = xlist';
psi0 = exp(i*(2*pi*ki/N)*xlist);

H = zeros(N, N);
for s= 1:(N-1)
    H(s,s+1) = -1;     H(s+1,s) = -1;
end
H(1,N) = -1;  H(N,1) = -1;
H(L+1, L+1) = U;

[VV,DD] = eig(H);
dd = diag(DD);
psi1 = VV'*psi0;
for s1 = 1:length(tlist)
    psi = VV*(exp(-i*tlist(s1)*dd).*psi1);   
    denlist(:, s1) = abs(psi).^2 + levitation *(s1-1);
end

h1 = figure;
plot(xlist, denlist)
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2)
% set(gca, 'position', [0.15  0.15  0.8  0.8] )
set(gca, 'fontsize', myfont)
% ylim([0 3])
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2 ,':')
% plot(dt*(0:Tmax), plist2, dt*(0:Tmax), plist_ana,'--')
xlabel('$ n $','fontsize',myfont,'Interpreter','latex');
ylabel('$|\Phi_n |^2$','fontsize',myfont,'Interpreter','latex');
% str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki),', n=',num2str(location));
% title(str,'fontsize',myfont)
% str = strcat('U=', num2str(U),'_N=',num2str(N),'_ki=',num2str(ki),'_n=',num2str(location),'.jpg');
% print(h1,'-djpeg',str)

str = strcat('snapshots.eps');
print(h1, '-depsc',str)