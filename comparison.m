% compare the ideal model and the realistic TBM;2016.may.08
clear all; close all; clc; tic; myfont = 22;

U = 1;
L = 100;   N = 2*L+1;
ki = 50;
M = N*2;                    % 2M+1 fictious levels
qi = 2*pi*ki/N;
deltaq = 2*pi/N;
deltaE = 2*sin(qi)*deltaq;
Elist = deltaE*(-M:M)';
g = U/N;
T = 2*pi/deltaE;
location = 11;
theta = 2*atan(g*T);
v = deltaE/ deltaq;
T1 = location/ v;
T2 = T - T1;

steps = 1000;
loop = 36;
dt = T/steps;
tlist = dt*(0:steps*loop);
amplist = zeros(1, length(tlist));
denlist = zeros(1, length(tlist));
amplist2 = zeros(1, length(tlist));
denlist2 = zeros(1, length(tlist));
amplist3 = zeros(1, length(tlist));
denlist3 = zeros(1, length(tlist));

xlist = -L:L;
xlist = xlist';
psi0 = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);
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
    time = tlist(s1);
    psi = VV*(exp(-i*time*dd).*psi1);   
    amplist3(s1) = psi(L+1 + location) ; 
    denlist3(s1) = abs(psi(L+1 + location ))^2;
end
denlist3 = N*denlist3;

Anlist = zeros(2*M + 1, 1);
floclist = zeros(1, 2*M + 1);
for s1 = -M : M
    floclist (s1 + M + 1) =  cos((qi + s1* deltaq) * location );
end

for s1 = 1: length(tlist)
    time = tlist(s1);
    pp = floor(time / T);
    time2 = time - pp*T;
    
    Anlist(M+1) = 1 - ( i*2*g*time2/(1+ i*g*T) );
    Anlist(1:M) = (2*g/(1+i*g*T))*( (exp(-i*time2*Elist(1:M))-1)./ Elist(1:M) );
    Anlist(M+2:2*M+1) = (2*g/(1+i*g*T))*( (exp(-i*time2*Elist(M+2:2*M+1))-1)./ Elist(M+2:2*M+1) );
    Anlist = exp(-i*theta*pp)*Anlist;
    
    amp = floclist*Anlist + i* sin(qi*location);
    %  amp = floclist*Anlist ;
    amplist(s1) = amp ;
    denlist(s1) = abs(amp)^2 ;
    
%     f1 = -( (deltaq*location - deltaE*time2) - floor((deltaq*location - deltaE*time2)/2/pi)*2*pi - pi) ;
%     f2 = -(deltaq*location - floor(deltaq*location /2/pi)*2*pi - pi);
%     f3 = - ( (-deltaq*location - deltaE*time2) - floor((-deltaq*location - deltaE*time2)/2/pi)*2*pi - pi) ;
%     f4 = - (-deltaq*location - floor(-deltaq*location /2/pi)*2*pi - pi);
%     amplist2(s1) = i* sin(qi*location)+exp(-i*theta*pp)* ( cos(qi*location)+ i*g/deltaE/(1+i*g*T)* ( exp(i*qi*location) * (-deltaE*time2 + (f1 - f2)) ...
%         + exp(-i*qi*location) * (-deltaE*time2 + (f3 - f4)) )) ;
    if time2 < T1 
        gammap = 0;
    else
        gammap = - 2*pi;
    end
        if time2 < T2 
        gammam = 0;
    else
        gammam = - 2*pi;
        end
    amplist2(s1) = i* sin(qi*location)+exp(-i*theta*pp)* ( cos(qi*location)+ i*g/deltaE/(1+i*g*T)* ( exp(i*qi*location) * gammap  ...
        + exp(-i*qi*location) * gammam )) ;
    denlist2(s1) = abs(amplist2(s1))^2;
end

h1 = figure;
plot(tlist/T, denlist, tlist/T, denlist2, tlist/T, denlist3,':' )
xlim([0 loop])
xlabel('t/T')
ylabel('$|\psi_n|^2$','fontsize',myfont,'Interpreter','latex');
str = strcat ('U=', num2str(U),', N=',num2str(N),', M=',num2str(M),', qi/\pi=',num2str(qi/pi),', n=',num2str(location));
title(str,'fontsize',myfont)
str = strcat('U=', num2str(U),'_N=',num2str(N),'_M=',num2str(M),'_qi2Pi=',num2str(qi/pi),'_n=',num2str(location),'.jpg');
print(h1,'-djpeg',str)

h2 = figure;
plot3(tlist/T, real(amplist), imag(amplist))
hold on 
plot3(tlist/T, real(amplist2), imag(amplist2))