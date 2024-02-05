% the idealized toy model;2016.may.08
clear all; close all; clc; tic; myfont = 22;

U = 0.5;
N = 200;                    % lattice size
M = 400;                   % 2M+1 fictious levels
qi = pi/2;
deltaE = 4*sin(qi)/N;
deltaq = 2*pi/N;
Elist = deltaE*(-M:M)';
g = U/N;
T = 2*pi/deltaE;
location = 22;
theta = 2*atan(g*T);

steps = 1000;
loop = 6;
dt = T/steps;
tlist = dt*(0:steps*loop);
amplist = zeros(1, length(tlist));
denlist = zeros(1, length(tlist));
amplist2 = zeros(1, length(tlist));
denlist2 = zeros(1, length(tlist));

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
    
    f1 = -0.5* ( (deltaq*location - deltaE*time2) - floor((deltaq*location - deltaE*time2)/2/pi)*2*pi - pi) ;
    f2 = -0.5* (deltaq*location - floor(deltaq*location /2/pi)*2*pi - pi);
    f3 = -0.5* ( (-deltaq*location - deltaE*time2) - floor((-deltaq*location - deltaE*time2)/2/pi)*2*pi - pi) ;
    f4 = -0.5* (-deltaq*location - floor(-deltaq*location /2/pi)*2*pi - pi);
    amplist2(s1) = i* sin(qi*location)+exp(-i*theta*pp)* ( cos(qi*location)+ g/deltaE/(1+i*g*T)* ( exp(i*qi*location) * (-i*deltaE*time2 +2*i* (f1 - f2)) ...
        + exp(-i*qi*location) * (-i*deltaE*time2 +2*i* (f3 - f4)) )) ;
    denlist2(s1) = abs(amplist2(s1))^2;
end

h1 = figure;
plot(tlist/T, denlist, tlist/T, denlist2 )
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
