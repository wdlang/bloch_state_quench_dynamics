% density at some point, a free fermion gas, 2016.sep.16
clear all; close all; clc; tic; 

pot = 1;
L = 120;
N = 2*L + 1;
location = 11;
xlist = (-L:L)';
qf = pi/2;
dt = 0.1;
Tmax = 650;
tlist = 0: dt : Tmax;
denlist = zeros(1, length(tlist));
overlap = zeros(1, length(tlist));
dq = 2*pi/N;
kf = floor(qf/dq);
klist = -kf: kf;
Nparticle = 2*kf + 1; 
blochs = zeros(N, length(klist));
for s1 = 1: length(klist)
    blochs (:,s1) = exp(i*2*pi*klist(s1)*xlist/N)/sqrt(N);
end
blochs2 = blochs';

H = zeros(N, N);
for s1 = 1: N
    H(s1, mod(s1, N)+1) = -1;
    H(mod(s1, N)+1, s1) = -1;
end
H(1,1) = pot;

[VV, DD] = eig(H);
dd = diag(DD);
waves0 = VV'*blochs;
for s1 = 1: length(tlist)
    time = tlist(s1);
    waves = VV*((exp(-i*time*dd)*ones(1,Nparticle)).*waves0);
    gram = blochs2* waves;
    overlap(s1) = abs (det (gram));
    
    denlist(s1) = sum(abs(waves(L+1+location, :)).^2);
end

h1= figure;
plot(tlist, overlap)

h2 = figure;
plot(tlist, denlist)