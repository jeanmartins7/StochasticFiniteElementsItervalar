function [dx] = func_viga(t,x)

global A m M3 Amp freq
%excitação externa harmônica


f = Amp*cos(2*pi*freq*t); %força harmônica

%indice da matriz onde será realizado a adição


%%
F = zeros(m,1);

F_ID = [zeros(m,1);inv(M3)*F];

dx = A*x + F_ID;
