function [dx] = func_viga(t,x)

global A m M3 Amp freq
%excita��o externa harm�nica


f = Amp*cos(2*pi*freq*t); %for�a harm�nica

%indice da matriz onde ser� realizado a adi��o


%%
F = zeros(m,1);

F_ID = [zeros(m,1);inv(M3)*F];

dx = A*x + F_ID;
