%%% PROGRAMA QUE TRANSFORMA ESPAÇO DE ESTADO PARA O DOMINIO DO TEMPO
%%% JEAN MARCEL PERES MARTINS
%%% ORIENTADOR: FABIAN ANDRES LARA MOLINA


function y_out=sys_model(ti,K3,M3)
global A Amp m freq

Amp = 0; %N
freq = 5; %Hz

C3 = 0.00062*K3;

[m,n]=size(M3);
A = [zeros(m) eye(m);
    -inv(M3)*K3 -inv(M3)*C3];

tspan = [1e-5:0.01:ti];
x0 = zeros(2*m,1);
x0((5),1) = -4.7756e-03;
x0((9),1) = -1.7366e-02;
x0((13),1) = -3.5166e-02;
x0((17),1) = -5.5570e-02;

[t,y] = ode15s(@func_viga, tspan,x0);

  y_out=y(end,(17));

end