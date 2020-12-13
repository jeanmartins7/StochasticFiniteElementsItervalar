%Start script
%%
%modeling Euler-Bernoulli beam in MFE
global A Amp m h L M3 NDE freq
%Parametros
N = 5;%Nelementos
E = 7e10;%modulo de young
Le = 5; %comprimento da barra
rho = 2700;%densidade da barra
h = 0.5;%comprimento
b = 0.5;%largura
A = b*h;%area
I = (b*h*(b^2 + h^2))/24;%inercia

Amp = 0; %N
freq = 5; %Hz

%Obtendo as matrizez de massa e rigidez
[K3,M3] = beam_ef_creating(N,rho,E,A,I,Le);
alfa=1e-2;
beta=1e-5;
C3 = 0.00062*K3;

[m,n]=size(M3);
A = [zeros(m) eye(m);
    -inv(M3)*K3 -inv(M3)*C3];

tspan = [1e-5:0.01:3];
x0 = zeros(2*m,1);
x0((5),1) = -4.7756e-03;
x0((9),1) = -1.7366e-02;
x0((13),1) = -3.5166e-02;
x0((17),1) = -5.5570e-02;


options=odeset('OutputFcn',@odeprog,'Events',@odeabort);
[t,S] = ode23(@func_viga, tspan,x0);
save S;
%%
load S
figure(1)
plot(t,S(:,(13)));
xlabel('time');
ylabel('dx(t)');
title('Deslocamento em X(t) na ponta da viga')

figure(2)
y = S(1,1:4:end);
plot(1:length(y),y,'-or')
hold on
plot([0 1],[0 y(1)],'-or')
title('Deflexão inicial da viga')

figure(3)
for ii = 1:length(t)
    
ani = subplot(1,1,1);
y = S(ii,1:4:end);


plot(1:length(y),y,'-or')
hold on
plot([0 1],[0 y(1)],'-or')

title('Simulação da viga')
axis(ani,'equal');
set(gca,'XLim',[0 5],'YLim',[-0.8 .80]);
str2 = ['Time elapsed: '  num2str(t(ii)) ' s'];
Time = text(2,1.5,str2);
pause(1e-100);
clf

   
end


