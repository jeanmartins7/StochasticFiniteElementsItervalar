%Start script

function intrval_solution
global A Amp m M3 freq N
clear all

clc

close all

global omega

%%
%Intervalos dos Parametros

%------------------------------

%  - rho

p_1_m=2700;

p_1_l=p_1_m-.05*p_1_m;

p_1_r=p_1_m+.05*p_1_m;

p_1_I=[p_1_l,p_1_r];
 

%Numero de intervalos

n_I=3;
    
 


%Options GA

options = gaoptimset('PopulationSize',n_I*10,'PopulationType','doubleVector');%,...

                    %'PopInitRange',[80.0e-3 180.0e-3 115 3 100.0e-3; 100.0e-3 200.0e-3 117 5 180.0e-3],...

                    %'PlotFcns',{@gaplot_GCI, @gaplotbestindiv, @gaplotdistance});

 

%--------------------------------------------------------------

% Analysis Intervalar

omega_f=1000;% Frequencia final

delta_omega=1;

omega_i=delta_omega;
iiii=1;
for omega=omega_i:delta_omega:omega_f           
    
        lb = [p_1_l];

        ub = [p_1_r];       
                

        %Interval Analysis for all Parameters

        [x_min,fval_min] = ga(@func1,n_I,[],[],[],[],lb,ub,[],options);

        [x_max,fval_max] = ga(@func1_max,n_I,[],[],[],[],lb,ub,[],options);

        y(iiii)=func1([p_1_m]);

        y_min(iiii)=fval_min;

        y_max(iiii)=-fval_max;

        p_min(iiii,:)=x_min;

        p_max(iiii,:)=x_max;

        

        iiii=iiii+1

end
%--------------------------------------------------------------

%--------------------------------------------------------------

 

%Plot of interval Analysis of all parameters

plot([delta_omega:delta_omega:omega_f],20*log10(y),'g');

hold on

plot([delta_omega:delta_omega:omega_f],20*log10(y_min),'b');

hold on

plot([delta_omega:delta_omega:omega_f],20*log10(y_max),'r');

title('Análise Intervalar','Interpreter','latex','fontsize',16)

ylabel('$|H(\omega)|$','Interpreter','latex','fontsize',22)

xlabel('$\omega$ [rad/s]','Interpreter','latex','fontsize',22)

box on

saveas(gcf,'interval_solution2gdl_massa2.tif');
saveas(gcf,'interval_solution2gdl_massa2.eps');


end

%-------------------------------------------------------------------------


 
 

 

%---------------------------------

 

function z=func1(GG) 
global A Amp m freq omega M3
    %modeling Euler-Bernoulli beam in MFE

    %Parametros
    N = 5;%Nelementos
    E = 7e10;%modulo de young
    Le = 5; %comprimento da barra
    h = 0.5;%comprimento
    b = 0.5;%largura
    Ae = b*h;%area
    I = (b*h*(b^2 + h^2))/24;%inercia
    

    %----------------------------------

    %Dominio da frequiencia
 

        rho = GG(1);
        [K3,M3] = beam_ef_creating(N,rho,E,Ae,I,Le);
        alfa=1e-2;
        beta=1e-5;
        C3=alfa*M3+beta*K3;
        zz = (K3-M3*omega^2 + i*omega*C3);
        H(:,:) = inv(zz);
        z=norm(abs(H(10,1)));
        


end

 

%---------------------------------

function z=func1_max(GG)
global A Amp m freq omega M3
  
    N = 5;%Nelementos
    E = 7e10;%modulo de young
    Le = 5; %comprimento da barra
    h = 0.5;%comprimento
    b = 0.5;%largura
    Ae = b*h;%area
    I = (b*h*(b^2 + h^2))/24;%inercia


    %----------------------------------

        rho = GG(1);
        [K3,M3] = beam_ef_creating(N,rho,E,Ae,I,Le);
        alfa=1e-2;
        beta=1e-5;
        C3=alfa*M3+beta*K3;
        zz = (K3-M3*omega^2+ i*omega*C3);
        H(:,:) = inv(zz);        
        z=norm(abs(H(10,1)));
        z=-z;
  
    %----------------------------------
end

%-----------------------------------