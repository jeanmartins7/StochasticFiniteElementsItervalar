function intrval_solution

clear all

clc

close all

global omega

%Intervalos dos Parametros

%------------------------------

 
%  - K 

p_1_m=1400;

p_1_l=1400-.1*p_1_m;

p_1_r=1400+.1*p_1_m;

p_1_I=[p_1_l,p_1_r];

 

 

%NTF 2  - C

p_2_m=50;%1%-4%2;

p_2_l=p_2_m-.3*p_2_m;%1%-12%2;

p_2_r=p_2_m+.3*p_2_m;%1%4%2;

p_2_I=[p_2_l,p_2_r];

 

%Interval - m

m=4.5;

p_3_m=4.5;%1%-4%2;

p_3_l=4.5-.1*p_3_m;%1%-12%2;

p_3_r=4.50+.1*p_3_m;%1%4%2;

p_3_I=[p_3_l,p_3_r];

 

 

%Numero de intervalos

n_I=3;

 


%Options GA

options = gaoptimset('PopulationSize',n_I*20,'PopulationType','doubleVector');%,...

                    %'PopInitRange',[80.0e-3 180.0e-3 115 3 100.0e-3; 100.0e-3 200.0e-3 117 5 180.0e-3],...

                    %'PlotFcns',{@gaplot_GCI, @gaplotbestindiv, @gaplotdistance});

 

%--------------------------------------------------------------

% Analysis Intervalar

omega_f=100;% Frequencia final

delta_omega=1;

omega_i=delta_omega;

iiii=omega_i/delta_omega;

for omega=omega_i:delta_omega:omega_f

    

        lb = [p_1_l p_2_m p_3_l];

        ub = [p_1_r p_2_m p_3_r];

 

        

        %Interval Analysis for all Parameters

        [x_min,fval_min] = ga(@func1,n_I,[],[],[],[],lb,ub,[],options);

        [x_max,fval_max] = ga(@func1_max,n_I,[],[],[],[],lb,ub,[],options);

        y(iiii)=func1([p_1_m p_2_m p_3_m]);

        y_min(iiii)=fval_min;

        y_max(iiii)=-fval_max;

        p_min(iiii,:)=x_min;

        p_max(iiii,:)=x_max;

        

        iiii=iiii+1;

end

%--------------------------------------------------------------

%--------------------------------------------------------------

 

%Plot of interval Analysis of all parameters

plot([delta_omega:delta_omega:omega_f],y,'g');

hold on

plot([delta_omega:delta_omega:omega_f],y_min,'b');

hold on

plot([delta_omega:delta_omega:omega_f],y_max,'r');

title('Análise Intervalar','Interpreter','latex','fontsize',16)

ylabel('$|H(\omega)|$','Interpreter','latex','fontsize',22)

xlabel('$\omega$ [rad/s]','Interpreter','latex','fontsize',22)

box on

%h=legend('$k$','$c$','$m$')

%set(h,'FontSize',22,'Interpreter','latex')

 

 

end

%-------------------------------------------------------------------------

 

 

 

 

%---------------------------------

 

function z=func1(GG)

 

global omega   

 

    k=GG(1);

    c=GG(2);

    m=GG(3);

    %----------------------------------

    %Dominio da frequiencia

    o=omega;

    k_d=-o^2*m+i*o*c+k; % ;;

    u=inv(k_d);

    z=norm(abs(u));

    %----------------------------------

    

%     z=round(z*1000)/1000;

end

 

%---------------------------------

function z=func1_max(GG)

 

global omega   

 

    k=GG(1);

    c=GG(2);

    m=GG(3);

    %----------------------------------

    %Dominio da frequiencia

    o=omega;

    k_d=-o^2*m+i*o*c+k; % 

    u=inv(k_d);

    z=norm(abs(u));

    %----------------------------------

%     z=round(z*1000)/1000;

    z=-z;

end

%-----------------------------------