clear all
close all
clc
% load EMF_2da_celda1
load EMF_YTec_1
h=60/3600;
nt=2e3;
%---- parametros del  modelo |  
% a = 0.765;b = 0.49;r1 = 0.0151;

% xo=[1.2667e-02   2.3000e-01   4.5000e-01];
% r1=xo(1);
% b=xo(2);
% a=xo(3);

a=1.4437e-02;   b=2.04658e-02;   r1=0.00132;
Q1 = 20;

A1=[1 0 ;1-exp(-h/b) exp(-h/b)];B1=[-h;(b-a)*(1-exp(-h/b))-h]/Q1;

a2=a; b2=b; r2=r1; Q2=Q1; EMF2=EMF;
a3=a; b3=b; r3=r1; Q3=Q1; EMF3=EMF;
a4=a; b4=b; r4=r1; Q4=Q1;EMF4=EMF;


A2=[1 0 ;1-exp(-h/b2) exp(-h/b2)];B2=[-h;(b2-a2)*(1-exp(-h/b2))-h]/Q2;
A3=[1 0 ;1-exp(-h/b3) exp(-h/b3)];B3=[-h;(b3-a3)*(1-exp(-h/b3))-h]/Q3;
A4=[1 0 ;1-exp(-h/b4) exp(-h/b4)];B4=[-h;(b4-a4)*(1-exp(-h/b4))-h]/Q4;



%---- CONDICIONES INICIALES  ----
V1(1)=3.34;
V2(1)=3.34;
V3(1)=3.34;
V4(1)=3.3;

ci1=interp1(EMF(:,1),EMF(:,2),V1(1),'linear','extrap');
ci2=interp1(EMF2(:,1),EMF2(:,2),V2(1),'linear','extrap');
ci3=interp1(EMF3(:,1),EMF3(:,2),V3(1),'linear','extrap');
ci4=interp1(EMF4(:,1),EMF4(:,2),V4(1),'linear','extrap');

Xbm1 = [ci1; ci1];
Xbm2 = [ci2; ci2];
Xbm3 = [ci3; ci3];
Xbm4 = [ci4; ci4];


I=0*Q1;%en Amp

%---------para ecualizar
% Rq1=20;Rq2=20;Rq3=20;Rq4=20;
Rq=ones(4,1); % 1 Ohm
Vce=0.5; %en volts

IR1=0; IR2=0; IR3=0; IR4=0; 
I1=I+IR1;
I2=I+IR2;
I3=I+IR3;
I4=I+IR4;



DT_on(1)=0;
DT_on(2)=0;
DT_on(3)=0;
DT_on(4)=0;
DT_off(1)=0;
DT_off(2)=0;
DT_off(3)=0;
DT_off(4)=0;
EQ(1)=0;
EQ(2)=0;
EQ(3)=0;
EQ(4)=0;
eps_min=0.02;
eps_max=0.03; 
eq=[];
numbat=4;
Vmin=2.6;
Vmax=3.35;
mx=[];

for i=2:nt-1
%      if EQ(1)==0 , IR1=0;  DT_off(1)=DT_off(1)+1;end
%      if EQ(2)==0 , IR2=0;  DT_off(2)=DT_off(2)+1; end
%      if EQ(3)==0 , IR3=0;  DT_off(3)=DT_off(3)+1;end
%      if EQ(4)==0 , IR4=0;  DT_off(4)=DT_off(4)+1; end
%      
%      if EQ(1)==1  && i > 3; IR1=V1(i-1)/Rq1; DT_on(1)=DT_on(1)+1; end 
%      if EQ(2)==1  && i > 3; IR2=V2(i-1)/Rq2; DT_on(2)=DT_on(2)+1;  end 
%      if EQ(3)==1  && i > 3; IR3=V3(i-1)/Rq3; DT_on(3)=DT_on(3)+1; end 
%      if EQ(4)==1  && i > 3; IR4=V4(i-1)/Rq4; DT_on(4)=DT_on(4)+1;  end 
    
    for k=1:numbat
        if EQ(k)==0 , IR(k)=0;  DT_off(k)=DT_off(k)+1;end
        if EQ(k)==1  && i > 3; IR(k)=(V(k)-Vce)/Rq(k); DT_on(k)=DT_on(k)+1; end
        Ies(k)=I+IR(k);
    end
    
%--------   MODELOS de las baterias  ---------------------  

%     I1(i)=I+IR1;
%     I2(i)=I+IR2;
%     I3(i)=I+IR3;
%     I4(i)=I+IR4;

    I1(i)=Ies(1);
    I2(i)=Ies(2);
    I3(i)=Ies(3);
    I4(i)=Ies(4);

%----- Modelo 1: ___________________    
    Xbm1 = A1*Xbm1 + B1*I1(i); % Xb
    V1(i) = interp1(EMF(:,2),EMF(:,1),Xbm1(2),'pchip','extrap') - I1(i)*r1;
    Smodelo(i) = Xbm1(1); % la columna 1 de Xbm es el Soc
    Xmodelo(i) = Xbm1(2); % la columna 2 de Xbm es el X

        
%----- Modelo 2: ___________________    
    Xbm2 = A2*Xbm2 + B2*I2(i);
    V2(i) = interp1(EMF2(:,2),EMF2(:,1),Xbm2(2),'pchip','extrap') - I2(i)*r2;
    Smodelo2(i) = Xbm2(1); % la columna 1 de Xbm es el Soc
    Xmodelo2(i) = Xbm2(2); % la columna 2 de Xbm es el X

    
%----- Modelo 3: ___________________    
    Xbm3 = A3*Xbm3 + B3*I3(i); % Xb
    V3(i) = interp1(EMF3(:,2),EMF3(:,1),Xbm3(2),'pchip','extrap') - I3(i)*r3;
    Smodelo3(i) = Xbm3(1); % la columna 1 de Xbm es el Soc
    Xmodelo3(i) = Xbm3(2); % la columna 2 de Xbm es el X

        
%----- Modelo 4: ___________________    
    Xbm4 = A4*Xbm4 + B4*I4(i);
    V4(i) = interp1(EMF4(:,2),EMF4(:,1),Xbm4(2),'pchip','extrap') - I4(i)*r4;
    Smodelo4(i) = Xbm4(1); % la columna 1 de Xbm es el Soc
    Xmodelo4(i) = Xbm4(2); % la columna 2 de Xbm es el X
    
    
    V=[V1(i) V2(i) V3(i) V4(i)];
%-----se invierte la corriente cuando pasa  los valores
    if ((min(V)< Vmin) || (max(V)>Vmax))
            I=-I;
    for k=1:numbat; EQ(k)=0;end
    end
%-----ecualizador-------------


    
    mv=3.7;
    for k=1:numbat
        if ((V(k) < mv)  && (EQ(k)==0)) % solo puede ser mínimo si no se está ecualizando
          mv = V(k);
        end
    end
    
    for k=1:numbat 
         if ((( mv -V(k)) > eps_min) && (EQ(k)==1)&& (DT_on(k)>5)) 
         EQ(k) = 0;
         DT_on(k)=0;

         end    
         if (((V(k) - mv) > eps_max) && (EQ(k)==0) && (DT_off(k)>10))
         EQ(k)= 1;
         DT_off(k)=0;
         end   
    end
    eq=[eq;EQ];
    mx=[mx;max(V-min(mv))];
end


%%
figure (1)
plot ([V1' V2' V3' V4']),legend('V1' ,'V2','V3' ,'V4')
ylim([Vmin Vmax])
% figure (2)
% plot (mx,'.','markersize',10),legend('maxima diferencia V-min(V)')
% figure(3),hist(mx,100)
% figure (4)
% plot (eq),legend('Eq1' ,'Eq2','Eq3' ,'Eq4')