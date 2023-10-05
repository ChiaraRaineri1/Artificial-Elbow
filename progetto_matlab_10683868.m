%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RAINERI CHIARA %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROGETTO MATLAB 2020/2021 MECCANICA APPLICATA ALLE MACCHINE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%%%INPUT
L15=70;          %[mm] 
L12=32;          %[mm] 
L24=34;          %[mm] 
L45=40;          %[mm] 
L56=17;          %[mm] 
L36=145;         %[mm] 
L16=87;          %[mm]
L13=200;         %[mm] 
p=2.1;           %[rad]
m=4.1;           %[rad] 
df=3.6416;       %[rad]

%Parametri dipendenti dal codice persona: 10683868 
% x=6 y=8
%L27=300+10*x+y 
%L2G=0.4*L27  
x=6;
y=8;
L27=300+10*x+y; %[mm]
L2G=0.4*L27;     %[mm]

%%
%%%%%%%%%%%% ANALISI CINEMATICA %%%%%%%%%%%%

%METODO CHIUSURA GLIFO
fprintf('RISOLUZIONE GLIFO\n\n');
% c + a = b
%VETTORI
a=L36;
b=L16;
c=L13;

%SISTEMA
% asse x: c*cos(p)+a*cos(alpha)=b*cos(beta);
% asse y: c*sin(p)+a*sin(alpha)=b*sin(beta);
%ANGOLI E GEOMETRIA
fprintf('Ricerca angoli con fsolve');
x1=[4.5 2.5];
F=@(x)[c*cos(p)+a*cos(x(1))-b*cos(x(2));
       c*sin(p)+a*sin(x(1))-b*sin(x(2))];
SOL=fsolve(F,x1);
alpha=SOL(1);
beta=SOL(2);
%DERIVATA 1^
%-a*alphap*sin(alpha)+b*betap*sin(beta)=-ap*cos(alpha);
% a*alphap*cos(alpha)-b*betap*cos(beta)=-ap*sin(alpha);
%DERIVATA 2^
%-a*alphapp*sin(alpha)+b*betapp*sin(beta)=-b*betap^2*cos(beta)-a*alphap^2*cos(alpha)+2*alphap*betap*sin(alpha)
% a*alphapp*cos(alpha)-b*betapp*cos(beta)=-b*betap^2*sin(beta)+a*alphap^2*sin(alpha)-2*alphap*betap*cos(alpha)

%%%%%%%%%%%% RISOLUZIONE %%%%%%%%%%%%

%Velocità motore 
Vm=4750;             %[giri/min]
ap=4750/60;          %[mm/s]

%VELOCITA'
fprintf('VELOCITA''\n')

%Matrice coefficienti
M=[-a*sin(alpha) +b*sin(beta);
    a*cos(alpha) -b*cos(beta)];
%Vettore termini noti 
N=[-ap*cos(alpha);
   -ap*sin(alpha)];
Y=M\N;
alphap=Y(1);
betap=Y(2);
disp(['Velocità angolare dell asta L36 è: ',num2str(alphap),'  rad/s']);
disp(['Velocità angolare dell asta L16 è: ',num2str(betap),'   rad/s']);
fprintf('Segno negativo: velocità in senso orario\n\n');

%ACCELERAZIONI
fprintf('ACCELERAZIONE\n');

%Vettore termini noti 
N=[-b*betap^2*cos(beta)+a*alphap^2*cos(alpha)+2*alphap*ap*sin(alpha);
   -b*betap^2*sin(beta)+a*alphap^2*sin(alpha)-2*alphap*ap*cos(alpha)];
Y=M\N;
alphapp=Y(1);
betapp=Y(2);
disp(['Accelerazione angolare dell asta L36 è: ',num2str(alphapp),'  rad/s^2']);
disp(['Accelerazione angolare dell asta L16 è: ',num2str(betapp),'   rad/s^2']);
fprintf('Segno negativo: velocità in senso orario\n\n\n');

%%
%METODO CHIUSURA QUADRILATERO
b1=L15;
e=L45;
f=L12;
d=L24;

%SISTEMA 
%b1*cos(gamma)+e*cos(eta)=f*cos(m)+d*cos(delta)
%b1*sin(gamma)+e*sin(eta)=f*sin(m)+d*sin(delta)

%ANGOLI E GEOMETRIA
fprintf('Ricerca angoli con fsolve');
x1=[5.128 2.7548];
F=@(x)[b1*cos(beta)+e*cos(x(1))-f*cos(m)-d*cos(x(2));
       b1*sin(beta)+e*sin(x(1))-f*sin(m)-d*sin(x(2))];
SOL=fsolve(F,x1);
gamma=SOL(1);
delta=SOL(2);

%DERIVATA 1^
%-e*gammap*sin(gamma)+d*deltap*sin(delta)= b1*beta*sin(beta)
% e*gammap*cos(gamma)-d*deltap*cos(delta)=-b1*beta*cos(beta)
%DERIVATA2 ^
%-e*gammapp*sin(gamma)+d*deltapp*sin(delta)= b1*betapp*sin(beta)+b1*betap^2*cos(beta)+e*gammap^2*cos(gamma)-d*deltap^2*cos(delta)
% e*gammapp*cos(gamma)-d*deltapp*cos(delta)=-b1*betapp*cos(beta)+b1*betap^2*sin(beta)+e*gammap^2*sin(gamma)-d*deltap^2*sin(delta)

%VELOCITA'
fprintf('VELOCITA''\n')

%Matrice coefficienti
M=[-e*sin(gamma)  d*sin(delta);
    e*cos(gamma) -d*cos(delta)];
%Vettore termini noti 
N=[+b1*betap*sin(beta);
   -b1*betap*cos(beta)];
Y=M\N;
gammap=Y(1);
deltap=Y(2);
disp(['Velocità angolare dell asta L45 è: ',num2str(gammap),'  rad/s']);
disp(['Velocità angolare dell asta L24 è: ',num2str(deltap),'  rad/s']);
fprintf('Segno negativo: velocità in senso orario\n\n');

%ACCELERAZIONI
fprintf('ACCELERAZIONE\n');

%Vettore termini noti 
N=[ b1*betapp*sin(beta)+b1*betap^2*cos(beta)+e*gammap^2*cos(gamma)-d*deltap^2*cos(delta);
   -b1*betapp*cos(beta)+b1*betap^2*sin(beta)+e*gammap^2*sin(gamma)-d*deltap^2*sin(delta)];
Y=M\N;
gammapp=Y(1);
deltapp=Y(2);
disp(['Accelerazione angolare dell asta L45 è: ',num2str(gammapp),'  rad/s^2']);
disp(['Accelerazione angolare dell asta L24 è: ',num2str(deltapp),'  rad/s^2']);
fprintf('Segno negativo: velocità in senso orario\n\n\n');

%%
%VELOCITA' E ACCELERAZIONI AVAMBRACCIO 
epsylon=(delta+df)-2*pi;
l2G=[L2G*cos(epsylon),L2G*sin(epsylon),0];
l2M=[L27*cos(epsylon),L27*sin(epsylon),0];
deltapv=[0 0 deltap];
deltappv=[0 0 deltapp];
VG= cross(deltapv,l2G);
VM= cross(deltapv,l2M);
fprintf('VELOCITA''AVAMBRACCIO E MANO\n');
disp(['Velocità baricentro avambraccio:',num2str(norm(VG)),' mm/s']);
disp(['Velocità mano:',num2str(norm(VM)),' mm/s']);
fprintf('\n');
%ACCELERAZIONE
AG=cross(deltappv,l2G)+cross(deltapv,cross(deltapv,l2G));
AM=cross(deltappv,l2M)+cross(deltapv,cross(deltapv,l2M));
fprintf('ACCELERAZIONE AVAMBRACCIO E MANO\n');
disp(['Accelerazione baricentro avambraccio:',num2str(norm(AG)),' mm/s^2']);
disp(['Accelerazione mano:',num2str(norm(AM)),' mm/s^2']);




%%
%SCHEMA cinematico 
x1=c*cos(pi-p);
y1=-c*sin(pi-p);
x2=x1-f*cos(m-pi);
y2=y1-f*sin(m-pi);
x3=0;
y3=0;
x4=x2-d*cos(pi-delta);
y4=y2+d*sin(pi-delta);
x5=x4-e*cos(2*pi-gamma);
y5=y4+e*sin(2*pi-gamma);
x6=a*sin(alpha-(3/2)*pi);
y6=-a*cos(alpha-(3/2)*pi);
x7=x2+l2M(1);
y7=y2+l2M(2);
xg=x2+L2G*cos(epsylon);
yg=y2+L2G*sin(epsylon);

figure(1);
title('SCHEMA SEMPLIFICATO BRACCIO MECCANICO');
grid on;
hold on;
l36=line([x3,x6],[y3,y6]);
set(l36,'color','r','linewidth',3);
l13=line([x3,x1],[y3,y1],'linestyle','--');
set(l13,'color','k','linewidth',1);
l16=line([x1,x6],[y1,y6]);
set(l16,'color','b','linewidth',4);
l12=line([x1,x2],[y1,y2],'linestyle','--');
set(l12,'color','k','linewidth',1);
l24=line([x2,x4],[y2,y4]);
set(l24,'color','y','linewidth',4);
l45=line([x4,x5],[y4,y5]);
set(l45,'color','g','linewidth',4);
l27=line([x2,x7],[y2,y7]);
set(l27,'color','m','linewidth',4);
scatter(x1,y1,'filled');
scatter(x2,y2,'filled');
scatter(x3,y3,'filled');
scatter(x4,y4,'filled');
scatter(x5,y5,'filled');
scatter(x6,y6,'filled');
scatter(x7,y7,'filled');
scatter(xg,yg,'filled');
text(x1+2,y1+2,'1:Cerniera fissa');
text(0+3,y3,'3:Cerniera fissa');
text(x6-10,y6+3,'6:Carrello');
text(x4-7,y4-4,'4');
text(x5+3,y5+4,'5');
text(x2,y2-6,'2:Cerniera fissa');
text(x7,y7-6,'7:Mano');
text(xg,yg,'G:Baricentro');
text(xg,yg+30,'Vg');
text(xg-5,yg-30,'Ag');
text(x7,y7+30,'Vm');
text(x7,y7-30,'Am');

gv=quiver(xg,yg,VG(1),VG(2),0.2);
ga=quiver(xg,yg,AG(1),AG(2),0.08);
%Rappresentazione accelerazione
mv=quiver(x7,y7,VM(1),VM(2),0.1);
ma=quiver(x7,y7,AM(1),AM(2),0.05);



%%
%%%%%%%%%%%% MOTO IN GRANDE %%%%%%%%%%%%

%Vettore tempo 
dt=0.005;
tmax=0.25; %[sec]
tempo=0:dt:tmax;

%INZIZIALIZZO VARIABILI E VETTORI
aptime=ap*(ones(size(tempo)));
atime=a*(ones(size(tempo)))+aptime.*tempo; %velocità costante

alphatime=zeros(size(tempo));
alphaptime=zeros(size(tempo));
alphapptime=zeros(size(tempo));

betatime=zeros(size(tempo));
betaptime=zeros(size(tempo));
betapptime=zeros(size(tempo));

gammatime=zeros(size(tempo));
gammaptime=zeros(size(tempo));
gammapptime=zeros(size(tempo));

deltatime=zeros(size(tempo));
deltaptime=zeros(size(tempo));
deltapptime=zeros(size(tempo));


epsylontime=zeros(size(tempo));
XMtime=zeros(size(tempo));
YMtime=zeros(size(tempo));
Vmtime=zeros(size(tempo));
Amtime=zeros(size(tempo));
%Istanti iniziali 
aptime(1)=ap;
atime(1)=a; %velocità costante

alphatime(1)=alpha;
alphaptime(1)=alphap;
alphapptime(1)=alphapp;

betatime(1)=beta;
betaptime(1)=betap;
betapptime(1)=betapp;

gammatime(1)=gamma;
gammaptime(1)=gammap;
gammapptime(1)=gammapp;

deltatime(1)=delta;
deltaptime(1)=deltap;
deltapptime(1)=deltapp;

YMtime(1)=y7;
XMtime(1)=x7;
Vmtime(1)=norm(VM);
Amtime(1)=norm(AM);
epsylontime(1)=epsylon;
%ciclo for
for i=2:length(tempo)
    
    %POSIZIONE-GLIFO
    x1=[alphatime(i-1) betatime(i-1)];
    F=@(x)[c*cos(p)+atime(i)*cos(x(1))-b*cos(x(2));
           c*sin(p)+atime(i)*sin(x(1))-b*sin(x(2))];
    SOL=fsolve(F,x1);
    alphatime(i)=SOL(1);
    betatime(i)=SOL(2);
    x1=[gammatime(i-1) deltatime(i-1)];
    F=@(x)[b1*cos(betatime(i))+e*cos(x(1))-f*cos(m)-d*cos(x(2));
           b1*sin(betatime(i))+e*sin(x(1))-f*sin(m)-d*sin(x(2))];
    SOL=fsolve(F,x1);
    gammatime(i)=SOL(1);
    deltatime(i)=SOL(2);
    epsylontime(i)=(deltatime(i)+df)-2*pi;
    
    XMtime(i)=x2+L27*cos(epsylontime(i));
    YMtime(i)=y2+L27*sin(epsylontime(i));
    
    %VELOCITA'GLIFO
    M=[-atime(i)*sin(alphatime(i)) +b*sin(betatime(i));
        atime(i)*cos(alphatime(i)) -b*cos(betatime(i))];
    %Vettore termini noti 
    N=[-aptime(i)*cos(alphatime(i));
       -aptime(i)*sin(alphatime(i))];
        Y=M\N;
    alphaptime(i)=Y(1);
    betaptime(i)=Y(2);
    %ACCELERAZIONE GLIFO
    N=[ -b*betaptime(i)^2*cos(betatime(i))+atime(i)*alphaptime(i)^2*cos(alphatime(i))+2*alphaptime(i)*aptime(i)*sin(alphatime(i));
        -b*betaptime(i)^2*sin(betatime(i))+atime(i)*alphaptime(i)^2*sin(alphatime(i))-2*alphaptime(i)*aptime(i)*cos(alphatime(i))];
    Y=M\N;
    alphapptime(i)=Y(1);
    betapptime(i)=Y(2);

    %VELOCITA' QUADRILATERO
    %Matrice coefficienti
    M=[-e*sin(gammatime(i))  d*sin(deltatime(i));
        e*cos(gammatime(i)) -d*cos(deltatime(i))];
    %Vettore termini noti 
    N=[+b1*betaptime(i)*sin(betatime(i));
       -b1*betaptime(i)*cos(betatime(i))];
    V=M\N;
    gammaptime(i)=V(1);
    deltaptime(i)=V(2);
    
    %ACCELERAZIONI
    M=[-e*sin(gammatime(i))  d*sin(deltatime(i));
        e*cos(gammatime(i)) -d*cos(deltatime(i))];
    N=[ -d*(deltaptime(i))^2*cos(deltatime(i))+b1*betapptime(i)*sin(betatime(i))+b1*(betaptime(i))^2*cos(betatime(i))+e*(gammaptime(i))^2*cos(gammatime(i));
        -d*(deltaptime(i))^2*sin(deltatime(i))-b1*betapptime(i)*cos(betatime(i))+b1*(betaptime(i))^2*sin(betatime(i))+e*(gammaptime(i))^2*sin(gammatime(i))];
    Y=M\N;
    gammapptime(i)=Y(1);
    deltapptime(i)=Y(2);

    Vmtime(i)=deltaptime(i)*L27;
    
    deltappvtime=[0 0 deltapptime(i)];
    deltapvtime=[0 0 deltaptime(i)];
    
    %Amtime(i)=sqrt((-deltapptime(i)*l2M(2)-deltaptime(i)^2*l2M(1))^2+(deltapptime(i)*l2M(1)-deltaptime(i)^2*l2M(2))^2);
    Amtime(i)=norm(cross(deltappvtime,l2M)+cross(deltapvtime,cross(deltapvtime,l2M)));

end





figure (2);
subplot(2,2,1);
hold on;
plot(tempo,XMtime);

title('SPOSTAMENTO ASSE X');
xlabel('tempo [sec]');
ylabel('Coordinata x [mm]');
grid on;
subplot(2,2,2)
hold on;

plot(tempo,YMtime);
title('SPOSTAMENTO ASSE Y');
xlabel('tempo [sec]');
ylabel('Coordinata y [mm]')
grid on;
subplot(2,2,3)
plot(tempo,Vmtime);
title('VELOCITA''MANO');
ylabel('mm/s');
xlabel('tempo [sec]');
grid on;
subplot(2,2,4)
plot(tempo,Amtime);
title('ACCELERAZIONE MANO');
ylabel('mm/s^2');
xlabel('tempo [sec]');
grid on;

%%
%GRAFICO IN MOVIMENTO 

%PUNTI FISSI
x1=c*cos(pi-p);
x1v=x1*(ones(size(tempo))); 
y1v=y1*(ones(size(tempo))); 
x2v=x2*(ones(size(tempo)));
y2v=y2*(ones(size(tempo)));
x3v=zeros(size(tempo));
y3v=zeros(size(tempo));

x4v=x2v-d*cos(pi-deltatime);
y4v=y2+d*sin(pi-deltatime);
x5v=x4v-e*cos(2*pi-gammatime);
y5v=y4v+e*sin(2*pi-gammatime);
x6v=atime.*sin(alphatime-(3/2)*pi);
y6v=-atime.*cos(alphatime-(3/2)*pi);
x7v=XMtime;
y7v=YMtime;


figure(3)
title('SCHEMA SEMPLIFICATO BRACCIO MECCANICO IN MOVIMENTO');
grid on;
for i=1:length(tempo)
    hold on;   %%dentro il ciclo for: tengo tutti i lati solo nell'istante i
    box on;
    cla
    %PUNTI FISSI
    l13z=line([x3v(i),x1v(i)],[y3v(i),y1v(i)],'linestyle','--');
    set(l13z,'color','k','linewidth',1);
    l12z=line([x1v(i),x2v(i)],[y1v(i),y2v(i)],'linestyle','--');
    set(l12z,'color','k','linewidth',1);
    scatter(x1v(i),y1v(i),'filled');
    scatter(x2v(i),y2v(i),'filled');
    scatter(x3v(i),y3v(i),'filled');
    %ASTE IN MOVIMENTO 
    l36z=line([x3v(i),x6v(i)],[y3v(i),y6v(i)]);
    set(l36z,'color','r','linewidth',3);
    scatter(x6v(i),y6v(i),'filled');
    l27z=line([x2v(i),x7v(i)],[y2v(i),y7v(i)]);
    set(l27z,'color','m','linewidth',4);
    scatter(x7v(i),y7v(i),'filled');
    l16z=line([x1v(i),x6v(i)],[y1v(i),y6v(i)]);
    set(l16z,'color','b','linewidth',4);
    l24z=line([x2v(i),x4v(i)],[y2v(i),y4v(i)]);
    set(l24z,'color','y','linewidth',4);
    l45z=line([x4v(i),x5v(i)],[y4v(i),y5v(i)]);
    set(l45z,'color','g','linewidth',4);
    scatter(x4v(i),y4v(i),'filled');
    scatter(x5v(i),y5v(i),'filled');
    scatter(x7v(i),y7v(i),'filled');
    pause(0.01)
    xlim([0 550])
    ylim([-270 50])
   
end

%%
%Azionamento motore
%determina la coppia che deve sviluppare il motore per mantenere il meccanismo in equilibrio.
%Bilancio di potenze 
VGm=VG*0.001;
VMm=VM*0.001;
AGm=AG*0.001;
AMm=AM*0.001;
Mmano=0.1+0.1*x+0.01*y; %[Kg]
Moggetto=x+y+1;  
Mma=Moggetto+Mmano;%[Kg]
Mav=1+0.1*x+0.01*y;
g=9.8; %[m/s^2]
G=[0 -g 0]; 
Jav=8e-3; %[kg*m^2];
omega=497.42; %[rad/s];

%Cm*omega+Mav*(G*VGm.')+Mma*(G*VMm.') => Sommatoria potenze 
%d(1/2Mav*VGm^2+1/2Jav*deltap^2+1/2Mma*VMm^2)/dt =>Derivata energia cinetica

Cm=(Mma*(VMm*AMm.')+Mav*(VGm*AGm.')+Jav*(deltapv*deltappv.')-Mav*(G*VGm.')-Mma*(G*VMm.'))/(omega);
F=(Mma*(VMm*AMm.')+Mav*(VGm*AGm.')+Jav*(deltapv*deltappv.')-Mav*(G*VGm.')-Mma*(G*VMm.'))/(ap*0.001);
disp(['Coppia necessaria all''equilibrio:',num2str(Cm),' [N*m]']);
disp(['Forza punto 6 all''equilibrio:',num2str(F),' [N]']);


%EQUILIBRIO 
M=[cos(alpha) 1                 0      0                  0       0 0  1               0;
   sin(alpha) 0                 1      0                  0       0 0  0               1;
   0         (abs(y1)-abs(y6)) (x1-x6) 0                  0       0 0  abs(y5)-abs(y6) (x5-x6);             
   0          0                 0      0                  0       1 0 -1               0;
   0          0                 0      0                  0       0 1  0               -1;
   0          0                 0      0                  0       0 0  abs(y4)-abs(y5) (x4-x5);
   0          0                 0     -1                  0       1 0  0                0;
   0          0                 0      0                  1       0 1  0                0;
   0          0                 0      -(abs(y2)-abs(y4)) (x2-x4) 0 0  0                0];

 N=[0;
    0;
    0;
    0;
    0;
    0;
   +Mav*abs(AGm(1))+Mma*abs(AMm(1));
   -Mav*g-Mma*g+Mav*abs(AGm(2))+Mma*abs(AMm(1));
   -Jav*deltapp-Mav*g*(xg-x4)-Mma*g*(x7-x4)+Mav*abs(AGm(2))*(xg-x4)+Mav*abs(AGm(1))*(abs(yg)-abs(y4))+Mma*abs(AMm(2))*(x7-x4)+Mma*abs(AMm(1))*(abs(y7)-abs(y4))];
    
   S=M\N; 
   
   disp(['Forza punto 6 all''equilibrio:',num2str(S(1)),' [N]']);
   disp(['Reazione H1 all''equilibrio:',num2str(S(2)),' [N]']);
   disp(['Reazione V1 all''equilibrio:',num2str(S(3)),' [N]']);
  
   
%%
%LEGGE C-omega del motore:
%Xr=Coppia
%Yr=Omega
%DICHIARO VETTORI 

Com=1.69; %[N/m]
omegas=586.43; %Velocità di sincronismo 
Coppiav=0:0.01:1.69;
omegav= Coppiav.*(-omegas/Com)+omegas;

%Coppia massima espressa 
omegamax=omega;
Coppiamax=omega*(-Com/omegas)+Com;

%Bilancio di potenze 
%Trascuro il rendimento e trascuro qualsiasi tipo di attrito 
% W1+Wdiss+W2=d(Ec)/dt;
% Wm=Coppiamax*omega > moto diretto 
%(Wm)*R+(W2-d(Ec2)/dt)=0 
%Massa oggetto
R=0.85;
Mo=((-Coppiamax*omega*R+Mav*(G*VGm.')+Mav*(VGm*AGm.')+Jav*(deltapv*deltappv.'))/((G*VMm.')-(VMm*AMm.')))-Mmano;
disp(['Massa dell''oggetto affinchè lavori in condizioni ideali :',num2str(Mo),' [Kg]']);

%Potenza erogata dal motore 
Wmv=(Coppiav.*omegav);
Wmm=Coppiamax*omega;

figure(4);
xlabel('Valori di coppia[N*m]');
ylabel('Velocità angolare[rad/s] e potenza [W]');
grid on;
hold on;
plot(Coppiav,omegav,'linewidth',2);
plot(Coppiav,Wmv,'linewidth',2,'color','r');
title('CURVA CARATTERISTICA COPPIA E POTENZA');
scatter(0,omegas,'filled');
scatter(Com,0,'filled');
scatter(Coppiamax,omegamax,'filled');
scatter(Coppiamax,Wmm,'filled');
p=line([Coppiamax,Coppiamax],[omegamax,0]);
set(p,'color','k','linestyle','--');
legend('Curva caratteristica ω','Curva caratteristica W','Velocità di sincronismo','Coppia allo spunto','Cmax,ωmax','potenza erogata C.I.');

%%
%VIBRAZIONI
fprintf('\n VIBRAZIONI \n');
%DATI IN INPUT
k=240;          %[N/m]
r=50;           %[N/m]
C1=0.001;       %[N*m]
O1=0.5;         %[rad/s]
C2=0.002;       %[N*m]
O2=1.2;         %[rad/s]
omegam=497.42;  %[rad/s]
omegag=deltapp;
vg=norm(VGm);   %modulo della velocità del baricentro in metri al secondo 
v7=norm(VMm);   %modulo della velocità della mano in metri al secondo 
pv=0.001/(2*pi);  %[m/rad]

%EQ legami cinematici
Tg=vg/omegam;
T7=v7/omegam;
Tr=omegag/omegam;

%Equivalenti 
Jeq=Mav*Tg^2+Mma*T7^2+Jav*Tr^2;
req=r*pv^2;
keq=k*pv^2;

%EQ DI MOTO (Coordinata libera at- rotazione motore)
%Forzanti:
%C1(t)=C1*e^(O1t), C1(t)=C2*e^(O2t)
%Eq di moto:
%Jeq*atpp+req*atp+keq*at=C1(t)+C2(t);

% FREQUENZA PROPRIA DEL SISTEMA (LIBERO, NON SMORZATO) 
%Jeq*atpp+kat=0
%at=A*e^(lt), atp=A*l*e^(l*t), atpp=A*l^2*e^(l*t)
%Jeq*l^2+k*l=0

omegalib=sqrt(keq/Jeq);
f=omegalib/(2*pi);
disp(['Pulsazione propria del sistema:',num2str(omegalib),' [rad/s]']);
disp(['Frequenza propria del sistema:',num2str(f),' [hz]']);
% SMORZAMENTO ADIMENSIONALE
rc=2*Jeq*omegalib; %Smorzamento critico 
h=req/(rc);  
disp(['Smorzamento adimensionale:',num2str(h),' [N]']);


%FUNZIONE DI TRASFERIMENTO C1 E C2 
%PRINCIPIO DI SOVRAPPOSIZIONE DEGLI EFFETTI 
%C1 
%at(t)=Ao*e^(i*O1*t);
%atp(t)=i*Ao*e^(i*O1*t);
%atpp(t)=-Ao*e^(i*O1*t);
%(at/C1)=1/(-Jeq*O1^2+i*req*O1+keq)
%Vettore omega1 
O1t=0:0.001:1.5;
%Funzione modulo della funzione di trasferimento 
fdt=1./(-Jeq*O1t.^2+1i*req.*O1t+keq);
modv=abs(fdt);
fasev=angle(fdt)*180/pi;

%punti alla pulsazione della coppia 1 e 2 
z1=1/(-Jeq*O1^2+1i*req*O1+keq); 
z2=1/(-Jeq*O2^2+1i*req*O2+keq);
mz1=abs(z1);
mz2=abs(z2);
fz1=angle(z1)*180/pi;
fz2=angle(z2)*180/pi;


%grafico fdt
figure(5);
subplot(2,1,1);
hold on;
plot(O1t,modv,'linewidth',2,'color','k');

title('FUNZIONE DI TRASFERIMENTO');
ylabel('MODULO FDT');
xlabel('PULSAZIONE [rad/s]');
scatter(O1,mz1,'filled');
scatter(O2,mz2,'filled');
grid on;
text(O2-0.2,mz2+1000000,'Zona sismografica');
subplot(2,1,2);
plot(O1t,fasev,'linewidth',2,'color','k');
hold on;
ylabel('FASE FDT');
xlabel('PULSAZIONE [rad/s]');
scatter(O1,fz1,'filled');
scatter(O2,fz2,'filled');
grid on;

%Risposta libera h<1
alfa=h*omegalib;
omegad=omegalib*sqrt(1-h^2); %pulsazione sistema smorzato
lambda1=-alfa+1i*omegad;
lambda2=-alfa-1i*omegad;

%condizioni iniziali atempo(0)=0, aptempo(0)=0
%0=X1l+X2l
%0=lambda1*X1l+lambda2*X2l
M=[1       1;
   lambda1 lambda2];
N=[0
   497.42 ];
R=M\N;
X1=R(1);
X2=R(2);

%%risposta a(t) 

ti=0:0.1:300;
alib=real(X1*exp(lambda1*ti)+X2*exp(lambda2*ti));
ac1=real(X1*exp(lambda1*ti)+X2*exp(lambda2*ti))+mz1*C1*cos(O1.*ti+fz1);
ac2=real(X1*exp(lambda1*ti)+X2*exp(lambda2*ti))+mz2*C2*cos(O2*ti+fz2);
acomp=real(X1*exp(lambda1*ti)+X2*exp(lambda2*ti))+mz1*C1*cos(O1.*ti+fz1)+mz2*C2*cos(O2*ti+fz2);

t=0:0.1:120;
c1atempo=mz1*C1*cos(O1*t+fz1);
c2atempo=mz2*C1*cos(O2*t+fz2);
%atempo=mz1*cos(O1.*t+fz1)+mz2*cos(O2*t+fz2);
atempo=c1atempo+c2atempo;

figure (6)
subplot(2,2,4)
plot(ti,acomp,'linewidth',1,'color','r'),grid on,title('RISPOSTA COMPLETA (C1+C2)');
subplot(2,2,1);
plot(ti,alib,'linewidth',1,'color','k'),grid on,title('RISPOSTA LIBERA');
subplot(2,2,3)
plot(ti,ac1,'linewidth',1,'color','b'),grid on,title('RISPOSTA LIBERA+C1');
subplot(2,2,2)
plot(ti,ac2,'linewidth',1,'color','g'),grid on,title('RISPOSTA LIBERA+C2');

figure(7)
title('RISPOSTE A REGIME');
hold on;
grid on;
box on;
plot(t,atempo,'color','k','linewidth',1);
plot(t,c1atempo,'color','r');
plot(t,c2atempo,'color','b');
xlabel('tempo [s]');
ylabel('a(t)[rad] a regime');
legend('Rispostaa regime C1+C2','Risposta a regmie C1','Risposta a regime C2');

grid on;
