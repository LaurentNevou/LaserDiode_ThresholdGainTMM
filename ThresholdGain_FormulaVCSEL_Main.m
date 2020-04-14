%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% lnev, 1 April 2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program computes the threshold gain in a VCSEL cavity defined by few 
% parameters (lambda0, na, nb, nc, LQW, N_DBRn and N_DBRp). It uses the formula 
% of the Transfer Matrix Methods aproach. While sweeping the Gain values, the 
% Transmission is computed. When the transmission diverges, it is where the gain
% equal the threshold gain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Analysis of multielement semiconductor lasers"
% K. J. Ebeling and L. A. Coldren
% Journal of Applied Physics 54, 2962 (1983); doi: 10.1063/1.332498
% https://doi.org/10.1063/1.332498
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Quantum Cascade Lasers", Oxford: Oxford University Press, 2013.
% Prof. Jerome Faist
% CHAPTER 10. MODE CONTROL,
% 10.2 Distributed feedback cavity
% 10.2.1 Multilayer approach, Bragg reflection condition, page 169
% https://www.amazon.com/Quantum-Cascade-Lasers-J%C3%A9r%C3%B4me-Faist/dp/0198528248
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emmanuel Rosencher, Optoelectronic, Cambridge Books Online
% Complement to Chapter 13
% 13.C Vertical cavity surface emitting lasers (VCSELs), page 671
% http://dx.doi.org/10.1017/CBO9780511754647.028
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cavity parameters

lambda0=1000e-9;            %% Central wavelength design [m]
na = 3;                     %% DBR refractive index-a, AlAs
nb = 3.6;                   %% DBR refractive index-b, GaAs
nc = 3.6;                   %% refractive index of the cavity, GaAs
lc = 2 * lambda0/(2*nc);    %% Lenght of the cavity [m]
LQW= 10e-9;                 %% quantum well thickness in which the gain will be [m]
N_DBRn=30;                  %% amount of DBR n-doped pairs
N_DBRp=20;                  %% amount of DBR p-doped pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_vec=linspace(800,1200,2000)*1e-9;            %% Wavelength [m]
lambda_vec=sort([lambda_vec lambda0]);              %% here, I make sure lambda0 is inside the vector lambda

Gain=[0:10:2000]*1e2;                               %% Gain [m-1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, I compute the transmission curve at Gain=0

for jj=1:length(lambda_vec)

    lambda=lambda_vec(jj);
    [T,R]=Transmission_VCSEL_f(lambda,0,lambda0,na,nb,nc,N_DBRn,N_DBRp,lc,LQW);
    
    Trans(jj,:)=T;
    Reflc(jj,:)=R;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In a VCSEL, the cavity is build such that only one longitudinal mode can live in the cavity
% Therefore, there is no need to scan in lambda! The mode is at lambda0

[T,R]=Transmission_VCSEL_f(lambda0,Gain,lambda0,na,nb,nc,N_DBRn,N_DBRp,lc,LQW);

idx_T = find( T==max(T) );
Gth   = Gain(idx_T);
LambdaGain=[lambda0 Gth max(T)];
display(strcat('lambda=',num2str(lambda0*1e6),'um ; ThGain=',num2str(Gth/100),'cm-1'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%X0fig=-3500; Y0fig=100;
X0fig=100; Y0fig=100;
Wfig=1000;Hfig=800;

figure('Name','Results','position',[X0fig Y0fig Wfig Hfig])
subplot(1,1,1)
hold on;grid on;
xscale=[lambda_vec(1) lambda_vec(end)]*1e6;
yscale1=[0 1];
yscale2=[0 LambdaGain(2)/100*1.5];

[AX,H1,H2] = plotyy(lambda_vec*1e6,Trans,LambdaGain(:,1)*1e6,LambdaGain(:,2)/100);

set(H1,'color','b','linewidth',1,'marker','none');
set(H2,'color','r','linestyle','none','marker','o');

set(AX(1),'ycolor','b','xlim',xscale,'ylim',yscale1,'ytick',[0:0.1:1],'fontsize',15);
set(AX(2),'ycolor','r','xlim',xscale,'ylim',yscale2,'ytick',0:200:LambdaGain(2)/10,'fontsize',15);

xlabel('lambda (um)')
ylabel(AX(1),'Transmission')
ylabel(AX(2),'Threshold Gain (cm-1)')
title(strcat('\lambda0=',num2str(lambda0*1e9),'nm; na=',num2str(na),'; nb=',num2str(nb),...
    '; nc=',num2str(nc),'; N-DBRn=',num2str(N_DBRn),'; N-DBRp=',num2str(N_DBRp),'; QWtick=',num2str(LQW*1e9),'nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,1,1,'fontsize',15)

semilogy(Gain/100,T,'g.-')
hold on; grid on;
ylim([1e-1 1e6])
semilogy(Gain/100,R,'b.-')

ylabel('Transmission','fontsize',15)
xlabel('Gain (cm-1)','fontsize',15)
title(strcat('\fontsize{15}ThresholdGain=',num2str(Gth/100),'cm-1 @\lambda=',num2str(lambda0*1e9),'nm'))
ylim([1e-1 1e6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%