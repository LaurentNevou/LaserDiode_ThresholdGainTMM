%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% lnev, 31 March 2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program computes the threshold gain in a DFB laser cavity defined by 4 
% parameters (n1, n2, Dn2, n3, L and LL). It uses the formula of the Transfer Matrix 
% Methods aproach. While sweeping the Gain values, the Transmission is computed.
% At a certain (lambda,gain) coupled values, the transmission diverges and it is
% where the gain equal the threshold gain.
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
% "Optoelectronic", Cambridge Books Online
% Prof. Emmanuel Rosencher,
% Complement to Chapter 13
% 13.A Distributed feedback (DFB) lasers, page 660
% http://dx.doi.org/10.1017/CBO9780511754647.028
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Optical Waves in Crystals"
% Prof. Amnon Yariv
% Chap 11, Guided Waves and Integrated Optics
% 11.4 Periodic Waveguide-Bragg Reflection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_vec=linspace(9.8,10.2,1000)*1e-6;            %% Wavelength [m]
%lambda_vec=linspace(938,942,1000)*1e-9;            %% Wavelength [m]

Gain=[0:0.1:50]*1e2;                                %% Gain [m-1]

Tmax=1e4; % Transmission value at which the algorithm stops

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cavity parameters

n1=1;                     %% refractive index of the air
n2=3.2;                   %% refractive index of the semiconductor (GaAs=3.2 @10um and GaAs=3.6 @1um)
Dn2=1.5e-2;               %% optical index variation in the DFB (na-nb)
n3=1;                     %% refractive index of the air (for HR coating, n3>500 and Tmax=1e3)
L=1e-3;                   %% Lenght of the cavity [meter]
lambda0=10e-6;            %% Central wavelength design [m]
LL=lambda0/(2*abs(n2));   %% DFB PERIOD thickness at lambda/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, I compute the transmission curve at Gain=0

for jj=1:length(lambda_vec)

    lambda=lambda_vec(jj);
    [T,R]=Transmission_DFB_f(lambda,0,L,LL,n1,n2,Dn2,n3);
    
    Trans(jj,:)=T;
    Reflc(jj,:)=R;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, I buid an algorithm in order to find the values (lambda,Gain) where T diverges

lambda  = lambda_vec(1);
jump    = lambda^2/2/n2/L;    % FabryPerot mode spacing
dlambda = jump/3;
Ttmp=0;
LambdaGain=[];TT=[];RR=[];

%figure

while lambda<lambda_vec(end)

cc=0;

while Ttmp<Tmax

    cc=cc+1;
    if cc>100
        lambda = lambda + lambda0 * 2/pi*Dn2/n2 / 2;
        dlambda = jump/100;
        cc=0;
    end
    
    [T,R]=Transmission_DFB_f(lambda,Gain,L,LL,n1,n2,Dn2,n3);
    Ttmp_new=max(T);
  
    if Ttmp_new>Ttmp
        lambda=lambda+dlambda;
        Ttmp=Ttmp_new;
    elseif Ttmp_new<Ttmp
        dlambda=-dlambda/2;
        lambda=lambda+dlambda;
        Ttmp=Ttmp_new;
    end

    % to plot the convergence  
%     plot(cc,lambda*1e6,'bo-')
%     hold on;grid on;
%     pause(0.001)

end
    %hold off
    idx_T = find( T==max(T) );
    Gth   = Gain(idx_T);
    LambdaGain=[LambdaGain;lambda Gth max(T)];
    TT=[TT T];
    RR=[RR R];
    display(strcat('lambda=',num2str(lambda*1e6),'um ; ThGain=',num2str(Gth/100),'cm-1'))

    jump    = lambda^2/2/n2/L;    % FabryPerot mode spacing
    dlambda = jump/100;
    lambda  = lambda + jump;
    Ttmp    = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Threshold gain formula

R1=((n1-n2)/(n1+n2))^2;                  % Reflection of the mirror facet
R2=((n2-n3)/(n2+n3))^2;                  % Reflection of the mirror facet
GthF=-1/2/L*log(R1*R2)/100;              % Threshold gain [cm-1]

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
yscale2=[0 20];

[AX,H1,H2] = plotyy(lambda_vec*1e6,Trans,LambdaGain(:,1)*1e6,LambdaGain(:,2)/100);

% for ii=1:length(LambdaGain(:,1))
%    plot([1 1]*LambdaGain(ii,1)*1e6,[0 2],'r-')
% end

set(H1,'color','b','linewidth',1,'marker','none');
set(H2,'color','r','linestyle','none','marker','o');

set(AX(1),'ycolor','b','xlim',xscale,'ylim',yscale1,'ytick',[0:0.1:1],'fontsize',15);
set(AX(2),'ycolor','r','xlim',xscale,'ylim',yscale2,'ytick',0:2:100,'fontsize',15);

xlabel('lambda (um)')
ylabel(AX(1),'Transmission')
ylabel(AX(2),'Threshold Gain (cm-1)')
title(strcat('L-cavity=',num2str(L*1e3),'mm; n1=',num2str(n1),'; n2=',num2str(n2),...
    '; n3=',num2str(n3),'; \Deltan2=',num2str(Dn2,'%.1e'),'; \Lambda=',num2str(LL*1e6,'%.2f'),'um'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,1,1,'fontsize',15)

semilogy(Gain/100,TT,'g.-')
hold on; grid on;
ylim([1e-1 1e6])
semilogy(Gain/100,RR,'b.-')
plot([1 1]*GthF,[1e-5 1e10],'k--')

ylabel('Transmission','fontsize',15)
xlabel('Gain (cm-1)','fontsize',15)

ylim([1e-1 1e6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%