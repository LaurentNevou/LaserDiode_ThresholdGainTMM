%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 14April2020, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the threhold gain of arbitrary cavity using the Matrix
% Transfer Method (TMM). The Gain is introduced inside the imaginary part of the
% optical index.
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

% Laurent s NOTES:
% Take care! In order to work for the VCSEL, the Electrical field MUST have a node
% at the QW. Also, for some reasons... the threshold gain computed with that method
% is exactly half of the one computed with the Formula

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz=10e-9;

lambda=1000*1e-9;            %% lambda could be a vector but it becomes very slow
%lambda=(800:1:1200)*1e-9;

%GGain=0;
%GGain=[-10:0.1:20]*1e2;     %% values for Fabry-Perot cavity [m-1]
GGain=[0:10:1500]*1e2;    %% values for VCSEL cavity [m-1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii =1:length(GGain)

Gain=GGain(ii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose between those 2 inputs
%input_GainFP                %% Gth < 20cm-1
input_GainVCSEL              %% Gth < 2000cm-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear z n nt t zz zv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z and the optical index n

t  = layer(:,1);
nt = layer(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(t)
  
  if j==1
    zz(1) = t(1);
    zv{1} = 0:dz:t(1); 
    z     = zv{1};
    n     = repmat((zv{j}'*0+1),[1 length(lambda)] ) .* repmat(nt(j,:),[length(zv{j}') 1]);
  else
    zz(j) = zz(end)+t(j);
    zv{j} = (zz(end-1)+dz):dz:zz(end);
    z     = [ z  zv{j} ];
    n     = [ n  ; repmat((zv{j}'*0+1),[1 length(lambda)] ) .* repmat(nt(j,:),[length(zv{j}') 1])  ];
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj=1:length(lambda)

  [AA,BB,psi] = TMM_f(zz,zv,nt(:,jj),nL,nR,lambda(jj));
  
  A(:,jj,ii)  = AA;
  B(:,jj,ii)  = BB;
  %PSI(:,jj,ii)= psi.'; %% sometimes, psi is longer than it should...
  PSI{jj,ii}= psi.';
  
end

  R(ii,:) = abs(B(1,:,ii)).^2;
  T(ii,:) = (nR/nL) * abs(A(end,:,ii)).^2 ;

end

idx=find(abs(lambda-lambda0)==min(abs(lambda-lambda0)));
idxT = find(T(:,idx)==max(T(:,idx)));
Gth = GGain(idxT);

display(strcat('lambda=',num2str(lambda0*1e9),'nm ; ThGain=',num2str(Gth/100),'cm-1'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X0fig=-3500; Y0fig=100;
X0fig=100; Y0fig=100;
Wfig=1000;Hfig=700;
FS=15;
LW=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(GGain)>1
    
    figure('Name','Results','position',[X0fig Y0fig Wfig Hfig])

    subplot(1,1,1,'fontsize',FS)

    semilogy(GGain*1e-2,R,'bo-')
    hold on;grid on;
    semilogy(GGain*1e-2,T,'go-')
    xlabel('Gain (cm-1)','fontsize',FS)
    ylabel('Transmission & Reflection','fontsize',FS)
    ylim([1e-1 1e6])
    title(strcat('\fontsize{15}ThresholdGain=',num2str(Gth/100),'cm-1 @\lambda=',num2str(lambda(idx)*1e9),'nm'))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(lambda)>1
    
    figure('Name','Results','position',[X0fig+100 Y0fig Wfig Hfig])
    subplot(1,1,1,'fontsize',FS)
    hold on;grid on;

    plot(lambda*1e9,T,'b-')
    xlabel('lambda (nm)','fontsize',FS)
    ylabel('Transmission','fontsize',FS)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Results','position',[X0fig+200 Y0fig Wfig Hfig])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,1,'fontsize',FS)
hold on;grid on;

plot(z*1e6,abs(n(:,idx)),'b')
xlabel('z (um)','fontsize',FS)
ylabel('Optical index','fontsize',FS)
xlim([z(1) z(end)]*1e6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,2,'fontsize',FS)
hold on;grid on;

plot(z*1e6,(abs(PSI{idx,idxT})).^2,'r')
xlabel('z (um)','fontsize',FS)
ylabel('|E|^2','fontsize',FS)
title(strcat('@\lambda0=',num2str(lambda(idx)*1e9),'nm'))
xlim([z(1) z(end)]*1e6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%