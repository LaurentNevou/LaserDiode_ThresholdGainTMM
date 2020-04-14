function[T,R]=Transmission_VCSEL_f(lambda,Gain,lambda0,na,nb,nc,N_DBRn,N_DBRp,lc,LQW)

c=2.99792458e8;             %% speed of light [m/s]
k=2*pi/lambda;              %% wave vector    [m-1]
w=2*pi*c/lambda;            %% pulsation      [s-1]

kk = Gain/2/w*c;            %% imaginary part of the optical refrative index
nGain = nc-1i*kk;           %% HERE IS THE MAJOR CHANGE!!! from PLUS (LOSSES) to MINUS (GAIN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DBR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

la=lambda0/(4*abs(na));     %% DBR layer-a thickness at lambda/4;
lb=lambda0/(4*abs(nb));     %% DBR layer-b thickness at lambda/4;

Pa = zeros(2,2,length(Gain));
Pa(1,1,:) = exp(+1i*k*na*la);
Pa(2,2,:) = exp(-1i*k*na*la);
Pb = zeros(2,2,length(Gain));
Pb(1,1,:) = exp(+1i*k*nb*lb);
Pb(2,2,:) = exp(-1i*k*nb*lb);

rab = (na-nb)./(na+nb);
rba = (nb-na)./(na+nb);

tab = 2*na./(na+nb);
tba = 2*nb./(na+nb);

Sab(1,1,:) = (1./tab) .* ones(1,length(Gain));
Sab(2,2,:) = (1./tab) .* ones(1,length(Gain));
Sab(1,2,:) = rab./tab;
Sab(2,1,:) = rab./tab;

Sba(1,1,:) = (1./tba) .* ones(1,length(Gain));
Sba(2,2,:) = (1./tba) .* ones(1,length(Gain));
Sba(1,2,:) = rba./tba;
Sba(2,1,:) = rba./tba;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fabry-Perot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pc = zeros(2,2,length(Gain));
Pc(1,1,:) = exp(+1i*k*nc*(lc-LQW)/2);
Pc(2,2,:) = exp(-1i*k*nc*(lc-LQW)/2);

Pgain = zeros(2,2,length(Gain));
Pgain(1,1,:) = exp(+1i*k*nGain*LQW);
Pgain(2,2,:) = exp(-1i*k*nGain*LQW);

% rac = (na-nc)./(na+nc);
% rca = (nc-na)./(na+nc);
% tac = 2*na./(na+nc);
% tca = 2*nc./(na+nc);

% Sac(1,1,:)=1./tac.*ones(1,length(Gain));
% Sac(2,2,:)=1./tac.*ones(1,length(Gain));
% Sac(1,2,:)=rac./tac;
% Sac(2,1,:)=rac./tac;

% Sca(1,1,:)=1./tca.*ones(1,length(Gain));
% Sca(2,2,:)=1./tca.*ones(1,length(Gain));
% Sca(1,2,:)=rca./tca;
% Sca(2,1,:)=rca./tca;

rbc = (nb-nc)./(nb+nc);
rcb = (nc-nb)./(nb+nc);
tbc = 2*nb./(nb+nc);
tcb = 2*nc./(nb+nc);

Sbc(1,1,:)=1./tbc.*ones(1,length(Gain));
Sbc(2,2,:)=1./tbc.*ones(1,length(Gain));
Sbc(1,2,:)=rbc./tbc;
Sbc(2,1,:)=rbc./tbc;

Scb(1,1,:)=1./tcb.*ones(1,length(Gain));
Scb(2,2,:)=1./tcb.*ones(1,length(Gain));
Scb(1,2,:)=rcb./tcb;
Scb(2,1,:)=rcb./tcb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj=1:length(Gain)
    % Unfortunately, I do not manage to do a multi-dimentionnal matrix product
    % without the loop. It would be much faster...
   
    DBRn(:,:,jj) = Pa(:,:,jj)*Sab(:,:,jj)*Pb(:,:,jj)*Sba(:,:,jj);
    DBRp(:,:,jj) = Pb(:,:,jj)*Sba(:,:,jj)*Pa(:,:,jj)*Sab(:,:,jj);
 
    %M(:,:,jj) =  DBRn(:,:,jj)^N_DBRn / Sba(:,:,jj) * Sbc(:,:,jj) * Pc(:,:,jj) * Scb(:,:,jj) * DBRp(:,:,jj)^N_DBRp / Sab(:,:,jj);
    M(:,:,jj) =  DBRn(:,:,jj)^N_DBRn / Sba(:,:,jj) * Sbc(:,:,jj) * Pc(:,:,jj) * Pgain(:,:,jj) * Pc(:,:,jj) * Scb(:,:,jj) * DBRp(:,:,jj)^N_DBRp / Sab(:,:,jj);
    
end

R = squeeze( abs( M(1,2,:)./M(2,2,:) ).^2 );
T = squeeze( na/na./abs(M(2,2,:)).^2 );

end
