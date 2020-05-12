function[T,R]=Transmission_DFB_f(lambda,Gain,L,LL,n1,n2,Dn2,n3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Optoelectronic", Cambridge Books Online
% Prof. Emmanuel Rosencher,
% Complement to Chapter 9
% 9.D Fabry–Perot cavities and Bragg reflectors, page 434
% http://dx.doi.org/10.1017/CBO9780511754647.028
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;             %% speed of light [m/s]
k=2*pi/lambda;              %% wave vector    [m-1]
w=2*pi*c/lambda;            %% pulsation      [s-1]
N=floor(L/LL);

kk=Gain/2/w*c;              %% imaginary part of the optical refrative index
n2a = n2     -1i*kk;  %% HERE IS THE MAJOR CHANGE!!! from PLUS (LOSSES) to MINUS (GAIN)
n2b = n2-Dn2 -1i*kk;  %% HERE IS THE MAJOR CHANGE!!! from PLUS (LOSSES) to MINUS (GAIN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DFB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = zeros(2,2,length(Gain));
P(1,1,:) = exp(+1i*k*n2a*LL/2);
P(2,2,:) = exp(-1i*k*n2a*LL/2);

r2ab = (n2a-n2b)./(n2a+n2b);
r2ba = (n2b-n2a)./(n2a+n2b);

t2ab = 2*n2a./(n2a+n2b);
t2ba = 2*n2b./(n2a+n2b);

S2ab(1,1,:) = (1./t2ab) .* ones(1,length(Gain));
S2ab(2,2,:) = (1./t2ab) .* ones(1,length(Gain));
S2ab(1,2,:) = r2ab./t2ab;
S2ab(2,1,:) = r2ab./t2ab;

S2ba(1,1,:) = (1./t2ba) .* ones(1,length(Gain));
S2ba(2,2,:) = (1./t2ba) .* ones(1,length(Gain));
S2ba(1,2,:) = r2ba./t2ba;
S2ba(2,1,:) = r2ba./t2ba;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fabry-Perot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r12 = (n1-n2a)./(n1+n2a);
r23 = (n2a-n3)./(n2a+n3);

t12 = 2*n1./(n1+n2a);
t23 = 2*n2a./(n2a+n3);

S12(1,1,:)=1./t12.*ones(1,length(Gain));
S12(2,2,:)=1./t12.*ones(1,length(Gain));
S12(1,2,:)=r12./t12;
S12(2,1,:)=r12./t12;

S23(1,1,:)=1./t23.*ones(1,length(Gain));
S23(2,2,:)=1./t23.*ones(1,length(Gain));
S23(1,2,:)=r23./t23;
S23(2,1,:)=r23./t23;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj=1:length(Gain)
    % Unfortunately, I do not manage to do a multi-dimentionnal matrix product
    % without the loop. It would be much faster...
    
    M(:,:,jj) = P(:,:,jj)*S2ba(:,:,jj)*P(:,:,jj)*S2ab(:,:,jj);
    M(:,:,jj) = S12(:,:,jj) * M(:,:,jj)^N * S23(:,:,jj);
end

R = squeeze( abs( M(1,2,:)./M(2,2,:) ).^2 );
T = squeeze( n3/n1./abs(M(2,2,:)).^2 );

end
