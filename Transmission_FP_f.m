function[T,R]=Transmission_FP_f(lambda,Gain,L,n1,n2,n3)

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
  
kk=Gain/2/w*c;              %% imaginary part of the optical refrative index
nn2=n2-1i*kk;  %% HERE IS THE MAJOR CHANGE!!! from PLUS (LOSSES) to MINUS (GAIN)

P = zeros(2,2,length(Gain));
P(1,1,:) = exp(+1i*k*nn2*L);
P(2,2,:) = exp(-1i*k*nn2*L);

r12 = (n1-nn2)./(n1+nn2);
r23 = (nn2-n3)./(nn2+n3);

t12 = 2*n1./(n1+nn2);
t23 = 2*nn2./(nn2+n3);

S12(1,1,:)=1./t12.*ones(1,length(Gain));
S12(2,2,:)=1./t12.*ones(1,length(Gain));
S12(1,2,:)=r12./t12;
S12(2,1,:)=r12./t12;

S23(1,1,:)=1./t23.*ones(1,length(Gain));
S23(2,2,:)=1./t23.*ones(1,length(Gain));
S23(1,2,:)=r23./t23;
S23(2,1,:)=r23./t23;

for jj=1:length(Gain)
    % Unfortunately, I do not manage to do a multi-dimentionnal matrix product
    % without the loop. It would be much faster...
    M(:,:,jj)=S12(:,:,jj)*P(:,:,jj)*S23(:,:,jj);

end

R = squeeze( abs( M(1,2,:)./M(2,2,:) ).^2 );
T = squeeze( n3/n1./abs(M(2,2,:)).^2 );
