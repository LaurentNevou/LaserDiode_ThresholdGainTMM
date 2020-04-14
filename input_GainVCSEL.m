%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% VCSEL cavity input for TMM model %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;             %% speed of light [m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda0 = 1000e-9;          %% Cavity Central wavelength
n1 = 3;                     %% AlAs
n2 = 3.6;                   %% GaAs
nc = 3.6;                   %% refractive index of the cavity
N_DBRn=30;                  %% amount of DBR n-doped pairs
N_DBRp=20;                  %% amount of DBR p-doped pairs
LQW = 10e-9;                %% quantum well thickness in which the gain will be [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=2*pi*c./lambda;           %% transformation of lambda in pulsation
kk=Gain/2./w*c;
nGain=nc-1i*kk;             %% HERE IS THE MAJOR CHANGE!!! from PLUS (LOSSES) to MINUS (GAIN)

% TAKE CARE!!! In order to have a node of the electrical field on the QW, 
% => XX must be odd  (1, 3, 5, ...) if n1<n2
% => XX must be even (2, 4, 6, ...) if n1>n2
XX=3;
lc= XX * lambda0/(2*nc);     
spacer=(lc-LQW)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l1=lambda0/(4*abs(n1));      %% DBR layer-1 thickness at lambda/4;
l2=lambda0/(4*abs(n2));      %% DBR layer-2 thickness at lambda/4;

DBR_n=[]; DBRn=[ l1 n1*ones(1,length(nGain)) ; l2 n2*ones(1,length(nGain)) ];
DBR_p=[]; DBRp=[ l2 n2*ones(1,length(nGain)) ; l1 n1*ones(1,length(nGain)) ];

for jj=1:N_DBRn
  DBR_n = [ DBR_n ; DBRn ];
end
for jj=1:N_DBRp
  DBR_p = [ DBR_p ; DBRp ];
end

layer=[ 

DBR_n 
spacer  nc*ones(1,length(nGain))
LQW     nGain 
spacer  nc*ones(1,length(nGain))
DBR_p 

];

nL=layer(1,2);
nR=layer(end,2);
