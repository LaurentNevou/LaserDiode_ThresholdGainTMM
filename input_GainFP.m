%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% FabryPerot cavity input for TMM model %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;             %% speed of light [m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda0 = 1000e-9;          %% Cavity Central wavelength [m]
Lc=1e-3;                    %% Cavity length [m]
nopt=3.55;                  %% Optical index of the ridge emitter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=2*pi*c./lambda;           %% transformation of lambda in pulsation
kk=Gain/2./w*c;
n3=nopt-1i*kk;              %% HERE IS THE MAJOR CHANGE!!! from PLUS (LOSSES) to MINUS (GAIN)

L0 = 1*lambda0./(2*abs(mean(n3)));   %% in order to get 1 FP mode @lambda0
N=Lc/L0;
L = round(N) * L0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

layer=[ 
1e-6 ones(1,length(n3))
L n3 
1e-6 ones(1,length(n3))
];

nL=layer(1,2);
nR=layer(end,2);
