function H=Ray_model(NT,NR,L)
% uncorrelated Rayleigh fading channel Model, CN(1,0)
%  Input :  L  : # of channel realization
%           NT    : number of transmitters
%           NR    : number of receivers
%  Output: H  : Channel vector


H = mean((randn(NT,NR,L)+1i*randn(NT,NR,L))/sqrt(2),3);
