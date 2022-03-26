function H=Ric_model(K_dB,NT,NR,L)
% Rician Channel Model
%   Input:
%       K_dB   : K factor [dB]
%       L      : # of channel realization
%       NT    : number of transmitters
%       NR    : number of receivers
%   Output:
%       h      : channel vector

K=10^(K_dB/10);
H = sqrt(K/(K+1)) + sqrt(1/(K+1))*Ray_model(NT,NR,L);