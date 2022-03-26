%% Angel Labrador
clc;
clear;

K_dB=3;    % K factor [dB] for the Rician Channel Model
SNR=20;     % SNR [dB]
Nm=4;       % Number of Antennas per UT
B=64;       % Number of Antennas at BS
Nr=32;      % Number of RIS elements
L=1000;     % Number of channel realization
UT=4;       % Number of Users
BW=10;      % Transmission Bandwidth [MHz]

% Data Generation
dsize=1000; % Data Size
% Generate QPSK-modulated data.
data = randi([0 3],dsize,UT);
M_order=64;
% modData = pskmod(data,M_order,pi/M_order);
modData = qammod(data,M_order);
% Encoder
ostbc = comm.OSTBCEncoder('NumTransmitAntennas',Nm,'SymbolRate',1/2);
txSig=cell(1,UT);
Q_m=cell(1,UT);
H_Am=cell(1,UT);
H_m=cell(1,UT);
W_m=cell(1,UT);
W_mn=cell(1,UT);
E_m=zeros(B,B,UT);
EW_m=zeros(Nm,Nm,UT);
P_m=zeros(1,UT);

% Power Constant elements
xi=0.3;       % efficiency of the transmit power amplifiers adopted at UT $m$
Pcm=10;                  % Power circutry each user [dbm]
P_BS=39;                % Power at BS [dbm]
P_RIS=5;                % Power of each RIS element [dbm]
%% Random Generate RIS Angle
theta = (2*pi)*rand(1,Nr); %Reflecting Angle RIS
Phi = diag(exp(1i*theta)); %phase shift Matrix for RIS, % Set Amp = 1
H_B=Ric_model(K_dB,B,Nr,L); % Channel RIS-BS

%% Data & Channel Generation for each UT
for m=1:UT
    txSig{m} = ostbc(modData(:,m));         % Transmited Signal
    Q_m{m}=cov(txSig{m});                   % Covariance Matrix for each UT
    H_Am{m}=Ric_model(K_dB,Nr,Nm,L);        % Channel UT-RIS
    H_m{m} = (H_B*Phi*H_Am{m});             % Channel Matrix per UT
    W_m{m} = pinv(H_m{m});                  % Precoding Matrix
    W_mn{m} = W_m{m}./vecnorm(W_m{m},2,1);  % Normalized
    E_m(:,:,m)=H_m{m}*Q_m{m}*H_m{m}';              % Inner Product
    EW_m(:,:,m)=W_mn{m}*E_m(:,:,m)*W_mn{m}';        % Inner Product with Precoding
    P_m(m)=xi*trace(Q_m{m});
end
I_B=eye(B);
I_BW=eye(Nm);

rho_SE=(mean(log(det(I_B+sum(E_m,3)/sqrt(1/10^(SNR/10))))));
rho_SEW=(mean(log(det(I_BW+sum(EW_m,3)/sqrt(1/10^(SNR/10))))));
P_total=sum(P_m)+UT*Pcm+P_BS+Nr*P_RIS;
rho_EE=BW*rho_SE/P_total;
rho_EEW=BW*rho_SEW/P_total;

