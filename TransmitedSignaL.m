function [txSig,Q_m] = TransmitedSignaL(UT,Nm,dsize,M_order)

% UT number of users
% Nm number of antennas
% dsize Data Size
% M_order QAM Modulation order

data = randi([0 3],dsize,UT);
% modData = pskmod(data,M_order,pi/M_order);
modData = qammod(data,M_order);
% Encoder
ostbc = comm.OSTBCEncoder('NumTransmitAntennas',Nm,'SymbolRate',1/2);
txSig=cell(1,UT);
Q_m=cell(1,UT);

    for m=1:UT
    txSig{m} = ostbc(modData(:,m));         % Transmited SignaL
    Q_m{m}=cov(txSig{m});                   % Covariance Matrix for each UT
    end
end