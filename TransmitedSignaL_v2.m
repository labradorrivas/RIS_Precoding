function [txSig,Q_m] = TransmitedSignaL_v2(UT,Nm,M_order,modtype)

% UT number of users
% Nm number of antennas
% M_order  Modulation order
% modtype 'PSK' or 'QAM'
pskModulator   = comm.PSKModulator(...
            'ModulationOrder',  2^M_order, ...
            'PhaseOffset',      0, ...
            'BitInput',         true);
modData=zeros(Nm,UT);        
switch modtype
       case 'PSK'
            data = randi([0 1],[Nm*M_order, UT]);
            for ut=1:UT
                modData(:,ut)=pskModulator(data(:,1));
            end
       case 'QAM'
            data = randi([0 3],[Nm*M_order, UT]);
            modData = qammod(data,M_order);
 end

txSig=modData;
Q_m=zeros(1,UT);

    for m=1:UT
    
    Q_m(m)=cov(txSig(:,m));             % Covariance Matrix for each UT
    end
end