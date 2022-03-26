K_dB=-10:5:20;    % K factor [dB] for the Rician Channel Model
B=8;       % Number of Antennas at BS
Nr=32;      % Number of RIS elements
L=100;     % Number of channel realization
% s = [randn(100,1); randn(100,1)+k];

for i=1:length(K_dB)
H=Ric_model(K_dB(i),B,Nr,L);
% histogram(abs(mean(H)),'Normalization','probability')
[f,x] = ksdensity(abs(mean(H)));
hold on
plot(x,f)
end