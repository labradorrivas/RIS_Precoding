%% Angel Labrador
clc;
clear all;

%% Plot IEEE Parameters
    set(0,'DefaultTextFontName','Times','DefaultTextFontSize',14,...
     'DefaultAxesFontName','Times','DefaultAxesFontSize',12,...
     'DefaultLineLineWidth',2,'DefaultLineMarkerSize',8)
    Color=[0    0  0;
            77  45  82;
            192 58  46;
            44  129 184;
            242 156 19;
            155 186 92;
            25  160 131;
            47  64  80]/255;
        
%% Scenario Data
iterations=100;
NL=100;       % # Realizations          
K_dB=3;     % K factor [dB] for the Rician Channel Model
SNR=0;      % SNR [dB]
Nm=[8];      % Number of Antennas per UT
B=8;       % Number of Antennas at BS
Nr=32;      % Number of RIS elements
L=1;      % Number of channel realization
UT=4;       % Number of Users
BW=100E6;     % Transmission Bandwidth [Hz]

% Power Constant elements
xi=1;                 % efficiency of the transmit power amplifiers adopted at UT $m$
Pcm=10;                 % Power circutry each user [dbm]
P_BS=39;                % Power at BS [dbm] 
    % a maximum output power limit of 38 dBm for medium range BSs, 24 dBm for local area BSs, and of 20 dBm for home BSs
P_RIS=5;                % Power of each RIS element [dbm]
P_max=(-10:5:35);      % Max power per UT [dbm]

% Data Generation
dsize=1;             % Data Size
M_order = 2;           % Modulation order
k = log2(M_order);     % Bits per symbol

%% 
time_wo=0;
time_svd=0;
time_zf=0;

for nm=1:length(Nm)
    
Wm_MAT=rand(Nm(nm),UT*B)+1i*rand(Nm(nm),UT*B);
Wm_MAT=Wm_MAT/norm(Wm_MAT,'fro');
rho_EE=zeros(length(P_max),4,NL,iterations);
rho_EEW=zeros(length(P_max),4,NL,iterations);
rho_SE=zeros(length(P_max),4,NL,iterations);
rho_SEW=zeros(length(P_max),4,NL,iterations);
rho_EEWZF=zeros(length(P_max),4,NL,iterations);
rho_SEWZF=zeros(length(P_max),4,NL,iterations);

W_norm=zeros(length(P_max),UT,NL,iterations);



    
    options.verbosity=0;
    options.stopfun = @mystopfun;
    options.maxiter = 100;
    options.nostalgia = 1.1;
    options.social = 1.1;


for nn=1:NL
    [txSig,Q_m]=TransmitedSignaL_v2(UT,Nm(nm),M_order,'PSK');
    diag_Q_m=diag(Q_m);
    %% Random Generate RIS Angle
theta = (2*pi)*rand(1,Nr);  % Reflecting Angle RIS % ones(1,Nr);%
Phi = (exp(1i*theta));      % Phase shift Matrix for RIS, % Set Amp = 1 ones(1,Nr);%
    H_Am=cell(1,UT);
    H_B=Ric_model(K_dB,B,Nr,L); % Channel RIS-BS
    for m=1:UT
    H_Am{m}=Ric_model(K_dB,Nr,Nm(nm),L);        % Channel UT-RIS
    end
for pp=1:length(P_max)
    Qm0={};
    Qm0{1,1} = ones(1,length(Q_m))*10^(P_max(pp)/10)/xi;
    Qm_init= Qm0{1,1};
    Qm0W=Qm0;
    Qm0WZF=Qm0;

    Phi0={};
    Phi0{1,1} = Phi';
    Phi0W=Phi0;
    Phi0WZF=Phi0;

    Wm0={};
    Wm0{1,1} = Wm_MAT/norm(Wm_MAT,'fro');
    Wm0W=Wm0;
    Wm0WZF=Wm0;

    %% Initial Random solution AND Max Power Transmit
    % Criterio w/o prec
    [~,rho_EE(pp,1,nn,:),rho_SE(pp,1,nn,:)]= rho_fun_PSO(Phi,Qm_init,2,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    % Criterio B
    [~,rho_EEW(pp,1,nn,:),rho_SEW(pp,1,nn,:)]= rho_fun_PSO(Phi,Qm_init,3,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    % Criterio prec-ZF
    [~,rho_EEWZF(pp,1,nn,:),rho_SEWZF(pp,1,nn,:)]= rho_fun_PSO(Phi,Qm_init,4,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    
for optz_iter=1:iterations


    %% Q_m Optimization
    tic
    % Criterio w/o prec
    problem1.M = positivefactory(1,UT); 
    problem1.cost  = @(x) (rho_fun_PSO(Phi,x,2,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS));
    [Q_m_new, ~, ~, ~] = pso(problem1, Qm0, options);
    Qm0=mat2cell(Q_m_new,1);
    [~,rho_EE(pp,2,nn,optz_iter),rho_SE(pp,2,nn,optz_iter)]=rho_fun_PSO(Phi,Q_m_new,2,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    time_wo=toc+time_wo;
    % Criterio B
    tic
    problem1W.M = positivefactory(1,UT);
    problem1W.cost  = @(x) rho_fun_PSO(Phi,x,3,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    [Q_m_newW, ~, ~, ~] = pso(problem1W, Qm0W, options);
    Qm0W=mat2cell(Q_m_newW,1);
    [~,rho_EEW(pp,2,nn,optz_iter),rho_SEW(pp,2,nn,optz_iter)]=rho_fun_PSO(Phi,Q_m_newW,3,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    time_svd=toc+time_svd;
    % Criterio prec-ZF
    tic
    problem1WZF.M = positivefactory(1,UT);
    problem1WZF.cost  = @(x) rho_fun_PSO(Phi,x,4,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    [Q_m_newWZF, ~, ~, ~] = pso(problem1WZF, Qm0WZF, options);
    Qm0WZF=mat2cell(Q_m_newWZF,1);
    [~,rho_EEWZF(pp,2,nn,optz_iter),rho_SEWZF(pp,2,nn,optz_iter)]=rho_fun_PSO(Phi,Q_m_newWZF,4,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    time_zf=toc+time_zf;
%%    Phi Optimization
    

    % Criterio w/o prec
    tic
    problem2.M = complexcirclefactory(Nr);
    problem2.cost  = @(x) rho_fun_PSO(x,Q_m_new,2,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);    
    [theta_new, ~, ~, ~] = pso(problem2, Phi0, options);
    Phi0=mat2cell(theta_new,Nr);
    [~,rho_EE(pp,3,nn,optz_iter),rho_SE(pp,3,nn,optz_iter)]=rho_fun_PSO(theta_new,Q_m_new,2,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);    
    time_wo=toc+time_wo;
    % Criterio B
    tic
    problem2W.M = complexcirclefactory(Nr);
    problem2W.cost  = @(x) rho_fun_PSO(x,Q_m_newW,3,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);    
    [theta_newW, ~, ~, ~] = pso(problem2W, Phi0W, options);
        Phi0W=mat2cell(theta_newW,Nr);
    [~,rho_EEW(pp,3,nn,optz_iter),rho_SEW(pp,3,nn,optz_iter)]=rho_fun_PSO(theta_newW,Q_m_newW,3,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);    
    time_svd=toc+time_svd;
    % Criterio prec-ZF
    tic
    problem2WZF.M = complexcirclefactory(Nr);
    problem2WZF.cost  = @(x) rho_fun_PSO(x,Q_m_newWZF,4,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);    
    [theta_newWZF, ~, ~, ~] = pso(problem2WZF, Phi0WZF, options);
        Phi0WZF=mat2cell(theta_newWZF,Nr);
    [~,rho_EEWZF(pp,3,nn,optz_iter),rho_SEWZF(pp,3,nn,optz_iter)]=rho_fun_PSO(theta_newWZF,Q_m_newWZF,4,Wm_MAT,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);    
    time_zf=toc+time_zf;
  %% W optz
    Wm0={};
    Wm0{1,1} = Wm_MAT/norm(Wm_MAT,'fro');
    
    % Criterio w/o prec
    tic
    problem3.M = spherecomplexfactory(Nm(nm),UT*B); % obliquecomplexfactory  euclideancomplexfactory
    problem3.cost  = @(x) rho_fun_PSO(theta_new,Q_m_new,1,x,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    [W_m_new, ~, ~, ~] = pso(problem3, Wm0, options);
    Wm0=mat2cell(W_m_new,B,Nr);
    W_m_cell=mat2cell(W_m_new,Nm(nm),ones(1,UT)*B);
    for uu=1:UT
        W_norm(pp,uu,nn,optz_iter)=norm(W_m_cell{uu},2);
    end
    [~,rho_EE(pp,4,nn,optz_iter),rho_SE(pp,4,nn,optz_iter)]=rho_fun_PSO(theta_new,Q_m_new,1,W_m_new,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    time_wo=toc+time_wo;
    % Criterio B
    tic
    problem3W.M = spherecomplexfactory(Nm(nm),UT*B); % obliquecomplexfactory  euclideancomplexfactory
    problem3W.cost  = @(x) rho_fun_PSO(theta_newW,Q_m_newW,1,x,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    [W_m_newW, ~, ~, ~] = pso(problem3W, Wm0W, options);
    Wm0W=mat2cell(W_m_newW,B,Nr);
    W_m_cell=mat2cell(W_m_newW,Nm(nm),ones(1,UT)*B);
    for uu=1:UT
        WW_norm(pp,uu,nn,optz_iter)=norm(W_m_cell{uu},2);
    end
    [~,rho_EEW(pp,4,nn,optz_iter),rho_SEW(pp,4,nn,optz_iter)]=rho_fun_PSO(theta_newW,Q_m_newW,1,W_m_newW,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    time_svd=toc+time_svd;
    % Criterio prec-ZF
    tic
    problem3WZF.M = spherecomplexfactory(Nm(nm),UT*B); % obliquecomplexfactory  euclideancomplexfactory
    problem3WZF.cost  = @(x) rho_fun_PSO(theta_newWZF,Q_m_newWZF,1,x,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    [W_m_newWZF, ~, ~, ~] = pso(problem3WZF, Wm0WZF, options);
    Wm0WZF=mat2cell(W_m_newWZF,B,Nr);
    W_m_cell=mat2cell(W_m_newWZF,Nm(nm),ones(1,UT)*B);
    for uu=1:UT
        WZF_norm(pp,uu,nn,optz_iter)=norm(W_m_cell{uu},2);
    end
    [~,rho_EEWZF(pp,4,nn,optz_iter),rho_SEWZF(pp,4,nn,optz_iter)]=rho_fun_PSO(theta_newWZF,Q_m_newWZF,1,W_m_newWZF,H_Am,H_B,'SE',SNR,B,Nr,BW,UT,Nm(nm),xi,P_max(pp),Pcm,P_BS,P_RIS);
    time_zf=toc+time_zf;
end

end
end

rho_EE_mean=mean(rho_EE(:,:,:,end),3);
rho_EEW_mean=mean(rho_EEW(:,:,:,end),3);
rho_SE_mean=mean(rho_SE(:,:,:,end),3);
rho_SEW_mean=mean(rho_SEW(:,:,:,end),3);
rho_EEWZF_mean=mean(rho_EEWZF(:,:,:,end),3);
rho_SEWZF_mean=mean(rho_SEWZF(:,:,:,end),3);
W_norm_max1=max(max(W_norm(:,:,:,end),[],3),[],2);
W_norm_max=W_norm_max1/max(W_norm_max1);

time_wo=time_wo/iterations/NL;
time_svd=time_svd/iterations/NL;
time_zf=time_zf/iterations/NL;

save(strcat('optz_EE_xi_',string(xi) ,'_',datestr(now,'yyyymmddHHMMSS')));
% 
% fig_ee=figure;
% fig_ee.Units='normalized';
% axes('NextPlot','replacechildren', 'ColorOrder',Color([2 3 4 5 6],:)); 
%     fig_ee=plot(P_max,rho_EE_mean(:,4),'-o',...
%         P_max,rho_EEW_mean(:,4),'-v',...
%         P_max,rho_EEWZF_mean(:,4),'-s',...
%         P_max,mean(rho_EE(:,1,:,1),3),'-p',...
%         P_max,mean(rho_EE(:,3,:,1),3),'->');
%     set(gca,'XMinorTick','on','YMinorTick','on');
%         set(fig_ee,'MarkerFaceColor','w');
%     grid on
% xlabel('$ P_{\textrm{max},m}$ (dBm)','Interpreter','latex');
% ylabel('EE  (bits/Joule)');
% % plot(rho_EEW_mean)
% % hold on
% % plot(rho_SEW_mean)
% legend('${\textsc{ee}^0_{\forall}(\vec{Q}^\iota_m,\vec{\Phi}^{\iota},\vec{W}^\iota)}$',...
%     '${\textsc{ee}^{SVD}_{\forall}(\vec{Q}^\iota_m,\vec{\Phi}^{\iota},\vec{W}^\iota)}$',...
%     '${\textsc{ee}^{ZF}_{\forall}(\vec{Q}^\iota_m,\vec{\Phi}^{\iota},\vec{W}^\iota)}$',...
%     '${\textsc{ee}^0_{\forall}(\vec{Q}^0_m,\vec{\Phi}^{0},\vec{W}^0)}$',...
%     '${\textsc{ee}^0_{\forall}(\vec{Q}^0_m,\vec{\Phi}^{\iota},\vec{W}^0)}$','Location','best','interpreter','latex');
% set(gcf,'renderer','Painters')
%  saveas(gcf,'EE_optz_EE_xi_03','epsc');
%  saveas(gcf,'EE_optz_EE_xi_03','fig');
% 
% 
% 
% 
% fig2=figure;
% fig2.Units='normalized';
% axes('NextPlot','replacechildren', 'ColorOrder',Color([2 2 2 2 3 3 3 3 4 4 4 4],:)); 
%     plot(P_max,rho_SE_mean(:,1),'--s',...
%         P_max,rho_SE_mean(:,2),'--d',...
%         P_max,rho_SE_mean(:,3),'--o',...
%         P_max,rho_SE_mean(:,4),'--x',...
%         P_max,rho_SEW_mean(:,1),':p',...
%         P_max,rho_SEW_mean(:,2),':>',...
%         P_max,rho_SEW_mean(:,3),':+',...
%         P_max,rho_SEW_mean(:,4),':v',...
%         P_max,rho_SEWZF_mean(:,1),'-.h',...
%         P_max,rho_SEWZF_mean(:,2),'-.<',...
%         P_max,rho_SEWZF_mean(:,3),'-.*',...
%         P_max,rho_SEWZF_mean(:,4),'-.^');
%     set(gca,'XMinorTick','on','YMinorTick','on');
%     grid on
% xlabel('$ P_{\textrm{max},m}$ (dBm)','Interpreter','latex');
% ylabel('SE  (bits/Joule)');
% % plot(rho_EEW_mean)
% % hold on
% % plot(rho_SEW_mean)
% legend('$\textsc{se}_{\vec{W}}(\textrm{Rndm})$ w/o prec','$\textsc{se}_{\vec{W}}(\vec{Q}_m)$ w/o prec','$\textsc{se}_{\vec{W}}(\vec{\Phi})$ w/o prec','$\textsc{se}_{\forall}(\vec{W})$ w/o prec',...
%     '$\textsc{se}_{\vec{W}}(\textrm{Rndm})$  prec-SVD','$\textsc{se}_{\vec{W}}(\vec{Q}_m)$ prec-SVD',...
%     '$\textsc{se}_{\vec{W}}(\vec{\Phi})$ prec-SVD','$\textsc{se}_{\forall}(\vec{W})$ prec-SVD',...
%     '$\textsc{se}_{\vec{W}}(\textrm{Rndm})$  prec-ZF','$\textsc{se}_{\vec{W}}(\vec{Q}_m)$  prec-ZF',...
%     '$\textsc{se}_{\vec{W}}(\vec{\Phi})$ prec-ZF','$\textsc{se}_{\forall}(\vec{W})$ prec-ZF','Location','best','interpreter','latex');
% set(gcf,'renderer','Painters')
% saveas(gcf,strcat('SE_optz_EE_xi_',string(xi)),'epsc');
% 
% 
% fig4=figure;
% fig4.Units='normalized';
% axes('NextPlot','replacechildren', 'ColorOrder',Color([1 2 3 4 5 6 7 8],:)); 
%     plot(rho_SE_mean(:,1),rho_EE_mean(:,1),'--s',...%        P_max,rho_SE_mean(:,2),'--x',...
%         rho_SE_mean(:,3),rho_EE_mean(:,3),'--o',...
%         rho_SE_mean(:,4),rho_EE_mean(:,4),'--x',...
%         rho_SEW_mean(:,1),rho_EEW_mean(:,1),':p',...%        P_max,rho_SEW_mean(:,2),':h',...
%         rho_SEW_mean(:,3),rho_EEW_mean(:,3),':+',...
%         rho_SEW_mean(:,4),rho_EEW_mean(:,4),':v',...
%         rho_SEWZF_mean(:,1),rho_EEWZF_mean(:,1),'-.h',...
%         rho_SEWZF_mean(:,3),rho_EEWZF_mean(:,3),'-.*',...
%         rho_SEWZF_mean(:,4),rho_EEWZF_mean(:,4),'-.^');
%     set(gca,'XMinorTick','on','YMinorTick','on');
%     grid on
% xlabel('SE  (bits/s/Hz)');
% ylabel('EE  (bits/Joule)');
% legend('Random w/o prec','Theta w/o prec','$\mathbf{W}_k$ w/o prec',...
%     'Random prec-SVD',...
%     'Theta prec-SVD','$\mathbf{W}_k$ prec-SVD',...
%     'Random prec-ZF',...
%     'Theta prec-ZF','$\mathbf{W}_k$ prec-ZF','Location','best','interpreter','latex');
% set(gcf,'renderer','Painters')
% saveas(gcf,strcat('SE_VS_EE_OPTZ_EE_xi_',string(xi)),'epsc');
% 
% 
% fig3=figure;
% fig3.Units='normalized';
% axes('NextPlot','replacechildren', 'ColorOrder',Color(2,:)); 
%     bar(P_max,W_norm_max);
%     set(gca,'XMinorTick','on','YMinorTick','on');
%     grid on
% xlabel('$ P_{\textrm{max},m}$ (dBm)','Interpreter','latex');
% ylabel('$\max{||\mathbf{W}_k||}$','interpreter','latex');
% set(gcf,'renderer','Painters')
% saveas(gcf,strcat('normW_optz_EE_xi_',string(xi)),'epsc'); % ,'_',datestr(now,'yyyymmddHHMMSS')
% 
% 
% rho_Nm_W(:,(2*nm-1):(2*nm))=[rho_SEW_mean(:,4) rho_EEW_mean(:,4)];
% 
  markers=["-s";"-o";"-x";"-p";"-+";"-v";"-h";"-*";":^";":>"];


figure
subplot(1,2,1)
% axes('NextPlot','replacechildren', 'ColorOrder',Color([2 2 2 2 3 3 3 3 4 4 4 4],:));
rand_loc=randi(iterations,1,length(P_max)) ;
for pm=1:length(P_max)
     value_iter=squeeze(mean(rho_EE(pm,4,:,:),3));
     plot(value_iter,markers(pm));
     text(rand_loc(pm),value_iter(rand_loc(pm)),"$ P_{\textrm{max},m}$ = " + string(P_max(pm))+" dBm",'interpreter','latex','FontSize',10,'BackgroundColor','w')
     hold on
 end
    hold off
    set(gca,'XMinorTick','on','YMinorTick','on');
    grid on
    %set(gca, 'YScale', 'log')
% legendStrings = "$ P_{\textrm{max},m}$ = " + string(P_max)+" dBm";
xlabel('Iteration #','Interpreter','latex');
ylabel(' ${\textsc{ee}^0_{\forall}(\vec{Q}^\iota_m,\vec{\Phi}^{\iota},\vec{W}^\iota)}$ (bits/Joule)','interpreter','latex');
% plot(rho_EEW_mean)
% hold on
% plot(rho_SEW_mean)
% legend(legendStrings,'Location','bestoutside','interpreter','latex');

subplot(1,2,2)
% axes('NextPlot','replacechildren', 'ColorOrder',Color([2 2 2 2 3 3 3 3 4 4 4 4],:)); 
 for pm=1:length(P_max)
     value_iter=squeeze(mean(rho_EEW(pm,4,:,:),3));
     plot(value_iter,markers(pm));
     text(rand_loc(pm),value_iter(rand_loc(pm)),"$ P_{\textrm{max},m}$ = " + string(P_max(pm))+" dBm",'interpreter','latex','FontSize',10,'BackgroundColor','w')
     hold on
 end
    hold off
    set(gca,'XMinorTick','on','YMinorTick','on');
    grid on
    %set(gca, 'YScale', 'log')
% legendStrings = "$ P_{\textrm{max},m}$ = " + string(P_max)+" dBm";
xlabel('Iteration #','Interpreter','latex');
ylabel(' ${\textsc{ee}^\textrm{SVD}_{\forall}(\vec{Q}^\iota_m,\vec{\Phi}^{\iota},\vec{W}^\iota)}$ (bits/Joule)','interpreter','latex');
% plot(rho_EEW_mean)
% hold on
% plot(rho_SEW_mean)
% legend(legendStrings,'Location','bestoutside','interpreter','latex');

set(gcf,'renderer','Painters')
saveas(gcf,'EE_optz_EE_periter','epsc');
saveas(gcf,'EE_optz_EE_periter','fig');
end



fig5=figure;
fig5.Units='normalized';
axes('NextPlot','replacechildren', 'ColorOrder',Color(1:length(Nm),:)); 
for nm=1:length(Nm)
    plot(P_max,rho_Nm_W(:,(2*nm)),markers(nm));
    hold on
end
set(gca,'XMinorTick','on','YMinorTick','on')
hold off
grid on
xlabel('$ P_{\textrm{max},m}$ (dBm)','Interpreter','latex');
ylabel('$\textsc{{ee}}^{\textsc{svd}}_{\forall}(\vec{Q}^\iota_m,\vec{\Phi}^{\iota},\vec{W}^\iota)$  (bits/Joule)','Interpreter','latex');
legend('$N_m=2$','$N_m=4$','$N_m=8$','$N_m=16$','$N_m=32$','Location','best','interpreter','latex');
set(gcf,'renderer','Painters')
saveas(gcf,strcat('EE_VS_Nm_optz_EE_xi_',string(xi)),'epsc'); % ,'_',datestr(now,'yyyymmddHHMMSS')
