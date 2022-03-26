function [rho_fo,rho_EE,rho_SE] = rho_fun_PSO(Phi,Q_m,Wm_source,Wm_MAT,H_Am,H_B,OPTZ,SNR,B,Nr,BW,UT,Nm,xi,P_max,Pcm,P_BS,P_RIS)

%       Wm_source
%         case 1 % random or external
%         case 2 % NP
%         case 3 % EB
%         case 4 % ZF
%         case 5 % MMSE
% H_Am=cell(1,UT);
H_m=cell(1,UT);
W_mn=cell(1,UT);
HW_m=cell(1,UT);
EW_m=cell(1,UT);
P_m=zeros(1,UT);
Penalty_W=zeros(1,UT);


% H_B=Ric_model(K_dB,B,Nr,L); % Channel RIS-BS

Phid=diag(Phi);

%% Data & Channel Generation for each UT
for m=1:UT
%     H_Am{m}=Ric_model(K_dB,Nr,Nm,L);        % Channel UT-RIS
    H_m{m} = (H_B*Phid*H_Am{m});             % Channel Matrix per UT   
    [U,~,V]=svd(H_m{m});    % SVD Decomposition
    switch Wm_source
        case 1 % random
            W_mn=mat2cell(Wm_MAT,Nm,ones(1,UT)*B);
            HW_m{m} = W_mn{m}*H_m{m};
        case 2 % NP
            HW_m{m} = H_m{m};
            W_mn{m}=0;
        case 3 % EB
            W_mn{m}=V/norm(V,'fro');
            HW_m{m} = U'*H_m{m}*W_mn{m};
        case 4 % ZF
            W_mn{m} = pinv(H_m{m}')/norm(pinv(H_m{m}'),'fro');
            HW_m{m} = W_mn{m}*H_m{m};
        case 5 % MMSE
            W_mn{m}=((1/sqrt(db2pow(-EbNoVec(n))/2)*eye(B)+H_m{m}*H_m{m}')^-1*H_m{m})';
            W_mn{m} = W_mn{m}/norm(W_mn{m},'fro');
            HW_m{m} = W_mn{m}*H_m{m}; 
    end

    EW_m{m}=HW_m{m}*Q_m(m)*HW_m{m}';            % Inner Product
    P_m(m)=xi*(real(trace(Q_m(m))));
    Penalty_W(m)=(norm(W_mn{m})-10^(P_max/10));
    Penalty_W4(m)=norm(W_mn{m});
    if Wm_source==2
        Penalty_W=Penalty_W*0;
    end
end
if Wm_source==1
    for uu=1:UT
        WW_norm(uu)=norm(W_mn{uu},2);
    end
   W_norm=mean(WW_norm);
else
    W_norm=0;
end

I_BW=eye(size(EW_m{1},1));
Penalty=(P_m-ones(1,UT)*10^(P_max/10));
rho_SE=abs(mean(log(det(I_BW+sum(cat(3,EW_m{:}),3)/sqrt(1/10^(SNR/10))))));
rho_SE_fo=-rho_SE+W_norm+1/2*sum(Penalty(Penalty>0).^2)+1/2*sum(Penalty_W(Penalty_W>0).^2);
% rho_SE_fo=W_norm+1/2*sum(Penalty(Penalty>0).^2)+1/2*sum(Penalty_W(Penalty_W>0).^2);

P_total=(sum(P_m))+UT*10^(Pcm/10)+10^(P_BS/10)+Nr*10^(P_RIS/10)/10;
rho_EE=abs(BW*rho_SE/P_total);
rho_EE_fo=-rho_EE+1/2*sum(Penalty(Penalty>0).^2)+100*sum(Penalty_W(Penalty_W>0).^2)+100*sum(Penalty_W4(Penalty_W4<0).^2);

switch OPTZ
       case 'EE'
            rho_fo = rho_EE_fo;
       case 'SE'
            rho_fo = rho_SE_fo;
end

end