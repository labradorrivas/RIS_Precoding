function [rho_fo,rho_EE,rho_SE] = rho_fun_PSO(H_Dm,H_B,Phi,H_Am,Q_m,CRITERIA,OPTZ,SNR,Nr,BW,UT,xi,P_max,Pcm,P_BS,P_RIS,varargin)

%       Wm_source
%         case 1 % NP
%         case 2 % EB
%         case 3 % ZF


EW_m=cell(1,UT);
P_m=zeros(1,UT);
Penalty_W=zeros(1,UT);
Penalty_W4=zeros(1,UT);
H_m=cell(1,UT);

for m=1:UT
       H_m{m}=(H_Dm{m}+H_B*diag(Phi)*H_Am{m});
end
if nargin==16
W_mn = precMAtrix(H_m,CRITERIA);
elseif nargin==17
    if CRITERIA~=3
     W_mn = mat2cell(varargin{1,1},size(H_Am{1},2),ones(1,UT)*size(H_Am{1},2));
    else
     W_mn = mat2cell(varargin{1,1},size(H_B,1),ones(1,UT)*size(H_Am{1},2));
    end
end

%% Data & Channel Generation for each UT
for m=1:UT
    if nargin==16
        [U,~,~]=svd(H_m{m});                     % SVD Decomposition
    switch CRITERIA
        case 1 % NP
            EW_m{m}=H_m{m}*Q_m(m)*H_m{m}';              % Inner Product
            P_m(m)=xi*(real(trace(Q_m(m))));
        case 2 % EB
            EW_m{m}=U'*H_m{m}*W_mn{m}*Q_m(m)*(U'*H_m{m}*W_mn{m})';     % Inner Product
            P_m(m)=xi*norm(W_mn{m});
        case 3 % ZF
            EW_m{m}=H_m{m}'*W_mn{m}*Q_m(m)*(H_m{m}'*W_mn{m})';             % Inner Product
            P_m(m)=xi*norm(W_mn{m});
    end
    elseif nargin==17
        if CRITERIA==3
            EW_m{m}=H_m{m}'*W_mn{m}*Q_m(m)*(H_m{m}'*W_mn{m})';             % Inner Product
        else
            EW_m{m}=H_m{m}*W_mn{m}*Q_m(m)*(H_m{m}*W_mn{m})';             % Inner Product
        end
            P_m(m)=xi*norm(W_mn{m});
    end
    

    
    Penalty_W(m)=(norm(W_mn{m})-10^(P_max/10));
    Penalty_W4(m)=norm(W_mn{m});
    if CRITERIA==1
        Penalty_W=Penalty_W*0;
    end
end


I_BW=eye(size(EW_m{1},1));
Penalty=(P_m-ones(1,UT)*10^(P_max/10));
rho_SE=abs(mean(log(det(I_BW+sum(cat(3,EW_m{:}),3)/sqrt(1/10^(SNR/10))))));
rho_SE_fo=-rho_SE+1/2*sum(Penalty(Penalty>0).^2)+1/2*sum(Penalty_W(Penalty_W>0).^2);

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