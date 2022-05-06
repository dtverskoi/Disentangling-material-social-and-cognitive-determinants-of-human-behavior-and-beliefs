%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum likelihood procedure for the alternative models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z,L] = MLEM(M,m0,p,Runs,options,eps)
% Dealing with muticollinearity
if size(M,2)>1
    [z,Vif,sValue,condInd,VarDecomp]=Collinearity(M);
    col=size(M,2);      % number of columns in M
    Q=cell(1,col);      % cell array keeps truck of independent and combined variables
    for k=1:col
        Q{k}=k;         % create cell array of indexes
    end
    
    while ~isempty(z) && col>1
        if length(z)==1
            Q(z)=[];
            M(:,z)=[];
        else
            if length(z)==2
                not_z=setdiff(1:col,z);    % not z
            elseif length(z)==3
                [~,not_z]=min(Vif);
                z=setdiff(z, not_z);       % adjust z
                not_z=setdiff(1:col,z);
            else
                [~,not_z]=min(Vif);
                z=setdiff(z,not_z);
            end
            comb_z=Q(z);
            comb=cat(1,comb_z{:});  % merge z in one cell
            Q(z)=[];
            Q=[Q comb];             % add to the end
            M=[M(:,not_z) mean(M(:,z),2)];
        end
        col=size(M,2);
        if col>1
            [z,Vif,sValue,condInd,VarDecomp]=Collinearity(M);
        end
    end
    [nrows,ncols] = cellfun(@size,Q);
    independent=cell2mat(Q(nrows==1));
    combined=cell2mat(Q(nrows>1))';
    not_dropped=union(independent,combined);
    dropped=setdiff(1:p, not_dropped);
    
    %%%%%%%%%%%%%% Estimation %%%%%%%%%%%%%%%%%%%%%%
    
    [~,K]=size(M); % K is the number of independent variables in regression; up to 4 corrrsponding to B1, B2, B3 and B4
    %CN(i,2)=K;
    init=Init(Runs,K+1);        % precomute initial conditions
    
    if K==1
        mloglik = @(w) min(-sum(log(normpdf((m0-w(1)*M(:,1))/w(2))/w(2))),10^10);
        A=[-1 0; 1 0; 0 -1];
        b=[0; 1; -eps];
    elseif K==2
        mloglik = @(w) min(-sum(log(normpdf((m0-w(1)*M(:,1)-w(2)*M(:,2))/w(3))/w(3))),10^10);
        A=[-1 0 0 ; 0 -1 0; 1 1 0; 0 0 -1];
        b=[0; 0; 1; -eps];
    elseif K==3
        mloglik = @(w) min(-sum(log(normpdf((m0-w(1)*M(:,1)-w(2)*M(:,2)-w(3)*M(:,3))/w(4))/w(4))),10^10);
        A=[-1 0 0 0; 0 -1 0 0 ; 0 0 -1 0; 1 1 1 0; 0 0 0 -1];
        b=[0; 0; 0; 1; -eps];
    else     % K=4
        mloglik = @(w) min(-sum(log(normpdf((m0-w(1)*M(:,1)-w(2)*M(:,2)-w(3)*M(:,3)-w(4)*M(:,4))/w(5))/w(5))),10^10);
        A=[-1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0; 1 1 1 1 0; 0 0 0 0 -1];
        b=[0; 0; 0; 0; 1; -eps];
    end
    Sol=nan(Runs,K+1);
    L_sol=nan(Runs,1);
    fsol=10^10;
    for run=1:Runs
        [sol, l_sol]=fmincon(mloglik,init(run,:),A,b,[],[],[],[],[],options);
        Sol(run,:)=sol;
        L_sol(run)=l_sol;
    end
    [L_best,ind_best]=min(L_sol);
    Sol_best=Sol(ind_best,:);
    L=L_best;                    % likelihood for individuals i
    Z=NaN(1,p+1);
    S_ind=length(independent);      % always first
    S_comb=length(combined);
    S_drop=length(dropped);
    E1=0;
    E2=0;
    if S_ind
        Z(1,independent)=Sol_best(1:S_ind);
        E1=sum(Z(independent));
    end
    if S_comb
        Z(1,combined)=Sol_best(S_ind+1:end-1)/S_comb; % unifor prior
        E2=sum(Z(combined));
    end
    if S_drop
        Z(1,dropped)=(1-E1-E2)/(2*S_drop);            % uniform prior; B0=B_droped
    end
    Z(1,end)=Sol_best(end);
else
    if size(M,2)==1
        init=Init(Runs,2);        % precomute initial conditions
        mloglik = @(w) min(-sum(log(normpdf((m0-w(1)*M(:,1))/w(2))/w(2))),10^10);
        A=[-1 0; 1 0; 0 -1];
        b=[0; 1; -eps];
        Sol=nan(Runs,2);
        L_sol=nan(Runs,1);
        fsol=10^10;
        for run=1:Runs
            [sol, l_sol]=fmincon(mloglik,init(run,:),A,b,[],[],[],[],[],options);
            Sol(run,:)=sol;
            L_sol(run)=l_sol;
        end
        [L_best,ind_best]=min(L_sol);
        Sol_best=Sol(ind_best,:);
        L=L_best;
        Z=Sol_best(:);
    else
        init=Init(Runs,1);        % precomute initial conditions
        mloglik = @(w) min(-sum(log(normpdf(m0/w(1))/w(1))),10^10);
        A=[-1];
        b=[-eps];
        Sol=nan(Runs,1);
        L_sol=nan(Runs,1);
        fsol=10^10;
        for run=1:Runs
            [sol, l_sol]=fmincon(mloglik,init(run,:),A,b,[],[],[],[],[],options);
            Sol(run,:)=sol;
            L_sol(run)=l_sol;
        end
        [L_best,ind_best]=min(L_sol);
        Sol_best=Sol(ind_best,:);
        L=L_best;
        Z=Sol_best(:);
    end
end
end

