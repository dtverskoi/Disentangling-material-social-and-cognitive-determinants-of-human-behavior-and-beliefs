%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% doit_clean.m estimates parameters of the utility function and beliefs dynamics with or without
% authority messaging. The functions Init.m and Collinearity.m are used to produce the above 
% estimates. Note that the function Collinearity.m uses the results of the function vif.m.

% The doit_clean.m code also includes: (1) estimates of the alternative models. To perform this
% analysis, the function MLEM.m is required; (2) bootstrap confidence intervals of the parameter 
% estimates; (3) analysis of gender differences; (4) effects of the differences in the model 
% parameters between individuals; (5) analysis of utility losses; (6) distributions of parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
rng('shuffle','twister')

%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%
N=150;          % number of individuals
T=35;           % number of time steps
Runs=10;         % number of different initial conditions in the optimization problem
eps=0.01;       % minimum st dev
e=0;
Equation=0:3;
Authority=0:0;

Debug=0;        % 0 or 1

ESTIMATION=1;
ALTERNATIVES=0;
COEF_GRAPHS=0;
LH_GRAPHS=0;
GENDER=0;
UTILITY=0;          % plot graphs of utility components
HISTOGRAMS=0;

save_graphs=0;

% Program
if ESTIMATION 
    Main=readmatrix('Main.csv');  %open data file
    % find parameters that we are interesting in: column 1 for an individual id;
    % 32 for x; 75 for xp; 70-74 for xe; and 82-86 for xn;
    % 33-37 for total contribution; 94 for session label;
    % 17 for active/inactive
    % 26,27,28, and 30 for missing data
    % 53 for material payoff
    Ind=[1 32 75 70 71 72 73 74 82 83 84 85 86 33 34 35 36 37 17 94 26 27 28 30 53];
    %  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
    
    D(1,:,:)=Main(1:N*T,Ind);
    D(2,:,:)=Main(N*T+1:end,Ind);
    save D.mat D
    
    Final=nan(2,4,7);           % all mean estimates
    C=cell(length(Equation),length(Authority));    
    options = optimoptions('fmincon','Algorithm','interior-point','Display','notify');
    resid=NaN(N,T,3,2);
    
    for equation=Equation
        if equation==0
            AFinal=NaN(10,7,2);   % Results
        end
        if equation~=0
            BFinal=NaN(5,5,2);   % Results
        end
        
        for authority=Authority
            
            % Useful
            n=nan(N,1);         % number of observations for each individual
            L=nan(N,1);         % likelihoods
            CN=nan(N,4);        % number of problematic variables in the design matirx and the actual number of regressors
            
            %Data
            if authority
                p=4;                    % number of parameters to estimate
            else
                p=3;                    % number of parameters to estimate
            end
            if equation
                p=p-1;
            end
            
            Par=nan(N,p+1);
            tic
            for i=1:N
                i
                I=find(D(authority+1,:,1)==i);                                                      % rows for individual i
                Di=squeeze(D(authority+1,I,:));                                                              % data matrix for individual i
                
                J1=find(Di(:,24)~=1);                                                   % not missed
                theta=min(30,(14-5*nanmean(Di(:,4:8),2)/12)/(1/6));                        % forward -looking best response (given expectations)
                
                J=find(Di(:,24)==1);                                                    % rounds with missing contibution x; all cases with missing data don't have x
                Di(J,:)=[];
                theta(J)=[];
                
                n(i)=size(Di,1);          % number of data points
                if n(i)>=30               % choose individuals with at least 30 observations
                    %%%%%%%%%%%%%%%%%% Data preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % extract data
                    x=Di(:,2);
                    y=Di(:,3);
                    tilde_y=mean(Di(:,9:13),2);
                    tilde_x=mean(Di(:,4:8),2);
                    X=mean(Di(1:end-1,14:18),2);
                    G=14*ones(n(i),1);
                    
                    if equation==1
                        m0=diff(y);
                        M=[x(1:end-1)       X  G(1:end-1)]-y(1:end-1);
                    elseif equation==2
                        m0=diff(tilde_y);
                        M=[y(1:end-1)       X  G(1:end-1)]-tilde_y(1:end-1);
                    elseif equation==3
                        m0=diff(tilde_x);
                        M=[tilde_y(1:end-1) X  G(1:end-1)]-tilde_x(1:end-1);
                    else
                        m0=x-theta;
                        M=[y tilde_y tilde_x G]-theta;
                    end
                    
                    % Perturbations
                    if ~authority
                        M(:,end)=[];
                    end
                    Md=M;
                    
                    % Dealing with muticollinearity
                    [z,Vif,sValue,condInd,VarDecomp]=Collinearity(M);
                    col=size(M,2);      % number of columns in M
                    Q=cell(1,col);      % cell array keeps track of independent and combined variables
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
                            Q=[Q comb];             % add at the end
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
                    
                    if Debug
                        i
                        Q
                        independent
                        combined
                        dropped
                    end
                    
                    %%%%%%%%%%%%%% Estimation %%%%%%%%%%%%%%%%%%%%%%
                    
                    [~,K]=size(M); % K is the number of independent variables in regression; up to 4 corrrsponding to B1, B2, B3 and B4
                    CN(i,2)=K;
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
                    fsol=10000;
                    for run=1:Runs
                        [sol, l_sol]=fmincon(mloglik,init(run,:),A,b,[],[],[],[],[],options);
                        Sol(run,:)=sol;
                        L_sol(run)=l_sol;
                    end
                    [L_best,ind_best]=min(L_sol);
                    Sol_best=Sol(ind_best,:);
                    L(i)=L_best;                    % likelihood for individuals i
                    
                    % Assign values to right parameters
                    S_ind=length(independent);
                    S_comb=length(combined);
                    S_drop=length(dropped);
                    E1=0;
                    E2=0;
                    if S_ind
                        Par(i,independent)=Sol_best(1:S_ind);
                        E1=sum(Par(i,independent));
                    end
                    if S_comb
                        Par(i,combined)=Sol_best(S_ind+1:end-1)/S_comb;  % unifor prior
                        E2=sum(Par(i,combined));
                    end
                    if S_drop
                        Par(i,dropped)=(1-E1-E2)/(2*S_drop);            % uniform prior; B0=B_droped
                    end
                    Par(i,end)=Sol_best(end);
                    if equation
                        Xpr=Par(i,1)*Md(:,1) + Par(i,2)*Md(:,2);
                        if authority
                            Xpr=Xpr+Par(i,3)*Md(:,3);
                        end
                        resid(i,1:n(i)-1,equation,authority+1)=m0-Xpr;
                    end
                end
            end
            toc
            writematrix(Par,['B5' num2str(equation+1) '.csv'])
            C{equation+1,authority+1}=Par;          % the whole set of parameters
            
            Used=find(n>=30);
            Estimate=nanmean(Par);
            BIC=2*nansum(L(Used))+(p+1)*log(sum(n(Used)))
            AIC=2*nansum(L(Used))+(p+1)*2
            
            if ~equation
                Final(authority+1, equation+1,1)=1-sum(Estimate(1:p),2);
            end
            
            Final(authority+1, equation+1, 2:size(Par,2))= Estimate(1:p);
            Final(authority+1, equation+1, end-1:end)= [Estimate(end) BIC];
            
            if ALTERNATIVES
                [AEstimate,ABIC,BEstimate,BBIC] = Alter(equation,authority,N,p,D,Used,Runs,options,eps);
                if ~equation
                    AFinal(:, 1, authority+1)=1-nansum(AEstimate(:,1:p),2);
                    AFinal(:, 2:size(Par,2), authority+1)= AEstimate(:,1:p);
                    AFinal(:, end-1:end, authority+1)= [AEstimate(:,end) ABIC'/10000];
                    AFinal(isnan(AFinal))=0;
                end
                    
                if equation
                    BFinal(:, 1:size(Par,2)-1, authority+1)= BEstimate(:,1:p);
                    BFinal(:, end-1:end, authority+1)= [BEstimate(:,end) BBIC'/10000];
                    BFinal(isnan(BFinal))=0;
                end
            end
            
        end
        save C.mat C          % save parameters
    end
    Omega=NaN(N,3,3,2);
    for authority=Authority
        for i=1:N
            Omega(i,:,:,authority+1)=cov(squeeze(resid(i,:,:,authority+1)),'omitrows');
        end
    end
end



if COEF_GRAPHS
    load('C.mat');              % load parameters
    f0=figure('Position',  [2000, 300,800, 400]);
    f1=figure('Position',  [2000, 300,800, 700]);
    
    for equation=Equation
        for authority=Authority
            Par=cell2mat(C(equation+1,authority+1));
            %%%%%%%%%%
            %Bootstrap
            %%%%%%%%%%
            Boot=1001;      % number of bootstrap samples
            if ~equation
                B=[1-sum(Par(:,1:end-1),2) Par(:,1:end-1)];
            else
                B=Par(:,1:end-1);
            end
            Coeff=NaN(Boot,size(B,2));
            for r=1:Boot
                % individual coefficients
                if r==1
                    IR=1:N;
                else
                    IR=randi(N,1,N);
                end
                % final estimates
                Coeff(r,:)=nanmean(B(IR,:),1);
            end
            Coef_sorted=sort(Coeff(2:Boot,:));
            CI=[Coeff(1,:); Coef_sorted((Boot-1)*0.025,:); Coef_sorted((Boot-1)*0.975,:)];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot the distribution of estimates and CI
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~equation
                set(0,'CurrentFigure',f0)
                subplot(1,2,authority+1)
            else
                set(0,'CurrentFigure',f1)
                subplot(3,2,2*equation-1+authority)
            end
            pars=size(CI,2);
            b=bar(CI(1,:));
            b.FaceColor = 'flat';
            set(gca,'FontSize',15);
            hold on
            er = errorbar(1:pars,CI(1,:),CI(1,:)-CI(2,:),CI(3,:)-CI(1,:),'CapSize',18);
            er.Color = [0.9 0 0];
            er.LineStyle = 'none';
            er.LineWidth=2;
            if ~equation
                b.CData(1,:) = [.52 0.73 .4];       %change color of the first bar
            end
            if ~authority
                if equation==1
                    xticklabels({'\alpha_1','\beta_1'})
                    ylabel('effects on $y$','Interpreter','Latex');
                    title('no messaging')
                elseif equation==2
                    xticklabels({'\alpha_2','\beta_2'})
                    ylabel('effects on $\tilde{y}$','Interpreter','Latex');
                elseif equation==3
                    xticklabels({'\alpha_3','\beta_3'})
                    ylabel('effects on $\tilde{x}$','Interpreter','Latex');
                else
                    xticklabels({'B_0','B_1','B_2','B_3','Interpreter','Latex'})
                    ylabel('effects on $x$','Interpreter','Latex');
                    title('no messaging')
                end
            else
                b.CData(pars,:) = [1 0 0];
                er5 = errorbar(pars,CI(1,pars),CI(1,pars)-CI(2,pars),CI(3,pars)-CI(1,pars),'CapSize',18);   %change color of the last
                if ~equation
                    b.CData(1,:) = [.52 0.73 .4];       %change color of the first bar
                end
                er5.Color = [0 0 0];
                er5.LineStyle = 'none';
                er5.LineWidth=2;
                if equation==1
                    xticklabels({'\alpha_1','\beta_1','\gamma_1'})
                    title('with messaging')
                elseif equation==2
                    xticklabels({'\alpha_2','\beta_2','\gamma_2'})
                elseif equation==3
                    xticklabels({'\alpha_3','\beta_3','\gamma_3'})
                else
                    xticklabels({'B_0','B_1','B_2','B_3','B_4'})
                    title('with messaging')
                end
            end
            xlim([0.5,pars+0.5])
            if ~equation
                ylim([0 0.5]);
            else
                ylim([0 0.8]);
            end
            hold off
        end
        
        if save_graphs
            set(0,'CurrentFigure',f0)
            print('utility_pars','-depsc2')
            set(0,'CurrentFigure',f1)
            print('beliefs_pars','-depsc2')
        end
    end
end


if GENDER
    MF=readmatrix('questionnaire.xlsx');
    gender=MF(:,20);
    gender=reshape(gender,N,2);
end

if LH_GRAPHS
    load('C.mat');              % load parameters
    load('D.mat');              % load trajectories
    F=14;
    for authority=Authority
        B=cell2mat(C(1,authority+1));
        B=[1-sum(B(:,1:end-1),2) B(:,1:end-1)];
        
        for par=1:size(B,2)+GENDER                 % B_0,...,
            if par<=size(B,2)
                low=find(B(:,par)<0.1);
                high=find(B(:,par)>=0.1);
                color1=[1 .7 .0];           % yellow
                color2=[0 .7 .0];           % green
                color_e=[0.9 0 0];          % red
            else
                low=find(gender(:,authority+1)==0);
                high=find(gender(:,authority+1)==1);
                color1=[0 1 1];             % blue
                color2=[1 0.5 0.8];         % pink
                color_e=[0 0 0];            % black
            end
            DD=squeeze(D(authority+1,:,:));
            Missing=find(DD(:,24)==1);
            DD(Missing,:)=nan;           % convert to nan
            U=[DD(:,2) DD(:,3) mean(DD(:,9:13),2) mean(DD(:,4:8),2) DD(:,25)]; % trajectories
            
            %%%%%%%%%%%%%%%%%%%%%%%%% Trajectories  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f1(authority+1,par)=figure('Position',  [2000, 300,800, 300]);
            for i=1:5
                subplot(1,5,i)
                Uc=reshape(U(:,i),N,T);
                Ulow=nanmean(Uc(low,:));
                Uhigh=nanmean(Uc(high,:));
                plot(1:T,Ulow,'color',color1,'LineWidth',3),
                hold on
                plot(1:T,Uhigh,'color',color2,'LineWidth',3);
                hold off
                xticks([15 30])
                if i==1
                    ylabel('x','interpreter','latex','FontSize',F);
                elseif i==2
                    ylabel('y','interpreter','latex','FontSize',F);
                elseif i==3
                    ylabel('$\tilde{y}$','interpreter','latex','FontSize',F);
                elseif i==4
                    ylabel('$\tilde{x}$','interpreter','latex','FontSize',F);
                else
                    ylabel('$\pi$','interpreter','latex','FontSize',F);
                end
                xlabel('round');
                xlim([0 35]);
                yticks(90:10:round(max(Uhigh)));
                if i<5
                    ylim([14 26]);
                    yticks([14 18 22 26]);
                end
            end
            
            if save_graphs
                print(f1(authority+1,par),['LH_traj_',num2str(authority),num2str(par-1)],'-depsc2');
            end
            %%%%%%%%%%%%%%%%%%%%%%%% Effects on x   %%%%%%%%%%%%%%%%%%%%%%%
            Boot=1000;
            S_low =nan(Boot,size(B,2));
            S_high=nan(Boot,size(B,2));
            for k=1:Boot
                S_low(k,:) =  nanmean(B(low(randi(length(low),length(low),1)),:));
                S_high(k,:)=  nanmean(B(high(randi(length(high),length(high),1)),:));
            end
            S_low=sort(S_low);
            S_high=sort(S_high);
            CI_low=abs(S_low([0.05*Boot 0.95*Boot],:)-mean(S_low));
            CI_high=abs(S_high([0.05*Boot 0.95*Boot],:)-mean(S_high));
            
            B_set=1:3:3*size(B,2);
            
            f2(authority+1,par)=figure;
            h1=bar(B_set,nanmean(B(low,:)),.25);
            hold on
            h2=bar(B_set+1,nanmean(B(high,:)),.25);
            er1=errorbar(B_set, nanmean(B(low,:)),CI_low(1,:),CI_low(2,:),'CapSize',18);
            er2=errorbar(B_set+1, nanmean(B(high,:)),CI_high(1,:),CI_high(2,:),'CapSize',18);
            er1.Color = color_e;
            er2.Color = color_e;
            h1.FaceColor = color1;
            h2.FaceColor = color2;
            hold off
            set(gca,'XTick',B_set+0.5)
            set(gca,'XTickLabel',{'B_0' 'B_1' 'B_2' 'B_3' 'B_4'})
            er1.LineStyle = 'none';
            er1.LineWidth=2;
            er2.LineStyle = 'none';
            er2.LineWidth=2;
            ylabel('effects on x','FontSize',14)
            set(gca,'YTick',0.0:0.2:0.6)
            xlim([0 B_set(end)+2])
            %                 ylim([0 0.6])
            if par<=size(B,2)
                title(['Differences in B',num2str(par-1)]);
                legend(strcat('low, n=', num2str(length(low))),strcat('high, n=', num2str(length(high))),'FontSize',20)
            else
                legend(strcat('males, n=', num2str(length(low))),strcat('females, n=', num2str(length(high))),'FontSize',20)
            end
            if save_graphs
                print(f2(authority+1,par),['LH_util_',num2str(authority),num2str(par-1)],'-depsc2');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%% Effects on y, tilde+y and tilde_x %%%%%%%%%%%
            
            B1=cell2mat(C(2,authority+1));
            B2=cell2mat(C(3,authority+1));
            B3=cell2mat(C(4,authority+1));
            
            BB=[B1(:,1:end-1) B2(:,1:end-1) B3(:,1:end-1)];
            
            S_low =nan(Boot,size(BB,2));
            S_high=nan(Boot,size(BB,2));
            for k=1:Boot
                S_low(k,:) =  nanmean(BB(low(randi(length(low),length(low),1)),:));
                S_high(k,:)=  nanmean(BB(high(randi(length(high),length(high),1)),:));
            end
            S_low=sort(S_low);
            S_high=sort(S_high);
            CI_low=abs(S_low([0.05*Boot 0.95*Boot],:)-mean(S_low));
            CI_high=abs(S_high([0.05*Boot 0.95*Boot],:)-mean(S_high));
            
            S_low=sort(S_low);
            S_high=sort(S_high);
            CI_low=abs(S_low([0.05*Boot 0.95*Boot],:)-mean(S_low));
            CI_high=abs(S_high([0.05*Boot 0.95*Boot],:)-mean(S_high));
            
            BB_set=1:3:3*size(BB,2);
            
            f3(authority+1,par)=figure;
            h1=bar(BB_set,nanmean(BB(low,:)),.25);
            hold on
            h2=bar(BB_set+1,nanmean(BB(high,:)),.25);
            er1=errorbar(BB_set, nanmean(BB(low,:)),CI_low(1,:),CI_low(2,:),'CapSize',18);
            er2=errorbar(BB_set+1, nanmean(BB(high,:)),CI_high(1,:),CI_high(2,:),'CapSize',18);
            hold off
            set(gca,'XTick',BB_set+0.5)
            if ~authority
                set(gca,'XTickLabel',{'\alpha_1' '\beta_1' '\alpha_2' '\beta_2' '\alpha_3' '\beta_3'})
            else
                set(gca,'XTickLabel',{'\alpha_1' '\beta_1' '\gamma_1' '\alpha_2' '\beta_2' '\gamma_2' '\alpha_3' '\beta_3' '\gamma_3'})
            end
            
            er1.LineStyle = 'none';
            er1.LineWidth=2;
            er2.LineStyle = 'none';
            er2.LineWidth=2;
            er1.Color = color_e;
            er2.Color = color_e;
            h1.FaceColor = color1;
            h2.FaceColor = color2;
            ylabel('effects on $y, \tilde{y}, \tilde{x}$','FontSize',14,'interpreter','latex')
            set(gca,'YTick',0.0:0.2:0.6)
            xlim([0 BB_set(end)+2])
            %             legend(strcat('low, n=', num2str(length(low))),strcat('high, n=', num2str(length(high))),'FontSize',20)
            if par<=size(B,2)
                title(['Differences in B',num2str(par-1)]);
            end
            if save_graphs
                print(f3(authority+1,par),['LH_belief_',num2str(authority),num2str(par-1)],'-depsc2');
            end
        end
    end
end

%%%%%%%%%%%%%%% Estimate A's and plot utility components %%%%%%%%%
if UTILITY
    load('C.mat');              % load parameters
    load('D.mat');              % load trajectories
    F=3;
    for authority=Authority
        
        DD=squeeze(D(authority+1,:,:));
        
        Par=cell2mat(C(1,authority+1));
        sumB=sum(Par(:,1:end-1),2);
        Bpar=[1-sumB Par(:,1:end-1)];
        Apar=nan(size(Bpar));
        D2=1/12;
        Apar(:,1)=(1-sumB)            ./(1-(1-D2)*sumB);
        Apar(:,2:end)=D2*Bpar(:,2:end)./(1-(1-D2)*sumB);
        
        %         if authority
        %             U1=nanmean(Apar);
        %         else
        %             U0=nanmean(Apar);
        %         end
        
        Missing=find(DD(:,24)==1);
        DD(Missing,:)=nan;           % convert to nan
        U=[DD(:,2) DD(:,3) mean(DD(:,9:13),2) mean(DD(:,4:8),2) mean(DD(:,14:18),2)]; % trajectories
        
        theta=min(30,(14-5*U(:,4)/12)/(1/6));
        Pi_lost=reshape(15*(theta+5*U(:,4))-0.5*(theta+5*U(:,4)).^2-theta,N,T);
        
        U0=Pi_lost.*Apar(:,1);
        U1=reshape((U(:,1)-U(:,2)).^2,N,T).*Apar(:,2);
        U2=reshape((U(:,1)-U(:,3)).^2,N,T).*Apar(:,3);
        U3=reshape((U(:,1)-U(:,4)).^2,N,T).*Apar(:,4);
        if authority
            U4=reshape((U(:,1)-14).^2,N,T).*Apar(:,5);
        end
        
        % payoofs
        Pi_real(authority+1,:)=nanmean(reshape(DD(:,25),N,T));
        Z=theta+5*U(:,4);
        v=theta./Z;
        W=30+v.*(15*Z-Z.^2/12)-theta;
        Pi_br(authority+1,:)=nanmean(reshape(W,N,T));
        
        figure
        plot( nanmean((Pi_br(authority+1,:)-Pi_real(authority+1,:)).*Apar(:,1)),'-o','color',[.52 0.73 .4],'LineWidth',F)
        hold on
        plot(nanmean(U1),'-o','color',[139 69 19]/255,'LineWidth',F)
        plot(nanmean(U2),'-o','LineWidth',F)
        plot(nanmean(U3),'-ob','LineWidth',F)
        if ~authority
            plot(nanmean(U1+U2+U3),'-k')
            legend('materail payoffs','cognitive dissonance','disapproval by peers','conformity w/ peers','total non-material','Location','NorthEast')
            title('no messaging')
        else
            plot(nanmean(U4),'-ob','LineWidth',F)
            plot(nanmean(U1+U2+U3+U4),'-k')
            legend('materail payoffs','cognitive dissonance','disapproval by peers','conformity w/ peers','conformity w/ authority','total non-material','Location','NorthEast')
            title('with messaging')
        end
        ylabel('loss of utility')
        xlabel('round')
        hold off
        set(gca,'XTick',0:10:30);
        set(gca,'YTick',0:10:50);
        ylim([0 40])
        if save_graphs
            print(['payoff_loss_',num2str(authority)],'-depsc2');
        end
    end
    figure
    plot(Pi_real(1,:),'-ob','LineWidth',F)
    hold on
    plot(Pi_br(1,:),'--b','LineWidth',F)
    plot(Pi_real(2,:),'-or','LineWidth',F)
    plot(Pi_br(2,:),'--r','LineWidth',F)
    xlabel('round')
    ylabel('payoffs')
    legend('$\pi_0$','$\pi_{BR,0}$','$\pi_1$','$\pi_{BR,1}$','Interpreter','Latex')
    hold off
    
    if save_graphs
        print('payoffs','-depsc2');
    end
end

if HISTOGRAMS
    HCa={'B_0','B_1','B_2','B_3','B_4'};
    HCb={'\alpha_1','\beta_1','\gamma_1'; '\alpha_2','\beta_2','\gamma_2'; '\alpha_3','\beta_3','\gamma_3'};
    
    equation=0;
    for authority=Authority
        figure
        set(gcf, 'Position',  [200, 300, 1000, 200]);
        H=C{1,authority+1};
        H=[1-nansum(H(:,1:end-1),2) H];
        for i=1:size(H,2)-1
            subplot(1,size(H,2)-1,i)
            histogram(H(:,i),15)
            set(gca,'FontSize',15)
            xlabel(HCa{i})
            xlim([0 1])
            ylim([0 70])
        end
        print(['hist',num2str(authority+1)],'-depsc2')
    end
    
    for authority=Authority
        figure
        if authority
            set(gcf, 'Position',  [200, 300, 600, 600]);
        else
            set(gcf, 'Position',  [200, 300, 400, 600]);
        end
        for equation=1:3
            H=C{equation+1,authority+1};
            for i=1:size(H,2)-1
                subplot(3,size(H,2)-1,(size(H,2)-1)*(equation-1)+i)
                histogram(H(:,i),15)
                set(gca,'FontSize',15)
                xlabel(HCb{equation,i})
                xlim([0 1])
                ylim([0 80])
            end
        end
        print(['histb',num2str(authority+1)],'-depsc2')
    end
end
