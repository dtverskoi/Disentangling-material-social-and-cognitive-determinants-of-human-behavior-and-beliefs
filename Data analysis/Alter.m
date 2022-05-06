%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The function Alter.m estimates all the alternative models. To perform this analysis, the function MLEM.m
% is required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AEstimate,ABIC,BEstimate,BBIC] = Alter(equation,authority,N,p,D,Used,Runs,options,eps)
if authority==0
    DD=squeeze(D(1,:,:));
else
    DD=squeeze(D(2,:,:));
end
D=DD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Alternative models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if equation==0
    BEstimate=[];
    BBIC=[];
    %Alternative models
    A=10;
    APar=NaN(N,A,p+1);
    La=NaN(N,A);
    tic
    for i=1:N
        I=find(D(:,1)==i);                                                      % rows for individual i
        Di=D(I,:);                                                              % data matrix for individual i
        
        J1=find(Di(:,24)~=1);                                                   % not missed
        th=min(30,(14-5*nanmean(Di(:,4:8),2)/12)/(1/6));                        % forward -looking best response (given expectations)
        
        J=find(Di(:,24)==1);                                                    % rounds with missing contibution x; all cases with missing data don't have x
        Di(J,:)=[];
        th(J)=[];
        
        n(i)=size(Di,1);          % number of data points
        if n(i)>=30               % choose individuals with at least 30 observations
            x=Di(:,2);
            y=Di(:,3);
            tilde_y=mean(Di(:,9:13),2);
            tilde_x=mean(Di(:,4:8),2);
            X=mean(Di(1:end-1,14:18),2);
            G=14*ones(n(i),1);
            m0=x-th;
            m0s=x(2:end);
            ths=min(30,(14-5*X/12)/(1/6)); 
            M=[y tilde_y tilde_x G]-th;
            if ~authority
                M(:,end)=[];
            end
            for a=1:A              % BR
                if a==1
                    m=m0;
                    M1=[];
                    p1=0;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,end)=Z;
                    La(i,a)=L1;
                end
                if a==2            % SL
                    m=m0-M(:,3);
                    M1=[];
                    p1=0;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,end)=Z;
                    APar(i,a,3)=1;
                    La(i,a)=L1;
                end
                if a==3            % BR+SL
                    m=m0;
                    M1=M(:,3);
                    p1=1;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,3)=Z(1);
                    APar(i,a,end)=Z(2);
                    La(i,a)=L1;
                end
                if a==4            % no cognitive dissonance
                    m=m0;
                    M1=M(:,2:end);
                    p1=p-1;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,2:end)=Z(1,:);
                    La(i,a)=L1;
                end
                if a==5            % no peer influence
                    m=m0;
                    M1=[M(:,1) M(:,3:end)];
                    p1=p-1;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,1)=Z(1);
                    APar(i,a,3:end)=Z(2:end);
                    La(i,a)=L1;
                end
                if a==6            % no conformity with peers
                    m=m0;
                    p1=p-1;
                    if authority
                        M1=[M(:,1:2) M(:,4:end)];
                    else
                        M1=M(:,1:2);
                    end
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,1:2)=Z(1:2);
                    APar(i,a,4:end)=Z(3:end);
                    La(i,a)=L1;
                end
                if a==7           % no conformity with the authority
                    if authority
                        m=m0;
                        M1=M(:,1:3);
                        p1=p-1;
                        [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                        APar(i,a,1:3)=Z(1:3);
                        APar(i,a,end)=Z(end);
                        La(i,a)=L1;
                    end
                end
                if a==9
                    m=m0s-ths;
                    M1=[];
                    p1=0;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,end)=Z;
                    La(i,a)=L1;
                end
                if a==10
                    m=m0s-X;
                    M1=[];
                    p1=0;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,a,end)=Z;
                    La(i,a)=L1;
                end
            end
        end
    end
    toc
    
    AEstimate=squeeze(nanmean(APar,1));
    pa=[1 1 2 p p p p p+1 1 1];
    ABIC=2*nansum(La(Used,:),1)+pa*log(sum(n(Used)));
    
    a=8;
    if a==8                       % homogeneity
        m0=[];
        M=[];
        for i=1:N
            I=find(D(:,1)==i);                                                      % rows for individual i
            Di=D(I,:);                                                              % data matrix for individual i
            J1=find(Di(:,24)~=1);                                                   % not missed
            th=min(30,(14-5*nanmean(Di(:,4:8),2)/12)/(1/6));                        % forward -looking best response (given expectations)
            J=find(Di(:,24)==1);                                                    % rounds with missing contibution x; all cases with missing data don't have x
            Di(J,:)=[];
            th(J)=[];
            
            n(i)=size(Di,1);          % number of data points
            if n(i)>=30               % choose individuals with at least 30 observations
                % extract data
                x=Di(:,2);
                y=Di(:,3);
                tilde_y=mean(Di(:,9:13),2);
                tilde_x=mean(Di(:,4:8),2);
                X=mean(Di(1:end-1,14:18),2);
                G=14*ones(n(i),1);
                m0=[m0; x-th];
                M8=[y tilde_y tilde_x G]-th;
                if ~authority
                    M8(:,end)=[];
                end
                M=[M; M8];
            end
        end
        
        p1=p;
        [Z,L1] = MLEM(M,m0,p1,Runs,options,eps);
        AEstimate(a,:)=Z(1,:);
        ABIC(a)=2*L1+pa(a)*log(sum(n(Used)));
    end  
end


if equation~=0
    AEstimate=[];
    ABIC=[];
    %Alternative models
    A=5;
    APar=NaN(N,A,p+1);
    La=NaN(N,A);
    tic
    
    for i=1:N
        i
        
        I=find(D(:,1)==i);                                                      % rows for individual i
        Di=D(I,:);                                                              % data matrix for individual i
        
        J1=find(Di(:,24)~=1);                                                   % not missed
        theta=min(30,(14-5*nanmean(Di(:,4:8),2)/12)/(1/6));                        % forward -looking best response (given expectations)
        
        J=find(Di(:,24)==1);                                                    % rounds with missing contibution x; all cases with missing data don't have x
        Di(J,:)=[];
        theta(J)=[];
        
        n(i)=size(Di,1);          % number of data points
        if n(i)>=30               % choose individuals with at least 30 observations
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
            end
            
            if ~authority
                M(:,end)=[];
            end
            
            for ab=1:A              % random noise
                if ab==1
                    m=m0;
                    M1=[];
                    p1=0;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,ab,end)=Z;
                    La(i,ab)=L1;
                end
                if ab==2            % no cognitive dissonance
                    m=m0;
                    M1=M(:,2:end);
                    p1=p-1;
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,ab,2:end)=Z(:);
                    La(i,ab)=L1;
                end
                if ab==3            % no conformity with peers
                    m=m0;
                    p1=p-1;
                    if authority
                        M1=[M(:,1) M(:,3:end)];
                    else
                        M1=M(:,1);
                    end
                    [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                    APar(i,ab,1)=Z(1);
                    APar(i,ab,3:end)=Z(2:end);
                    La(i,ab)=L1;
                end
                if ab==4           % no conformity with the authority
                    if authority
                        m=m0;
                        M1=M(:,1:2);
                        p1=p-1;
                        [Z,L1] = MLEM(M1,m,p1,Runs,options,eps);
                        APar(i,ab,1:2)=Z(1:2);
                        APar(i,ab,end)=Z(end);
                        La(i,ab)=L1;
                    end
                end
            end
        end
    end
    
    BEstimate=squeeze(nanmean(APar,1));
    pb=[1 p p p p+1];
    BBIC=2*nansum(La(Used,:),1)+pb*log(sum(n(Used)));
    
    ab=5;
    if ab==5                      % homogeneity
        m0=[];
        M=[];
        for i=1:N
            I=find(D(:,1)==i);                                                      % rows for individual i
            Di=D(I,:);                                                              % data matrix for individual i
            
            J1=find(Di(:,24)~=1);                                                   % not missed
            theta=min(30,(14-5*nanmean(Di(:,4:8),2)/12)/(1/6));                        % forward -looking best response (given expectations)
            
            J=find(Di(:,24)==1);                                                    % rounds with missing contibution x; all cases with missing data don't have x
            Di(J,:)=[];
            theta(J)=[];
            
            n(i)=size(Di,1);          % number of data points
            if n(i)>=30               % choose individuals with at least 30 observations
                x=Di(:,2);
                y=Di(:,3);
                tilde_y=mean(Di(:,9:13),2);
                tilde_x=mean(Di(:,4:8),2);
                X=mean(Di(1:end-1,14:18),2);
                G=14*ones(n(i),1);
                
                if equation==1
                    m0=[m0; diff(y)];
                    M5=[x(1:end-1)       X  G(1:end-1)]-y(1:end-1);
                elseif equation==2
                    m0=[m0; diff(tilde_y)];
                    M5=[y(1:end-1)       X  G(1:end-1)]-tilde_y(1:end-1);
                elseif equation==3
                    m0=[m0; diff(tilde_x)];
                    M5=[tilde_y(1:end-1) X  G(1:end-1)]-tilde_x(1:end-1);
                end
                
                if ~authority
                    M5(:,end)=[];
                end
                M=[M; M5];
            end
        end
        p1=p;
        [Z,L1] = MLEM(M,m0,p1,Runs,options,eps);
        BEstimate(ab,:)=Z(1,:);
        BBIC(ab)=2*L1+pb(ab)*log(sum(n(Used)));
    end
end
end

