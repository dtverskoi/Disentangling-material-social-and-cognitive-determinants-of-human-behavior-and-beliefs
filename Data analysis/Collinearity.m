%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function Collinearity.m performs multicollinearity detection. It uses the results of the
% vif.m function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [z,Vif,sValue,condInd,VarDecomp]=Collinearity(M)
% identifies problematic variables in design matrix M
% uses threshold of 30 for condition number, 0.5 for varDecomp, and 10 for vif

[sValue,condInd,VarDecomp]=collintest(M,'display','off');
ind=condInd>30;                                    
K_belsley=find(sum(double(VarDecomp(ind,:)>0.5),1));    % indices of problematic variables according to collintest
Vif=vif(M);
% K_vif=[find(Vif>10) find(isnan(Vif))];                  % indices of problematic variables according to vif
K_vif=find(Vif>10);
if sum(double(isnan(Vif))) || sum(double(isinf(Vif)))
    z=sort(K_belsley);
else
    z=sort(union(K_belsley, K_vif));  
end
