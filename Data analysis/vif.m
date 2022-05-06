function V=vif(X)

[n,p]=size(X);
V=zeros(1,p);
if p>=2
    for i=1:p
        pred=setdiff(1:p,i);
        md1=fitlm(X(:,pred),X(:,i),'Intercept',false);
        R2=md1.Rsquared.Ordinary;
        V(i)=1/(1-R2);
    end
else
    if p==1
        md1=fitlm(zeros(size(X(:,1))),X(:,1),'Intercept',false);
        R2=md1.Rsquared.Ordinary;
        V(1)=1/(1-R2);
    end
end