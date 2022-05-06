%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions Init.m produces initial conditions for the maximum likelihood 
% procedure used in the doit_clean.m and MLEM.m files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z=Init(Runs,d)
% d>=3 is the number of dimensions for the vector of initial conditions
% d-1 components come from a broken stick distribution, d-th component fom
% a uniform distribution
if d==1
    z=[5*rand(Runs,1)+1];
else
    if d==2
        z=[rand(Runs,1) 5*rand(Runs,1)+1];
    else
        f0=sort(rand(Runs,d-1),2);
        f1=[f0(:,1) diff(f0,1,2) 1-f0(:,end)];
        jj=randperm(d);
        f=f1(:,jj);
        z=[f(:,1:d-1) 5*rand(Runs,1)+1];
        %z=[f(:,1) diff(f,1,2) 1-f(:,end) 5*rand(Runs,1)+1];
    end
end

end
