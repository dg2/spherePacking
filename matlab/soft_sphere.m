%function [R,alpha]=soft_sphere(K,C);
% Soft version of the one-class SVM / support estimation using
% hyperspheres.

function [R,alpha]=soft_sphere(K,C);

N=length(K);
if (C<1/N)
    fprintf('C is too small, changing to 1.1/N\n');
    C=1.1/N;
end
% Configure the optimization
f=-diag(K)';
options=optimset('Display','off');
%A=[];
%b=[];
A=-1*eye(N);
b=zeros(N,1);
THRESH=C/1000;

Aeq=ones(1,N);
beq=1;
%LB=0;
%UB=inf;
LB=[];
UB=ones(N,1)*C;
% Solve the optimization problem
warning off
if (is_octave)
    [alpha,W]=qp([], 2*K,f',Aeq,beq,LB,UB,[],A,b);
else
    [alpha,W,flag]=quadprog(2*K,f,A,b,Aeq,beq,LB,UB,[],options);
end
        
warning on

% Radius of the spere
% The radius can be obtained using any point i such that 0<alpha(i)<C
i=find(alpha>0+THRESH & alpha<C-THRESH,1);
if (isempty(i))
    fprintf('All support vectors are outside the ball -> Can not directly determine radius\n');
    % Obtain the distances from the different non-SV points to the center
    idx_try=setdiff(1:N,find(alpha>(C-THRESH)));
%    Kxx=K(idx_try,:);
    D=sqrt(diag(K)-2*K*alpha+alpha'*K*alpha);
    D=D(idx_try);
    R=max(D);
else
    R=sqrt(K(i,i)-2*K(i,:)*alpha+alpha'*K*alpha);
end
