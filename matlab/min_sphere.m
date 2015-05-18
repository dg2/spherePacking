% function [R,alpha]=min_sphere(K);
% Find the radius R and center (written as a support vector expansion with 
% coefficients alpha) of the smallest hypersphere containing the points 
% defined by the kernel matrix K.

% Dario Garcia, 2010

function [R,alpha]=min_sphere(K);

N=length(K);
TOL=1e-6;
% Configure the optimization
f=-diag(K)';
options=optimset('Display','off');
A=-1*eye(N);
b=zeros(N,1);

Aeq=ones(1,N);
beq=1;
LB=[];
UB=[];
% Solve the optimization problem
warning off
% 
% We want to minimise x^T*K*x + f^T*x s.t. Aeq * x = beq, A * x <= b 
% Minimisers traditionally work with 0.5 * x^t * K * ..., so that's
% why we multiply 2*K.
if (is_octave)
    [alpha,W]=qp([], 2*K,f',Aeq,beq,LB,UB,[],A,b);
else
    [alpha,W]=quadprog(2*K,f,A,b,Aeq,beq,LB,UB,[],options);
end

warning on
% Radius of the spere
R=sqrt(-W);
if ( abs(sum(alpha)-1) > TOL)
    fprintf('Optimization problem!');
    keyboard
end
