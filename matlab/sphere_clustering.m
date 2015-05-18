% function [labels,spheres,spheres_init]=sphere_clustering(X,NUM_CLUSTERS,C,kernel_kind,width)
% Sphere-based Clustering for sets of vectors
%   Uses a naive, quadprog based, solution of the minimal sphere problem.
% Input arguments:
% X -> Data as cell array of DxNi matrices
% NUM_CLUSTERS
% C -> Penalty parameter for soft-sphere packing (if hard packing is
% desired, C=[])
% kernel_kind, width -> Kernel parameters

% Darío García, 01/06/2010

function [labels,spheres,spheres_init,D]=sphere_clustering(X,NUM_CLUSTERS,C,kernel_kind,width)


if (nargin<3)
    C=[];
end
if (nargin<4)
    kernel_kind='linear';
end
if (nargin<5)
    width=[];
end

N=length(X);
R=cell(N,1);
alpha=R;
SVs=R;
THRESH_COEF=0.1;
K=cell(N,1);

% Learn the individual spheres
for i=1:N
    % Pack the sequence
    tmp=struct();
    Ker=calckernel(kernel_kind,width,X{i}');
    if (isempty(C))
        [R{i},alpha_tmp]=min_sphere(Ker);
    else        
        [R{i},alpha_tmp]=soft_sphere(Ker,C);        
    end
    THRESH=THRESH_COEF/length(alpha_tmp);
    indexes=find(alpha_tmp>THRESH);
    SVs{i}=X{i}(:,indexes);
    alpha{i}=alpha_tmp(indexes);
    K{i}=Ker(indexes,indexes);
end
spheres=struct('R',R,'SVs',SVs,'alpha',alpha);
spheres_init=spheres;

% Distances between centers
D=zeros(N);
for i=1:N-1
    alpha_i=spheres(i).alpha;
    for j=i+1:N
        alpha_j=spheres(j).alpha;
        Kcross=calckernel(kernel_kind,width,spheres(i).SVs',spheres(j).SVs');
        D(i,j)=sqrt(alpha_i'*K{i}*alpha_i+alpha_j'*K{j}*alpha_j-2*alpha_j'*Kcross*alpha_i);
    end
end
D=D+D';
done=0;

labels=[1:N]';
NUM_SPHERES=N;
R=[spheres(:).R];
tmp=R'*ones(1,N);
tmp=(tmp+tmp')/2;


% Hierarchical sphere fusion

while (~done)
    % We want to build the smallest sphere encompassing two existing
    % spheres
    
    R=[spheres(:).R];
    tmp=R'*ones(1,NUM_SPHERES);
    tmp=(tmp+tmp');
    RADIUS1=R'*ones(1,NUM_SPHERES);
    RADIUS2=RADIUS1';
    MAX_DIAM=max(2*RADIUS1,2*RADIUS2);
    J=max(D+tmp+1e3*eye(size(D)),MAX_DIAM); % Increase the elements in the diagonal so that they are not the minimum
    [val,ind]=min(J(:));
    [i1,i2]=ind2sub(size(J),ind);

    
    % Create the data for the new sphere
    R_new=val/2;
    
    tmp=struct();    
    tmp.SVs=[spheres(i1).SVs spheres(i2).SVs];
    
    a_t=[-1*spheres(i1).alpha;spheres(i2).alpha;];
    a_t=a_t/D(i1,i2);
    tmp.alpha=[spheres(i1).alpha;zeros(length(spheres(i2).alpha),1)]+a_t*(val/2-spheres(i1).R);
    tmp.R=R_new;

            
    %%
    % Update distance matrix
    %%    
    idx_right=setdiff(1:NUM_SPHERES,[i1 i2]);
    NUM_SPHERES=NUM_SPHERES-1;
    D=D(idx_right,idx_right)/2; % So we can later on do D+D'    
    % Labels for the new sphere
    idx_new=(find(labels==i1 | labels==i2));
    for i=1:length(idx_right)
        labels(find(labels==idx_right(i)))=i;
    end
    labels(idx_new)=NUM_SPHERES;
    
    % Update the structure holding the active spheres
    spheres=spheres(idx_right);           
    spheres=[spheres;tmp];   

    K_i=calckernel(kernel_kind,width,spheres(end).SVs');
    alpha_i=spheres(end).alpha;
    for j=1:NUM_SPHERES-1
        K_j=calckernel(kernel_kind,width,spheres(j).SVs');
        alpha_j=spheres(j).alpha;
        Kcross=calckernel(kernel_kind,width,spheres(end).SVs',spheres(j).SVs');
        D(NUM_SPHERES,j)=sqrt(alpha_i'*K_i*alpha_i+alpha_j'*K_j*alpha_j-2*alpha_j'*Kcross*alpha_i);
    end
    D(NUM_SPHERES,NUM_SPHERES)=0;
    D=D+D';
    
    % End condition 
    if (NUM_SPHERES==NUM_CLUSTERS)
        done=1;
    end
end

