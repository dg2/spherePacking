% Sphere-packing Clustering for sets of vectors
% In this demo we show the evolution of the sphere packing algorithm for
% clustering sets of vectors, i.e. the hierarchical merging of spheres
% representing the support of the sets. 


% Darío García, 31/05/2010

% Generate some data 
mu{1}=[0 0]';
mu{2}=[3 3]';
N=10;
labels_real=(rand(N,1)>0.5)+1;
X=cell(N,1);
R=X;
alpha=X;
SVs=X;
center=X;
L=50;
THRESH=1/(10*L);
NUM_CLUSTERS=2;
COLORS=['r','g','b','c','m','y','k'];
C=[];
%C=0.2;

for i=1:N
    X{i}=randn(2,L)+repmat(mu{labels_real(i)},1,L);
    % Pack the sequence
    tmp=struct();
    Ker=calckernel('linear',[],X{i}');
    if (isempty(C))
        [R{i},alpha_tmp]=min_sphere(Ker);
    else
        [R{i},alpha_tmp]=soft_sphere(Ker,C);
    end
    indexes=find(alpha_tmp>THRESH);
    SVs{i}=X{i}(:,indexes);
    alpha{i}=alpha_tmp(indexes);
    center{i}=sum(repmat(alpha_tmp',2,1).*X{i},2);
end
spheres=struct('R',R,'SVs',SVs,'alpha',alpha,'center',center);
spheres_init=spheres;

D=zeros(N);
for i=1:N-1
    K_i=calckernel('linear',[],spheres(i).SVs');
    alpha_i=spheres(i).alpha;
    for j=i+1:N
        K_j=calckernel('linear',[],spheres(j).SVs');
        alpha_j=spheres(j).alpha;
        Kcross=calckernel('linear',[],spheres(i).SVs',spheres(j).SVs');
        D(i,j)=sqrt(alpha_i'*K_i*alpha_i+alpha_j'*K_j*alpha_j-2*alpha_j'*Kcross*alpha_i);
        if (abs(D(i,j)-norm(spheres(i).center-spheres(j).center))>1e-5)
            error('Distances are not preserved: ¿Threshold is too high?');
        end
    end
end
D=D+D';
done=0;
NUM_SPHERES=N;
R=[spheres(:).R];
tmp=R'*ones(1,N);
tmp=(tmp+tmp')/2;

% Initial plot
figure;
hold on
for i=1:N
    scatter(X{i}(1,:),X{i}(2,:),5,COLORS(mod(i,length(COLORS))+1));
    circle(spheres(i).center,spheres(i).R,100,sprintf('%c--',COLORS(mod(i,length(COLORS))+1)));
end
idx_spheres=num2cell(1:N);
labels=[1:N]';
hold off

% Main loop
while (~done)
    %% Find the two closests spheres
    R=[spheres(:).R];
    tmp=R'*ones(1,NUM_SPHERES);
    tmp=(tmp+tmp');
    RADIUS1=R'*ones(1,NUM_SPHERES);
    RADIUS2=RADIUS1';
    MAX_RAD=max(2*RADIUS1,2*RADIUS2);
    J=max(D+tmp+1e8*eye(size(D)),MAX_RAD); % Increase the elements in the diagonal so that they are not the minimum
    [val,ind]=min(J(:));
    [i1,i2]=ind2sub(size(J),ind);
    fprintf('Fusing spheres %i and %i\n',i1,i2);

    
    %% Find the center and radius of the encompassing sphere
    vec=(spheres(i2).center-spheres(i1).center)/D(i1,i2);
    c_new=spheres(i1).center+vec*(val/2-spheres(i1).R);   
    R_new=val/2;
    
    tmp=struct();    
    tmp.SVs=[spheres(i1).SVs spheres(i2).SVs];
    
    a_t=[-1*spheres(i1).alpha;spheres(i2).alpha;];
    a_t=a_t/D(i1,i2);
    tmp.alpha=[spheres(i1).alpha;zeros(length(spheres(i2).alpha),1)]+a_t*(val/2-spheres(i1).R);    
    tmp.R=R_new;
    tmp.center=c_new;
  
       
    % A good bound on the slack ratio (and simple to calculate) involves turning the kernel matrix into a distance matrix
    % and look for the furthest apart points. All the required kernel
    % matrices are already calculated (or need to be calculated somewhere),
    % despite in this code we may be recalculating them for simplicity
    Ker_tmp=calckernel('linear',[],[spheres(i1).SVs spheres(i2).SVs]');
    D_tmp=sqrt(kernel_to_dist(Ker_tmp));
    bound_svs=0.5*max(D_tmp(:));
  
    %%
    % Update distance matrix
    %%    
    idx_right=setdiff(1:NUM_SPHERES,[i1 i2]);
    NUM_SPHERES=NUM_SPHERES-1;
    D=D(idx_right,idx_right)/2; % So we can later on do D+D'    
    % Labels for the new sphere
    idx_new=(find(labels==i1 | labels==i2));
    for i=1:length(idx_right)
        labels(labels==idx_right(i))=i;
    end
    labels(idx_new)=NUM_SPHERES;
    
    % Update the structure holding the active spheres
    spheres=spheres(idx_right);           
    spheres=[spheres;tmp];   


    K_i=calckernel('linear',[],spheres(end).SVs');
    alpha_i=spheres(end).alpha;
    for j=1:NUM_SPHERES-1
        K_j=calckernel('linear',[],spheres(j).SVs');
        alpha_j=spheres(j).alpha;
        Kcross=calckernel('linear',[],spheres(end).SVs',spheres(j).SVs');
        D(NUM_SPHERES,j)=sqrt(alpha_i'*K_i*alpha_i+alpha_j'*K_j*alpha_j-2*alpha_j'*Kcross*alpha_i);
        if (abs(D(NUM_SPHERES,j)-norm(spheres(NUM_SPHERES).center-spheres(j).center))>1e-5)
            keyboard
        end
    end
    D(NUM_SPHERES,NUM_SPHERES)=0;
    D=D+D';

    
    % Plot the results
    clf
    hold on    
    for i=1:NUM_SPHERES-1    
        tmp=[X{find(labels==i)}];
        scatter(tmp(1,:),tmp(2,:),5,COLORS(mod(i,length(COLORS))+1));
        circle(spheres(i).center,spheres(i).R,100,sprintf('%c--',COLORS(mod(i,length(COLORS))+1)));
        plot(spheres(i).SVs(1,:),spheres(i).SVs(2,:),sprintf('s%c',COLORS(mod(i,length(COLORS))+1)));
    end
    tmp=[X{find(labels==NUM_SPHERES)}];
    scatter(tmp(1,:),tmp(2,:),20,'r');
    circle(spheres(NUM_SPHERES).center,spheres(NUM_SPHERES).R,100,'r');
    plot(spheres(NUM_SPHERES).SVs(1,:),spheres(NUM_SPHERES).SVs(2,:),'sr');

     hold off
    
    % Compare the "encompassing" sphere with the minimum spheres containing
    % just the data vectors in the corresponding sequences
    data=[X{labels==NUM_SPHERES}];
    Ker=calckernel('linear',[],data');
    
    if (~isempty(C))
        [R_tmp,alpha_tmp]=soft_sphere(Ker,C);
    else
        [R_tmp,alpha_tmp]=min_sphere(Ker);
    end
    if (R_tmp==0)
        keyboard;
    end
    fprintf('Radius of the new sphere: %f \n', R_new);
    fprintf('Radius of the minimal sphere: %f => Slack ratio %.4f \n',R_tmp,R_new/R_tmp);

    
    % Compare also with the bounds
    % Data-dependent bound: 
    fprintf('SVs-based bound: %f => Bound on the slack ratio %.4f \n',bound_svs,R_new/bound_svs);
       
    if (NUM_SPHERES==NUM_CLUSTERS)
        done=1;
    end

    keyboard;
    
end
err=cluster_error(labels,labels_real);
fprintf('Error: %.3f\n',err);
