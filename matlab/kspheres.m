function [assign,R,SVs]=kspheres(Ker,C);

MAX_ITER=1000;
iter=0;
N=length(Ker);
converged=0;
SVs=cell(C,1);

% Initial labeling: Choose some points as centroids and assign by minimal
% distance
tmp=randperm(N);
idx_init=tmp(1:C);
for i=1:C
   % keyboard
    D(:,i)=sqrt(diag(Ker)-2*Ker(:,idx_init(i))+Ker(idx_init(i),idx_init(i)));    
    [dummy,idx_init(i+1)]=max(D(:,i));
end

% Assign
[dummy,assign]=min(D,[],2);

% Main loop
while (~converged && iter<MAX_ITER)
    iter=iter+1;
    assign_old=assign;
    % Obtain the centres of the minimal spheres and the distances
    for c=1:C
        idx=find(assign==c);
        Kc=Ker(idx,idx);
        [R(c,iter),alpha{c}]=min_sphere(Kc);
        % We can threshold alpha to get a sparse vector
        D(:,c)=sqrt(diag(Ker)-2*Ker(:,idx)*alpha{c}+alpha{c}'*Ker(idx,idx)*alpha{c});    
        SVs{c}=idx(alpha{c}>THRESH);        
    end  
    % Check if the centers are inside the other balls (check that distance
    % between centers is larger than sum of radiuses)
    
    
    % Re-assign
    [dummy,assign]=min(D,[],2);
    fprintf('Iter. %i -> Sum of radiuses: %f, \t Mean distortion: %f\n',iter,sum(R(:,iter)),mean(dummy));
    
    % TODO: Why are we converging but having different radius sums? Looks
    % like sometimes the balls overlap ... Should we force the center of
    % one ball outside the other? We should have some additional constraint
    
    % Check for convergency
    if (all(assign==assign_old))
        fprintf('Convergence!\n');
        converged=1;
    end
    

end