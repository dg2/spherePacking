function [X, labels_real]=generate_data(N, L)
    mu{1}=[0 0]';
    mu{2}=[3 3]';
    labels_real=(rand(N,1)>0.5)+1;
    X=cell(N,1);
    R=X;
    alpha=X;
    SVs=X;
    center=X;
    THRESH=1/(10*L);
    NUM_CLUSTERS=2;
    COLORS=['r','g','b','c','m','y','k'];
    C=[];
    %C=0.2;

    for i=1:N
        X{i}=randn(2,L)+repmat(mu{labels_real(i)},1,L);
    end
