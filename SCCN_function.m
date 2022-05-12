%load('/Users/bright1/Dropbox/research1-btw/Realdata/De0_Feb07/WWR.mat');
function [A_ID, Aindex, Alist, B_ID, Bindex, Blist]=SCCN(W0,cutoff0, SA, SB, r, lambda, num_skips, kmeans_iter, fig)
    %%%% This function applies the SCCN algorithm to extract sub-area pairs
    %%%% from a ROI pair. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs:
    %%%%%  W0:      an n-by-m vFC inference matrix, where each entry is a test -log(p)
    %%%%%           value from the raw data. 
    %%%%%  cutoff0: the pre-selected cutoff for W0 to first exclude most
    %%%%%           non0informative false-positive edges
    %%%%%  SA, SB: the infrastructure graphs that store the spatial-contiguity
    %%%%%           constraint info 
    %%%%%  r: the cutoff used in the objective function 
    %%%%%  lambda: the turning parameter in the objective function
    %%%%%  num_skips: the number of skips between two iterations in Algorihtm 1 due to
    %%%%%  ultra-high dimension of brain imaging data. 
    %%%%%  kmeans_iter: numbers of iterations set for the k-means clustering
    %%%%%  fig: fig=1 if we want to view the output performance of different cluster sizes; fig=0 otherwise
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs:
    
    %%%%%  A_ID:     the cluster index of each cluster in ROI A in a power descending
    %%%%%           order. i.e. A_ID(1) will be the cluster index of voxels in ROI A coming from the most
    %%%%%           concentrated sub-area pair
    %%%%%  Aindex:    the cluster index of every voxel in ROI A
    %%%%%  Alist:   the reordered voxel index, voxel in the same cluster are
    %%%%%           permuted together in such way: [find(Aindx==A_ID(1))
    %%%%%           find(Aindx==A_ID(2)) ... find(Aindx==A_ID(K))]
    %%%%%  B_ID, Bindex, Blist are defined similarly as above. 
        
    %% Data Processing 
    nrow=size(W0,1);
    ncol=size(W0,2);
    
    W=W0; W(W<cutoff0)=0; %perform screening on W0
    %figure;imagesc(W)
    
    WA=W*W.'.*SA;  WB=W.'*W.*SB; %incorporate spatial info W
        
    %% Spectral Clustering 
    degsA = sum(WA, 2);
    DA    = sparse(1:size(WA, 1), 1:size(WA, 2), degsA);
    LA=DA-WA; %Laplacian matrix for ROI A
    
    degsB = sum(WB, 2);
    DB    = sparse(1:size(WB, 1), 1:size(WB, 2), degsB);
    LB=DB-WB; %Laplacian matrix for ROI B
    
    
    %%%%% eigen decomposition on LA and LB
    [Ua, Ev_a] = eigs(LA,50,'smallestreal'); 
    %figure;hist(diag(Ev_a),50);
    %figure;plot(diag(Ev_a),'*')
    
    %UB
    [Ub, Ev_b] = eigs(LB,50,'smallestreal');
    %figure;hist(diag(Ev_b),50);
    %figure;plot(diag(Ev_b),'*')
    
  
    %%%%% Determine the optimal numbers of clusters Ka, Kb, for ROI A and ROI B,
    %%%%% respectively 
    Cqual=[]; %the output values of the objective function for each (Ka, Kb) pair
    for Ka=1:num_skips:nrow 
         Ca=kmeans(Ua,Ka,'Replicates',kmeans_iter); 
        for Kb=1:num_skips:ncol 
            Cb=kmeans(Ub,Kb,'Replicates',kmeans_iter);
            
            output=[];          
            for i=1:Ka         
                Alist=find(Ca==i); %the i-th partion for ROI A             
                for j=1:Kb
                    Blist=find(Cb==j); %the j-th partion for ROI B                 
                    Wsub= W (Alist,Blist ); %submatrix of the sub-area pair (Ui,Vj)   
                    supraWsub=sum(Wsub(find(Wsub>r)));  
                    ab= length(Alist) * length (Blist); %size of the sub-area pair                      
                    output(i,j)=( supraWsub )^lambda * ( supraWsub / ab )^(2-lambda); %equivalent to the objective function                   
                end                    
            end                          
            Cqual(Ka,Kb)= sum (sum (output)) ;
        end
    end 

    if fig==1
    figure; surf(Cqual(1:num_skips:nrow,1:num_skips:ncol));
    end
     
    %%%find the optimal Ka, Kb that gives the maximum output value:
    K=find(Cqual == max(max(Cqual)) );
    K = K(1); % in case there are multiple maximizers
    [Ka,Kb]=ind2sub([size(Cqual,1),size(Cqual,2)],K);

    %% Plug in the optimal Ka, Kb, and find the cluster ID for each voxel 
    Ca1=kmeans(Ua,Ka,'Replicates',kmeans_iter);
    Cb1=kmeans(Ub,Kb,'Replicates',kmeans_iter);

    output1=[];          
    for i=1:Ka         
        Alist=find(Ca1==i); %the i-th partion for ROI A             
        for j=1:Kb
            Blist=find(Cb1==j); %the j-th partion for ROI B                 
            Wsub= W (Alist,Blist ); %submatrix of the sub-area pair (Ui,Vj)   
            supraWsub=sum(Wsub(find(Wsub>r)));  
            ab= length(Alist) * length (Blist); %size of the sub-area pair                      
            output1(i,j)=( supraWsub )^lambda * ( supraWsub / ab )^(2-lambda); %equivalent to the objective function                   
        end                    
    end 

    diagscore=output1;
    diagscore(isnan(diagscore) |  isinf(diagscore) )=0;
    %max(max(diagscore))

    %% Results: 
    a=diagscore;
    [R,C] = ndgrid(1:size(a,1),1:size(a,2));
    [b,idx] = sort(a(:),'descend'); % b: sorted output values
    voxel_ID=[R(idx),C(idx)]; % the corresonding voxel (i,j) position to the sorted values

    IndA=voxel_ID(:,1);
    [bb,i,j]=unique(IndA, 'first');
    A_ID=IndA(sort(i));  %the cluster index of each cluster in a power descending order for ROI A
    Alist=[]; % the reordered voxel index in ROI A
    Ca=Aindex; %the cluster index of each voxel in ROI A 
    for i=1:length(A_ID)
        Alist=[Alist;find(Aindex==A_ID(i))];
    end
    
    IndB=voxel_ID(:,2);
    [bb,i,j]=unique(IndB, 'first');
    B_ID=IndB(sort(i)); %the cluster index of each cluster in a power descending order for ROI B
    Cb=Bindex; %the cluster index of each voxel in ROI B
    Blist=[]; % the reordered voxel index in ROI B
    for i=1:length(B_ID)
        Blist=[Blist;find(Bindex==B_ID(i))];
    end

    %Worder= W(Alist,Blist');
    %figure;
    % imagesc(Worder ); colorbar 
    % colormap jet

end 


    




