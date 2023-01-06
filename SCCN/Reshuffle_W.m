function [A_ID, Aindex, Alist, B_ID, Bindex, Blist]=Reshuffle_W(Ka,Kb,Uc_idx,Ud_idx) 
    %%%% This function resuffles the input inference matrix W based on the
    %%%% extracted network structures by SCCN, where the detected densely
    %%%% altered sub-networks are pushed to the top for clearer
    %%%% visualization.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs:
    %%%%%  W:      an n-by-m connectivity inference matrix. Specifically, W= -log(p_ij), where p_ij is the p-value to test if the
    %%%%%          voxel-pair (i,j) is associated with a clinical/behaviral condition of interest
    %%%%%  S_A, S_B: the infrastructure graphs that store the spatial-contiguity
    %%%%%           constraint info 
    %%%%%  r: the cutoff used in the objective function 
    %%%%%  lambda: the turning parameter in the objective function
    %%%%%  num_skips: the number of skips between iterations to efficiently process
    %%%%%  ultra-high dimension data. 
    %%%%%  kmeans_iter: numbers of iterations set for k-means clustering
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