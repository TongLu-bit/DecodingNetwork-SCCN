function [A_ID, Alist, B_ID, Blist]=Reshuffle_W(Ka,Kb,A_idx,B_idx,W, r, lambda) 
    %%%% This function resuffles the input inference matrix W based on the
    %%%% extracted network structures by SCCN, where the detected densely
    %%%% altered sub-networks are pushed to the top for clearer
    %%%% visualization.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs:
    %%%%% Ka: number of clusters (i.e., sub-areas) partitioned for ROI A by SCCN
    %%%%% Kb: number of clusters (i.e., sub-areas) partitioned for ROI B by SCCN
    %%%%% A_idx: the node membership of ROI A
    %%%%% B_idx: the node membership of ROI B
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs:
    
    %%%%%  A_ID:     the cluster index of each cluster in ROI A in a power descending
    %%%%%           order. i.e. A_ID(1) will be the cluster index of voxels in ROI A coming from the most
    %%%%%           concentrated sub-area pair
    %%%%%  Alist:   the reordered voxel index, voxel in the same cluster are
    %%%%%           permuted together in such way: [find(Aindx==A_ID(1))
    %%%%%           find(Aindx==A_ID(2)) ... find(Aindx==A_ID(K))]
    %%%%%  B_ID, Blist are defined similarly as above. 
    
    
    % Calculate the output of objective function for the cluster numbers, Ka, Kb, determined by
    % function 'SCCN_alg' 
    output=[];          
    for i=1:Ka         
        Alist=find(A_idx==i); %the i-th partion for ROI A             
        for j=1:Kb
            Blist=find(B_idx==j); %the j-th partion for ROI B                 
            Wsub= W (Alist,Blist ); %submatrix of the sub-area pair (Ui,Vj)   
            supraWsub=sum(Wsub(find(Wsub>r)));  
            ab= length(Alist) * length (Blist); %size of the sub-area pair                      
            output(i,j)=( supraWsub )^lambda * ( supraWsub / ab )^(2-lambda); %equivalent to the objective function                   
        end                    
    end 

    diagscore=output;
    diagscore(isnan(diagscore) |  isinf(diagscore) )=0;
    %max(max(diagscore))

    %% Reshuffle the inference matrix W in order of sub-area-pair density: 
    [R,C] = ndgrid(1:size(diagscore,1),1:size(diagscore,2)); %create a rectangular grid 
    [sorted_score,idx] = sort(diagscore(:),'descend');
    voxel_ID=[R(idx),C(idx)]; % the corresonding voxel-pair (i,j) position based on the sorted scores

    IndA=voxel_ID(:,1);
    [sorted_A,i,j]=unique(IndA, 'first');
    A_ID=IndA(sort(i));  %the cluster index of each cluster in a power descending order for ROI A
    Alist=[]; % the reordered voxel index for ROI A
    for i=1:length(A_ID)
        Alist=[Alist;find(A_idx==A_ID(i))];
    end
    
    IndB=voxel_ID(:,2);
    [bb,i,j]=unique(IndB, 'first');
    B_ID=IndB(sort(i)); %the cluster index of each cluster in a power descending order for ROI B
    Blist=[]; % the reordered voxel index for ROI B
    for i=1:length(B_ID)
        Blist=[Blist;find(B_idx==B_ID(i))];
    end

    %Visualize reshuffled W:
    Worder= W(Alist,Blist');
    figure;
    imagesc(Worder ); colorbar; colormap jet;
    set(gca, 'clim', [0 15]);
    title('Reshuffled W - After SCCN');    
end 