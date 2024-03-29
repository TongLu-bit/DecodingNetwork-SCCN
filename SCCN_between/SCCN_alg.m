function [Ka, Kb, A_idx,B_idx]=SCCN_alg(W, S_A, S_B, r, lambda, num_skips, kmeans_iter, fig)
    %%%% This function applies the SCCN algorithm to yield true sub-area
    %%%% pair network structure from a ROI pair. 
    
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
    
    %%%%% Ka: number of clusters (i.e., sub-areas) partitioned for ROI A by SCCN
    %%%%% Kb: number of clusters (i.e., sub-areas) partitioned for ROI B by SCCN
    %%%%% A_idx: the node membership of ROI A
    %%%%% B_idx: the node membership of ROI B
        
    %% Incorporate spatial-contiguity constraint into W
    nrow=size(W,1);ncol=size(W,2);
    WA=W*W.'.*S_A;  WB=W.'*W.*S_B; %
        
    %% Ratio-cut partitioning WA and WB
    % Calculate the Laplacian matrix for ROI A and B
    degsA = sum(WA, 2);
    DA    = sparse(1:size(WA, 1), 1:size(WA, 2), degsA);
    LA=DA-WA; 
    
    degsB = sum(WB, 2);
    DB    = sparse(1:size(WB, 1), 1:size(WB, 2), degsB);
    LB=DB-WB; 
    
    %%%%% Eigen decomposition on Laplacian matrix LA and LB
    [Ua, Ev_a] = eigs(LA,50,'smallestreal'); 
    %figure;hist(diag(Ev_a),50);
    %figure;plot(diag(Ev_a),'*')
    
    [Ub, Ev_b] = eigs(LB,50,'smallestreal');
    %figure;hist(diag(Ev_b),50);
    %figure;plot(diag(Ev_b),'*')
    
  
    %%%%% Determine the optimal numbers of clusters Ka, Kb, for ROI A and ROI B respectively 
    Cqual=[]; %the output values of the objective function for each (Ka, Kb) pair
    for Ka=1:num_skips:nrow 
         disp('Display Ka(# of cluters) to show iteration progress:');
         Ka  
         Ca=kmeans(Ua,Ka,'Replicates',kmeans_iter);   %cluters for nodes in ROI A
        for Kb=1:num_skips:ncol 
            Cb=kmeans(Ub,Kb,'Replicates',kmeans_iter);  %cluters for nodes in ROI B
            
            output=[];          
            for i=1:Ka         
                Alist=find(Ca==i); %the i-th partion for ROI A             
                for j=1:Kb
                    Blist=find(Cb==j); %the j-th partion for ROI B                 
                    Wsub= W (Alist,Blist ); %submatrix of sub-area pair (U_c,V_d)   
                    supraWsub=sum(Wsub(find(Wsub>r)));  
                    ab= length(Alist) * length (Blist); %size of sub-area pair (U_c,V_d)                        
                    output(i,j)=( supraWsub )^lambda * ( supraWsub / ab )^(2-lambda); %equivalent to the objective function                   
                end                    
            end                          
            Cqual(Ka,Kb)= sum (sum (output)) ;
        end
    end 

    if fig==1
    figure; surf(Cqual(1:num_skips:nrow,1:num_skips:ncol));
    title('Outputs of objective function given different values of (Ka, Kb)')
    xlabel('Ka (# clusters for ROI A)')
    ylabel('Kb (# clusters for ROI B)')
    zlabel('Outputs of objective function')
    end
     
    %%%%% find the optimal Ka, Kb that gives the maximum output value:
    K=find(Cqual == max(max(Cqual)) );
    K = K(1); % in case there are multiple maximizers
    [Ka,Kb]=ind2sub([size(Cqual,1),size(Cqual,2)],K); %optimal numbers of clusters 

    %%%%% Sub-area network structure
    A_idx=kmeans(Ua,Ka,'Replicates',kmeans_iter); %node memberships of ROI A
    B_idx=kmeans(Ub,Kb,'Replicates',kmeans_iter);  %node memberships of ROI B
end 