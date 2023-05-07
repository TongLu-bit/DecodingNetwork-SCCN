function [W_output, CID, Clist, C, z1]=SCCN_within(W, SA, threshold, lambda, num_iter, skip_step)
%%%% This function is for detecting covariated-related subnetworks


%%% Inputs:
%%%%%  W:   a 1 by n vector of t test p-values from the raw data. If
%%%%%           data contains k nodes, n=k*(k-1)/2

%%%%%  SA:  infrastructure graph for the ROI, say ROI A.  Each entry S_ii' in SA is a spatial-adjacency indicator variable between nodes i and i' in ROI A, 
%%%%%       where S_ii0= 1 if Euclidean distance d_ii'<= epsilon; S_ii'=0, ow.
%%%%%       e.g., in a 2D space, epsilon=sqrt(2); in a 3D space, epsilon=sqrt(3)
%%%%%  (toy example to construct SA in a 2D space:
% % % N=30; %suppose a 30*30 grid space
% % % indx_1d_A=[1:N^2];
% % % indx2d=reshape(indx_1d_A,N,N); % 30*30
% % % [indxI,indxJ]=ind2sub(size(indx2d),indx_1d_A);
% % % IJ=[indxI',indxJ'];
% % % A_1d=pdist(IJ,'chebychev');  % 900choose2 pairs
% % % A_2d=squareform(A_1d);  %900*900
% % % A_2d_adj=A_2d;              %only adjacent vovels have distance=1
% % % A_2d_adj(A_2d_adj>1.5)=0;    %all other distance=0
% % % A_spa=A_2d_adj;  %900*900 matrix A
% % % figure;imagesc(A_spa);


%%%%%  threshold:      A threshold on the entries of W. We only do the clustering on the
%%%%%           significant edges.
%%%%% lambda:          turning parameter. to balance size and density of
%%%%% sub-networks
%%%%% num_iter:         number of iterations to repeat SCCN_within sub-networkk
%%%%% detection operations in order to find optimal cluster number K
%%%%% skip_step :       number of steps to skip in interations


%%% Outputs:
%%%%%  Cindx:    the cluster/subnetwork index of every non-isolated node
%%%%%  CID:     the cluster/subnetwork index of every cluster in a power descending
%%%%%           order. i.e. CID(1) will be the cluster index of the most
%%%%%           concentrated cluster
%%%%%  Clist:   the reordered node index, nodes in the same cluster are
%%%%%           permuted together in such way: [find(Cindx==CID(1))
%%%%%           find(Cindx==CID(2)) ... find(Cindx==CID(K))]
%%%%%  C:        cluster information
%%%%%  z1:       indexes of nodes that pass screening



%% Preprocessing of data
warning('off');
%nlogp = -log(P)
% W=squareform(nlogp);
W1=W;
W(W1<threshold)=0;%Threshold on the p-values

%figure;imagesc(W)
z1=find(sum(W)>0); %Exclude the isolated nodes

% if(isempty(z1))
%     % if after the screening the matrix is all zero
%     Cindx = ones(1,size(W1,1));
%     CID=1;
%     Clist = 1:size(W1,1);
% else

W=W(z1,z1);

SA_new=SA(z1,z1);

WA=W.*SA_new;  
%figure; imagesc(WA);ax=gca;ax.FontSize=18;colorbar; colormap jet;

%Obtain Laplacian Matrix
degsA = sum(WA, 2);
DA    = sparse(1:size(WA, 1), 1:size(WA, 2), degsA);
LA=DA-WA;


%% Determine the number of clusters K
Mk=[];
Qual=[];
lenW=length(find(WA>0))/2 % # of sig. pts in Wa  e.g, lenW=137
diff   = eps;




[Ua, Ev_a] = eigs(LA,40,'smallestreal');  
%figure;hist(diag(Ev_a),80);
% figure;plot(diag(Ev_a),'*');
% diag(Ev_a);


for m=1:num_iter
    Prp_net=[];
    for K=1:skip_step:size(WA,1)
        K
%         try
%             [Ua, Ev_a] = eigs(LA,K, diff);
%             %returns k eigenvalues based on the value of sigma. 
%             %(The eigenvalues closest to the number sigma.)
%         catch
%             [Ua, Ev_a] = eigs(LA,K);
%         end

       
        C=kmeans(Ua,K,'Replicates',5);
        %Specify iter*2 replicates to help find a lower, local minimum.
       
       %'Replicates' ? Number of times to repeat clustering 
       %using new initial cluster centroid positions
       
        indx=[]; %indx
        A_net=[];% #edge in the net
        net_V=[];% size of each cluster
        C_net=[];
        for k=1:K
            indx=[indx;find(C==k)];
            net_V(k)=length(find(C==k));
            WC=WA(find(C==k), find(C==k));
            C_net(k) = sum(WC(find(WC>0)))/2;   % sum Wij (all sig. pts in cluster k)
            A_net(k)=(net_V(k)*(net_V(k)-1))/2;
        end
        Prp_net(K)=sum(C_net)^lambda *  ( sum(C_net) / sum(A_net) )^(2-lambda);
    end
        K = find(Prp_net == max(Prp_net)); %Prp_net?????index
    %K = K(1);%In case several k's give the same Prp_net value
    Mk(m,:)=[K max(Prp_net)]; %Mk????m?
    
    Qual(:,m)=Prp_net;
end

 figure;plot(Qual,'x');
 title('Optimzation of Objective Function','FontSize',20,'FontWeight','bold','Color','k');
xlabel("Number of clusters",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("Output value of Objective Function",'FontSize',20,'FontWeight','bold','Color','k');


K=Mk(find(Mk(:,2)==max(Mk(:,2))),1);
K=K(1)  %451 for 900*900

    
%% Find the cluster ID for each of the nodes

try
    [Ua, Ev_a] = eigs(LA,K, diff);
catch
   [Ua, Ev_a] = eigs(LA,K);
end

C=kmeans(Ua,K,'Replicates',8);
indx=[]; 
A_net=[];
net_V=[];
C_net=[];

for k=1:K
    indx=[indx;find(C==k)];
    net_V(k)=length(find(C==k));
    WC=WA(find(C==k), find(C==k));
%     C_net(k)=length(find(WC>0))/2;
    C_net(k) = sum(WC(find(WC>0)))/2;
    A_net(k)=(net_V(k)*(net_V(k)-1))/2;
end


%Pard_diagscore=(C_net).^2 ./(A_net)/lenW;  

Pard_diagscore=(C_net).^lambda * (C_net/(A_net))^(2-lambda);  
%Pard_diagscore=sum(C_net)^1.2 *  ( sum(C_net) / sum(A_net) )^(2-1.2);

%isnan(A) returns a logical array containing 1 (true) 
%where the elements of A are NaN
Pard_diagscore(isnan(Pard_diagscore))=0;  %change NaN to 0

[Pard_diagscore_sort,Pard_diagscore_sortID]=sort(Pard_diagscore,'descend');
   
Pard_inx_imporance=[];
for i=1:K
   Pard_inx_imporance=[  Pard_inx_imporance; find(C==Pard_diagscore_sortID(i))];
end

Cindx = 1:size(W,1); 
%z1=find(sum(W)>0); %give the lables of the non-isolated nodes

Cindx(z1)=C; 

%setdiff(A,B): Find the values in A that are not in B.
Cindx(setdiff(1:size(W,1),z1))=-1;
CID=Pard_diagscore_sortID;
Clist = z1(Pard_inx_imporance);
Clist = [Clist setdiff(1:size(W,1),z1)];
W_output=W1(Clist,Clist);

figure;
imagesc(W_output);colormap jet;colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';
title('Community subnetwork in heatmap','FontSize',20,'FontWeight','bold','Color','k');
xlabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');

end
