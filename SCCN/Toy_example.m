% This is a top example to implement SCCN on a synthetic data. 

%%%%%%%%%%%%%%% Step 1. Load, construct, and process input data %%%%%%%%%%%%%%%
% 1.1. Load W0, the inference connectivity matrix between region A and B
% Details: 1. 900 voxels in A and 1600 voxels in B, where Region A is a
% 30*30 square, Region B is a 40*40 square
%          2. W0= -log(p_ij), where p_ij is the p-value to test if the
%         voxel-pair (i,j) is associated with a clinical/behaviral
%         condition of interest
%          3. Voxel i is from region A; voxel j is from region B

load('W0.mat');
figure; imagesc(W0);colormap jet;colormap;
title('W0 - raw data');

% 1.2. Perform screening on W0 
% Details: 1. The post-screened inference matrix W can effectively exclude 
%             most non-informative false-positive edges
%          2. We choose p=0.05 as the screening thereshold
%          3. The real threshold for W0= -log(p_ij) is -log(0.05)~=3
cutoff=3;
W=W0;W(W<cutoff)=0;
figure; imagesc(W);colormap jet;colormap;
title('W - Before SCCN');

%1.3. Construct S_A, S_B, the infrastructure graph
% Details: 1. S_ii'=1 if d_ii' <=epsilon; S_ii'=0 otherwise
%          2. Since the similation data is in 2D, epsilon set to be sqrt(2)
%          3. S_A, S_B prescribe the spatial-contiguity constraints 

%Construct S_A
N_A=sqrt(size(W,1)); %side length of Region A 
indx_A=reshape(1:N_A^2,N_A,N_A); 
[indx_A_x,indx_A_y]=ind2sub(size(indx_A),1:N_A^2);
IJ_A=[indx_A_x',indx_A_y']; %2D coordinates of each node in A 

A_1d=pdist(IJ_A,'chebychev');  %pairwise chebychev distance between all 900 voxels in A
A_2d=squareform(A_1d);  %900*900 distance matrix
S_A=A_2d;              
S_A(S_A>sqrt(2))=0;    %set distance between non-adjacent voxels to 0
figure;imagesc(S_A);

%Construct S_B
N_B=sqrt(size(W,2)); %side length of Region B 
indx_B=reshape(1:N_B^2,N_B,N_B); 
[indx_B_x,indx_B_y]=ind2sub(size(indx_B),1:N_B^2);
IJ_B=[indx_B_x',indx_B_y']; %2D coordinates of each node in B

B_1d=pdist(IJ_B,'chebychev');  %pairwise chebychev distance between all 900 voxels in A
B_2d=squareform(B_1d);  %1600*1600 distance matrix
S_B=B_2d;              
S_B(S_B>sqrt(2))=0;    %set distance between non-adjacent voxels to 0
figure;imagesc(S_B);


%%%%%%%%%%%%%%% Step 2. Implement SCCN with input {W, S_A, S_B} %%%%%%%%%%%%%%%
%set parameters based on prior domain knowledge 
r=3.5;  %threshold for W, set p=0.03 <=>-log(p)~=3.5
lambda=1.4;  %tuning paramter in objective function
num_skips=350; %the number of skips between iterations 
kmeans_iter=3; %numbers of iterations set for k-means clustering
fig=1 ;% 1 if we want to view the output performance of different cluster sizes; 0 otherwise

%Now call function SCCN_alg, see details described in the file "SCCN_alg.m"
[Ka, Kb, A_idx,B_idx]=SCCN_alg(W, S_A, S_B, r, lambda, num_skips, kmeans_iter, fig);

% Resuffles the input inference matrix W based on the extracted network
% structures by the SCCN_alg function above
[A_ID, Alist, B_ID, Blist]=Reshuffle_W(Ka, Kb, A_idx, B_idx, W, r, lambda); 
