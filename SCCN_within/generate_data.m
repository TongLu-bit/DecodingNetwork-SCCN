function [W, threshold_GT, case_mtx, ctrl_mtx] = generate_data(Da, Da_idx1, Da_idx2, case_num, ctrl_num, mu0, mu1, cohensd,show_figure)
%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%
%%% Da: a matrix describe the spatial location of covariate-related subnetworks
%%% Da_idx1: node index of dense subnetwork1
%%% Da_idx2: node index of dense subnetwork2

%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%
%%% W: Voxel-wise Inference Matrix
%%% threshold_GT: optimal threshold for p
%%% case_mtx: n*p matrix for case, where n=sample size, p=# FCs
%%% ctrl_mtx: n*p matrix for control, where n=sample size, p=# FCs

%%% Tong Lu, May 2023


%% initialization
N=size(Da,1)^2;
N_edge = N*(N-1)/2;

% an N*N matrix storing node-wise FC for one subject
mat_fc=zeros(N,N);
mat_fc(Da_idx1,Da_idx1)=mat_fc(Da_idx1,Da_idx1)+1;
mat_fc(Da_idx2,Da_idx2)=mat_fc(Da_idx2,Da_idx2)+1;
%figure;imagesc(mat_fc);
for i = 1:900
    mat_fc(i,i) = 0;
end
case_idx = find(squareform(mat_fc));
ctrl_idx = 1:N_edge;ctrl_idx(case_idx) = [];

case_mtx = zeros(case_num, N_edge);
sigmma0=(mu1-mu0)/cohensd;
sigmma1=sigmma0;

case_mtx(:,case_idx) = randn(case_num, length(case_idx))*sigmma1 + mu1; %randn, N(0,1)
case_mtx(:,ctrl_idx) = randn(case_num, length(ctrl_idx))*sigmma0 + mu0;  %randn: N(0,1)

% define controlling subjects' FC matrices
ctrl_mtx = randn(ctrl_num, N_edge)*sigmma0 + mu0;

%% obtain Voxel-wise Inference Matrix
[~,p_vec]=ttest2(case_mtx,ctrl_mtx);
W=squareform(-log(p_vec)); 

if show_figure==1
    figure; imagesc(W);ax=gca;ax.FontSize=18;colorbar; colormap jet;
    xlabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');
    ylabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');
    title("Input W",'FontSize',45,'FontWeight','bold','Color','k');
end


%% find a good cut for p-value
threshold_vec = [0.001 0.005 0.01 0.05 0.1];
f1score_vec = zeros(length(threshold_vec),1);
for i = 1:length(threshold_vec)
    target = ones(1,N_edge);
    target(ctrl_idx)=0;
    output = p_vec <= threshold_vec(i);
    [~,cm,~,~] = confusion(target,output);
    f1score_vec(i) = cm(2,2)/(cm(2,2) + 0.5*(cm(1,2) + cm(2,1)));
end

[~,thresh_idx] = max(f1score_vec);
threshold_GT=threshold_vec(thresh_idx);
% sig_edge = sum(p_vec <= threshold_GT);

end