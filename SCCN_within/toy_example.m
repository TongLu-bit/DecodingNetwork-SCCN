%%%%% This is a toy example which provides: 
% (1) simulations to generate a unipartite graph G with two dense subgraphs G1 and G2, both of
%%%%% which contain nodes with spatially contiguity(SC); 
% (2) implementation of SCCN-within (with SC integrated) on the simulated G to reveal latent topological
%%%%% structures (i.e., identify G0 (null), G1, G2); 
% (3) implementation of SCCN-within (without SC integrated) on the simulated G to reveal latent topological
%%%%% structures (i.e., identify G0 (null), G1, G2); 
% (4) network-based permutation test on extracted sub-networks for significance
% (5) plot detected subgraphs that pass the permutation test
% (6) performance evaluation of selected sub-networks in terms of TPR and FPR 

clc;clear;close all


%% (1) Simulations of G with two dense subgraphs G1 and G2

% construct G, ground truth
Da=zeros(30,30);
Da(3:12, 3:11)=1; %13*14=182 nodes
Da(10:25, 18:25)=2; %22*11=242 nodes

Da_show=Da;
temp_idx=find(Da_show==2);
Da_show(temp_idx)=1;
figure; imagesc(Da_show);ax=gca;ax.FontSize=18;ax.FontWeight='bold'; colormap summer;
title("ROI A, |R|=900 nodes",'FontSize',45,'FontWeight','bold','Color','k');

% add noise
rate_in=0.3;rate_out=0.3; %noise ratios
show_figure=1;
[Da_idx1,Da_idx2,Da_0]  = pick_case_idx(Da,rate_in, rate_out,show_figure);

 
% initialization
case_num=100;ctrl_num=100;
mu0=0;mu1=1;
cohensd=1.8;
show_figure=0;

% generate data
[W_input, threshold_GT, case_mtx, ctrl_mtx] = generate_data(Da, Da_idx1, Da_idx2, case_num, ctrl_num, mu0, mu1, cohensd,show_figure);

figure; imagesc(W_input);ax=gca;ax.FontSize=18;colorbar; colormap jet;
xlabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');
title("Input W",'FontSize',45,'FontWeight','bold','Color','k');

%% (2) implementation of **SCCN-within** (with SC integrated) on G

% construct infrastructure graph SA (spatial contiguity)
N1=size(Da,1);
indx_1d_A=[1:N1^2];
indx2d=reshape(indx_1d_A,N1,N1); % 30*30
[indxI,indxJ]=ind2sub(size(indx2d),indx_1d_A);
IJ=[indxI',indxJ'];

A_1d=pdist(IJ,'chebychev');  % 900choose2 pairs
A_2d=squareform(A_1d);  %900*900
A_2d_adj=A_2d;              %only adjacent vovels have distance=1
A_2d_adj(A_2d_adj>1.5)=0;    %all other distance=0
A_spa=A_2d_adj;  %900*900 matrix A

figure;imagesc(A_spa);


% initialization of parameters in proposed obj. function (see variable details in
% the function SCCN_within.m
SA=A_spa;
p0=0.008;
sum_cutoff=0;
lambda=1.5;
num_iter=1;
skip_step=20;
show_figure=1;show_progress=1;

[CID_SCCN, W_output, Clist_SCCN]=SCCN_within(W_input, SA, p0, sum_cutoff, lambda, num_iter,skip_step,show_figure,show_progress);


figure;
imagesc(W_output);colormap jet;colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';
title('Output W: SCCN with SC','FontSize',40,'FontWeight','bold','Color','k');
xlabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');



%% (3)implementation of **SCCN-within** (without SC integrated) on G
p0=threshold_GT;
num_iter=1
skip_step=50
show_progress=1;

tic
[CID_SICERS,W_SICERS, Clist_SICERS]=SICERS_skip(W_input,p0,num_iter,skip_step,show_progress);
elapsed_SICERS = toc;

figure;
imagesc(W_SICERS);colormap jet;colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';
title('Output W: SCCN without SC','FontSize',40,'FontWeight','bold','Color','k');
xlabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');
ylabel("Voxels in ROI A",'FontSize',20,'FontWeight','bold','Color','k');

%% (4) Newwork-based permutation test
Wt_orig = [case_mtx;ctrl_mtx];
Wp_orig = W_input;

Clist=Clist_SCCN;
CID=CID_SCCN;
r_vec = [0.1 0.05 0.01 0.005 0.001]; %threshold for p-values


T_orig =  custom_statistic(Wp_orig, Clist, CID, r_vec);

M=50; %number of permutation
T_vec = zeros(M,1);
show_figure=0;show_progress=0;skip_step=50;
%node_density = @(C) sum(sum(C))/size(C,1);
for m=1:M 
    disp(['current iteration m= ',num2str(m)])
    Wt_c=Wt_orig(randperm(case_num+ctrl_num),:);
    Wt1_c=Wt_c(1:case_num,:);
    Wt2_c=Wt_c(case_num+1:end,:);
    [~,ptemp]=ttest2(Wt1_c,Wt2_c);
    
    Wp=squareform(-log(ptemp));
   
    [CID, Wp_output, Clist]=SCCN_within(Wp, SA, p0, sum_cutoff, lambda, num_iter,skip_step,show_figure,show_progress);  

    [~,T_vec(m)] = custom_statistic(Wp, Clist, CID, r_vec);    
end

CID=CID_SCCN; %%% remember to update!
P_value = zeros(length(T_orig),1); %multiple testing p-value
for i = 1:length(T_orig)
    if CID(i) <= 5
        P_value(i) = 1;
    else
        P_value(i) = sum((T_vec - T_orig(i)) >0 )/M;
    end
end

%% (5) plot detected subgraphs that pass the permutation test
CID=CID_SCCN;
Clist=Clist_SCCN;%%% remember to update!


cur_list_vec=[];
up = length(find(P_value<0.05)); 
for i = 1:up
    % load cluster node list
    if i == 1
        cur_list = Clist(1:CID(1));
    else
        cur_list = Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i));
    end
    cur_list_vec=[cur_list_vec cur_list];
end 
 

Da_found=zeros(30,30);Da_found(cur_list_vec)=1;
figure;
imagesc(Da_found);colormap parula;colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';
title('Recovered ROI A','FontSize',45,'FontWeight','bold','Color','k');

%% (6) performance evaluation of selected sub-networks in terms of TPR and FPR 
% we calculate the edge-wise inference of the subgraphs
% the confusion matrix may look like
%                         (i,j) in Gc hat       |  (i,j) not in Gc hat
% ---------------------------------------------------------------------
% (i,j) in Gc        |      TP                  |    FN
% ---------------------------------------------------------------------
% (i,j) not in Gc    |      FP                  |    TN
% ---------------------------------------------------------------------
Da=zeros(30,30);
Da(3:12, 3:11)=1; %13*14=182 nodes
Da(10:25, 18:25)=2; %22*11=242 nodes
Da_idx10=find(Da==1);
Da_idx20=find(Da==2);

true_pos_nodes=[Da_idx10;Da_idx20];
true_neg_nodes=setdiff(1:900,true_pos_nodes);

tpr_vec=[]; fpr_vec=[];
show_figure=0;show_progress=0;  Wt_orig = [case_mtx;ctrl_mtx];skip_step=20;
r_vec = [0.1 0.05 0.01 0.005 0.001]; %threshold for p-values
%figure;imagesc(W);colormap jet;colorbar;ax=gca;ax.FontSize=18;ax.FontWeight='bold';

for d=1:500 %simulate 500 datasets and compute the mean (s.d) of TPR and FPR

    text = ['Current simulated dataset: d=',num2str(d)];
    disp(text)

    [W, ~, ~, ~] = generate_data(Da, Da_idx1, Da_idx2, case_num, ctrl_num, mu0, mu1, cohensd,show_figure);
    
    [CID_SCCN, W_output, Clist_SCCN]=SCCN_within(W_input, SA, p0, sum_cutoff, lambda, num_iter,skip_step,show_figure,show_progress);
    T_orig =  custom_statistic(W, Clist_SCCN, CID_SCCN, r_vec);    
    
%     [CID_SICERS,W_SICERS, Clist_SICERS]=SICERS_skip(W,p0,num_iter,skip_step,show_progress);
%      T_orig =  custom_statistic(W, Clist_SICERS, CID_SICERS, r_vec);    

    CID=CID_SCCN; %%% remember to update
    % need to compute T_vec from step (4)
    P_value = zeros(length(T_orig),1);
    for i = 1:length(T_orig)
        if CID(i) <= 5
            P_value(i) = 1;
        else
            P_value(i) = sum((T_vec - T_orig(i)) >-1 )/M;
        end
    end
    
    % obtain indexes of nodes in Gc hat
    cur_list_vec=[];
    Clist=Clist_SCCN;%%% remember to update
   
    up = length(find(P_value<0.05)); 
    for i = 1:up
        % load cluster node list
        if i == 1
            cur_list = Clist(1:CID(1));
        else
            cur_list = Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i));
        end
        cur_list_vec=[cur_list_vec cur_list];
    end 

    predict_pos_nodes=cur_list_vec;

    % TPR: TP/(TP+FN)
    tp_nodes=intersect(predict_pos_nodes,true_pos_nodes,'stable');
    tpr=length(tp_nodes)/length(true_pos_nodes); 
    tpr_vec=[tpr_vec tpr];
    
    
    %% FPR: FP/(TN+FP)
    fp_nodes=intersect(predict_pos_nodes,true_neg_nodes,'stable');
    fpr=length(fp_nodes)/length(true_neg_nodes);
    fpr_vec=[fpr_vec fpr];
    
     
end 

display(['mean of TPR is ',num2str(mean(tpr_vec))] )
display(['std of TPR is ',num2str(std(tpr_vec))] )

display(['mean of FPR is ',num2str(mean(fpr_vec))] )
display(['std of FPR is ',num2str(std(fpr_vec))] )






