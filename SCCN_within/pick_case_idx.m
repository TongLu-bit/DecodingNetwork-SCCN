function [Da_idx1,Da_idx2,Da_0]  = pick_case_idx(Da,rate_in, rate_out,show_figure)
%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%
%%% Da: a matrix describe the spatial location of covariate-related subnetworks
%%% rate_in: proportion of noise added outside of covariate-related subnetworks
%%% rate_out: proportion of anti-significant nodes excluded from covariate-related subnetworks
%%% show_figure: =1 if one wants to visualize the simulated spatial location with noise 

%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%
%%% Da_idx1: node index of dense subnetwork1
%%% Da_idx2: node index of dense subnetwork2
%%% Da_0: Da after noise added

%%% Tong Lu, May 2023

%% Constract graph with noise 
    full_length=size(Da,1);
    
    % Add noise around dense sub-areas 
    rowA1=randi([1 full_length],1,rate_in*full_length); %generate a 1*50 vector, containing values from discrete Uniform(1,30)
    colA1=randi([1 full_length],1,rate_in*full_length);
    for i=1:rate_in*full_length
        Da(rowA1(i),colA1(i))=1;
    end
    
    rowA2=randi([1 full_length],1,rate_in*full_length);
    colA2=randi([1 full_length],1,rate_in*full_length);
    for i=1:rate_in*full_length
        Da(rowA2(i),colA2(i))=2;
    end
    
    
    % De-highlight significant points from dense sub-areas 
    rowA01=randi([1 full_length],1,rate_out*full_length);
    colA01=randi([1 full_length],1,rate_out*full_length);
    for i=1:rate_out*full_length
        Da(rowA01(i),colA01(i))=0;
    end
    
    rowA02=randi([1 full_length],1,rate_out*full_length);
    colA02=randi([1 full_length],1,rate_out*full_length);
    for i=1:rate_out*full_length
        Da(rowA02(i),colA02(i))=0;
    end


%% Get the index of covariate-related nodes for cases
    Da_idx1=find(Da==1);
    Da_idx2=find(Da==2);

    Da_0=Da;
    Da_0(Da_idx2)=1;
    if show_figure==1
        figure; imagesc(Da_0);ax=gca;ax.FontSize=18;ax.FontWeight='bold'; colormap summer;
        title("ROI A, |R|=900 nodes",'FontSize',45,'FontWeight','bold','Color','k');
    end 
end