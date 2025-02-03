function [triggered, significance] = ecg_bna_get_triggered_parameters(site,triggers, width_in_samples,realD,cfg)
significance=struct;
% ecg_bna_get_triggered_split_shuffled - gets the time-frequency
% pow,phs spectrogram and phsBP for a specified time window around all shuffled Rpeak onset for a single site for
% given trials (usually trials belonging to a condition) in a session
%
% USAGE:
%	triggered = ecg_bna_get_triggered_split_shuffled( site_lfp,
%	cond_ecg, state, ecg_bna_cfg )
%
% INPUTS:
%		site_lfp            - struct containing LFP data for all trials of
%		a session from a single site
%       cond_ecg         - ecg data of trials whose time freq spectrogram
%       are to be obtaied
%       state               - a cell array specifying time window during
%       which LFP tfs should be obtained
% OUTPUTS:
%		triggered   - struct containing LFP tfs for all given
%		trials
%
% REQUIRES:	lfp_tfa_baseline_normalization ???
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_Rpeak_evoked_LFP,
% ecg_bna_get_shuffled_Rpeak_evoked_LFP




triggered = struct;
w_samples=width_in_samples(1):width_in_samples(2);

%% get rid of some events to get a matrix (!!!!)
n_shuffles=size(triggers,1);
NNN=triggers==0;
AAA=all(NNN,1);
TT=triggers(:,~AAA);
NNN=TT==0;
nRpeaks=sum(~NNN,2);

%% reducing rpeaks - this can lead to 0 remaining shuffled Rpeaks quite easily...
nRpeaks=sum(~NNN,2);
minvalid=min(nRpeaks);
nrows=size(NNN,2);
todrawfrom=1:nrows;
for sh = 1: n_shuffles
    F=find(NNN(sh,:));
    todraw=todrawfrom;
    todraw(F)=[];
    randindex=rand(numel(todraw),1);
    [~, randindex]=sort(randindex);
    ndrawn=nrows-numel(F)-minvalid;
    toremove(sh,:)=   [F todraw(randindex(1:ndrawn))]+(sh-1)*nrows;
end

trigs=TT';
trigs(toremove)=[];
%trigs=reshape(trigs,[n_shuffles,minvalid]); %% check dimensions
trigs=reshape(trigs,[minvalid,n_shuffles]); %% check dimensions
trigs=trigs';

% %% alternative: remove indexes where more than half are invalid
% 
% trigs=zeros(n_shuffles,max(nRpeaks));
% 
% [~,sortix]=sort(-nRpeaks);
% for sh = sortix'
%     trigs(sh,1:nRpeaks(sh))=TT(sh,~TT(sh,:)==0);
% end
%     
% toremove=sum(trigs==0,1)>size(trigs,1)/2;
% TT(:,toremove)=[];
% for tr = 1:size(TT,2)
%     F=TT(:,tr)==0;
%     F
%     
%     todraw=todrawfrom;
%     todraw(F)=[];
%     randindex=rand(numel(todraw),1);
%     [~, randindex]=sort(randindex);
%     ndrawn=nrows-numel(F)-minvalid;
%     toremove(sh,:)=   [F todraw(randindex(1:ndrawn))]+(sh-1)*nrows;
% end
% 
% %%%

tfs=[site.tfs];
%% check pre-allocation time saving (?)
% 
    FN_in={'pow','powbp','lfp','pha','pha','phabp'};
    FN={'pow','powbp','lfp','pha','itpc','itpcbp'};

    for f=1:numel(FN_in)  
        fn=FN{f};      
        fni=FN_in{f};
        S2=size(tfs.(fni),1);
        S3=numel(w_samples);
        
        triggered.(fn).mean        = zeros(1,S2,S3);
        triggered.(fn).std         = zeros(1,S2,S3);
        triggered.(fn).conf95      = zeros(2,S2,S3);
    end

%% NOW: loop through window samples to compute shuffled parameters all at once!

for s = 1: numel(w_samples)
    w=w_samples(s);
    t=trigs+w;
    triggered.ntriggers=minvalid;
    
    for f=1:numel(FN_in)
        fn=FN{f};
        fni=FN_in{f};
        
        AA=      reshape(tfs.(fni)(:,t),size(tfs.(fni),1),size(t,1),size(t,2));
        if ismember(fn,{'itpc','itpcbp'})
            BB=abs(mean(AA,3));
        elseif ismember(fn,{'pha'})
            BB=angle(mean(AA,3));
        else
            BB=mean(AA,3);
        end
        
        triggered.(fn).mean(1,:,s)        = mean(BB,2);
        triggered.(fn).std(1,:,s)         = std(BB,0,2);
        
        %% check dimensions
        triggered.(fn).conf95(1,:,s)      = prctile(BB,97.5,2);
        triggered.(fn).conf95(2,:,s)      = prctile(BB,2.5,2);
        %triggered.(fn).complete(:,:,s)    = BB';
        complete.(fn)(:,:,s)    = BB';
    end
end

if nargin>3 %% surrogate data,add signficance test here !
    pthreshold=0.05;
    clusterthreshold=0.20;
    for f=1:numel(FN_in)
        fn=FN{f};
        dat=complete.(fn);
        s_ix=1:size(dat,1);
        
        
        %tThreshold = abs(tinv(pthreshold, numel(s_ix)-1));   
        
        tThreshold = (1-pthreshold)/2;   
        % find (for each surrogate) the max cluster significance (somehow i
        % feel this should be done PER frequency, hmmmm....)
        for s=s_ix
            dat_s=squeeze(dat(s,:,:));
            est_mean=squeeze(mean(dat(s_ix~=s,:,:),1));
            est_var=squeeze(var(dat(s_ix~=s,:,:),1));
            %T=(dat_s-est_mean)./sqrt(est_var); %% we actually care about pos/neg here
            
            T=squeeze(sum(dat(s_ix~=s,:,:)<repmat(dat(s,:,:),[size(dat,1)-1,1,1]),1));
            T=T/(size(dat,1)-1)-0.5;
            
            sigpos=T>tThreshold;
            signeg=T<-tThreshold;
            
            CCpos = bwconncomp(sigpos);
            CCneg = bwconncomp(signeg);
            sumTpos = cellfun(@(x) sum(T(x)),CCpos.PixelIdxList);
            sumTneg = cellfun(@(x) sum(T(x)),CCneg.PixelIdxList);
            T_max(s) = max([sumTpos 0]);
            T_min(s) = min([sumTneg 0]);
        end
        
        %% now find significant clusters in the real data        
        T_max_thr=prctile(T_max,100*(1-clusterthreshold));   % sum of     
        T_min_thr=prctile(T_min,100*clusterthreshold);    % sum of    
        est_mean=squeeze(mean(dat,1));
        est_var=squeeze(var(dat,1));
        
        %T=(squeeze(realD.(fn).mean)-est_mean)./sqrt(est_var); 
        
        T=squeeze(sum(dat<repmat(realD.(fn).mean,[size(dat,1),1,1]),1));
        T=T/(size(dat,1)-1)-0.5;
        
        sigpos=T>tThreshold;        
        signeg=T<-tThreshold;
        
        CCpos = bwconncomp(sigpos);
        CCneg = bwconncomp(signeg);
        sigCpos = cellfun(@(x) sum(T(x))>=T_max_thr,CCpos.PixelIdxList);
        sigCneg = cellfun(@(x) sum(T(x))<=T_min_thr,CCneg.PixelIdxList);
        sigpixpos=vertcat(CCpos.PixelIdxList{sigCpos});
        sigpixneg=vertcat(CCneg.PixelIdxList{sigCneg});
        significance.(fn)=zeros(size(realD.(fn).mean));
        significance.(fn)(sigpixpos)=true;
        significance.(fn)(sigpixneg)=true;
    end
   
    if false
        for f=1:numel(FN_in)
            fn=FN{f};
            
            null_mean=squeeze(triggered.(fn).mean);
            null_std=squeeze(triggered.(fn).std);
            actual_tf=squeeze(realD.(fn).mean);
            alpha=0.05;
            cluster_thresh=0.05;
            
            [~, ~, ~, significance_map] = tfClusterPermTestWithMeanStd(actual_tf, null_mean, null_std, alpha, cluster_thresh);
            
            significance.(fn)=zeros(size(realD.(fn).mean));
            significance.(fn)(1:numel(significance_map))=significance_map;
        end
    end
end


end

function [clusters, cluster_stats, p_values, significance_map] = tfClusterPermTestWithMeanStd(actual_tf, null_mean, null_std, alpha, cluster_thresh)
% Time-frequency cluster-based permutation test using mean and std of permuted data
%
% The function is renamed to tfClusterPermTestWithMeanStd to reflect the new input parameters.
% Instead of taking the full 3D permuted data matrix, the function now takes the 2D mean and standard deviation of the permuted data as input.
% The null distribution of cluster-level statistics is generated using a parametric approach, where random t-maps are created based on the provided null mean and standard deviation.
% The p-values are computed by comparing the actual cluster statistics to the null distribution.
%
% Inputs:
% actual_tf: 2D matrix [time x frequency] of actual signal
% null_mean: 2D matrix [time x frequency] of mean of permuted data
% null_std: 2D matrix [time x frequency] of standard deviation of permuted data
% alpha: significance level (e.g., 0.05)
% cluster_thresh: threshold for initial clustering (e.g., 0.05)
%
% Outputs:
% clusters: cell array of cluster indices (combining positive and negative clusters)
% cluster_stats: statistics for each cluster
% p_values: p-values for each cluster
% significance_map: binary map of significant clusters

% Get dimensions
[n_time, n_freq] = size(actual_tf);

% Compute t-statistics at each point
t_map = (actual_tf - null_mean) ./ (null_std + eps);

% Threshold the t-statistics
pos_thresh_map = t_map > tinv(1-cluster_thresh/2, inf);
neg_thresh_map = t_map < -tinv(1-cluster_thresh/2, inf);

% Find clusters in thresholded maps
pos_connected_clusters = bwconncomp(pos_thresh_map, 8);  % 8-connected neighborhoods
neg_connected_clusters = bwconncomp(neg_thresh_map, 8);

% Initialize outputs
n_pos_clusters = pos_connected_clusters.NumObjects;
n_neg_clusters = neg_connected_clusters.NumObjects;
pos_clusters = pos_connected_clusters.PixelIdxList;
neg_clusters = neg_connected_clusters.PixelIdxList;

% Combine clusters
n_total_clusters = n_pos_clusters + n_neg_clusters;
clusters = cell(1, n_total_clusters);
cluster_stats = zeros(1, n_total_clusters);
p_values = ones(1, n_total_clusters);

% Copy positive clusters and their statistics
for i = 1:n_pos_clusters
    clusters{i} = pos_clusters{i};
    cluster_stats(i) = sum(t_map(pos_clusters{i}));
end

% Copy negative clusters and their statistics
for i = 1:n_neg_clusters
    clusters{i + n_pos_clusters} = neg_clusters{i};
    cluster_stats(i + n_pos_clusters) = sum(t_map(neg_clusters{i}));
end

significance_map = false(n_time, n_freq);

if n_total_clusters == 0
    return;
end

% Build null distribution of cluster statistics using parametric approach
num_perms = 10000;
% Generate random maps for all permutations at once
perm_maps = randn(n_time, n_freq, num_perms) .* repmat(null_std,[1,1,num_perms]) + repmat(null_mean,[1,1,num_perms]);
pos_null_cluster_stats = zeros(1, num_perms);
neg_null_cluster_stats = zeros(1, num_perms);

% Process each permutation
for i = 1:num_perms
    curr_map = perm_maps(:,:,i);
    
    % Find positive clusters
    pos_perm_thresh_map = curr_map > tinv(1-cluster_thresh/2, inf);
    pos_perm_clusters = bwconncomp(pos_perm_thresh_map, 8);
    
    if pos_perm_clusters.NumObjects > 0
        pos_perm_cluster_stats = zeros(1, pos_perm_clusters.NumObjects);
        for c = 1:pos_perm_clusters.NumObjects
            pos_perm_cluster_stats(c) = sum(curr_map(pos_perm_clusters.PixelIdxList{c}));
        end
        pos_null_cluster_stats(i) = max(pos_perm_cluster_stats);
    end
    
    % Find negative clusters
    neg_perm_thresh_map = curr_map < -tinv(1-cluster_thresh/2, inf);
    neg_perm_clusters = bwconncomp(neg_perm_thresh_map, 8);
    
    if neg_perm_clusters.NumObjects > 0
        neg_perm_cluster_stats = zeros(1, neg_perm_clusters.NumObjects);
        for c = 1:neg_perm_clusters.NumObjects
            neg_perm_cluster_stats(c) = sum(curr_map(neg_perm_clusters.PixelIdxList{c}));
        end
        neg_null_cluster_stats(i) = min(neg_perm_cluster_stats);
    end
end

% Compute p-values
for i = 1:n_pos_clusters
    p_values(i) = mean(pos_null_cluster_stats >= cluster_stats(i));
end
for i = 1:n_neg_clusters
    p_values(i + n_pos_clusters) = mean(neg_null_cluster_stats <= cluster_stats(i + n_pos_clusters));
end

% Create significance map
for i = 1:n_total_clusters
    if p_values(i) < alpha
        significance_map(clusters{i}) = true;
    end
end
end

