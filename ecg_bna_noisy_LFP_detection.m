function [concat_raw, noisy_trials_lfp_mean , noisy_trials_lfp_zscore] = ecg_bna_noisy_LFP_detection(concat_raw)
%% Settings for detection of noisy parts - partially taken from LFP pipeline
% configuration for lfp noise rejection
cfg_noise = [];
% whether or not to apply noise rejection
% Set to 0 to accept all trials
% Set to 1 to run the noise trial detection methods
cfg_noise.detect = 1;
% combination of methods to be used - future use
% currently all methods are used together
cfg_noise.methods = {'amp', 'std', 'diff' , 'zscore'};
% threshold for lfp raw amplitude (number of std deviations from mean)
cfg_noise.amp_thr = 6; %6
% number of consecutive samples beyond threshold to be considered for marking
% a noisy trial
cfg_noise.amp_N = 10;
% no of standard deviations of trial LFP w.r.t LFP std of all trials
cfg_noise.std_thr = 4;
% threshold for lfp derivative (number of std deviations from mean)
cfg_noise.diff_thr = 6;
% number of consecutive samples beyond threshold to be considered for marking
% a noisy trial
cfg_noise.diff_N = 10;
cfg_noise.zscore_thr = 4;
% % threshold for lfp power in number of standard deviations from mean
% cfg_noise.pow_thr = 4;
% % whether single trials should be plotted
% cfg_noise.plottrials = 0;
%%
concat_site_diff_lfp = [nan diff(concat_raw)];

% now get the mean and std of LFP for each site
site_lfp_mean = mean(concat_raw);
site_lfp_std = std(concat_raw);
% get lfp raw threshold
lfp_raw_minbound = site_lfp_mean - cfg_noise.amp_thr * site_lfp_std;
lfp_raw_maxbound = site_lfp_mean + cfg_noise.amp_thr * site_lfp_std;

% concatenate all derivatives to a 1d array
arr_concat_diff_lfp = concat_site_diff_lfp;
% now get the mean and std of LFP for each site
lfp_diff_mean = nanmean(arr_concat_diff_lfp);
lfp_diff_std = nanstd(arr_concat_diff_lfp);
% get lfp derivative threshold
lfp_diff_minbound = lfp_diff_mean - cfg_noise.diff_thr * lfp_diff_std;
lfp_diff_maxbound = lfp_diff_mean + cfg_noise.diff_thr * lfp_diff_std;
% arrays to store the trails marked as noisy by each method
noisy_trials_lfp_mean   = zeros(1,length(concat_raw));
noisy_trials_lfp_std    = zeros(1,length(concat_raw));
noisy_trials_lfp_diff   = zeros(1,length(concat_raw));
noisy_trials_lfp_zscore    = zeros(1,length(concat_raw));
noisy_trials_lfp_pow    = zeros(1,length(concat_raw));

% now find the noisy parts

concat_raw_zscored = zscore(concat_raw);
% a LFP data is noisy if any LFP datapoint in the trial is beyond
% threshold for N consecutive samples
for d = 1:length(concat_raw)
    if concat_raw(d) < lfp_raw_minbound || concat_raw(d) > lfp_raw_maxbound
        noisy_trials_lfp_mean(d)=1;
%         concat_raw(d)= nan;
    end
    if concat_raw_zscored(d)> cfg_noise.zscore_thr
        noisy_trials_lfp_zscore(d)=1;
%         concat_raw(d)= nan;
    end
end
% plot(noisy_trials_lfp_zscore),hold on, plot(noisy_trials_lfp_mean)
x = find(noisy_trials_lfp_mean);
y = find(noisy_trials_lfp_zscore);
if isempty(intersect(x,y))
    z = sort([x,y]);
elseif ~isempty(intersect(x,y))
    z = sort(unique([x,y]));
end

mask = ones(1,length(concat_raw));
mask(z) = 0;
mask = find(mask);
trial_lfp_std = std(concat_raw(mask));

idx_shift = [1 , find(diff(z)>1)+1 , length(z)];
idx_shift = unique(idx_shift);
if ~isempty(z) && length(z)>1
    for ii = 1:length(diff(idx_shift))
        tmp = [idx_shift(ii), idx_shift(ii+1)];
        lng = diff(tmp);
        if (z(tmp(1))-1 > 0)
            subs_mean = mean([concat_raw(z(tmp(1))-1)  , concat_raw(z(tmp(2)-1)+1)]);
            subs_std = std([concat_raw(z(tmp(1))-1)  ,  concat_raw(z(tmp(2)-1)+1)]);
%             concat_raw(z(tmp(1):tmp(2)-1)) = subs_mean.*ones(1,lng);
            %% Linear Interpolation:
            x = [z(tmp(1))-1 , z(tmp(2)-1)+1];
            v = [concat_raw(z(tmp(1))-1)  , concat_raw(z(tmp(2)-1)+1)];
            xq = z(tmp(1))-1 : z(tmp(2)-1)+1 ;
            vq = interp1(x, v, xq);
            
            concat_raw(z(tmp(1):tmp(2)-1)) = vq(2:end-1);
        else
            subs_mean = concat_raw(z(tmp(2)-1)+1);
%             concat_raw(z(tmp(1):tmp(2)-1)) = subs_mean.*ones(1,lng);
            %% Linear Interpolation:
            x = [z(tmp(1)) , z(tmp(2)-1)+1];
            v = [concat_raw(z(tmp(1)))  , concat_raw(z(tmp(2)-1)+1)];
            xq = z(tmp(1)) : z(tmp(2)-1)+1 ;
            vq = interp1(x, v, xq);
            
            concat_raw(xq) = vq;
        end
    end
end
% 


% % a trial is noisy if for more than N consecutive samples,
% % derivative of LFP is greater or less than mean +/- 2*std of
% % LFP derivative of all trials
%
% for d = 1:length(concat_site_diff_lfp)
%     if concat_site_diff_lfp(d) > lfp_diff_maxbound || concat_site_diff_lfp(d) < lfp_diff_minbound
%         noisy_trials_lfp_diff(d)=1;
%     end
% end
%
% % a trial is noisy if the spectral power for 50% of frequency bins at any time bin
% % is greater than mean power + 2*std of all trials at that time bin
% pow_mean_f = nanmean(concat_site_lfp_pow{t}, 3);
% pow_std_f = nanstd(concat_site_lfp_pow{t}, 0, 3);
% for tbin = 1:size(concat_site_lfp_pow{t}, 3)
%     pow_t = concat_site_lfp_pow{t}(:,:,tbin);
%     noisyfbins = sum(pow_t > pow_mean_f + cfg_noise.pow_thr*pow_std_f);
%     if noisyfbins / size(concat_site_lfp_pow{t}, 2) > 0.5
%         site_lfp.trials(t).noisy = 1;
%         noisy_trials_lfp_pow(t) = 1;
%         break;
%     end
% end
end