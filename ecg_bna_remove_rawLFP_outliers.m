function allSitesData = ecg_bna_remove_rawLFP_outliers(sitesdir,sitefiles,cfg,ts_original)
ts=round(cfg.lfp.timestep/ts_original);
allSitesData = [];
% extreme_z_threshold = 3;
% noise_perc_bin = 0.05;
temp_allSites = struct;
nSites = length(sitefiles);
for s = 1: nSites
    load([sitesdir,filesep,sitefiles(s).name],'sites');
    temp_allSites(s).site = sites;
    temp_allSites(s).name = sitefiles(s).name;
    temp_allSites(s).inconsistant_ch = 0;
    temp_allSites(s).noisy_site = 0;
    %     temp_allSites(s).noise_perc_bin = noise_perc_bin;
    %     temp_allSites(s).extm_zval = extreme_z_threshold;
end

targets = {'VPL_R', 'VPL_L', 'dPul_R', 'dPul_L','MD_L','MD_R'};
allSites_targets = arrayfun(@(x) (x.site.target), temp_allSites, 'UniformOutput', false);
inval_targets = ~ismember(allSites_targets,targets);
temp_allSites(inval_targets) = [];

trials_block = arrayfun(@(x) unique(x.site.block), temp_allSites, 'UniformOutput', false);
blocks_with_LFP = unique([trials_block{:}]);

%sampling_table = table;
dummy = 1;
for s = 1: length(temp_allSites)
    samples_past=0;
    % you are not resampling here :P
    %samples_past_resampled=0; 
    
    for b = 1:length(blocks_with_LFP)
        B = blocks_with_LFP(b);
        if sum(ismember(unique(temp_allSites(s).site.block),B))>0
            sampling_table.site(dummy)  = s;
            sampling_table.block(dummy) = B;
            
            b_samples=[temp_allSites(s).site.block ==B];
            %bs=samples_past_resampled+1;
            %be=floor(sum(temp_allSites(s).site.LFP_samples(b_samples))/ts)+samples_past_resampled;
            bs_original=samples_past+1;
            be_original=sum(temp_allSites(s).site.LFP_samples(b_samples))+samples_past;
            
            sampling_table.start_sample(dummy)  = bs_original;
            sampling_table.end_sample(dummy)    = be_original;
            
            samples_past=be_original;
            %samples_past_resampled=be;
            dummy = dummy+1;
            %     tmp_block_idx(b,:) = arrayfun(@(x) ismember(unique(x.site.block),blocks_with_LFP(b)), temp_allSites, 'UniformOutput', false);
        end
    end
end


lfp_data = arrayfun(@(x) x.site.LFP, temp_allSites, 'UniformOutput', false);
lfp_sizes = cellfun(@(x) size(x, 2), arrayfun(@(s) s.site.LFP, temp_allSites, 'UniformOutput', false));
% common_size = mode(lfp_sizes);
%
% inconsistent_sites = find(lfp_sizes ~= common_size);
% lfp_data(inconsistent_sites) = [];
% lfp_tagets(inconsistent_sites) = [];
% if ~isempty(inconsistent_sites)
%     [temp_allSites(inconsistent_sites).inconsistant_ch] = deal(1);
% end

% idea was not to check across all sites
all_LFP_concat = horzcat(lfp_data{:});
all_sites_25perc = prctile(abs(all_LFP_concat),25);
all_sites_75perc = prctile(abs(all_LFP_concat),75);

% concat_raw_zscored = zscore(concat_raw, 0, 1);
% extreme_counts = sum(abs(concat_raw_zscored) > extreme_z_threshold, 2);
% noisy_channels = find(extreme_counts > noise_perc_bin * size(concat_raw, 2)); % More than 5% of timepoints exceed threshold
noisy_site = [];
for s = 1:length(temp_allSites)
    
% but rather across all sites other than the one you are currently checking 
site_idx=1:numel(lfp_data);
all_LFP_concat = horzcat(lfp_data{site_idx~=s});
all_sites_25perc = prctile(abs(all_LFP_concat),25);
all_sites_75perc = prctile(abs(all_LFP_concat),75);
    
    tmp_site_25perc = prctile(abs(lfp_data{s}),25);
    tmp_site_75perc = prctile(abs(lfp_data{s}),75);
    %yline(tmp_site_25perc,'b--','LineWidth',3)
    if (tmp_site_75perc < all_sites_25perc)||(tmp_site_25perc > all_sites_75perc)
        noisy_site = [noisy_site, s];
    end
end

if ~isempty(noisy_site)
    [temp_allSites(noisy_site).noisy_site] = deal(1);
end
sites_to_remove=[temp_allSites.noisy_site]==1;
allSitesData = temp_allSites(~sites_to_remove);
allSitesData = rmfield(allSitesData,{'inconsistant_ch','noisy_site'});



%% sites that are removed dont need to be fixed, correct?
%% if so, remove them from sampling table immediately
%% alternatively, you can do noise detection BEFORE creating the sampling table, so that you do not run into this problem in the first place
to_remove_from_sampling_table=ismember(sampling_table.site,find(sites_to_remove)); % please check this line, I'm not sure it does what it supposed to
sampling_table.site(to_remove_from_sampling_table)  = [];
sampling_table.block(to_remove_from_sampling_table) = [];
sampling_table.start_sample(to_remove_from_sampling_table)  = [];
sampling_table.end_sample(to_remove_from_sampling_table)    = [];

% clear temp_allSites, because it's huge and should not be needed any more ?
clear temp_allSites

% trials_block = arrayfun(@(x) unique(x.site.block), temp_allSites, 'UniformOutput', false);
% blocks_with_LFP = unique([trials_block{:}]);
lfp_tagets = arrayfun(@(x) x.site.target, allSitesData, 'UniformOutput', false);
site_names = arrayfun(@(x) x.name, allSitesData, 'UniformOutput', false); %% are these site names?
target_list = unique(lfp_tagets);

%% okay here the point was to average across all sites in one hemisphere, not in one target only..

% i have never used contains, bu ti guess i would do something like this: 
% sites_hemisphere is supposed to be 1 for left hemispheres and 2 for
% right hemispheres
sites_hemisphere=contains(lfp_tagets,'_L')+2*contains(lfp_tagets,'_R'); %look for the _ as well, because f.e. VPL contains L 
hemisphere_list = unique(sites_hemisphere);

% if (sum(contains(target_list,'L'))>0 && sum(contains(target_list,'L'))<length(target_list))...
%         || (sum(contains(target_list,'R'))>0 && sum(contains(target_list,'R'))<length(target_list))
    
    for h = 1: length(hemisphere_list)
        %target = target_list{tr};
        hemisphere_sites = find(sites_hemisphere==hemisphere_list(h));
        
        
        for b=1:numel(blocks_with_LFP)
            B=blocks_with_LFP(b);
            %% thing is i dont trust this as it is im afraid (sampling table and temp_allSites have very different indexing i believe?)
            %% -> hence the removing noise sites from sampling table beforehand
            
            %% i tend to create the indexes i will need inside this loop at the start so it doesnt get confusing
            block_idx=ismember(sampling_table.block,B) & ismember(sampling_table.site,hemisphere_sites);
            site_with_blockB_idx = sampling_table.site(block_idx);
            %site_with_blockB_idx = site_with_blockB_idx([temp_allSites(site_with_blockB_idx).noisy_site]==0);
            block_lfp_concat = [];
            for s = 1: length(site_with_blockB_idx)
                %% i believe its easier to track like this (just a matter of taste i guess)
                idx=block_idx & sampling_table.site==site_with_blockB_idx(s);
                start_sample = sampling_table.start_sample(idx);
                end_sample   = sampling_table.end_sample(idx);
                block_lfp_concat = [block_lfp_concat; allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample)];
            end
            block_lfp_concat_mean = squeeze(mean(block_lfp_concat,1));
            
            for s = 1: length(site_with_blockB_idx)
                %% i believe its easier to track like this (just a matter of taste i guess)
                idx=block_idx & sampling_table.site==site_with_blockB_idx(s);
                start_sample = sampling_table.start_sample(idx);
                end_sample   = sampling_table.end_sample(idx);
                allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample) = allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample) - block_lfp_concat_mean;
            end
        end
    end
    
% else
%     for b=1:numel(blocks_with_LFP)
%         B=blocks_with_LFP(b);
%         site_with_blockB_idx = sampling_table.site(ismember(sampling_table.block,B));
%         site_with_blockB_idx = site_with_blockB_idx([temp_allSites(site_with_blockB_idx).noisy_site]==0);
%         block_lfp_concat = [];
%         for s = 1: length(site_with_blockB_idx)
%             start_sample = sampling_table.start_sample(ismember(sampling_table.block,B)& sampling_table.site==site_with_blockB_idx(s));
%             end_sample   = sampling_table.end_sample(ismember(sampling_table.block,B)& sampling_table.site==site_with_blockB_idx(s));
%             block_lfp_concat = [block_lfp_concat; temp_allSites(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample)];
%         end
%         block_lfp_concat_mean = squeeze(mean(block_lfp_concat,1));
%         
%         for s = 1: length(site_with_blockB_idx)
%             start_sample = sampling_table.start_sample(ismember(sampling_table.block,B)& sampling_table.site==site_with_blockB_idx(s));
%             end_sample   = sampling_table.end_sample(ismember(sampling_table.block,B)& sampling_table.site==site_with_blockB_idx(s));
%             allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample) = allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample) - block_lfp_concat_mean;
%         end
%     end
%     
% end


% if (sum(contains(target_list,'L'))>0 && sum(contains(target_list,'L'))<length(target_list))...
%         || (sum(contains(target_list,'R'))>0 && sum(contains(target_list,'R'))<length(target_list))
%
%     concat_site_lfp_mean = nan(length(target_list),size(concat_raw,2));
%
%     for tr = 1: length(target_list)
%         target = target_list{tr};
%         site_idx = (ismember(lfp_tagets, target));
%         concat_site_lfp_mean(tr,:) = squeeze(mean(concat_raw(site_idx,:),1));
%
%         site_to_use = find(~([temp_allSites.inconsistant_ch] ==1 ...
%             | [temp_allSites.noisy_site] ==1 )...
%             & ismember(arrayfun(@(x) x.site.target, temp_allSites, 'UniformOutput', false), target));
%
%         for s = 1: length(site_to_use)
%             allSitesData(end+1).site = temp_allSites(site_to_use(s)).site;
%             allSitesData(end).site.LFP = allSitesData(s).site.LFP - concat_site_lfp_mean(tr,:);
%             allSitesData(end).name = temp_allSites(site_to_use(s)).name;
%         end
%     end
% else
%     concat_site_lfp_mean =  squeeze(mean(concat_raw,1));
%
%     site_to_use = find(~([temp_allSites.inconsistant_ch] ==1 | [temp_allSites.noisy_site] ==1 ));
%
%     for s = 1: length(site_to_use)
%         allSitesData(end+1).site = temp_allSites(site_to_use(s)).site;
%         allSitesData(end).site.LFP = allSitesData(s).site.LFP - concat_site_lfp_mean;
%         allSitesData(end).name = temp_allSites(site_to_use(s)).name;
%     end
%
% end
% end
