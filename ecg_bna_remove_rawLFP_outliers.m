function allSitesData = ecg_bna_remove_rawLFP_outliers(sitesdir,sitefiles)

allSitesData = [];
extreme_z_threshold = 3;
noise_perc_bin = 0.05;
temp_allSites = struct;
nSites = length(sitefiles);
for s = 1: nSites
    load([sitesdir,filesep,sitefiles(s).name],'sites');
    temp_allSites(s).site = sites;
    temp_allSites(s).name = sitefiles(s).name;
    temp_allSites(s).inconsistant_ch = 0;
    temp_allSites(s).noisy_ch = 0;
    temp_allSites(s).noise_perc_bin = noise_perc_bin;
    temp_allSites(s).extm_zval = extreme_z_threshold;
end
% trials_block = [temp_allSites(s).site.block];
% blocks_with_LFP = unique(trials_block);
lfp_tagets = arrayfun(@(x) x.site.target, temp_allSites, 'UniformOutput', false);
lfp_data = arrayfun(@(x) x.site.LFP, temp_allSites, 'UniformOutput', false);
lfp_sizes = cellfun(@(x) size(x, 2), arrayfun(@(s) s.site.LFP, temp_allSites, 'UniformOutput', false));
common_size = mode(lfp_sizes);
inconsistent_sites = find(lfp_sizes ~= common_size);
lfp_data(inconsistent_sites) = [];
lfp_tagets(inconsistent_sites) = [];
temp_allSites(inconsistent_sites).inconsistant_ch = 1;


concat_raw = vertcat(lfp_data{:});
concat_raw_zscored = zscore(concat_raw, 0, 1);

extreme_counts = sum(abs(concat_raw_zscored) > extreme_z_threshold, 2);
noisy_channels = find(extreme_counts > noise_perc_bin * size(concat_raw, 2)); % More than 5% of timepoints exceed threshold
if ~isempty(noisy_channels)
    temp_allSites(noisy_channels).noisy_ch = 1;
end

target_list = unique(lfp_tagets);
if (sum(contains(target_list,'L'))>0 && sum(contains(target_list,'L'))<length(target_list))...
        || (sum(contains(target_list,'R'))>0 && sum(contains(target_list,'R'))<length(target_list))
    for tr = 1: length(target_list)
        target = target_list{tr};
        site_idx = find((ismember(lfp_tagets, target)));
        concat_site_lfp_mean(tr,:) = squeeze(mean(concat_raw(site_idx,:),1));
        
        site_to_use = find(~([temp_allSites.inconsistant_ch] ==1 ...
            | [temp_allSites.noisy_ch] ==1 )...
            & ismember(arrayfun(@(x) x.site.target, temp_allSites, 'UniformOutput', false), target));
        
        for s = 1: length(site_to_use)
            allSitesData(end+1).site = temp_allSites(site_to_use(s)).site;
            allSitesData(end).site.LFP = allSitesData(s).site.LFP - concat_site_lfp_mean(tr,:);
            allSitesData(end).name = temp_allSites(site_to_use(s)).name;
        end
    end
else
    concat_site_lfp_mean =  squeeze(mean(concat_raw,1));
    
    site_to_use = find(~([temp_allSites.inconsistant_ch] ==1 | [temp_allSites.noisy_ch] ==1 ));
    
    for s = 1: length(site_to_use)
        allSitesData(end+1).site = temp_allSites(site_to_use(s)).site;
        allSitesData(end).site.LFP = allSitesData(s).site.LFP - concat_site_lfp_mean;
        allSitesData(end).name = temp_allSites(site_to_use(s)).name;
    end
    
end

end