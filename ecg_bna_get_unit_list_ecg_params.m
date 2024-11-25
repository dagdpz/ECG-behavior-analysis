function ecg_bna_get_unit_list_ecg_params(cfg)

basepath_to_save=[cfg.unit_lists '_ECG'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

list_of_lists = dir([cfg.unit_lists filesep '*.mat']);

% list_of_lists_600             = dir([cfg.unit_lists filesep '*_600.mat']); % find files after exclusion
% list_of_lists_noCB            = dir([cfg.unit_lists filesep '*noCB*.mat']);
% list_of_lists_withCB          = dir([cfg.unit_lists filesep '*withCB.mat']);
% list_of_lists_CB_corr         = dir([cfg.unit_lists filesep '*noCB_corr.mat']);
% list_of_lists_CB_corr_ccs     = dir([cfg.unit_lists filesep '*noCB_corr_ccs.mat']);
% list_of_lists_CB_excl         = dir([cfg.unit_lists filesep '*_excluded*.mat']);
% list_of_lists_before          = dir([cfg.unit_lists filesep '*_before_exclusion.mat']);
% list_of_lists_highAmp         = dir([cfg.unit_lists filesep '*_high_amplitude.mat']);
% list_of_lists_lowAmp          = dir([cfg.unit_lists filesep '*_low_amplitude.mat']);
% list_of_lists_lowAmp_ccs_any  = dir([cfg.unit_lists filesep '*_low_amplitude_ccs_any.mat']);
% list_of_lists_lowAmp_ccs_both = dir([cfg.unit_lists filesep '*_low_amplitude_ccs_both.mat']);

% lists2drop = ismember({list_of_lists.name}, {list_of_lists_600.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_noCB.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_withCB.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_CB_corr.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_CB_corr_ccs.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_CB_excl.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_before.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_highAmp.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_lowAmp.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_lowAmp_ccs_any.name}) | ...
%     ismember({list_of_lists.name}, {list_of_lists_lowAmp_ccs_both.name});
% 
% list_of_lists = list_of_lists(~lists2drop);

for listNum = 1:length(list_of_lists)
    
    filename = [list_of_lists(listNum).folder filesep list_of_lists(listNum).name];
    
    % stable for at least one task block OR one rest block
    load(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres')
    
    % figure out number of R-peaks
    dt = ecg_bna_load_variables(cfg,unit_ids, 'per_unit_-0.25-0.25s', 'Output', {'NrEvents'}, 0);
    
    enoughRpeaks = zeros(length(cfg.condition), length(unit_ids));
    for conNum = 1:length(cfg.condition)
        
        L = cfg.condition(conNum).name;
        
        enoughRpeaks(conNum,:) = dt.(L).NrEvents > cfg.spk.unit_exclusion.nCardiacCycles;
        
    end
    
    if strfind(list_of_lists(listNum).name, 'stable') & isempty(strfind(list_of_lists(listNum).name, 'excluded'))
        ids_enough_Rpeaks = all(enoughRpeaks);
        ids_both = 3*ones(length(unit_ids),1);
    elseif strfind(list_of_lists(listNum).name, 'selected') & isempty(strfind(list_of_lists(listNum).name, 'excluded'))
        ids_enough_Rpeaks = any(enoughRpeaks);
        load(filename, 'ids_both')
    else
        continue % don't create ECG-related dataset for this unit list
    end
    
    clear dt
    
    % update ids_both according to the sufficiency of R-peaks collected - I
    % have to update cases where there were not enough R-peaks
    ids_both_tmp = zeros(size(enoughRpeaks));
    ids_both_tmp(1,enoughRpeaks(1,:)==1) = 1; % for rest
    ids_both_tmp(2,enoughRpeaks(2,:)==1) = 2; % for task
    
    ids_both_ecg_based_sum = sum(ids_both_tmp); % come up with the same condition encoding
    
    % construct the resulting ids_both that considers also the number of
    % R-peaks
    ids_both_ecg_based = zeros(1,size(enoughRpeaks,2));
    % 1. If ids_both == ids_both_tmp, we're fine, take indices over
    ids_both_ecg_based(eq(ids_both,ids_both_ecg_based_sum')) = ids_both_ecg_based_sum(eq(ids_both,ids_both_ecg_based_sum'));
    % 2. If spike-based is 3 AND ecg-based is lower, take only one
    % condition
    ids_to_take = ids_both == 3 & ids_both_ecg_based_sum' < 3;
    ids_both_ecg_based(ids_to_take) = ids_both_ecg_based_sum(ids_to_take);
    % There are no other suitable conditions
    ids_both = ids_both_ecg_based;
    
    clear enoughRpeaks
    
    dt = ecg_bna_load_variables(cfg, unit_ids, 'cardioballistic', 'data', {'distance2thr', 'AMP_MI', 'cc_PSTH_feature', 'pp_PSTH_feature'}, 0);
    enough_bins = zeros(length(cfg.condition), length(unit_ids));
    for conNum = 1:length(cfg.condition)
        L = cfg.condition(conNum).name;
        
        enough_bins(conNum,:) = dt.(L).distance2thr > 1;
        
    end
%     ids_enough_bins = all(enough_bins);
%     sum(ids_enough_bins)
    
    % Considerations of exclusion of units by fit significance and the
    % goodness-of-fit parameters:
    % 1. There are fits with p < 0.05 that I consider significant: some of
    % them are true positive, some of them are not.
    % 2. I want to take (mostly) true positives and drop most false
    % positives using a threshold for R-squared
    % 3. I will consider significant those fits that have p-value <
    % 0.05 and explain considerable variance of the data (R-sq. > 0.3)
    
%     [~, h_AMP_MI] = bonf_holm([out_rest.AMP_MI(:,2), out_task.AMP_MI(:,2)]); % p < 0.05 corresponds to significant cosine fit implying a significant cardioballistic effect
    
    h_AMP_MI = zeros(length(cfg.condition), length(unit_ids));
    high_Rsq = zeros(length(cfg.condition), length(unit_ids));
    AMP_pp   = zeros(length(cfg.condition), length(unit_ids));
    AMP_cc   = zeros(length(cfg.condition), length(unit_ids));
    for conNum = 1:length(cfg.condition)
        
        L = cfg.condition(conNum).name;
        
        h_AMP_MI(conNum,:) = dt.(L).AMP_MI(:,2) < 0.01;
        high_Rsq(conNum,:) = dt.(L).AMP_MI(:,4) > 0.3; % to drop false positives
        AMP_pp(conNum,:)   = dt.(L).pp_PSTH_feature(1,:); % take the first feature - AMP
        AMP_cc(conNum,:)   = dt.(L).cc_PSTH_feature(1,:);
        
    end
    
    % p < 0.05 - significant correlation between phase PSTH and phase
    % dynamic of AMP feature
    AMP_pp_neg = AMP_pp;
    AMP_pp_neg(AMP_cc >= 0) = NaN;
    [~, h_AMP_cc_neg] = bonf_holm(AMP_pp_neg);
    
    % p < 0.05 - significant MI according to Mosher's procedure
    
    % choose units with positive and significant correlation coefficients
    AMP_pp_pos              = AMP_pp;
    AMP_pp_pos(AMP_cc <= 0) = NaN;
    h_AMP_cc_pos            = AMP_pp_pos < 0.05;
%     [~, h_AMP_cc_pos] = bonf_holm(AMP_pp_pos);
    
%     with_CBE = any(h_AMP_MI' & high_Rsq');
    with_CBE = any(h_AMP_MI);
    
%     no_cardioballistic_effect = ~any(h_AMP_MI' & high_Rsq'); % if a unit is affected in at least one condition we'll get rid of it
    no_cardioballistic_effect = ~any(h_AMP_MI);
    
    no_cardioballistic_effect_corrected = all(~h_AMP_MI | enough_bins);
    
    no_cardioballistic_effect_corrected_ccs = all(~h_AMP_MI | enough_bins | (h_AMP_cc_neg & AMP_cc < 0));
    
    high_amplitude_only = all(enough_bins);
    
    low_amplitude_only  = ~any(enough_bins);
    
    low_amplitude_ccs_any  = ~any(enough_bins) & any(h_AMP_cc_pos);
    
    low_amplitude_ccs_both = ~any(enough_bins) & all(h_AMP_cc_pos);
    
    with_CBE_ccs_any  = any(h_AMP_cc_pos);
    
    with_CBE_ccs_both = all(h_AMP_cc_pos);
    
    % find unit ids and targets with enough R-peak counts
    unit_ids_600    = unit_ids(ids_enough_Rpeaks);
    targets_600     = targets(ids_enough_Rpeaks);
    sites_600       = sites(ids_enough_Rpeaks);
    grid_x_600      = grid_x(ids_enough_Rpeaks);
    grid_y_600      = grid_y(ids_enough_Rpeaks);
    depths_600      = depths(ids_enough_Rpeaks);
    hemispheres_600 = hemispheres(ids_enough_Rpeaks);
    ids_both_600    = ids_both(ids_enough_Rpeaks);
    % without cardioballistic effect
    unit_ids_noCB    = unit_ids(ids_enough_Rpeaks & no_cardioballistic_effect);
    targets_noCB     = targets(ids_enough_Rpeaks & no_cardioballistic_effect);
    sites_noCB       = sites(ids_enough_Rpeaks & no_cardioballistic_effect);
    grid_x_noCB      = grid_x(ids_enough_Rpeaks & no_cardioballistic_effect);
    grid_y_noCB      = grid_y(ids_enough_Rpeaks & no_cardioballistic_effect);
    depths_noCB      = depths(ids_enough_Rpeaks & no_cardioballistic_effect);
    hemispheres_noCB = hemispheres(ids_enough_Rpeaks & no_cardioballistic_effect);
    ids_both_noCB    = ids_both(ids_enough_Rpeaks & no_cardioballistic_effect);
    % with cardioballistic effect but still having enough R-peaks for the
    % analysis
    unit_ids_withCB    = unit_ids(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    targets_withCB     = targets(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    sites_withCB       = sites(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    grid_x_withCB      = grid_x(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    grid_y_withCB      = grid_y(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    depths_withCB      = depths(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    hemispheres_withCB = hemispheres(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    ids_both_withCB    = ids_both(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    % without cardioballistic effect and high spike amplitude
    unit_ids_noCB_corr    = unit_ids(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    targets_noCB_corr     = targets(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    sites_noCB_corr       = sites(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    grid_x_noCB_corr      = grid_x(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    grid_y_noCB_corr      = grid_y(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    depths_noCB_corr      = depths(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    hemispheres_noCB_corr = hemispheres(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    ids_both_noCB_corr    = ids_both(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    % without cardioballistic effect, high spike amplitude, and with
    % non-significant correlation coefficient between phase PSTH and AMP 
    % phase dynamic
    unit_ids_noCB_corr_ccs    = unit_ids(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    targets_noCB_corr_ccs     = targets(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    sites_noCB_corr_ccs       = sites(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    grid_x_noCB_corr_ccs      = grid_x(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    grid_y_noCB_corr_ccs      = grid_y(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    depths_noCB_corr_ccs      = depths(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    hemispheres_noCB_corr_ccs = hemispheres(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    ids_both_noCB_corr_ccs    = ids_both(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    % don't use CBE analysis, choose units with high spike amplitude
    unit_ids_high_amplitude    = unit_ids(ids_enough_Rpeaks & high_amplitude_only);
    targets_high_amplitude     = targets(ids_enough_Rpeaks & high_amplitude_only);
    sites_high_amplitude       = sites(ids_enough_Rpeaks & high_amplitude_only);
    grid_x_high_amplitude      = grid_x(ids_enough_Rpeaks & high_amplitude_only);
    grid_y_high_amplitude      = grid_y(ids_enough_Rpeaks & high_amplitude_only);
    depths_high_amplitude      = depths(ids_enough_Rpeaks & high_amplitude_only);
    hemispheres_high_amplitude = hemispheres(ids_enough_Rpeaks & high_amplitude_only);
    ids_both_high_amplitude    = ids_both(ids_enough_Rpeaks & high_amplitude_only);
    % choose units with amplitude next to the threshold
    unit_ids_low_amplitude    = unit_ids(ids_enough_Rpeaks & low_amplitude_only);
    targets_low_amplitude     = targets(ids_enough_Rpeaks & low_amplitude_only);
    sites_low_amplitude       = sites(ids_enough_Rpeaks & low_amplitude_only);
    grid_x_low_amplitude      = grid_x(ids_enough_Rpeaks & low_amplitude_only);
    grid_y_low_amplitude      = grid_y(ids_enough_Rpeaks & low_amplitude_only);
    depths_low_amplitude      = depths(ids_enough_Rpeaks & low_amplitude_only);
    hemispheres_low_amplitude = hemispheres(ids_enough_Rpeaks & low_amplitude_only);
    ids_both_low_amplitude    = ids_both(ids_enough_Rpeaks & low_amplitude_only);
    % [any] low amplitude and significant cc between phase dynamics
    unit_ids_low_amplitude_ccs_any    = unit_ids(ids_enough_Rpeaks & low_amplitude_ccs_any);
    targets_low_amplitude_ccs_any     = targets(ids_enough_Rpeaks & low_amplitude_ccs_any);
    sites_low_amplitude_ccs_any       = sites(ids_enough_Rpeaks & low_amplitude_ccs_any);
    grid_x_low_amplitude_ccs_any      = grid_x(ids_enough_Rpeaks & low_amplitude_ccs_any);
    grid_y_low_amplitude_ccs_any      = grid_y(ids_enough_Rpeaks & low_amplitude_ccs_any);
    depths_low_amplitude_ccs_any      = depths(ids_enough_Rpeaks & low_amplitude_ccs_any);
    hemispheres_low_amplitude_ccs_any = hemispheres(ids_enough_Rpeaks & low_amplitude_ccs_any);
    ids_both_low_amplitude_ccs_any    = ids_both(ids_enough_Rpeaks & low_amplitude_ccs_any);
    % [both] low amplitude and significant cc between phase dynamics
    unit_ids_low_amplitude_ccs_both    = unit_ids(ids_enough_Rpeaks & low_amplitude_ccs_both);
    targets_low_amplitude_ccs_both     = targets(ids_enough_Rpeaks & low_amplitude_ccs_both);
    sites_low_amplitude_ccs_both       = sites(ids_enough_Rpeaks & low_amplitude_ccs_both);
    grid_x_low_amplitude_ccs_both      = grid_x(ids_enough_Rpeaks & low_amplitude_ccs_both);
    grid_y_low_amplitude_ccs_both      = grid_y(ids_enough_Rpeaks & low_amplitude_ccs_both);
    depths_low_amplitude_ccs_both      = depths(ids_enough_Rpeaks & low_amplitude_ccs_both);
    hemispheres_low_amplitude_ccs_both = hemispheres(ids_enough_Rpeaks & low_amplitude_ccs_both);
    ids_both_low_amplitude_ccs_both    = ids_both(ids_enough_Rpeaks & low_amplitude_ccs_both);
    % selected units - without low amplitude and pos. ccs (any)
    unit_ids_noLow_amplitude_ccs_any    = unit_ids(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    targets_noLow_amplitude_ccs_any     = targets(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    sites_noLow_amplitude_ccs_any       = sites(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    grid_x_noLow_amplitude_ccs_any      = grid_x(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    grid_y_noLow_amplitude_ccs_any      = grid_y(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    depths_noLow_amplitude_ccs_any      = depths(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    hemispheres_noLow_amplitude_ccs_any = hemispheres(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    ids_both_noLow_amplitude_ccs_any    = ids_both(ids_enough_Rpeaks & ~low_amplitude_ccs_any);
    % selected units - without low amplitude and pos. ccs (both)
    unit_ids_noLow_amplitude_ccs_both    = unit_ids(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    targets_noLow_amplitude_ccs_both     = targets(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    sites_noLow_amplitude_ccs_both       = sites(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    grid_x_noLow_amplitude_ccs_both      = grid_x(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    grid_y_noLow_amplitude_ccs_both      = grid_y(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    depths_noLow_amplitude_ccs_both      = depths(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    hemispheres_noLow_amplitude_ccs_both = hemispheres(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    ids_both_noLow_amplitude_ccs_both    = ids_both(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    % units with CBE - positively correlated with neuronal phase dynamic in
    % any condition
    unit_ids_with_CBE_ccs_any    = unit_ids(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    targets_with_CBE_ccs_any     = targets(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    sites_with_CBE_ccs_any       = sites(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    grid_x_with_CBE_ccs_any      = grid_x(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    grid_y_with_CBE_ccs_any      = grid_y(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    depths_with_CBE_ccs_any      = depths(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    hemispheres_with_CBE_ccs_any = hemispheres(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    ids_both_with_CBE_ccs_any    = ids_both(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    % units with CBE - positively correlated with neuronal phase dynamic in
    % both conditions
    unit_ids_with_CBE_ccs_both    = unit_ids(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    targets_with_CBE_ccs_both     = targets(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    sites_with_CBE_ccs_both       = sites(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    grid_x_with_CBE_ccs_both      = grid_x(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    grid_y_with_CBE_ccs_both      = grid_y(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    depths_with_CBE_ccs_both      = depths(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    hemispheres_with_CBE_ccs_both = hemispheres(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    ids_both_with_CBE_ccs_both    = ids_both(ids_enough_Rpeaks & ~low_amplitude_ccs_both);
    
    % save unit list for 600
    unit_ids    = unit_ids_600;
    targets     = targets_600;
    sites       = sites_600;
    grid_x      = grid_x_600;
    grid_y      = grid_y_600;
    depths      = depths_600;
    hemispheres = hemispheres_600;
    ids_both    = ids_both_600;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_600.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % noCB
    unit_ids    = unit_ids_noCB;
    targets     = targets_noCB;
    sites       = sites_noCB;
    grid_x      = grid_x_noCB;
    grid_y      = grid_y_noCB;
    depths      = depths_noCB;
    hemispheres = hemispheres_noCB;
    ids_both    = ids_both_noCB;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_noCB.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % withCB
    unit_ids    = unit_ids_withCB;
    targets     = targets_withCB;
    sites       = sites_withCB;
    grid_x      = grid_x_withCB;
    grid_y      = grid_y_withCB;
    depths      = depths_withCB;
    hemispheres = hemispheres_withCB;
    ids_both    = ids_both_withCB;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_withCB.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % noCB + amp criterion
    unit_ids    = unit_ids_noCB_corr;
    targets     = targets_noCB_corr;
    sites       = sites_noCB_corr;
    grid_x      = grid_x_noCB_corr;
    grid_y      = grid_y_noCB_corr;
    depths      = depths_noCB_corr;
    hemispheres = hemispheres_noCB_corr;
    ids_both    = ids_both_noCB_corr;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_noCB_corr.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % noCB + amp corterion + phase PSTH vs. feature dynamic correlation
    unit_ids    = unit_ids_noCB_corr_ccs;
    targets     = targets_noCB_corr_ccs;
    sites       = sites_noCB_corr_ccs;
    grid_x      = grid_x_noCB_corr_ccs;
    grid_y      = grid_y_noCB_corr_ccs;
    depths      = depths_noCB_corr_ccs;
    hemispheres = hemispheres_noCB_corr_ccs;
    ids_both    = ids_both_noCB_corr_ccs;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_noCB_corr_ccs.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % only high amplitude
    unit_ids    = unit_ids_high_amplitude;
    targets     = targets_high_amplitude;
    sites       = sites_high_amplitude;
    grid_x      = grid_x_high_amplitude;
    grid_y      = grid_y_high_amplitude;
    depths      = depths_high_amplitude;
    hemispheres = hemispheres_high_amplitude;
    ids_both    = ids_both_high_amplitude;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_high_amplitude.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % only low amplitude
    unit_ids    = unit_ids_low_amplitude;
    targets     = targets_low_amplitude;
    sites       = sites_low_amplitude;
    grid_x      = grid_x_low_amplitude;
    grid_y      = grid_y_low_amplitude;
    depths      = depths_low_amplitude;
    hemispheres = hemispheres_low_amplitude;
    ids_both    = ids_both_low_amplitude;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_low_amplitude.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % [any] low amplitude and significant cc between phase dynamics
    unit_ids    = unit_ids_low_amplitude_ccs_any;
    targets     = targets_low_amplitude_ccs_any;
    sites       = sites_low_amplitude_ccs_any;
    grid_x      = grid_x_low_amplitude_ccs_any;
    grid_y      = grid_y_low_amplitude_ccs_any;
    depths      = depths_low_amplitude_ccs_any;
    hemispheres = hemispheres_low_amplitude_ccs_any;
    ids_both    = ids_both_low_amplitude_ccs_any;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_low_amplitude_ccs_any.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % [both] low amplitude and significant cc between phase dynamics
    unit_ids    = unit_ids_low_amplitude_ccs_both;
    targets     = targets_low_amplitude_ccs_both;
    sites       = sites_low_amplitude_ccs_both;
    grid_x      = grid_x_low_amplitude_ccs_both;
    grid_y      = grid_y_low_amplitude_ccs_both;
    depths      = depths_low_amplitude_ccs_both;
    hemispheres = hemispheres_low_amplitude_ccs_both;
    ids_both    = ids_both_low_amplitude_ccs_both;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_low_amplitude_ccs_both.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % [any] no low amplitude and no pos. significant cc between phase dynamics
    unit_ids    = unit_ids_noLow_amplitude_ccs_any;
    targets     = targets_noLow_amplitude_ccs_any;
    sites       = sites_noLow_amplitude_ccs_any;
    grid_x      = grid_x_noLow_amplitude_ccs_any;
    grid_y      = grid_y_noLow_amplitude_ccs_any;
    depths      = depths_noLow_amplitude_ccs_any;
    hemispheres = hemispheres_noLow_amplitude_ccs_any;
    ids_both    = ids_both_noLow_amplitude_ccs_any;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_noLow_amplitude_ccs_any.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % [both] no low amplitude and no pos. significant cc between phase dynamics
    unit_ids    = unit_ids_noLow_amplitude_ccs_both;
    targets     = targets_noLow_amplitude_ccs_both;
    sites       = sites_noLow_amplitude_ccs_both;
    grid_x      = grid_x_noLow_amplitude_ccs_both;
    grid_y      = grid_y_noLow_amplitude_ccs_both;
    depths      = depths_noLow_amplitude_ccs_both;
    hemispheres = hemispheres_noLow_amplitude_ccs_both;
    ids_both    = ids_both_noLow_amplitude_ccs_both;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_noLow_amplitude_ccs_both.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % [any] with CBE based on correlation
    unit_ids    = unit_ids_with_CBE_ccs_any;
    targets     = targets_with_CBE_ccs_any;
    sites       = sites_with_CBE_ccs_any;
    grid_x      = grid_x_with_CBE_ccs_any;
    grid_y      = grid_y_with_CBE_ccs_any;
    depths      = depths_with_CBE_ccs_any;
    hemispheres = hemispheres_with_CBE_ccs_any;
    ids_both    = ids_both_with_CBE_ccs_any;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_with_CBE_ccs_any.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % [both] with CBE based on correlation
    unit_ids    = unit_ids_with_CBE_ccs_both;
    targets     = targets_with_CBE_ccs_both;
    sites       = sites_with_CBE_ccs_both;
    grid_x      = grid_x_with_CBE_ccs_both;
    grid_y      = grid_y_with_CBE_ccs_both;
    depths      = depths_with_CBE_ccs_both;
    hemispheres = hemispheres_with_CBE_ccs_both;
    ids_both    = ids_both_with_CBE_ccs_both;
    save([basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_with_CBE_ccs_both.mat'], 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    
    % 4. save the table with unit counts
    % 600
    unique_areas = unique(targets_600);
    unit_counts = count_units_by_target(targets_600, unique_areas);
    % noCB
    unit_counts_noCB = count_units_by_target(targets_noCB, unique_areas);
    % withCB
    unit_counts_withCB = count_units_by_target(targets_withCB, unique_areas);
    % noCB + amp
    unit_counts_noCB_corr = count_units_by_target(targets_noCB_corr, unique_areas);
    % noCB + amp + high/low cc PSTH vs. feature dynamic
    unit_counts_noCB_corr_ccs = count_units_by_target(targets_noCB_corr_ccs, unique_areas);
    % only high amplitude
    unit_counts_highAMP = count_units_by_target(targets_high_amplitude, unique_areas);
    % only low amplitude
    unit_counts_low_amplitude = count_units_by_target(targets_low_amplitude, unique_areas);
    % [any] low amplitude + sig amp
    unit_counts_low_amplitude_ccs_any = count_units_by_target(targets_low_amplitude_ccs_any, unique_areas);
    % [both] low amplitude + sig amp
    unit_counts_low_amplitude_ccs_both = count_units_by_target(targets_low_amplitude_ccs_both, unique_areas);
    % [any] no low amplitude and no pos. significant cc between phase dynamics
    unit_counts_noLow_amplitude_ccs_any = count_units_by_target(targets_noLow_amplitude_ccs_any, unique_areas);
    % [both] no low amplitude and no pos. significant cc between phase dynamics
    unit_counts_noLow_amplitude_ccs_both = count_units_by_target(targets_noLow_amplitude_ccs_both, unique_areas);
    % [any] with CBE based on correlation
    unit_counts_with_CBE_ccs_any = count_units_by_target(targets_noLow_amplitude_ccs_both, unique_areas);
    % [both] with CBE based on correlation
    unit_counts_with_CBE_ccs_both = count_units_by_target(targets_noLow_amplitude_ccs_both, unique_areas);
    T = ...
        table(unique_areas, unit_counts, unit_counts_noCB, ...
        unit_counts_withCB, unit_counts_noCB_corr, ...
        unit_counts_noCB_corr_ccs, unit_counts_highAMP, ...
        unit_counts_low_amplitude, unit_counts_low_amplitude_ccs_any, ...
        unit_counts_low_amplitude_ccs_both, ...
        unit_counts_noLow_amplitude_ccs_any, unit_counts_noLow_amplitude_ccs_both, ...
        unit_counts_with_CBE_ccs_any, unit_counts_with_CBE_ccs_both);
    writetable(T, [basepath_to_save filesep list_of_lists(listNum).name(1:end-4) '_600.xls'])
    %     clear T unique_areas ic unit_counts filename unit_ids_after_exclusion targets_after_exclusion
    
end

end

function unit_counts = count_units_by_target(targets, unique_areas)
[unqTargets, ~, ic] = unique(targets);
unit_counts = accumarray(ic,1);
Lia = ismember(unique_areas, unqTargets);
unit_counts_tmp = double(Lia);
unit_counts_tmp(Lia) = unit_counts;
unit_counts = unit_counts_tmp;
end
