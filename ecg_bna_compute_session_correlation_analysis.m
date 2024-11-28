function ecg_bna_compute_session_correlation_analysis(trials,population,Rpeaks,cfg)
% Function to perform the correlation analysis between RR durations and
% corresponding FR within those cycles.
% This analysis had two important considerations. To be able to compare the
% results of shifting data related to one another between task and rest, we
% made spike and ECG data continuous and omitted invalid R-peaks of invalid
% spike periods. Took all trials, not only 

basepath_to_save=[cfg.SPK_root_results_fldr filesep 'correlation_analysis'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

%% load list of selected units - won't process all the crap
if cfg.spk.apply_exclusion_criteria
    
    unit_list = load([cfg.SPK_root_results_fldr filesep 'unit_lists' filesep cfg.spk.unit_list]);
    unitList = unique(unit_list.unit_ids);
    
    % figure out which units take from this session
    selected_this_session = ismember({population.unit_ID}, unitList);
    population = population(selected_this_session);
    
end

Rblocks=[Rpeaks.block];

for unitNum = 1:length(population)
    
    %% get to processing
    disp(['Processing unit ' num2str(unitNum) ' out of ' num2str(length(population))])
    
    pop=population(unitNum);
    
    T=ph_get_unit_trials(pop,trials);
    
    %% Make sure we only take overlapping blocks
    blocks_unit=unique([pop.block]);
    blocks=intersect(blocks_unit,Rblocks);
    
    %% preallocate 'data' structure
    data.unitId            = pop.unit_ID;
    data.target            = pop.target;
    data.channel           = pop.channel;
    data.unit              = pop.block_unit{2,1};
    data.quantSNR          = pop.avg_SNR;
    data.Single_rating     = pop.avg_single_rating;
    data.stability_rating  = pop.avg_stability;
    data.thresholds_microV = single([NaN; NaN; NaN; NaN]);
    data.FR                = single(mean(pop.FR_average));
    data.cc_lag_list       = cfg.correlation.lag_list;
    data.criteria          = pop.criteria;
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        data.(L).valid_RRinterval_starts = single(nan(1,1));
        data.(L).valid_RRinterval_ends   = single(nan(1,1));
        data.(L).AT_one_stream           = single(nan(1,1));
        data.(L).timeRRstart             = single(nan(1,1));
        data.(L).FRbyRR_Hz               = single(nan(1,1));
        data.(L).cycleDurations_s        = single(nan(1,1));
        data.(L).HR_Hz                   = single(nan(1,1));
        data.(L).n_cycles                = single(nan(1,1));
        data.(L).pearson_r               = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).pearson_p               = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).permuted_p              = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).FR_n_cycles             = single(nan(1,1));
        data.(L).FR_pearson_r            = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).FR_pearson_p            = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).FR_permuted_p           = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).RR_n_cycles             = single(nan(1,1));
        data.(L).RR_pearson_r            = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).RR_pearson_p            = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).RR_permuted_p           = single(nan(length(cfg.correlation.lag_list), 1));
        
        data.(L).Hz_n_cycles             = single(nan(1,1));
        data.(L).Hz_pearson_r            = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).Hz_pearson_p            = single(nan(length(cfg.correlation.lag_list), 1));
        data.(L).Hz_permuted_p           = single(nan(length(cfg.correlation.lag_list), 1));
        
    end
    
    %%
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        
        % 1. take all the trials
        % 2. figure out spike and RR data for those trials
        % 3. figure out indices for invalid (wrong block, wrong condition, not completed)
        % 4. replace data belonging to wrong trials with NANs (both for spike and RR)
        % 5. replace data from RR intervals belonging to ITIs with NANs
        
        % What to do with inter-block intervals? - As one doesn't know how
        % many RR-intervals are there between blocks, let's replace those
        % with 12 RR slots - those will be nans and the corresponding FRs
        % as well
        
        %% get condition AND valid block trials only
        % create special trial configuration
        CT = ecg_bna_get_condition_trials(T, cfg.condition(c));
        tr=ismember([T.block],blocks) & CT;
        if sum(tr)<=1 || (~isfield(Rpeaks, 'RPEAK_ts_insp') && cfg.process_Rpeaks_inhalation_exhalation) ...
                      || (~isfield(Rpeaks, 'RPEAK_ts_exp')  && cfg.process_Rpeaks_inhalation_exhalation) % do calculations only if number of trials > 1
            continue
        end
        popcell=pop.trial(tr);
        trcell=T(tr);
        
        % 0. Prepare data variables
        % 0.1. Figure out house-keeping variables
        states_onset               = {trcell.states_onset};
        states                     = {trcell.states};
        TDT_ECG1_t0_from_rec_start = {trcell.TDT_ECG1_t0_from_rec_start};
        block_nums                 = {trcell.block};
        state2_times               = cellfun(@(x,y) x(y == 2), states_onset, states, 'Uniformoutput', false); % trial starts = state 2
        state90_times              = cellfun(@(x,y) x(y == 90), states_onset, states, 'Uniformoutput', false); % trial ends = state 90
%         state98_times              = cellfun(@(x,y) x(y == 98), states_onset, states, 'Uniformoutput', false); % trial ends = state 98
        % 0.2. Compute RR-intervals
        valid_RRinterval_ends      = single([Rpeaks.(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
        if any(isnan([Rpeaks.(['RPEAK_ts' cfg.condition(c).Rpeak_field])]))
            thisisinteresing=1
        end
        %valid_RRinterval_ends      = valid_RRinterval_ends(~isnan(valid_RRinterval_ends)); %% somehow, many of those can be NaN ?????
        valid_RRinterval_starts    = single(valid_RRinterval_ends - [Rpeaks.(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
        % 1. Figure out RR-intervals lying within trials
        trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
        trial_ends90_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
        %         trial_ends98_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state98_times, TDT_ECG1_t0_from_rec_start, block_nums);
        
        RR_within_trial90_idx = false(length(valid_RRinterval_starts),1)';
        %         RR_within_trial98_idx = false(length(valid_RRinterval_starts),1);
        for RRnum = 1:length(RR_within_trial90_idx)
            % indexing for FR variable
            RR_within_trial90_idx(RRnum) = ...
                any(valid_RRinterval_starts(RRnum)>trial_starts_one_stream & ...
                valid_RRinterval_ends(RRnum)<trial_ends90_one_stream);
            
            %             % indexing for RR variable
            %             RR_within_trial98_idx(RRnum) = ...
            %                 any(valid_RRinterval_starts(RRnum)>trial_starts_one_stream & ...
            %                 valid_RRinterval_ends(RRnum)<trial_ends98_one_stream);
        end
              
        % 2. Prepare spiking data
        % 2.1. take arrival times and the corresponding waveforms
        AT = {popcell.arrival_times};
        WF = {popcell.waveforms};
        % 2.2. choose only those that happen after MP-state 1
        idx_after_state1 = cellfun(@(x,y) x>y, AT, state2_times, 'Uniformoutput', false);
        % 2.3. add TDT_ECG1_t0_from_rec_start and Rpeak block offset to spike times
        AT_one_stream_cell = cellfun(@(x,y,z,a) x(y)+z+Rpeaks([Rpeaks.block] == a).offset, AT, idx_after_state1, TDT_ECG1_t0_from_rec_start, block_nums, 'Uniformoutput', false);
        WF_one_stream_cell = cellfun(@(x,y) x(y,:), WF, idx_after_state1, 'Uniformoutput', false);
        % 2.4. merge all spike times and waveforms together
        AT_one_stream = cat(1, AT_one_stream_cell{:});
        WF_one_stream = cat(1, WF_one_stream_cell{:});
        if size(WF_one_stream,1) < 10
            continue
        end
        
        % check the number of spikes left after computing phases
        if length(valid_RRinterval_starts) < 3
            continue
        end
        
        % put data in the output structure - this is to create rasters as
        % in Kim et al., 2019
        data.(L).valid_RRinterval_starts = valid_RRinterval_starts(RR_within_trial90_idx);
        data.(L).valid_RRinterval_ends   = valid_RRinterval_ends(RR_within_trial90_idx);
        data.(L).AT_one_stream           = AT_one_stream;
        
        % % figure out how many R peaks I need to cover long intervals between R-peaks
        % find long intervals between R-peaks
        % long_RRs = diff(valid_RRinterval_starts) > 0.6; % hard coded value for a long RR-interval
        %divMed   = floor(diff(valid_RRinterval_starts)/0.4); % 0.4 - median RR for Magnus
        %modRR    = mod(diff(valid_RRinterval_starts),0.4);
        
        % find intervals between R-peaks that are long enough and decide on
        % the number of additional R-peaks for replacement
        % ids          = find(long_RRs);
        % n_add_Rpeaks = divMed(long_RRs) + floor(modRR(long_RRs)/0.25) -1;
        
        
        n_add_Rpeaks = round(diff(valid_RRinterval_starts)/0.4) -1; %% these are cycle starts! so the number of cycles -1 (!)
        ids          = find(n_add_Rpeaks);
        n_add_Rpeaks = n_add_Rpeaks(n_add_Rpeaks~=0);
        
        % loop through long intervals to fill them with NaNs backwards to
        % keep index of former elements unchanged
        for intRRnum = length(ids):-1:1
            valid_RRinterval_starts = [valid_RRinterval_starts(1:ids(intRRnum)) nan(1, n_add_Rpeaks(intRRnum)) valid_RRinterval_starts(ids(intRRnum)+1:end)];
            valid_RRinterval_ends = [valid_RRinterval_ends(1:ids(intRRnum)) nan(1, n_add_Rpeaks(intRRnum)) valid_RRinterval_ends(ids(intRRnum)+1:end)];
            RR_within_trial90_idx = [RR_within_trial90_idx(1:ids(intRRnum)) false(1, n_add_Rpeaks(intRRnum)) RR_within_trial90_idx(ids(intRRnum)+1:end)];
        end
        
        % II. Compute unit firing rate per RR-interval
        data.(L).timeRRstart           = valid_RRinterval_starts;
        [FRbyRR_Hz, cycleDurations_s]  = ...
            computeFRperCycle(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream);
        FRbyRR_Hz(~RR_within_trial90_idx)        = NaN; 
        HR_Hz   = 1/cycleDurations_s;
        HR_Hz   = HR_Hz';
        % option 2
%         FRbyRR_Hz(~RR_within_trial90_idx)        = NaN;
%         cycleDurations_s(~RR_within_trial98_idx) = NaN;
        % option 1
%         FRbyRR_Hz(~RR_within_trial90_idx)        = NaN;
%         cycleDurations_s(~RR_within_trial90_idx) = NaN;
        data.(L).FRbyRR_Hz        = FRbyRR_Hz;
        data.(L).cycleDurations_s = cycleDurations_s;
        data.(L).HR_Hz            = HR_Hz;
        
        % compute correlation with different lag        
        [data.(L).n_cycles, data.(L).pearson_r, data.(L).pearson_p, data.(L).permuted_p] = ...
            compute_correlation_by_lag(cfg,FRbyRR_Hz,cycleDurations_s);
        
        % autocorrelation FRbyRR_Hz
        [data.(L).FR_n_cycles, data.(L).FR_pearson_r, data.(L).FR_pearson_p, data.(L).FR_permuted_p] = ...
            compute_correlation_by_lag(cfg,FRbyRR_Hz,FRbyRR_Hz);
        
        % autocorrelation cycleDurations_s
        [data.(L).RR_n_cycles, data.(L).RR_pearson_r, data.(L).RR_pearson_p, data.(L).RR_permuted_p] = ...
            compute_correlation_by_lag(cfg,cycleDurations_s,cycleDurations_s);
        
        
        % compute correlation with different lag between FR in Hz and HR in
        % Hz
        [data.(L).Hz_n_cycles, data.(L).Hz_pearson_r, data.(L).Hz_pearson_p, data.(L).Hz_permuted_p] = ...
            compute_correlation_by_lag(cfg, FRbyRR_Hz, HR_Hz);
        
    end
    save([basepath_to_save filesep data.unitId '_' data.target '_correlation.mat'], 'data', '-v7.3')
    clear data

end
end

function [n_cycles, pearson_r, pearson_p, permuted_p] = compute_correlation_by_lag(cfg,FRbyRR_Hz,cycleDurations_s)
n_permutations = cfg.correlation.n_permutations;
tail           = cfg.correlation.tail;
alpha_level    = cfg.correlation.alpha_level;
stat           = cfg.correlation.stat;
reports        = cfg.correlation.reports;
seed_state     = cfg.correlation.seed_state;

Min_overlap=cfg.spk.unit_exclusion.nCardiacCycles/2;

% Initialize output structure to avoid broadcasting issues
lag_list        = cfg.correlation.lag_list;
nLags           = length(lag_list);
n_cycles   = zeros(1, nLags);
pearson_r  = zeros(1, nLags);
pearson_p  = zeros(1, nLags);
permuted_p = zeros(1, nLags);

parfor lagNum = 1:nLags
    curr_lag = lag_list(lagNum);
    lag_abs=abs(curr_lag);
    lag_sign=sign(curr_lag);
    if lag_sign<0
        % If lag is negative, FR follows RR
        fr_hz = FRbyRR_Hz(1+lag_abs:end);
        rr_s  = cycleDurations_s(1:end-lag_abs);
    else
        % If lag is positive, RR follows FR
        fr_hz = FRbyRR_Hz(1:end-lag_abs);
        rr_s  = cycleDurations_s(1+lag_abs:end);
    end
    
    % compute correlation coefficient
    valid = ~isnan(fr_hz) & ~isnan(rr_s);
    [temp_r, temp_p] = corrcoef(fr_hz(valid), rr_s(valid));
    n_cycles(lagNum)  = sum(valid);
    if sum(valid)>Min_overlap
        pearson_r(lagNum) = temp_r(2,1);
        pearson_p(lagNum) = temp_p(2,1);
        
        % i dont think this makes sense at all unfortunately... the
        % multicomparison correction should be due to multiple lags, no?
        % NO! skewdness of both parameters to correlate causes a chance
        % level correlation different from 0 !!!
        permuted_p(lagNum) = mult_comp_perm_corr(fr_hz(valid), rr_s(valid), n_permutations, tail, alpha_level, stat, reports, seed_state);
    else
        
        pearson_r(lagNum) = NaN;
        pearson_p(lagNum) = NaN;
        permuted_p(lagNum)= NaN;
    end
end

end

function [FRbyRR_Hz, cycleDurations] = computeFRperCycle(interval_starts, interval_ends, eventTimes)
% make inputs vertical vectors
interval_starts = interval_starts';
interval_ends = interval_ends';
eventTimes = eventTimes';

% Calculate ECG cycle durations
cycleDurations = interval_ends - interval_starts;

% cycleNums = arrayfun(@(x) sum(x > interval_starts & x < interval_ends), eventTimes, 'UniformOutput', false);
spikeCounts = arrayfun(@(x,y) sum(eventTimes > x & eventTimes < y), interval_starts, interval_ends);
FRbyRR_Hz = spikeCounts ./ cycleDurations;

% figure,
% scatter(cycleDurations, FRbyRR_Hz)
end