function ecg_bna_compute_session_phase_response_analysis(trials,population,Rpeaks,cfg)
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
        
        % if there are fewer than 10 spikes, don't bother and proceed to
        % other units
        if size(WF_one_stream,1) < 10
            continue
        end
        
        % if there are fewer than 3 RR intervals, don't bother and proceed 
        % to other units
        if length(valid_RRinterval_starts) < 3
            continue
        end
        
        %% phase response analysis
        
        % 1. compute time from the previous heartbeat to the current spike
        RR_starts_90 = valid_RRinterval_starts(RR_within_trial90_idx);
        RR_ends_90   = valid_RRinterval_ends(RR_within_trial90_idx);
        
        % find index of previous heartbeat for each spike
        tmp      = AT_one_stream > RR_starts_90;
        tmp_cell = mat2cell(tmp,ones(length(AT_one_stream),1),length(RR_starts_90));
        RR_ind   = cellfun(@(x) find(x,1,'last'), tmp_cell,'UniformOutput',false);
        empty_id = cellfun(@isempty,RR_ind);
        RR_ind   = [RR_ind{~empty_id}];
        
        % time from the previous heartbeat to the current spike
        deltaT = 1000 * (AT_one_stream(~empty_id)' - RR_starts_90(RR_ind)); % ms
        
        % 2. cardiac cycle change
        RR_durations = RR_ends_90(RR_ind) - RR_starts_90(RR_ind);
        deltaRR = 1000 * [diff() 0]; % ms
%         deltaRR = deltaRR(1:end-1);
        deltaRR(deltaRR == 0) = NaN;
        
        % put data in the output structure - this is to create rasters as
        % in Kim et al., 2019
        data.(L).valid_RRinterval_starts = valid_RRinterval_starts(RR_within_trial90_idx);
        data.(L).valid_RRinterval_ends   = valid_RRinterval_ends(RR_within_trial90_idx);
        data.(L).AT_one_stream           = AT_one_stream;
        
        

end
end
