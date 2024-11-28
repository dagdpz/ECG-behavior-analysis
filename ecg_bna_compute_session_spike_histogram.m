function Output=ecg_bna_compute_session_spike_histogram(trials,population,Rpeaks,sessions_info,cfg)
Sanity_check=0; % ECG triggered ECG, turn off since typically there is no ECG data in the spike format

for numTiming = 1:length(cfg.analyse_states)

    curr_analyse_states = cfg.analyse_states{numTiming};
    BINS=(curr_analyse_states{1,3}:cfg.time.PSTH_binwidth:curr_analyse_states{1,4})*1000;
    offset_blocks_Rpeak=[Rpeaks.offset];
    Rblocks=[Rpeaks.block];

    % for pseudocircular fitting
    bin_ids = 1+curr_analyse_states{5}:length(BINS)-curr_analyse_states{6};

    if cfg.spk.apply_exclusion_criteria
        load([cfg.unit_lists filesep cfg.spk.unit_list], 'unit_ids')
        
        pop_units  = {population.unit_ID};
        inc_ids    = ismember(pop_units, unit_ids);
        population = population(inc_ids);
        
    end
    
    for u=1:numel(population)
        tic
        pop=population(u);
        disp(['Processing ' pop.unit_ID])
        T=ph_get_unit_trials(pop,trials);
        
        %% Make sure we only take overlapping blocks
        blocks_unit=unique([pop.block]);
        blocks=intersect(blocks_unit,Rblocks);
        b=ismember(Rblocks,blocks);
        
        % preallocate 'Output' structure
        for c=1:numel(cfg.condition)
            L=cfg.condition(c).name;
            % unit id and recording quality parameters
            Output.unit_ID               = pop.unit_ID;
            Output.target                = pop.target;
            Output.quantSNR              = pop.avg_SNR;
            Output.Single_rating         = pop.avg_single_rating;
            Output.stability_rating      = pop.avg_stability;
            Output.criteria              = pop.criteria;
            % ecg-related variables
            Output.(L).SD                = single(nan(1, length(BINS)));
            Output.(L).SD_STD            = single(nan(1, length(BINS)));
            Output.(L).SD_SEM            = single(nan(1, length(BINS)));
            Output.(L).SDCL              = single(nan(1, length(BINS)));
            Output.(L).SDCu              = single(nan(1, length(BINS)));
            Output.(L).SDP               = single(nan(1, length(BINS)));
            Output.(L).SDPCL             = single(nan(1, length(BINS)));
            Output.(L).SDPCu             = single(nan(1, length(BINS)));
            Output.(L).sig_all           = single(zeros(1, length(BINS)));
            Output.(L).sig               = single(zeros(1, length(BINS)));
            Output.(L).sig_FR_diff       = single(nan(1));
            Output.(L).sig_time          = single(nan(1));
            Output.(L).sig_n_bins        = single(zeros(1));
            Output.(L).sig_sign          = single(zeros(1));
            Output.(L).NrTrials          = single(nan(1));
            Output.(L).NrEvents          = single(nan(1));
            Output.(L).FR                = single(nan(1));
            Output.(L).raster            = single(nan(1));
            Output.(L).Rts               = single(nan(1)); % RR ends
            Output.(L).Rds               = single(nan(1)); % RR durations
            Output.(L).Rds_perm          = single(nan(1));
            Output.(L).SDsubstractedSDP            = single(nan(1, length(BINS)));
            Output.(L).SDsubstractedSDP_normalized = single(nan(1, length(BINS)));
            Output.(L).FR_ModIndex_SubtrSDP        = single(nan(1));
            Output.(L).FR_ModIndex_PcS             = single(nan(1));
            
            Output.bin_centers                     = BINS(bin_ids);
            Output.phase_bin_centers               = 2 * pi * BINS(bin_ids) / range(BINS(bin_ids));
                                                     
            Output.(L).lowIBI.SD          = single(nan(1,length(BINS)));
            Output.(L).lowIBI.SD_STD      = single(nan(1,length(BINS)));
            Output.(L).lowIBI.SD_SEM      = single(nan(1,length(BINS)));
            Output.(L).lowIBI.SDP         = single(nan(1,length(BINS)));
            Output.(L).lowIBI.SDPCL       = single(nan(1,length(BINS)));
            Output.(L).lowIBI.SDPCu       = single(nan(1,length(BINS)));
            Output.(L).lowIBI.sig_all     = single(zeros(1,length(BINS)));
            Output.(L).lowIBI.sig         = single(zeros(1,length(BINS)));
            Output.(L).lowIBI.sig_FR_diff = single(nan(1));
            Output.(L).lowIBI.sig_time    = single(nan(1));
            Output.(L).lowIBI.sig_n_bins  = single(zeros(1));
            Output.(L).lowIBI.sig_sign    = single(zeros(1));
            Output.(L).lowIBI.NrEvents    = single(nan(1));
            Output.(L).lowIBI.SDsubstractedSDP            = single(nan(1,length(BINS)));
            Output.(L).lowIBI.SDsubstractedSDP_normalized = single(nan(1,length(BINS)));
            Output.(L).lowIBI.FR_ModIndex_SubtrSDP        = single(nan(1));
            Output.(L).lowIBI.FR_ModIndex_PcS             = single(nan(1));
            
            Output.(L).highIBI.SD          = single(nan(1,length(BINS)));
            Output.(L).highIBI.SD_STD      = single(nan(1,length(BINS)));
            Output.(L).highIBI.SD_SEM      = single(nan(1,length(BINS)));
            Output.(L).highIBI.SDP         = single(nan(1,length(BINS)));
            Output.(L).highIBI.SDPCL       = single(nan(1,length(BINS)));
            Output.(L).highIBI.SDPCu       = single(nan(1,length(BINS)));
            Output.(L).highIBI.sig_all     = single(zeros(1,length(BINS)));
            Output.(L).highIBI.sig         = single(zeros(1,length(BINS)));
            Output.(L).highIBI.sig_FR_diff = single(nan(1));
            Output.(L).highIBI.sig_time    = single(nan(1));
            Output.(L).highIBI.sig_n_bins  = single(zeros(1));
            Output.(L).highIBI.sig_sign    = single(zeros(1));
            Output.(L).highIBI.NrEvents    = single(nan(1));
            Output.(L).highIBI.SDsubstractedSDP            = single(nan(1,length(BINS)));
            Output.(L).highIBI.SDsubstractedSDP_normalized = single(nan(1,length(BINS)));
            Output.(L).highIBI.FR_ModIndex_SubtrSDP        = single(nan(1));
            Output.(L).highIBI.FR_ModIndex_PcS             = single(nan(1));
            
%             Output.(L).linear.yfit                 = single(nan(1,length(bin_ids)));
%             Output.(L).linear.coefs                = single([NaN NaN]);
%             Output.(L).linear.rsquared             = single(NaN);
%             Output.(L).linear.adjrsquared          = single(NaN);
%             Output.(L).linear.sse                  = single(NaN);
%             Output.(L).linear.dfe                  = single(NaN);
%             Output.(L).linear.rmse                 = single(NaN);
%             Output.(L).linear.pvalue               = single([NaN; NaN]);
%             Output.(L).linear.aic                  = single(NaN);
%             Output.(L).linear.bic                  = single(NaN);
%             
%             Output.(L).cosine.average              = single(nan(1,length(bin_ids)));
%             Output.(L).cosine.startPoint           = single([NaN NaN NaN]);
%             Output.(L).cosine.yfit                 = single(nan(1,length(bin_ids)));
%             Output.(L).cosine.coefs                = single([NaN NaN NaN]);
%             Output.(L).cosine.rsquared             = single(NaN);
%             Output.(L).cosine.adjrsquared          = single(NaN);
%             Output.(L).cosine.sse                  = single(NaN);
%             Output.(L).cosine.dfe                  = single(NaN);
%             Output.(L).cosine.rmse                 = single(NaN);
%             Output.(L).cosine.pvalue               = single(NaN);
%             Output.(L).cosine.aic                  = single(NaN);
%             Output.(L).cosine.bic                  = single(NaN);
%             
%             Output.(L).vonMisesPos.average         = single(nan(1,length(bin_ids)));
%             Output.(L).vonMisesPos.startPoint      = single([NaN NaN NaN NaN]);
%             Output.(L).vonMisesPos.yfit            = single(nan(1,length(bin_ids)));
%             Output.(L).vonMisesPos.coefs           = single([NaN NaN NaN NaN]);
%             Output.(L).vonMisesPos.rsquared        = single(NaN);
%             Output.(L).vonMisesPos.adjrsquared     = single(NaN);
%             Output.(L).vonMisesPos.sse             = single(NaN);
%             Output.(L).vonMisesPos.dfe             = single(NaN);
%             Output.(L).vonMisesPos.rmse            = single(NaN);
%             Output.(L).vonMisesPos.pvalue          = single(NaN);
%             Output.(L).vonMisesPos.aic             = single(NaN);
%             Output.(L).vonMisesPos.bic             = single(NaN);
%             
%             Output.(L).vonMisesNeg.average         = single(nan(1,length(bin_ids)));
%             Output.(L).vonMisesNeg.startPoint      = single([NaN NaN NaN NaN]);
%             Output.(L).vonMisesNeg.yfit            = single(nan(1,length(bin_ids)));
%             Output.(L).vonMisesNeg.coefs           = single([NaN NaN NaN NaN]);
%             Output.(L).vonMisesNeg.rsquared        = single(NaN);
%             Output.(L).vonMisesNeg.adjrsquared     = single(NaN);
%             Output.(L).vonMisesNeg.sse             = single(NaN);
%             Output.(L).vonMisesNeg.dfe             = single(NaN);
%             Output.(L).vonMisesNeg.rmse            = single(NaN);
%             Output.(L).vonMisesNeg.pvalue          = single(NaN);
%             Output.(L).vonMisesNeg.aic             = single(NaN);
%             Output.(L).vonMisesNeg.bic             = single(NaN);
            
            % median split - IBI low
            Output.(L).IBI_median                  = single(NaN);
            
            Output.(L).lowIBI.SD                   = single(nan(1, length(BINS)));
            Output.(L).highIBI.SD                  = single(nan(1, length(BINS)));
            
            Output.(L).lowIBI_timeRRstart             = single(NaN);
            Output.(L).lowIBI_timeRRend               = single(NaN);
            Output.(L).lowIBI_starts                  = single(NaN);
            Output.(L).lowIBI_meanHR_bpm              = single(NaN);
            Output.(L).lowIBI_medianHR_bpm            = single(NaN);
            Output.(L).lowIBI_stdHR_bpm               = single(NaN);
            Output.(L).lowIBI_SDNN_ms                 = single(NaN);
            
            % median split - IBI high
            Output.(L).highIBI_timeRRstart            = single(NaN);
            Output.(L).highIBI_timeRRend              = single(NaN);
            Output.(L).highIBI_starts                 = single(NaN);
            Output.(L).highIBI_meanHR_bpm             = single(NaN);
            Output.(L).highIBI_medianHR_bpm           = single(NaN);
            Output.(L).highIBI_stdHR_bpm              = single(NaN);
            Output.(L).highIBI_SDNN_ms                = single(NaN);
            
%             Output.(L).lowIBI_linear.yfit             = single(nan(1,length(bin_ids)));
%             Output.(L).lowIBI_linear.coefs            = single([NaN NaN]);
%             Output.(L).lowIBI_linear.rsquared         = single(NaN);
%             Output.(L).lowIBI_linear.adjrsquared      = single(NaN);
%             Output.(L).lowIBI_linear.sse              = single(NaN);
%             Output.(L).lowIBI_linear.dfe              = single(NaN);
%             Output.(L).lowIBI_linear.rmse             = single(NaN);
%             Output.(L).lowIBI_linear.pvalue           = single([NaN; NaN]);
%             Output.(L).lowIBI_linear.aic                = single(NaN);
%             Output.(L).lowIBI_linear.bic                = single(NaN);
%             
%             Output.(L).lowIBI_cosine.average          = single(nan(1,length(bin_ids)));
%             Output.(L).lowIBI_cosine.startPoint       = single([NaN NaN NaN]);
%             Output.(L).lowIBI_cosine.yfit             = single(nan(1,length(bin_ids)));
%             Output.(L).lowIBI_cosine.coefs            = single([NaN NaN NaN]);
%             Output.(L).lowIBI_cosine.rsquared         = single(NaN);
%             Output.(L).lowIBI_cosine.adjrsquared      = single(NaN);
%             Output.(L).lowIBI_cosine.sse              = single(NaN);
%             Output.(L).lowIBI_cosine.dfe              = single(NaN);
%             Output.(L).lowIBI_cosine.rmse             = single(NaN);
%             Output.(L).lowIBI_cosine.pvalue           = single(NaN);
%             Output.(L).lowIBI_cosine.aic                = single(NaN);
%             Output.(L).lowIBI_cosine.bic                = single(NaN);
%             
%             Output.(L).lowIBI_vonMisesPos.average     = single(nan(1,length(bin_ids)));
%             Output.(L).lowIBI_vonMisesPos.startPoint  = single([NaN NaN NaN NaN]);
%             Output.(L).lowIBI_vonMisesPos.yfit        = single(nan(1,length(bin_ids)));
%             Output.(L).lowIBI_vonMisesPos.coefs       = single([NaN NaN NaN NaN]);
%             Output.(L).lowIBI_vonMisesPos.rsquared    = single(NaN);
%             Output.(L).lowIBI_vonMisesPos.adjrsquared = single(NaN);
%             Output.(L).lowIBI_vonMisesPos.sse         = single(NaN);
%             Output.(L).lowIBI_vonMisesPos.dfe         = single(NaN);
%             Output.(L).lowIBI_vonMisesPos.rmse        = single(NaN);
%             Output.(L).lowIBI_vonMisesPos.pvalue      = single(NaN);
%             Output.(L).lowIBI_vonMisesPos.aic                = single(NaN);
%             Output.(L).lowIBI_vonMisesPos.bic                = single(NaN);
%             
%             Output.(L).lowIBI_vonMisesNeg.average     = single(nan(1,length(bin_ids)));
%             Output.(L).lowIBI_vonMisesNeg.startPoint  = single([NaN NaN NaN NaN]);
%             Output.(L).lowIBI_vonMisesNeg.yfit        = single(nan(1,length(bin_ids)));
%             Output.(L).lowIBI_vonMisesNeg.coefs       = single([NaN NaN NaN NaN]);
%             Output.(L).lowIBI_vonMisesNeg.rsquared    = single(NaN);
%             Output.(L).lowIBI_vonMisesNeg.adjrsquared = single(NaN);
%             Output.(L).lowIBI_vonMisesNeg.sse         = single(NaN);
%             Output.(L).lowIBI_vonMisesNeg.dfe         = single(NaN);
%             Output.(L).lowIBI_vonMisesNeg.rmse        = single(NaN);
%             Output.(L).lowIBI_vonMisesNeg.pvalue      = single(NaN);
%             Output.(L).lowIBI_vonMisesNeg.aic         = single(NaN);
%             Output.(L).lowIBI_vonMisesNeg.bic         = single(NaN);
%             
%             Output.(L).highIBI_linear.yfit            = single(nan(1,length(bin_ids)));
%             Output.(L).highIBI_linear.coefs           = single([NaN NaN]);
%             Output.(L).highIBI_linear.rsquared        = single(NaN);
%             Output.(L).highIBI_linear.adjrsquared     = single(NaN);
%             Output.(L).highIBI_linear.sse             = single(NaN);
%             Output.(L).highIBI_linear.dfe             = single(NaN);
%             Output.(L).highIBI_linear.rmse            = single(NaN);
%             Output.(L).highIBI_linear.pvalue          = single([NaN; NaN]);
%             Output.(L).highIBI_linear.aic             = single(NaN);
%             Output.(L).highIBI_linear.bic             = single(NaN);
%             
%             Output.(L).highIBI_cosine.average         = single(nan(1,length(bin_ids)));
%             Output.(L).highIBI_cosine.startPoint      = single([NaN NaN NaN]);
%             Output.(L).highIBI_cosine.yfit            = single(nan(1,length(bin_ids)));
%             Output.(L).highIBI_cosine.coefs           = single([NaN NaN NaN]);
%             Output.(L).highIBI_cosine.rsquared        = single(NaN);
%             Output.(L).highIBI_cosine.adjrsquared     = single(NaN);
%             Output.(L).highIBI_cosine.sse             = single(NaN);
%             Output.(L).highIBI_cosine.dfe             = single(NaN);
%             Output.(L).highIBI_cosine.rmse            = single(NaN);
%             Output.(L).highIBI_cosine.pvalue          = single(NaN);
%             Output.(L).highIBI_cosine.aic             = single(NaN);
%             Output.(L).highIBI_cosine.bic             = single(NaN);
%             
%             Output.(L).highIBI_vonMisesPos.average    = single(nan(1,length(bin_ids)));
%             Output.(L).highIBI_vonMisesPos.startPoint = single([NaN NaN NaN NaN]);
%             Output.(L).highIBI_vonMisesPos.yfit       = single(nan(1,length(bin_ids)));
%             Output.(L).highIBI_vonMisesPos.coefs      = single([NaN NaN NaN NaN]);
%             Output.(L).highIBI_vonMisesPos.rsquared   = single(NaN);
%             Output.(L).highIBI_vonMisesPos.adjrsquared= single(NaN);
%             Output.(L).highIBI_vonMisesPos.sse        = single(NaN);
%             Output.(L).highIBI_vonMisesPos.dfe        = single(NaN);
%             Output.(L).highIBI_vonMisesPos.rmse       = single(NaN);
%             Output.(L).highIBI_vonMisesPos.pvalue     = single(NaN);
%             Output.(L).highIBI_vonMisesPos.aic        = single(NaN);
%             Output.(L).highIBI_vonMisesPos.bic        = single(NaN);
%             
%             Output.(L).highIBI_vonMisesNeg.average    = single(nan(1,length(bin_ids)));
%             Output.(L).highIBI_vonMisesNeg.startPoint = single([NaN NaN NaN NaN]);
%             Output.(L).highIBI_vonMisesNeg.yfit       = single(nan(1,length(bin_ids)));
%             Output.(L).highIBI_vonMisesNeg.coefs      = single([NaN NaN NaN NaN]);
%             Output.(L).highIBI_vonMisesNeg.rsquared   = single(NaN);
%             Output.(L).highIBI_vonMisesNeg.adjrsquared= single(NaN);
%             Output.(L).highIBI_vonMisesNeg.sse        = single(NaN);
%             Output.(L).highIBI_vonMisesNeg.dfe        = single(NaN);
%             Output.(L).highIBI_vonMisesNeg.rmse       = single(NaN);
%             Output.(L).highIBI_vonMisesNeg.pvalue     = single(NaN);
%             Output.(L).highIBI_vonMisesNeg.aic        = single(NaN);
%             Output.(L).highIBI_vonMisesNeg.bic        = single(NaN);
            
        end
        
        for c=1:numel(cfg.condition)
            L=cfg.condition(c).name;
            %% get condition AND valid block trials only
            CT = ecg_bna_get_condition_trials(T, cfg.condition(c));
            tr=ismember([T.block],blocks) & CT;
            popcell=num2cell(pop.trial(tr));
            trcell=num2cell(T(tr));
            
            % add trial onset time to each spike so its basically one stream again
            % also, make sure spikes aren't counted twice (because previous trial is appended in beginning;
            % removing overlapping spikes here            % add trial onset time         % add block separator
            arrival_times=cellfun(@(x,y) y.arrival_times(y.arrival_times>x.states_onset(x.states==2)) + x.TDT_ECG1_t0_from_rec_start+offset_blocks_Rpeak(Rblocks==x.block),trcell,popcell,'uniformoutput',false);
            trial_onsets=cellfun(@(x) x.TDT_ECG1_t0_from_rec_start+offset_blocks_Rpeak(Rblocks==x.block),trcell);
            trial_ends=cellfun(@(x) x.states_onset(x.states==90)+x.TDT_ECG1_t0_from_rec_start+offset_blocks_Rpeak(Rblocks==x.block),trcell,'uniformoutput',false); % no clue why this needs to be nonuniformoutput, it did work earlier so this is confusing...
            trial_ends=[trial_ends{:}];
            
            if numel(trial_onsets)<=1 || (~isfield(Rpeaks, 'RPEAK_ts_insp') && cfg.process_Rpeaks_inhalation_exhalation) || (~isfield(Rpeaks, 'RPEAK_ts_exp') && cfg.process_Rpeaks_inhalation_exhalation)
                continue; % out(1).nrblock_combinedFiles might be empty! and even if there is 1 trial we're not processing
            end
            
            %% compute spike density as one continuous vector across all concatenated trials (hmmm there migth be a problem with interleaved trial types here)
            AT=vertcat(arrival_times{:});
            AT(AT>trial_ends(end))=[];
            [SD_all_trials,RAST,PSTH_time,SD_1ms,RAST_1ms,PSTH_time_1ms]=ecg_bna_spike_density(AT,trial_onsets,trial_ends,cfg.time);
            
            %% retrieve R-peaks and their reshuffles
            [RPEAK_ts,RPEAK_ts_perm,RPEAK_dur,RPEAK_dur_perm] = retrieve_Rpeaks_and_reshuffles(Rpeaks,b,cfg,c);
            
            %% define which parts of the continous PSTH are during a trial
            
            during_trial_index     = ecg_bna_define_during_trial_index(trial_onsets, trial_ends, PSTH_time, curr_analyse_states, cfg.time.PSTH_binwidth);
            during_trial_index_1ms = ecg_bna_define_during_trial_index(trial_onsets, trial_ends, PSTH_time_1ms, curr_analyse_states, 0.001); % for rasters with 1-ms bins
            
            realPSTHs     = compute_PSTH(RPEAK_ts,RPEAK_dur,RAST,SD_all_trials,PSTH_time,during_trial_index,curr_analyse_states,cfg.time.PSTH_binwidth);
            realPSTHs_1ms = compute_PSTH(RPEAK_ts,RPEAK_dur,RAST_1ms,SD_1ms,PSTH_time_1ms,during_trial_index_1ms,curr_analyse_states,0.001); % for rasters with 1-ms bins
            shuffledPSTH  = compute_PSTH(RPEAK_ts_perm,RPEAK_dur_perm,RAST,SD_all_trials,PSTH_time,during_trial_index,curr_analyse_states,cfg.time.PSTH_binwidth);
            SD            = do_statistics(realPSTHs,shuffledPSTH,BINS,cfg.time.significance_window{numTiming});
            
            Output.(L).SD                          = SD.SD_mean ;
            Output.(L).SD_STD                      = SD.SD_STD;
            Output.(L).SD_SEM                      = SD.SD_SEM ;
            Output.(L).SDCL                        = SD.SDconf(1,:);
            Output.(L).SDCu                        = SD.SDconf(2,:);
            Output.(L).SDP                         = SD.SDPmean ;
            Output.(L).SDPCL                       = SD.SDPconf(1,:) ;
            Output.(L).SDPCu                       = SD.SDPconf(2,:) ;
            Output.(L).sig_all                     = SD.sig_all;
            Output.(L).sig                         = SD.sig;
            Output.(L).sig_FR_diff                 = SD.sig_FR_diff;
            Output.(L).sig_time                    = SD.sig_time;
            Output.(L).sig_n_bins                  = SD.sig_n_bins;
            Output.(L).sig_sign                    = SD.sig_sign;
            Output.(L).NrTrials                    = sum(tr);
            Output.(L).NrEvents                    = realPSTHs.n_events;
            Output.(L).FR                          = mean(SD_all_trials); %% not too sure this was the intended one...
            Output.(L).raster                      = logical(realPSTHs_1ms.raster); % logical replaces all numbers >0 with 1 and reduces memory load
            %Output.(L).Rts                        = single(realPSTHs.RTs{1}); % unless we need this, dont save it!
            Output.(L).Rds                         = hist(realPSTHs.RDs{1},cfg.time.histbins); % put RR durations to plot those in the histograms later
            Output.(L).Rds_perm                    = hist([shuffledPSTH.RDs{:}],cfg.time.histbins);
            Output.(L).SDsubstractedSDP            = Output.(L).SD - Output.(L).SDP; % spikes/s, difference between mean and jittered data
            Output.(L).SDsubstractedSDP_normalized = Output.(L).SDsubstractedSDP ./ Output.(L).SDP *100; % percent signal change
            Output.(L).FR_ModIndex_SubtrSDP        = max(Output.(L).SDsubstractedSDP) - min(Output.(L).SDsubstractedSDP); % difference between max and min FR
            Output.(L).FR_ModIndex_PcS             = max(Output.(L).SDsubstractedSDP_normalized) - min(Output.(L).SDsubstractedSDP_normalized); % difference between max and min % signal change
            
            % implement median split to heart-cycle durations
            Output.(L).IBI_median           = median(realPSTHs.RDs{1});
            
            % reshuffle based on the computed median
            % lowIBI - set up settings
            cfg.time.IBI       = 1;
            cfg.time.IBI_low   = 1;
            cfg.time.IBI_high  = 0;
            cfg.time.IBI_thrsh = Output.(L).IBI_median * ones(1, length(cfg.condition));
            Rpeaks_lowIBI      = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info,cfg.time);
            for XXX=1:numel(Rpeaks_lowIBI)
                Rpeaks_lowIBI(XXX).RPEAK_ts=Rpeaks_lowIBI(XXX).RPEAK_ts-Rpeaks_lowIBI(XXX).offset+Rpeaks(XXX).offset;
                Rpeaks_lowIBI(XXX).shuffled_ts=Rpeaks_lowIBI(XXX).shuffled_ts-Rpeaks_lowIBI(XXX).offset+Rpeaks(XXX).offset;
                Rpeaks_lowIBI(XXX).offset=Rpeaks(XXX).offset;
            end
            [RPEAK_ts_lowIBI,RPEAK_ts_perm_lowIBI,RPEAK_dur_lowIBI,RPEAK_dur_perm_lowIBI] = ...
                retrieve_Rpeaks_and_reshuffles(Rpeaks_lowIBI, b, cfg ,c);
            
            % highIBI - set up settings
            cfg.time.IBI       = 1;
            cfg.time.IBI_low   = 0;
            cfg.time.IBI_high  = 1;
            cfg.time.IBI_thrsh = Output.(L).IBI_median * ones(1, length(cfg.condition));
            Rpeaks_highIBI      = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info,cfg.time);
            for XXX=1:numel(Rpeaks_highIBI)
                Rpeaks_highIBI(XXX).RPEAK_ts=Rpeaks_highIBI(XXX).RPEAK_ts-Rpeaks_highIBI(XXX).offset+Rpeaks(XXX).offset;
                Rpeaks_highIBI(XXX).shuffled_ts=Rpeaks_highIBI(XXX).shuffled_ts-Rpeaks_highIBI(XXX).offset+Rpeaks(XXX).offset;
                Rpeaks_highIBI(XXX).offset=Rpeaks(XXX).offset;
            end
            [RPEAK_ts_highIBI,RPEAK_ts_perm_highIBI,RPEAK_dur_highIBI,RPEAK_dur_perm_highIBI] = ...
                retrieve_Rpeaks_and_reshuffles(Rpeaks_highIBI, b, cfg, c);
            
            % !!! IMPORTANT !!! erase settings related to median split
            % reshuffles and re-initiate them for the next unit
            cfg.time = rmfield(cfg.time, {'IBI', 'IBI_low', 'IBI_high', 'IBI_thrsh'});
            
            lowIBIids  = realPSTHs.RDs{1} < Output.(L).IBI_median;
            highIBIids = realPSTHs.RDs{1} > Output.(L).IBI_median;
            
            % check the number of spikes left after computing phases
            if realPSTHs.n_events < 100 || sum(lowIBIids) < 3 || sum(highIBIids) < 3
                continue
            end
            
%             Output.(L).linear               = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(:, bin_ids)', realPSTHs.n_events, 'linear');
%             Output.(L).cosine               = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(:, bin_ids)', realPSTHs.n_events, 'cosine');
%             Output.(L).cosine.coefs(2)      = Output.(L).cosine.coefs(2)+min(Output.phase_bin_centers);
%             Output.(L).vonMisesPos          = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(:, bin_ids)', realPSTHs.n_events, 'vonMises', 1);
%             Output.(L).vonMisesPos.coefs(4) = Output.(L).vonMisesPos.coefs(4)+min(Output.phase_bin_centers);
%             Output.(L).vonMisesNeg          = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(:, bin_ids)', realPSTHs.n_events, 'vonMises', -1);
%             Output.(L).vonMisesNeg.coefs(4) = Output.(L).vonMisesNeg.coefs(4)+min(Output.phase_bin_centers);
           
            
            % compute time-domain analysis for lowIBI
            realPSTHs_lowIBI     = compute_PSTH(RPEAK_ts_lowIBI, RPEAK_dur_lowIBI, RAST, SD_all_trials, PSTH_time, during_trial_index, curr_analyse_states, cfg.time.PSTH_binwidth);
            shuffledPSTH_lowIBI  = compute_PSTH(RPEAK_ts_perm_lowIBI,RPEAK_dur_perm_lowIBI,RAST,SD_all_trials,PSTH_time,during_trial_index,curr_analyse_states,cfg.time.PSTH_binwidth);
            SD_lowIBI            = do_statistics(realPSTHs_lowIBI,shuffledPSTH_lowIBI,BINS,cfg.time.significance_window{numTiming});
            
            Output.(L).lowIBI.SD          = SD_lowIBI.SD_mean;
            Output.(L).lowIBI.SD_STD      = SD_lowIBI.SD_STD;
            Output.(L).lowIBI.SD_SEM      = SD_lowIBI.SD_SEM ;
            Output.(L).lowIBI.SDCL        = SD_lowIBI.SDconf(1,:);
            Output.(L).lowIBI.SDCu        = SD_lowIBI.SDconf(2,:);
            Output.(L).lowIBI.SDP         = SD_lowIBI.SDPmean ;
            Output.(L).lowIBI.SDPCL       = SD_lowIBI.SDPconf(1,:) ;
            Output.(L).lowIBI.SDPCu       = SD_lowIBI.SDPconf(2,:) ;
            Output.(L).lowIBI.sig_all     = SD_lowIBI.sig_all;
            Output.(L).lowIBI.sig         = SD_lowIBI.sig;
            Output.(L).lowIBI.sig_FR_diff = SD_lowIBI.sig_FR_diff;
            Output.(L).lowIBI.sig_time    = SD_lowIBI.sig_time;
            Output.(L).lowIBI.sig_n_bins  = SD_lowIBI.sig_n_bins;
            Output.(L).lowIBI.sig_sign    = SD_lowIBI.sig_sign;
            Output.(L).lowIBI.NrEvents    = realPSTHs_lowIBI.n_events;
            Output.(L).lowIBI.SDsubstractedSDP            = Output.(L).lowIBI.SD - Output.(L).lowIBI.SDP; % spikes/s, difference between mean and jittered data
            Output.(L).lowIBI.SDsubstractedSDP_normalized = Output.(L).lowIBI.SDsubstractedSDP ./ Output.(L).lowIBI.SDP *100; % percent signal change
            Output.(L).lowIBI.FR_ModIndex_SubtrSDP        = max(Output.(L).lowIBI.SDsubstractedSDP) - min(Output.(L).lowIBI.SDsubstractedSDP); % difference between max and min FR
            Output.(L).lowIBI.FR_ModIndex_PcS             = max(Output.(L).lowIBI.SDsubstractedSDP_normalized) - min(Output.(L).lowIBI.SDsubstractedSDP_normalized); % difference between max and min % signal change
        
            % compute time-domain analysis for highIBI
            realPSTHs_highIBI     = compute_PSTH(RPEAK_ts_highIBI,RPEAK_dur_highIBI,RAST,SD_all_trials,PSTH_time,during_trial_index,curr_analyse_states,cfg.time.PSTH_binwidth);
            shuffledPSTH_highIBI  = compute_PSTH(RPEAK_ts_perm_highIBI,RPEAK_dur_perm_highIBI,RAST,SD_all_trials,PSTH_time,during_trial_index,curr_analyse_states,cfg.time.PSTH_binwidth);
            SD_highIBI            = do_statistics(realPSTHs_highIBI,shuffledPSTH_highIBI,BINS,cfg.time.significance_window{numTiming});
            
            Output.(L).highIBI.SD          = SD_highIBI.SD_mean;
            Output.(L).highIBI.SD_STD      = SD_highIBI.SD_STD;
            Output.(L).highIBI.SD_SEM      = SD_highIBI.SD_SEM ;
            Output.(L).highIBI.SDCL        = SD_highIBI.SDconf(1,:);
            Output.(L).highIBI.SDCu        = SD_highIBI.SDconf(2,:);
            Output.(L).highIBI.SDP         = SD_highIBI.SDPmean ;
            Output.(L).highIBI.SDPCL       = SD_highIBI.SDPconf(1,:) ;
            Output.(L).highIBI.SDPCu       = SD_highIBI.SDPconf(2,:) ;
            Output.(L).highIBI.sig_all     = SD_highIBI.sig_all;
            Output.(L).highIBI.sig         = SD_highIBI.sig;
            Output.(L).highIBI.sig_FR_diff = SD_highIBI.sig_FR_diff;
            Output.(L).highIBI.sig_time    = SD_highIBI.sig_time;
            Output.(L).highIBI.sig_n_bins  = SD_highIBI.sig_n_bins;
            Output.(L).highIBI.sig_sign    = SD_highIBI.sig_sign;
            Output.(L).highIBI.NrEvents    = realPSTHs_highIBI.n_events;
            Output.(L).highIBI.SDsubstractedSDP            = Output.(L).highIBI.SD - Output.(L).highIBI.SDP; % spikes/s, difference between mean and jittered data
            Output.(L).highIBI.SDsubstractedSDP_normalized = Output.(L).highIBI.SDsubstractedSDP ./ Output.(L).highIBI.SDP *100; % percent signal change
            Output.(L).highIBI.FR_ModIndex_SubtrSDP        = max(Output.(L).highIBI.SDsubstractedSDP) - min(Output.(L).highIBI.SDsubstractedSDP); % difference between max and min FR
            Output.(L).highIBI.FR_ModIndex_PcS             = max(Output.(L).highIBI.SDsubstractedSDP_normalized) - min(Output.(L).highIBI.SDsubstractedSDP_normalized); % difference between max and min % signal change
            
            % low IBI
%             Output.(L).lowIBI_linear               = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(lowIBIids, bin_ids)', sum(lowIBIids), 'linear');
            
%             Output.(L).lowIBI_linear_smoothed      = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.mean(bin_ids)', sum(lowIBIids), 'linear');
            
%             Output.(L).lowIBI_cosine               = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(lowIBIids, bin_ids)', sum(lowIBIids), 'cosine');
%             Output.(L).lowIBI_cosine.coefs(2)      = Output.(L).lowIBI_cosine.coefs(2)+min(Output.phase_bin_centers);
            
            
            
%             Output.(L).lowIBI_vonMisesPos          = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(lowIBIids, bin_ids)', sum(lowIBIids), 'vonMises', 1);
%             Output.(L).lowIBI_vonMisesPos.coefs(4) = Output.(L).lowIBI_vonMisesPos.coefs(4)+min(Output.phase_bin_centers);
            
%             Output.(L).lowIBI_vonMisesNeg          = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(lowIBIids, bin_ids)', sum(lowIBIids), 'vonMises', -1);
%             Output.(L).lowIBI_vonMisesNeg.coefs(4) = Output.(L).lowIBI_vonMisesNeg.coefs(4)+min(Output.phase_bin_centers);
            
            % high IBI
%             Output.(L).highIBI_linear               = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(highIBIids, bin_ids)', sum(highIBIids), 'linear');
            
%             Output.(L).highIBI_cosine               = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(highIBIids, bin_ids)', sum(highIBIids), 'cosine');
%             Output.(L).highIBI_cosine.coefs(2)      = Output.(L).highIBI_cosine.coefs(2)+min(Output.phase_bin_centers);
            
%             Output.(L).highIBI_vonMisesPos          = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(highIBIids, bin_ids)', sum(highIBIids), 'vonMises', 1);
%             Output.(L).highIBI_vonMisesPos.coefs(4) = Output.(L).highIBI_vonMisesPos.coefs(4)+min(Output.phase_bin_centers);
            
%             Output.(L).highIBI_vonMisesNeg          = ecg_bna_fit_neuronal_data(cfg, Output.phase_bin_centers-min(Output.phase_bin_centers), realPSTHs.raster(highIBIids, bin_ids)', sum(highIBIids), 'vonMises', -1);
%             Output.(L).highIBI_vonMisesNeg.coefs(4) = Output.(L).highIBI_vonMisesNeg.coefs(4)+min(Output.phase_bin_centers);
            
            clear realPSTHs shuffledPSTH SD
            
            %% The part following here is internal sanity check and should be turned off in general since there typically is no ECG data in the spike format
            if Sanity_check %% this needs to be fixed as well, this might be incorrect after least update...
                ECG_data=[pop.trial(tr).TDT_ECG1];
                ECG_time=cellfun(@(x) [1/x.TDT_ECG1_SR:1/x.TDT_ECG1_SR:numel(x.TDT_ECG1)/x.TDT_ECG1_SR]+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart+offset_blocks_Rpeak(Rblocks==x.block),trcell,'uniformoutput',false);
                CG_time=[ECG_time{:}];
                tt=0;
                clear ECG_to_plot
                for t=1:numel(RPEAK_ts)
                    ECg_t_idx=CG_time>RPEAK_ts(t)-0.5 & CG_time<RPEAK_ts(t)+0.5;
                    if round(sum(ECg_t_idx)-trcell{1}.TDT_ECG1_SR)~=0
                        continue
                    end
                    tt=tt+1;
                    ECG_to_plot(tt,:)=ECG_data(ECg_t_idx);
                end
                
                figure
                filename=[unit_ID '_block_' num2str(block) '_ECG_average'];
                lineProps={'color','b','linewidth',1};
                shadedErrorBar(1:size(ECG_to_plot,2),mean(ECG_to_plot,1),sterr(ECG_to_plot,1),lineProps,1);
                export_fig([cfg.per_session_folder, filesep, filename], '-pdf','-transparent') % pdf by run
                close(gcf);
                
                figure
                filename=[unit_ID '_block_' num2str(block) '_ECG_per5trials'];
                plot(ECG_to_plot(1:5:end,:)');
                export_fig([cfg.per_session_folder, filesep, filename], '-pdf','-transparent') % pdf by run
                close(gcf);
            end
        end
        %% save output
        output_folder = [cfg.per_session_folder '_' num2str(curr_analyse_states{3}) '-' num2str(curr_analyse_states{4}) 's'];
        if ~exist(output_folder,'dir')
            mkdir(output_folder)
        end
        save([output_folder, filesep, pop.unit_ID, '_', pop.target],'Output')
        clear Output
        toc
    end
end
end

function during_trial_index = ecg_bna_define_during_trial_index(trial_onsets, trial_ends, PSTH_time, curr_analyse_states, PSTH_binwidth)

trial_onset_samples=ceil((trial_onsets-PSTH_time(1))/PSTH_binwidth);
trial_ends_samples=floor((trial_ends-PSTH_time(1))/PSTH_binwidth);
trial_onset_samples(trial_onset_samples==0)=1;
during_trial_index=false(size(PSTH_time));
drop_samples = trial_onset_samples < 1 | trial_ends_samples < 1; % in very rare cases samples have negative values, drop those
trial_onset_samples = trial_onset_samples(~drop_samples);
trial_ends_samples = trial_ends_samples(~drop_samples);

bins_window_start_relative_to_trigger=curr_analyse_states{3}/PSTH_binwidth;
bins_window_end_relative_to_trigger=curr_analyse_states{4}/PSTH_binwidth;
for t=1:numel(trial_onset_samples)
    during_trial_index(trial_onset_samples(t)-bins_window_start_relative_to_trigger:trial_ends_samples(t)-bins_window_end_relative_to_trigger)=true;
end

end

function out=compute_PSTH(RPEAK_ts,RPEAK_dur,RAST,SD,PSTH_time,during_trial_index,curr_analyse_states,PSTH_binwidth)
RPEAK_samples=round((RPEAK_ts-PSTH_time(1))/PSTH_binwidth);
bins_before=round(curr_analyse_states{3}/PSTH_binwidth);
bins_after=round(curr_analyse_states{4}/PSTH_binwidth);
bins=bins_before:bins_after;

%% remove samples that would land outside
rpeaks2exclude = ...
    RPEAK_samples>=numel(SD)-bins_after | RPEAK_samples<=-bins_before;

RPEAK_ts(rpeaks2exclude)       = NaN;
RPEAK_samples(rpeaks2exclude)  = NaN;
RPEAK_dur(rpeaks2exclude)      = NaN;

%% preallocate "out" structure
out.raster     = NaN;
out.PSTH       = NaN;
out.mean       = nan(1,numel(bins));
out.std        = nan(1,numel(bins));
out.SEM        = nan(1,numel(bins));
out.n_events   = NaN;
out.RTs        = {};
out.RDs        = {};

%% loop through rows of RPEAK_samples: 1 row for real, nReshuffles rows of reshuffled data
for p=1:size(RPEAK_samples,1)
    RT=RPEAK_ts(p,~isnan(RPEAK_samples(p,:)));      % take time-stamps
    RS=RPEAK_samples(p,~isnan(RPEAK_samples(p,:))); % take samples
    RD=RPEAK_dur(p,~isnan(RPEAK_samples(p,:)));     % take durations
    
    within_trial_idx = during_trial_index(RS);
    
    RT=RT(within_trial_idx);
    RS=RS(within_trial_idx);
    RD=RD(within_trial_idx);
    
    idx_by_trial     = bsxfun(@plus,RS',bins); % bsxfun produces a RPEAK-(nBins-1) matrix with samples taken from SDF for each R-peak
    if size(RPEAK_samples,1) == 1
        out.raster     = RAST(idx_by_trial);
    end
    out.PSTH = SD(idx_by_trial);
    out.mean(p,:)=mean(out.PSTH, 1);
    out.std(p,:) = std(out.PSTH, [], 1);
    out.SEM(p,:)=sterr(out.PSTH);
    out.n_events(p)=numel(RS);
    out.RTs{p}=RT;
    out.RDs{p}=RD;
end
end

function [SD,RAST,PSTH_time,SD_1ms,RAST_1ms,PSTH_time_1ms]=ecg_bna_spike_density(AT,trial_onsets,trial_ends,cfg)
%% make SD across all trials appended (no average)!
switch cfg.kernel_type
    case 'gaussian'
        Kernel     = normpdf(-5*cfg.gaussian_kernel:cfg.PSTH_binwidth:5*cfg.gaussian_kernel,0,cfg.gaussian_kernel);
    case 'box'
        n_bins=round(2*cfg.gaussian_kernel/cfg.PSTH_binwidth);
        Kernel=ones(1,n_bins)/n_bins/cfg.PSTH_binwidth; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
end
PSTH_time     = trial_onsets(1):cfg.PSTH_binwidth:trial_ends(end);
PSTH_time_1ms = trial_onsets(1):0.001:trial_ends(end);
RAST          = hist(AT,PSTH_time);
RAST_1ms      = hist(AT,PSTH_time_1ms);
SD            = conv(RAST,Kernel,'same');
SD_1ms        = conv(RAST_1ms,Kernel,'same');
end

function [SD,BINS]=do_statistics(Real,Shuffled,BINS,significance_window)
definition='max';

SD_mean=Real.mean;
SDconf(1,:)=abs(prctile(Real.PSTH,2.5,1)-SD_mean);
SDconf(2,:)=abs(prctile(Real.PSTH,97.5,1)-SD_mean);

SDPmean=nanmean(Shuffled.mean,1);
SDPconf(1,:)=abs(prctile(Shuffled.mean,2.5,1)-SDPmean);
SDPconf(2,:)=abs(prctile(Shuffled.mean,97.5,1)-SDPmean);

SD.SD_mean=Real.mean;
SD.SD_STD=Real.std;
SD.SD_SEM=Real.SEM;
SD.SDPmean=SDPmean;
SD.SDPconf=SDPconf;
SD.SDconf = SDconf;

%% signficance
sig_to_check=BINS>significance_window(1)*1000 & BINS<significance_window(2)*1000;
pos_diff=SD_mean-(SDPmean+SDPconf(2,:));
neg_diff=SD_mean-(SDPmean-SDPconf(1,:));
sig_above=pos_diff>0&sig_to_check;
sig_below=neg_diff<0&sig_to_check;

sig_idx_above_start=find(diff([0 sig_above 0])>0);
sig_idx_above_end=find(diff([0 sig_above 0])<0);
sig_idx_below_start=find(diff([0 sig_below 0])>0);
sig_idx_below_end=find(diff([0 sig_below 0])<0);

% find longest period of significance
[ma,m_ia]=max(sig_idx_above_end-sig_idx_above_start);
[mb,m_ib]=max(sig_idx_below_end-sig_idx_below_start);
% find maximum deviation from surrogates
[max_pos_diff,max_idx]=max(pos_diff.*sig_to_check);
[max_neg_diff,min_idx]=min(neg_diff.*sig_to_check);

m_imax=find(sig_idx_above_start<=max_idx & sig_idx_above_end>=max_idx);
m_imin=find(sig_idx_below_start<=min_idx & sig_idx_below_end>=min_idx);

if strcmp(definition,'max')
    if abs(max_pos_diff)>abs(max_neg_diff)
        sig_sign=1;
        m_i=m_imax;
        sig_start_end=[sig_idx_above_start(m_i):sig_idx_above_end(m_i)-1];
    else
        sig_sign=-1;
        m_i=m_imin;
        sig_start_end=[sig_idx_below_start(m_i):sig_idx_below_end(m_i)-1];
    end
elseif strcmp(definition,'dur')
    if ma>mb
        sig_sign=1;
        m_i=m_ia;
        sig_start_end=[sig_idx_above_start(m_i):sig_idx_above_end(m_i)-1];
    else
        sig_sign=-1;
        m_i=m_ib;
        sig_start_end=[sig_idx_below_start(m_i):sig_idx_below_end(m_i)-1];
    end
end

m=numel(sig_start_end);

if m==0
    sig_sign=0;
    t_start_end=[NaN NaN];
    maxidx=1;
    max_FR_diff=NaN;
else
    [max_FR_diff, maxidx]=max(abs(SD_mean(sig_start_end)-SDPmean(sig_start_end)));
    t_start_end=BINS(sig_start_end);
end
sig_period=BINS>=t_start_end(1) & BINS<=t_start_end(end) ;
max_time=t_start_end(maxidx);

SD.sig_all=sig_above-sig_below;
SD.sig=(sig_above&sig_period)-(sig_below&sig_period);
SD.sig_FR_diff=max_FR_diff;
SD.sig_time=max_time;
SD.sig_n_bins=m;
SD.sig_sign=sig_sign;
end

function out = select2D_cat_nans(input_array,idx)

cell_input_array = mat2cell(input_array, ones(size(input_array,1),1), size(input_array,2));
cell_idx         = mat2cell(idx, ones(size(idx,1),1), size(idx,2));
cell_selected    = cellfun(@(x,y) x(y), cell_input_array, cell_idx, 'UniformOutput', false);
len_selected     = cellfun(@length, cell_selected, 'UniformOutput', true);
max_len          = max(len_selected);
elem2add         = max_len - len_selected;
elem2add         = mat2cell(elem2add, ones(length(elem2add),1), 1);
nan2add          = cellfun(@(x) nan(1,x), elem2add, 'UniformOutput', false);
out              = cellfun(@(x,y) cat(2,x,y), cell_selected, nan2add, 'UniformOutput', false);
out              = cat(1,out{:});

end

function [RPEAK_ts,RPEAK_ts_perm,RPEAK_dur,RPEAK_dur_perm] = retrieve_Rpeaks_and_reshuffles(Rpeaks,b,cfg,c)

if contains(cfg.condition(c).Rpeak_field, 'extrasystole')
    RPEAK_ts  = [Rpeaks(b).RPEAK_ts_extrasystole];
    RPEAK_dur = [Rpeaks(b).RPEAK_dur_extrasystole];
else
    RPEAK_ts=[Rpeaks(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])];
    RPEAK_ts_perm=[Rpeaks(b).(['shuffled_ts' cfg.condition(c).Rpeak_field])];
    RPEAK_dur = [Rpeaks(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])];
    RPEAK_dur_perm = [Rpeaks(b).(['shuffled_dur' cfg.condition(c).Rpeak_field])];
end

end
