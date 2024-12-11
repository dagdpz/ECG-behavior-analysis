function ecg_bna_compute_session_ECG_related_spikePhase(trials,population,Rpeaks,sessions_info,cfg)

basepath_to_save=[cfg.SPK_root_results_fldr filesep 'cardioballistic'];
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
    b=ismember(Rblocks,blocks);
    
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
    data.criteria          = pop.criteria;
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        data.(L).meanHR_bpm                     = single(NaN);
        data.(L).medianHR_bpm                   = single(NaN);
        data.(L).stdHR_bpm                      = single(NaN);
        data.(L).SDNN_ms                        = single(NaN);
        data.(L).spike_phases_radians           = single(NaN);
        data.(L).spike_phases_histogram         = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).spike_phases_histogram2        = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).spike_phases_histogram_smoothed= single(nan(1,cfg.phase.N_phase_bins));
        data.(L).histogram_MI                   = single(NaN);
        data.(L).histogram_p                    = single(NaN);
        data.(L).histogram_phase                = single(NaN);
        data.(L).rsquared                       = single(NaN);
        
        data.(L).lowIBI_spike_phases_histogram_smoothed  = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI_spike_phases_histogram_smoothed = single(nan(1,cfg.phase.N_phase_bins));
        
        data.(L).SD                = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).SD_STD            = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).SD_SEM            = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).SDP               = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).SDPCL             = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).SDPCu             = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).sig_all           = single(zeros(1,cfg.phase.N_phase_bins));
        data.(L).sig               = single(zeros(1,cfg.phase.N_phase_bins));
        data.(L).sig_FR_diff       = single(nan(1));
        data.(L).sig_time          = single(nan(1));
        data.(L).sig_n_bins        = single(zeros(1));
        data.(L).sig_sign          = single(zeros(1));
        data.(L).SDsubstractedSDP            = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).SDsubstractedSDP_normalized = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).FR_ModIndex_SubtrSDP        = single(nan(1));
        data.(L).FR_ModIndex_PcS             = single(nan(1));
        
        data.(L).lowIBI.SD          = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.SD_STD      = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.SD_SEM      = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.SDP         = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.SDPCL       = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.SDPCu       = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.sig_all     = single(zeros(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.sig         = single(zeros(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.sig_FR_diff = single(nan(1));
        data.(L).lowIBI.sig_time    = single(nan(1));
        data.(L).lowIBI.sig_n_bins  = single(zeros(1));
        data.(L).lowIBI.sig_sign    = single(zeros(1));
        data.(L).lowIBI.SDsubstractedSDP            = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.SDsubstractedSDP_normalized = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).lowIBI.FR_ModIndex_SubtrSDP        = single(nan(1));
        data.(L).lowIBI.FR_ModIndex_PcS             = single(nan(1));
        
        data.(L).highIBI.SD          = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.SD_STD      = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.SD_SEM      = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.SDP         = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.SDPCL       = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.SDPCu       = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.sig_all     = single(zeros(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.sig         = single(zeros(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.sig_FR_diff = single(nan(1));
        data.(L).highIBI.sig_time    = single(nan(1));
        data.(L).highIBI.sig_n_bins  = single(zeros(1));
        data.(L).highIBI.sig_sign    = single(zeros(1));
        data.(L).highIBI.SDsubstractedSDP            = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.SDsubstractedSDP_normalized = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).highIBI.FR_ModIndex_SubtrSDP        = single(nan(1));
        data.(L).highIBI.FR_ModIndex_PcS             = single(nan(1));
        
%         data.(L).linear.yfit                    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).linear.coefs                   = single([NaN NaN]);
%         data.(L).linear.rsquared                = single(NaN);
%         data.(L).linear.adjrsquared             = single(NaN);
%         data.(L).linear.sse                     = single(NaN);
%         data.(L).linear.dfe                     = single(NaN);
%         data.(L).linear.rmse                    = single(NaN);
%         data.(L).linear.pvalue                  = single([NaN; NaN]);
%         data.(L).linear.aic                     = single(NaN);
%         data.(L).linear.bic                     = single(NaN);
        
%         data.(L).linear_smoothed.yfit                    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).linear_smoothed.coefs                   = single([NaN NaN]);
%         data.(L).linear_smoothed.rsquared                = single(NaN);
%         data.(L).linear_smoothed.adjrsquared             = single(NaN);
%         data.(L).linear_smoothed.sse                     = single(NaN);
%         data.(L).linear_smoothed.dfe                     = single(NaN);
%         data.(L).linear_smoothed.rmse                    = single(NaN);
%         data.(L).linear_smoothed.pvalue                  = single([NaN; NaN]);
%         data.(L).linear_smoothed.aic                     = single(NaN);
%         data.(L).linear_smoothed.bic                     = single(NaN);
        
%         data.(L).cosine.average                 = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).cosine.startPoint              = single([NaN NaN NaN]);
%         data.(L).cosine.yfit                    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).cosine.coefs                   = single([NaN NaN NaN]);
%         data.(L).cosine.rsquared                = single(NaN);
%         data.(L).cosine.adjrsquared             = single(NaN);
%         data.(L).cosine.sse                     = single(NaN);
%         data.(L).cosine.dfe                     = single(NaN);
%         data.(L).cosine.rmse                    = single(NaN);
%         data.(L).cosine.pvalue                  = single(NaN);
%         data.(L).cosine.aic                     = single(NaN);
%         data.(L).cosine.bic                     = single(NaN);
        
%         data.(L).cosine_smoothed.average                 = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).cosine_smoothed.startPoint              = single([NaN NaN NaN]);
%         data.(L).cosine_smoothed.yfit                    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).cosine_smoothed.coefs                   = single([NaN NaN NaN]);
%         data.(L).cosine_smoothed.rsquared                = single(NaN);
%         data.(L).cosine_smoothed.adjrsquared             = single(NaN);
%         data.(L).cosine_smoothed.sse                     = single(NaN);
%         data.(L).cosine_smoothed.dfe                     = single(NaN);
%         data.(L).cosine_smoothed.rmse                    = single(NaN);
%         data.(L).cosine_smoothed.pvalue                  = single(NaN);
%         data.(L).cosine_smoothed.aic                     = single(NaN);
%         data.(L).cosine_smoothed.bic                     = single(NaN);
        
%         data.(L).vonMisesPos.average            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesPos.startPoint         = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesPos.yfit               = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesPos.coefs              = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesPos.rsquared           = single(NaN);
%         data.(L).vonMisesPos.adjrsquared        = single(NaN);
%         data.(L).vonMisesPos.sse                = single(NaN);
%         data.(L).vonMisesPos.dfe                = single(NaN);
%         data.(L).vonMisesPos.rmse               = single(NaN);
%         data.(L).vonMisesPos.pvalue             = single(NaN);
%         data.(L).vonMisesPos.aic                = single(NaN);
%         data.(L).vonMisesPos.bic                = single(NaN);
        
%         data.(L).vonMisesPos_smoothed.average            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesPos_smoothed.startPoint         = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesPos_smoothed.yfit               = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesPos_smoothed.coefs              = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesPos_smoothed.rsquared           = single(NaN);
%         data.(L).vonMisesPos_smoothed.adjrsquared        = single(NaN);
%         data.(L).vonMisesPos_smoothed.sse                = single(NaN);
%         data.(L).vonMisesPos_smoothed.dfe                = single(NaN);
%         data.(L).vonMisesPos_smoothed.rmse               = single(NaN);
%         data.(L).vonMisesPos_smoothed.pvalue             = single(NaN);
%         data.(L).vonMisesPos_smoothed.aic                = single(NaN);
%         data.(L).vonMisesPos_smoothed.bic                = single(NaN);
        
%         data.(L).vonMisesNeg.average            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesNeg.startPoint         = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesNeg.yfit               = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesNeg.coefs              = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesNeg.rsquared           = single(NaN);
%         data.(L).vonMisesNeg.adjrsquared        = single(NaN);
%         data.(L).vonMisesNeg.sse                = single(NaN);
%         data.(L).vonMisesNeg.dfe                = single(NaN);
%         data.(L).vonMisesNeg.rmse               = single(NaN);
%         data.(L).vonMisesNeg.pvalue             = single(NaN);
%         data.(L).vonMisesNeg.aic                = single(NaN);
%         data.(L).vonMisesNeg.bic                = single(NaN);
        
%         data.(L).vonMisesNeg_smoothed.average            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesNeg_smoothed.startPoint         = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesNeg_smoothed.yfit               = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).vonMisesNeg_smoothed.coefs              = single([NaN NaN NaN NaN]);
%         data.(L).vonMisesNeg_smoothed.rsquared           = single(NaN);
%         data.(L).vonMisesNeg_smoothed.adjrsquared        = single(NaN);
%         data.(L).vonMisesNeg_smoothed.sse                = single(NaN);
%         data.(L).vonMisesNeg_smoothed.dfe                = single(NaN);
%         data.(L).vonMisesNeg_smoothed.rmse               = single(NaN);
%         data.(L).vonMisesNeg_smoothed.pvalue             = single(NaN);
%         data.(L).vonMisesNeg_smoothed.aic                = single(NaN);
%         data.(L).vonMisesNeg_smoothed.bic                = single(NaN);
        
        % median split - IBI low
        data.(L).IBI_median                     = single(NaN);
        
        data.(L).lowIBI_timeRRstart             = single(NaN);
        data.(L).lowIBI_timeRRend               = single(NaN);
        data.(L).lowIBI_starts                  = single(NaN);
        data.(L).lowIBI_meanHR_bpm              = single(NaN);
        data.(L).lowIBI_medianHR_bpm            = single(NaN);
        data.(L).lowIBI_stdHR_bpm               = single(NaN);
        data.(L).lowIBI_SDNN_ms                 = single(NaN);
        
%         data.(L).lowIBI_linear.yfit             = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_linear.coefs            = single([NaN NaN]);
%         data.(L).lowIBI_linear.rsquared         = single(NaN);
%         data.(L).lowIBI_linear.adjrsquared      = single(NaN);
%         data.(L).lowIBI_linear.sse              = single(NaN);
%         data.(L).lowIBI_linear.dfe              = single(NaN);
%         data.(L).lowIBI_linear.rmse             = single(NaN);
%         data.(L).lowIBI_linear.pvalue           = single([NaN; NaN]);
%         data.(L).lowIBI_linear.aic              = single(NaN);
%         data.(L).lowIBI_linear.bic              = single(NaN);
        
%         data.(L).lowIBI_linear_smoothed.yfit             = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_linear_smoothed.coefs            = single([NaN NaN]);
%         data.(L).lowIBI_linear_smoothed.rsquared         = single(NaN);
%         data.(L).lowIBI_linear_smoothed.adjrsquared      = single(NaN);
%         data.(L).lowIBI_linear_smoothed.sse              = single(NaN);
%         data.(L).lowIBI_linear_smoothed.dfe              = single(NaN);
%         data.(L).lowIBI_linear_smoothed.rmse             = single(NaN);
%         data.(L).lowIBI_linear_smoothed.pvalue           = single([NaN; NaN]);
%         data.(L).lowIBI_linear_smoothed.aic              = single(NaN);
%         data.(L).lowIBI_linear_smoothed.bic              = single(NaN);
        
%         data.(L).lowIBI_cosine.average          = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_cosine.startPoint       = single([NaN NaN NaN]);
%         data.(L).lowIBI_cosine.yfit             = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_cosine.coefs            = single([NaN NaN NaN]);
%         data.(L).lowIBI_cosine.rsquared         = single(NaN);
%         data.(L).lowIBI_cosine.adjrsquared      = single(NaN);
%         data.(L).lowIBI_cosine.sse              = single(NaN);
%         data.(L).lowIBI_cosine.dfe              = single(NaN);
%         data.(L).lowIBI_cosine.rmse             = single(NaN);
%         data.(L).lowIBI_cosine.pvalue           = single(NaN);
%         data.(L).lowIBI_cosine.aic              = single(NaN);
%         data.(L).lowIBI_cosine.bic              = single(NaN);
        
%         data.(L).lowIBI_cosine_smoothed.average          = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_cosine_smoothed.startPoint       = single([NaN NaN NaN]);
%         data.(L).lowIBI_cosine_smoothed.yfit             = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_cosine_smoothed.coefs            = single([NaN NaN NaN]);
%         data.(L).lowIBI_cosine_smoothed.rsquared         = single(NaN);
%         data.(L).lowIBI_cosine_smoothed.adjrsquared      = single(NaN);
%         data.(L).lowIBI_cosine_smoothed.sse              = single(NaN);
%         data.(L).lowIBI_cosine_smoothed.dfe              = single(NaN);
%         data.(L).lowIBI_cosine_smoothed.rmse             = single(NaN);
%         data.(L).lowIBI_cosine_smoothed.pvalue           = single(NaN);
%         data.(L).lowIBI_cosine_smoothed.aic              = single(NaN);
%         data.(L).lowIBI_cosine_smoothed.bic              = single(NaN);
        
%         data.(L).lowIBI_vonMisesPos.average     = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesPos.startPoint  = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesPos.yfit        = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesPos.coefs       = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesPos.rsquared    = single(NaN);
%         data.(L).lowIBI_vonMisesPos.adjrsquared = single(NaN);
%         data.(L).lowIBI_vonMisesPos.sse         = single(NaN);
%         data.(L).lowIBI_vonMisesPos.dfe         = single(NaN);
%         data.(L).lowIBI_vonMisesPos.rmse        = single(NaN);
%         data.(L).lowIBI_vonMisesPos.pvalue      = single(NaN);
%         data.(L).lowIBI_vonMisesPos.aic         = single(NaN);
%         data.(L).lowIBI_vonMisesPos.bic         = single(NaN);
        
%         data.(L).lowIBI_vonMisesPos_smoothed.average     = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesPos_smoothed.startPoint  = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesPos_smoothed.yfit        = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesPos_smoothed.coefs       = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesPos_smoothed.rsquared    = single(NaN);
%         data.(L).lowIBI_vonMisesPos_smoothed.adjrsquared = single(NaN);
%         data.(L).lowIBI_vonMisesPos_smoothed.sse         = single(NaN);
%         data.(L).lowIBI_vonMisesPos_smoothed.dfe         = single(NaN);
%         data.(L).lowIBI_vonMisesPos_smoothed.rmse        = single(NaN);
%         data.(L).lowIBI_vonMisesPos_smoothed.pvalue      = single(NaN);
%         data.(L).lowIBI_vonMisesPos_smoothed.aic         = single(NaN);
%         data.(L).lowIBI_vonMisesPos_smoothed.bic         = single(NaN);
        
%         data.(L).lowIBI_vonMisesNeg.average     = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesNeg.startPoint  = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesNeg.yfit        = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesNeg.coefs       = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesNeg.rsquared    = single(NaN);
%         data.(L).lowIBI_vonMisesNeg.adjrsquared = single(NaN);
%         data.(L).lowIBI_vonMisesNeg.sse         = single(NaN);
%         data.(L).lowIBI_vonMisesNeg.dfe         = single(NaN);
%         data.(L).lowIBI_vonMisesNeg.rmse        = single(NaN);
%         data.(L).lowIBI_vonMisesNeg.pvalue      = single(NaN);
%         data.(L).lowIBI_vonMisesNeg.aic         = single(NaN);
%         data.(L).lowIBI_vonMisesNeg.bic         = single(NaN);
        
%         data.(L).lowIBI_vonMisesNeg_smoothed.average     = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesNeg_smoothed.startPoint  = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesNeg_smoothed.yfit        = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).lowIBI_vonMisesNeg_smoothed.coefs       = single([NaN NaN NaN NaN]);
%         data.(L).lowIBI_vonMisesNeg_smoothed.rsquared    = single(NaN);
%         data.(L).lowIBI_vonMisesNeg_smoothed.adjrsquared = single(NaN);
%         data.(L).lowIBI_vonMisesNeg_smoothed.sse         = single(NaN);
%         data.(L).lowIBI_vonMisesNeg_smoothed.dfe         = single(NaN);
%         data.(L).lowIBI_vonMisesNeg_smoothed.rmse        = single(NaN);
%         data.(L).lowIBI_vonMisesNeg_smoothed.pvalue      = single(NaN);
%         data.(L).lowIBI_vonMisesNeg_smoothed.aic         = single(NaN);
%         data.(L).lowIBI_vonMisesNeg_smoothed.bic         = single(NaN);
        
        % median split - IBI high
        data.(L).highIBI_timeRRstart            = single(NaN);
        data.(L).highIBI_timeRRend              = single(NaN);
        data.(L).highIBI_starts                 = single(NaN);
        data.(L).highIBI_meanHR_bpm             = single(NaN);
        data.(L).highIBI_medianHR_bpm           = single(NaN);
        data.(L).highIBI_stdHR_bpm              = single(NaN);
        data.(L).highIBI_SDNN_ms                = single(NaN);
        
%         data.(L).highIBI_linear.yfit            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_linear.coefs           = single([NaN NaN]);
%         data.(L).highIBI_linear.rsquared        = single(NaN);
%         data.(L).highIBI_linear.adjrsquared     = single(NaN);
%         data.(L).highIBI_linear.sse             = single(NaN);
%         data.(L).highIBI_linear.dfe             = single(NaN);
%         data.(L).highIBI_linear.rmse            = single(NaN);
%         data.(L).highIBI_linear.pvalue          = single([NaN; NaN]);
%         data.(L).highIBI_linear.aic             = single(NaN);
%         data.(L).highIBI_linear.bic             = single(NaN);
        
%         data.(L).highIBI_linear_smoothed.yfit            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_linear_smoothed.coefs           = single([NaN NaN]);
%         data.(L).highIBI_linear_smoothed.rsquared        = single(NaN);
%         data.(L).highIBI_linear_smoothed.adjrsquared     = single(NaN);
%         data.(L).highIBI_linear_smoothed.sse             = single(NaN);
%         data.(L).highIBI_linear_smoothed.dfe             = single(NaN);
%         data.(L).highIBI_linear_smoothed.rmse            = single(NaN);
%         data.(L).highIBI_linear_smoothed.pvalue          = single([NaN; NaN]);
%         data.(L).highIBI_linear_smoothed.aic             = single(NaN);
%         data.(L).highIBI_linear_smoothed.bic             = single(NaN);
        
%         data.(L).highIBI_cosine.average         = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_cosine.startPoint      = single([NaN NaN NaN]);
%         data.(L).highIBI_cosine.yfit            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_cosine.coefs           = single([NaN NaN NaN]);
%         data.(L).highIBI_cosine.rsquared        = single(NaN);
%         data.(L).highIBI_cosine.adjrsquared     = single(NaN);
%         data.(L).highIBI_cosine.sse             = single(NaN);
%         data.(L).highIBI_cosine.dfe             = single(NaN);
%         data.(L).highIBI_cosine.rmse            = single(NaN);
%         data.(L).highIBI_cosine.pvalue          = single(NaN);
%         data.(L).highIBI_cosine.aic             = single(NaN);
%         data.(L).highIBI_cosine.bic             = single(NaN);
        
%         data.(L).highIBI_cosine_smoothed.average         = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_cosine_smoothed.startPoint      = single([NaN NaN NaN]);
%         data.(L).highIBI_cosine_smoothed.yfit            = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_cosine_smoothed.coefs           = single([NaN NaN NaN]);
%         data.(L).highIBI_cosine_smoothed.rsquared        = single(NaN);
%         data.(L).highIBI_cosine_smoothed.adjrsquared     = single(NaN);
%         data.(L).highIBI_cosine_smoothed.sse             = single(NaN);
%         data.(L).highIBI_cosine_smoothed.dfe             = single(NaN);
%         data.(L).highIBI_cosine_smoothed.rmse            = single(NaN);
%         data.(L).highIBI_cosine_smoothed.pvalue          = single(NaN);
%         data.(L).highIBI_cosine_smoothed.aic             = single(NaN);
%         data.(L).highIBI_cosine_smoothed.bic             = single(NaN);
        
%         data.(L).highIBI_vonMisesPos.average    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesPos.startPoint = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesPos.yfit       = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesPos.coefs      = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesPos.rsquared   = single(NaN);
%         data.(L).highIBI_vonMisesPos.adjrsquared= single(NaN);
%         data.(L).highIBI_vonMisesPos.sse        = single(NaN);
%         data.(L).highIBI_vonMisesPos.dfe        = single(NaN);
%         data.(L).highIBI_vonMisesPos.rmse       = single(NaN);
%         data.(L).highIBI_vonMisesPos.pvalue     = single(NaN);
%         data.(L).highIBI_vonMisesPos.aic        = single(NaN);
%         data.(L).highIBI_vonMisesPos.bic        = single(NaN);
        
%         data.(L).highIBI_vonMisesPos_smoothed.average    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesPos_smoothed.startPoint = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesPos_smoothed.yfit       = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesPos_smoothed.coefs      = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesPos_smoothed.rsquared   = single(NaN);
%         data.(L).highIBI_vonMisesPos_smoothed.adjrsquared= single(NaN);
%         data.(L).highIBI_vonMisesPos_smoothed.sse        = single(NaN);
%         data.(L).highIBI_vonMisesPos_smoothed.dfe        = single(NaN);
%         data.(L).highIBI_vonMisesPos_smoothed.rmse       = single(NaN);
%         data.(L).highIBI_vonMisesPos_smoothed.pvalue     = single(NaN);
%         data.(L).highIBI_vonMisesPos_smoothed.aic        = single(NaN);
%         data.(L).highIBI_vonMisesPos_smoothed.bic        = single(NaN);
        
%         data.(L).highIBI_vonMisesNeg.average    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesNeg.startPoint = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesNeg.yfit       = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesNeg.coefs      = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesNeg.rsquared   = single(NaN);
%         data.(L).highIBI_vonMisesNeg.adjrsquared= single(NaN);
%         data.(L).highIBI_vonMisesNeg.sse        = single(NaN);
%         data.(L).highIBI_vonMisesNeg.dfe        = single(NaN);
%         data.(L).highIBI_vonMisesNeg.rmse       = single(NaN);
%         data.(L).highIBI_vonMisesNeg.pvalue     = single(NaN);
%         data.(L).highIBI_vonMisesNeg.aic        = single(NaN);
%         data.(L).highIBI_vonMisesNeg.bic        = single(NaN);
        
%         data.(L).highIBI_vonMisesNeg_smoothed.average    = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesNeg_smoothed.startPoint = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesNeg_smoothed.yfit       = single(nan(1,cfg.phase.N_phase_bins));
%         data.(L).highIBI_vonMisesNeg_smoothed.coefs      = single([NaN NaN NaN NaN]);
%         data.(L).highIBI_vonMisesNeg_smoothed.rsquared   = single(NaN);
%         data.(L).highIBI_vonMisesNeg_smoothed.adjrsquared= single(NaN);
%         data.(L).highIBI_vonMisesNeg_smoothed.sse        = single(NaN);
%         data.(L).highIBI_vonMisesNeg_smoothed.dfe        = single(NaN);
%         data.(L).highIBI_vonMisesNeg_smoothed.rmse       = single(NaN);
%         data.(L).highIBI_vonMisesNeg_smoothed.pvalue     = single(NaN);
%         data.(L).highIBI_vonMisesNeg_smoothed.aic        = single(NaN);
%         data.(L).highIBI_vonMisesNeg_smoothed.bic        = single(NaN);
        
        data.(L).waveforms_microvolts           = single(nan(1,32)); % 32 - number of points per spike
        data.(L).waveforms_upsampled_microvolts = single(nan(1,128)); % 128 - number of points per upsampled spike 4*32
        data.(L).waveforms_byBin_microvolts     = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).AMP_microV                     = single(NaN);
        data.(L).AMP_voltageBinned              = single(nan(1,500));
        data.(L).AMP_voltageBins                = single(nan(1,500));
        data.(L).HW_ms                          = single(NaN);
        data.(L).TPW_ms                         = single(NaN);
        data.(L).REP_ms                         = single(NaN);
        data.(L).AMP_microV_byBin               = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).HW_ms_byBin                    = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).TPW_ms_byBin                   = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).REP_ms_byBin                   = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).AMP_reshuffled_avg             = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).AMP_lowerPrctile_2_5           = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).AMP_upperPrctile_97_5          = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).HW_reshuffled_avg              = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).HW_lowerPrctile_2_5            = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).HW_upperPrctile_97_5           = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).TPW_reshuffled_avg             = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).TPW_lowerPrctile_2_5           = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).TPW_upperPrctile_97_5          = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).REP_reshuffled_avg             = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).REP_lowerPrctile_2_5           = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).REP_upperPrctile_97_5          = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).AMP_MI                         = single(nan(1,5));
        data.(L).HW_MI                          = single(nan(1,5));
        data.(L).TPW_MI                         = single(nan(1,5));
        data.(L).REP_MI                         = single(nan(1,5));
        data.(L).AMP_MI_binned                  = single(nan(1,5));
        data.(L).HW_MI_binned                   = single(nan(1,5));
        data.(L).TPW_MI_binned                  = single(nan(1,5));
        data.(L).REP_MI_binned                  = single(nan(1,5));
        data.(L).AMP_microV_byBin_smoothed      = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).HW_ms_byBin_smoothed           = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).TPW_ms_byBin_smoothed          = single(nan(1,cfg.phase.N_phase_bins));
        data.(L).REP_ms_byBin_smoothed          = single(nan(1,cfg.phase.N_phase_bins));
        %         data.(L).allCorr                        = single(nan(1,6));
        %         data.(L).allLinMod                      = single(nan(6,2));
        data.(L).AMP_max_consec_bins            = single(nan(1,1));
        data.(L).AMP_modulation_index           = single(nan(1,1));
        data.(L).HW_max_consec_bins             = single(nan(1,1));
        data.(L).HW_modulation_index            = single(nan(1,1));
        data.(L).TPW_max_consec_bins            = single(nan(1,1));
        data.(L).TPW_modulation_index           = single(nan(1,1));
        data.(L).REP_max_consec_bins            = single(nan(1,1));
        data.(L).REP_modulation_index           = single(nan(1,1));
        data.(L).distance2thr                   = single(nan(1,1));
        data.(L).cc_PSTH_feature                = single(nan(1,4));
        data.(L).pp_PSTH_feature                = single(nan(1,4));
        data.(L).pperm_PSTH_feature             = single(nan(1,4)); % permuted p for the correlation between each feature and phase PSTH
        data.(L).timeRRstart                    = single(nan(1,1));
    end
    
    % find the corresponding WC file and load it
    chNum = data.channel;
    blkNum = unique([pop.block]);
    WCfile = ph_figure_out_waveclus_file_by_channel_and_blocks(chNum, blkNum, cfg.Input_WC);
    WC = load(WCfile, 'thr');
    if length(WC.thr) == 4
        data.thresholds_microV = 10^6 * WC.thr;
        data.thresholds_microV(3:4) = -1*data.thresholds_microV(3:4);
    end
    clear WC
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        
        %% get condition AND valid block trials only
        CT = ecg_bna_get_condition_trials(T, cfg.condition(c));
        tr=ismember([T.block],blocks) & CT;
        if sum(tr)<=1 || (~isfield(Rpeaks, 'RPEAK_ts_insp') && cfg.process_Rpeaks_inhalation_exhalation) || (~isfield(Rpeaks, 'RPEAK_ts_exp') && cfg.process_Rpeaks_inhalation_exhalation) % do calculations only if number of trials > 1
            continue
        end
        popcell=pop.trial(tr);
        trcell=T(tr);
        
        % 0. Prepare data variables
        states_onset               = {trcell.states_onset};
        states                     = {trcell.states};
        TDT_ECG1_t0_from_rec_start = {trcell.TDT_ECG1_t0_from_rec_start};
        block_nums                 = {trcell.block};
        state2_times               = cellfun(@(x,y) x(y == 2), states_onset, states, 'Uniformoutput', false); % trial starts = state 2
        state90_times              = cellfun(@(x,y) x(y == 90), states_onset, states, 'Uniformoutput', false); % trial ends = state 90
        
        % compute RR-intervals
        valid_RRinterval_ends      = single([Rpeaks(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
        valid_RRinterval_starts    = single(valid_RRinterval_ends - [Rpeaks(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
        % compute shuffled RR-intervals
        shuffled_RRinterval_ends   = single([Rpeaks(b).(['shuffled_ts' cfg.condition(c).Rpeak_field])]);
        shuffled_RRinterval_starts = single(shuffled_RRinterval_ends - [Rpeaks(b).(['shuffled_dur' cfg.condition(c).Rpeak_field])]);
        
        shuffled_RRinterval_ends   = mat2cell(shuffled_RRinterval_ends, ones(size(shuffled_RRinterval_ends,1),1), size(shuffled_RRinterval_ends,2));
        shuffled_RRinterval_starts = mat2cell(shuffled_RRinterval_starts, ones(size(shuffled_RRinterval_starts,1),1), size(shuffled_RRinterval_starts,2));
        
        % 0. figure out RR-intervals lying within trials
        trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
        trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
        
        % Create a logical array by vectorized comparison
        RR_within_trial_idx = any(bsxfun(@gt,valid_RRinterval_starts',trial_starts_one_stream) & ...
            bsxfun(@lt,valid_RRinterval_ends',trial_ends_one_stream), 2);
        
        % Use logical indexing to filter valid RR intervals
        valid_RRinterval_starts = valid_RRinterval_starts(RR_within_trial_idx);
        valid_RRinterval_ends   = valid_RRinterval_ends(RR_within_trial_idx);
        
        % reshuffled - get rid of RRs beyond the current set of trials
        for shuffNum = 1:length(shuffled_RRinterval_starts)
            shuffledRR_within_trial_idx = any(bsxfun(@gt,shuffled_RRinterval_starts{shuffNum}',trial_starts_one_stream) & ...
                bsxfun(@lt,shuffled_RRinterval_ends{shuffNum}',trial_ends_one_stream), 2);
            
            shuffled_RRinterval_ends{shuffNum}   = shuffled_RRinterval_ends{shuffNum}(shuffledRR_within_trial_idx); 
            shuffled_RRinterval_starts{shuffNum} = shuffled_RRinterval_starts{shuffNum}(shuffledRR_within_trial_idx);
        end
        
        % compute parameters of heart activity
        data.(L).meanHR_bpm   = mean(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).medianHR_bpm = median(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).stdHR_bpm    = std(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).SDNN_ms      = std(1000 * (valid_RRinterval_ends - valid_RRinterval_starts));
        
        % implement median split to heart-cycle durations
        RRs   = valid_RRinterval_ends - valid_RRinterval_starts;
        M_IBI = median(valid_RRinterval_ends - valid_RRinterval_starts);
        
        data.(L).IBI_median = M_IBI;
        
        % reshuffle based on the computed median
        %% lowIBI - set up settings
        cfg.time.IBI       = 1;
        cfg.time.IBI_low   = 1;
        cfg.time.IBI_high  = 0;
        cfg.time.IBI_thrsh = data.(L).IBI_median* ones(1, length(cfg.condition));
        Rpeaks_lowIBI      = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info,cfg.time);
        for XXX=1:numel(Rpeaks_lowIBI)
            Rpeaks_lowIBI(XXX).RPEAK_ts=Rpeaks_lowIBI(XXX).RPEAK_ts-Rpeaks_lowIBI(XXX).offset+Rpeaks(XXX).offset;
            Rpeaks_lowIBI(XXX).shuffled_ts=Rpeaks_lowIBI(XXX).shuffled_ts-Rpeaks_lowIBI(XXX).offset+Rpeaks(XXX).offset;
            Rpeaks_lowIBI(XXX).offset=Rpeaks(XXX).offset;
        end
        
        %% compute within trial lowIBI RR intervals for real and shuffled data
        % compute RR-intervals
        lowIBI_valid_RRinterval_ends      = single([Rpeaks_lowIBI(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
        lowIBI_valid_RRinterval_starts    = single(lowIBI_valid_RRinterval_ends - [Rpeaks_lowIBI(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
        % compute shuffled RR-intervals
        lowIBI_shuffled_RRinterval_ends   = single([Rpeaks_lowIBI(b).(['shuffled_ts' cfg.condition(c).Rpeak_field])]);
        lowIBI_shuffled_RRinterval_starts = single(lowIBI_shuffled_RRinterval_ends - [Rpeaks_lowIBI(b).(['shuffled_dur' cfg.condition(c).Rpeak_field])]);
        
        lowIBI_shuffled_RRinterval_ends   = mat2cell(lowIBI_shuffled_RRinterval_ends, ones(size(lowIBI_shuffled_RRinterval_ends,1),1), size(lowIBI_shuffled_RRinterval_ends,2));
        lowIBI_shuffled_RRinterval_starts = mat2cell(lowIBI_shuffled_RRinterval_starts, ones(size(lowIBI_shuffled_RRinterval_starts,1),1), size(lowIBI_shuffled_RRinterval_starts,2));
        
%         % 0. figure out RR-intervals lying within trials
%         lowIBI_trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks_lowIBI([Rpeaks_lowIBI.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
%         lowIBI_trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks_lowIBI([Rpeaks_lowIBI.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
        
        % Create a logical array by vectorized comparison
        lowIBI_RR_within_trial_idx = any(bsxfun(@gt,lowIBI_valid_RRinterval_starts',trial_starts_one_stream) & ...
                bsxfun(@lt,lowIBI_valid_RRinterval_ends',trial_ends_one_stream), 2);
        
        % Use logical indexing to filter valid RR intervals
        lowIBI_valid_RRinterval_starts = lowIBI_valid_RRinterval_starts(lowIBI_RR_within_trial_idx);
        lowIBI_valid_RRinterval_ends   = lowIBI_valid_RRinterval_ends(lowIBI_RR_within_trial_idx);
        
        % reshuffled - get rid of RRs beyond the current set of trials
        for shuffNum = 1:length(lowIBI_shuffled_RRinterval_starts)
            lowIBI_shuffledRR_within_trial_idx = any(bsxfun(@gt,lowIBI_shuffled_RRinterval_starts{shuffNum}',trial_starts_one_stream) & ...
                bsxfun(@lt,lowIBI_shuffled_RRinterval_ends{shuffNum}',trial_ends_one_stream), 2);
            
            lowIBI_shuffled_RRinterval_ends{shuffNum}   = lowIBI_shuffled_RRinterval_ends{shuffNum}(lowIBI_shuffledRR_within_trial_idx); 
            lowIBI_shuffled_RRinterval_starts{shuffNum} = lowIBI_shuffled_RRinterval_starts{shuffNum}(lowIBI_shuffledRR_within_trial_idx);
        end
        
        lowIBI_RRs   = lowIBI_valid_RRinterval_ends - lowIBI_valid_RRinterval_starts;
        
        %% highIBI - set up settings
        cfg.time.IBI       = 1;
        cfg.time.IBI_low   = 0;
        cfg.time.IBI_high  = 1;
        cfg.time.IBI_thrsh = data.(L).IBI_median* ones(1, length(cfg.condition));
        Rpeaks_highIBI      = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info,cfg.time);
        for XXX=1:numel(Rpeaks_highIBI)
            Rpeaks_highIBI(XXX).RPEAK_ts=Rpeaks_highIBI(XXX).RPEAK_ts-Rpeaks_highIBI(XXX).offset+Rpeaks(XXX).offset;
            Rpeaks_highIBI(XXX).shuffled_ts=Rpeaks_highIBI(XXX).shuffled_ts-Rpeaks_highIBI(XXX).offset+Rpeaks(XXX).offset;
            Rpeaks_highIBI(XXX).offset=Rpeaks(XXX).offset;
        end
        
        % !!! IMPORTANT !!! erase settings related to median split
        % reshuffles and re-initiate them for the next unit
        cfg.time = rmfield(cfg.time, {'IBI', 'IBI_low', 'IBI_high', 'IBI_thrsh'});
        
        %% compute within trial highIBI RR intervals for real and shuffled data
        % compute RR-intervals
        highIBI_valid_RRinterval_ends      = single([Rpeaks_highIBI(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
        highIBI_valid_RRinterval_starts    = single(highIBI_valid_RRinterval_ends - [Rpeaks_highIBI(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
        % compute shuffled RR-intervals
        highIBI_shuffled_RRinterval_ends   = single([Rpeaks_highIBI(b).(['shuffled_ts' cfg.condition(c).Rpeak_field])]);
        highIBI_shuffled_RRinterval_starts = single(highIBI_shuffled_RRinterval_ends - [Rpeaks_highIBI(b).(['shuffled_dur' cfg.condition(c).Rpeak_field])]);
        
        highIBI_shuffled_RRinterval_ends   = mat2cell(highIBI_shuffled_RRinterval_ends, ones(size(highIBI_shuffled_RRinterval_ends,1),1), size(highIBI_shuffled_RRinterval_ends,2));
        highIBI_shuffled_RRinterval_starts = mat2cell(highIBI_shuffled_RRinterval_starts, ones(size(highIBI_shuffled_RRinterval_starts,1),1), size(highIBI_shuffled_RRinterval_starts,2));
        
%         % 0. figure out RR-intervals lying within trials
%         trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks_highIBI([Rpeaks_highIBI.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
%         trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks_highIBI([Rpeaks_highIBI.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
        
        % Create a logical array by vectorized comparison
        highIBI_RR_within_trial_idx = any(bsxfun(@gt,highIBI_valid_RRinterval_starts',trial_starts_one_stream) & ...
                bsxfun(@lt,highIBI_valid_RRinterval_ends',trial_ends_one_stream), 2);
        
        % Use logical indexing to filter valid RR intervals
        highIBI_valid_RRinterval_starts = highIBI_valid_RRinterval_starts(highIBI_RR_within_trial_idx);
        highIBI_valid_RRinterval_ends   = highIBI_valid_RRinterval_ends(highIBI_RR_within_trial_idx);
        
        % reshuffled - get rid of RRs beyond the current set of trials
        for shuffNum = 1:length(highIBI_shuffled_RRinterval_starts)
            highIBI_shuffledRR_within_trial_idx = any(bsxfun(@gt,highIBI_shuffled_RRinterval_starts{shuffNum}',trial_starts_one_stream) & ...
                bsxfun(@lt,highIBI_shuffled_RRinterval_ends{shuffNum}',trial_ends_one_stream), 2);
            
            highIBI_shuffled_RRinterval_ends{shuffNum}   = highIBI_shuffled_RRinterval_ends{shuffNum}(highIBI_shuffledRR_within_trial_idx); 
            highIBI_shuffled_RRinterval_starts{shuffNum} = highIBI_shuffled_RRinterval_starts{shuffNum}(highIBI_shuffledRR_within_trial_idx);
        end
        
        highIBI_RRs   = highIBI_valid_RRinterval_ends - highIBI_valid_RRinterval_starts;
        
%         lowIBIids  = realPSTHs.RDs{1} < Output.(L).IBI_median;
%         highIBIids = realPSTHs.RDs{1} > Output.(L).IBI_median;
        
        data.(L).lowIBI_timeRRstart   = lowIBI_valid_RRinterval_starts;
        data.(L).lowIBI_timeRRend     = lowIBI_valid_RRinterval_ends;
        data.(L).lowIBI_starts        = lowIBI_RRs;
        data.(L).lowIBI_meanHR_bpm    = mean(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_medianHR_bpm  = median(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_stdHR_bpm     = std(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_SDNN_ms       = std(1000 * data.(L).lowIBI_starts);
        
        data.(L).highIBI_timeRRstart  = highIBI_valid_RRinterval_starts;
        data.(L).highIBI_timeRRend    = highIBI_valid_RRinterval_ends;
        data.(L).highIBI_starts       = highIBI_RRs;
        data.(L).highIBI_meanHR_bpm   = mean(60 ./ data.(L).highIBI_starts);
        data.(L).highIBI_medianHR_bpm = median(60 ./ data.(L).highIBI_starts);
        data.(L).highIBI_stdHR_bpm    = std(60 ./ data.(L).highIBI_starts);
        data.(L).highIBI_SDNN_ms      = std(1000 * data.(L).highIBI_starts);
        
        % 1. take arrival times and the corresponding waveforms
        AT = {popcell.arrival_times};
        WF = {popcell.waveforms};
        % 2. choose only those that happen after MP-state 1
        idx_after_state1 = cellfun(@(x,y) x>y, AT, state2_times, 'Uniformoutput', false);
        % 3. add TDT_ECG1_t0_from_rec_start and Rpeak block offset to spike times
        AT_one_stream_cell = cellfun(@(x,y,z,a) x(y)+z+Rpeaks([Rpeaks.block] == a).offset, AT, idx_after_state1, TDT_ECG1_t0_from_rec_start, block_nums, 'Uniformoutput', false);
        WF_one_stream_cell = cellfun(@(x,y) x(y,:), WF, idx_after_state1, 'Uniformoutput', false);
        % 4. merge all spike times and waveforms together
        AT_one_stream = cat(1, AT_one_stream_cell{:});
        WF_one_stream = cat(1, WF_one_stream_cell{:});
        if size(WF_one_stream,1) < 10
            continue
        end
        % 5. calculate heart cycle phase where individual spikes ended up
        [eventPhases, eventsTaken, cycleNums_withSpikes] = ...
            DAG_eventPhase(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream); % all data
        
        [lowIBI_eventPhases, ~, lowIBI_cycleNums_withSpikes] = ...
            DAG_eventPhase(data.(L).lowIBI_timeRRstart, data.(L).lowIBI_timeRRend, AT_one_stream);
        
        [highIBI_eventPhases, ~, highIBI_cycleNums_withSpikes] = ...
            DAG_eventPhase(data.(L).highIBI_timeRRstart, data.(L).highIBI_timeRRend, AT_one_stream);
        
        % check the number of spikes left after computing phases
        if length(eventPhases) < 3 || length(lowIBI_eventPhases) < 3 || length(highIBI_eventPhases) < 3 || ...
                length(valid_RRinterval_starts) < 100 || length(data.(L).lowIBI_timeRRstart) < 3 || length(data.(L).highIBI_timeRRstart) < 3
            continue
        end
        
        % for reshuffled data - calculate heart cycle phase where individual spikes ended up
        shuffled_SDF         = nan(cfg.time.n_permutations,cfg.phase.N_phase_bins);
        shuffled_SDF_lowIBI  = nan(cfg.time.n_permutations,cfg.phase.N_phase_bins);
        shuffled_SDF_highIBI = nan(cfg.time.n_permutations,cfg.phase.N_phase_bins);
        
        tic
%         medianIBI = data.(L).IBI_median;
        phase_bin_centers = cfg.phase.phase_bin_centers;
        parfor shuffNum = 1:length(shuffled_RRinterval_starts)
            % all data
            [shuffled_eventPhases, ~, shuffled_cycleNums_withSpikes] = ...
                DAG_eventPhase(shuffled_RRinterval_starts{shuffNum}, shuffled_RRinterval_ends{shuffNum}, AT_one_stream);
            
            % all data
            shuffled_histogram = ...
                hist3([shuffled_eventPhases, shuffled_cycleNums_withSpikes], ...
                'ctrs', {phase_bin_centers 1:length(shuffled_RRinterval_starts{shuffNum})});
            
            shuffled_SDF(shuffNum,:) = mean(average_smooth_data(shuffled_histogram,cfg));
            
            % compute median IBI
%             currRRids  = unique(shuffled_cycleNums_withSpikes);
%             currRRdurs = shuffled_RRinterval_ends{shuffNum}(currRRids) - shuffled_RRinterval_starts{shuffNum}(currRRids);
            
            % low IBI
%             lowIBI_ids = currRRdurs < medianIBI;
            [lowIBI_shuffled_eventPhases, ~, lowIBI_shuffled_cycleNums_withSpikes] = ...
                DAG_eventPhase(lowIBI_shuffled_RRinterval_starts{shuffNum}, lowIBI_shuffled_RRinterval_ends{shuffNum}, AT_one_stream);
            
            % high IBI
%             highIBI_ids = currRRdurs > medianIBI;
            [highIBI_shuffled_eventPhases, ~, highIBI_shuffled_cycleNums_withSpikes] = ...
                DAG_eventPhase(highIBI_shuffled_RRinterval_starts{shuffNum}, highIBI_shuffled_RRinterval_ends{shuffNum}, AT_one_stream);
            
            % low IBI
            shuffled_histogram_lowIBI = ...
                hist3([lowIBI_shuffled_eventPhases, lowIBI_shuffled_cycleNums_withSpikes], ...
                'ctrs', {phase_bin_centers 1:length(lowIBI_RRs)});
            
            shuffled_SDF_lowIBI(shuffNum,:) = mean(average_smooth_data(shuffled_histogram_lowIBI,cfg));
            
            % high IBI
            shuffled_histogram_highIBI = ...
                hist3([highIBI_shuffled_eventPhases, highIBI_shuffled_cycleNums_withSpikes], ...
                'ctrs', {phase_bin_centers 1:length(highIBI_RRs)});
            
            shuffled_SDF_highIBI(shuffNum,:) = mean(average_smooth_data(shuffled_histogram_highIBI,cfg));
            
        end
        toc
        
        % 6. Put results for real data together
        data.(L).spike_phases_radians    = eventPhases;
        data.(L).spike_phases_histogram  = hist(data.(L).spike_phases_radians, cfg.phase.phase_bin_centers) / length(valid_RRinterval_starts) / length(cfg.phase.phase_bin_centers); % compute overal phase histogram
        data.(L).spike_phases_histogram2 = hist3([eventPhases, cycleNums_withSpikes], 'ctrs', {cfg.phase.phase_bin_centers 1:length(valid_RRinterval_starts)});
        
        data.(L).lowIBI_spike_phases_radians          = lowIBI_eventPhases;
        data.(L).lowIBI_spike_phases_histogram        = hist(data.(L).lowIBI_spike_phases_radians, cfg.phase.phase_bin_centers);
        data.(L).lowIBI_spike_phases_histogram2       = hist3([lowIBI_eventPhases, lowIBI_cycleNums_withSpikes], 'ctrs', {cfg.phase.phase_bin_centers 1:length(data.(L).lowIBI_timeRRstart)});
        
        data.(L).highIBI_spike_phases_radians         = highIBI_eventPhases;
        data.(L).highIBI_spike_phases_histogram       = hist(data.(L).highIBI_spike_phases_radians, cfg.phase.phase_bin_centers);
        data.(L).highIBI_spike_phases_histogram2      = hist3([highIBI_eventPhases, highIBI_cycleNums_withSpikes], 'ctrs', {cfg.phase.phase_bin_centers 1:length(data.(L).highIBI_timeRRstart)});
        
        % all data
        real_data = average_smooth_data(data.(L).spike_phases_histogram2,cfg);
        SD        = ecg_bna_do_statistics(real_data,shuffled_SDF,1:cfg.phase.N_phase_bins);
        
        data.(L).SD          = SD.SD_mean;
        data.(L).SD_STD      = SD.SD_STD;
        data.(L).SD_SEM      = SD.SD_SEM ;
        data.(L).SDP         = SD.SDPmean ;
        data.(L).SDPCL       = SD.SDPconf(1,:) ;
        data.(L).SDPCu       = SD.SDPconf(2,:) ;
        data.(L).sig_all     = SD.sig_all;
        data.(L).sig         = SD.sig;
        data.(L).sig_FR_diff = SD.sig_FR_diff;
        data.(L).sig_time    = SD.sig_time;
        data.(L).sig_n_bins  = SD.sig_n_bins;
        data.(L).sig_sign    = SD.sig_sign;
        data.(L).SDsubstractedSDP            = data.(L).SD - data.(L).SDP; % spikes/s, difference between mean and jittered data
        data.(L).SDsubstractedSDP_normalized = data.(L).SDsubstractedSDP ./ data.(L).SDP *100; % percent signal change
        data.(L).FR_ModIndex_SubtrSDP        = max(data.(L).SDsubstractedSDP) - min(data.(L).SDsubstractedSDP); % difference between max and min FR
        data.(L).FR_ModIndex_PcS             = max(data.(L).SDsubstractedSDP_normalized) - min(data.(L).SDsubstractedSDP_normalized); % difference between max and min % signal change
        
        real_data_lowIBI = average_smooth_data(data.(L).lowIBI_spike_phases_histogram2,cfg);
        SD_lowIBI        = ecg_bna_do_statistics(real_data_lowIBI,shuffled_SDF_lowIBI,1:cfg.phase.N_phase_bins);
        
        data.(L).lowIBI.SD          = SD_lowIBI.SD_mean;
        data.(L).lowIBI.SD_STD      = SD_lowIBI.SD_STD;
        data.(L).lowIBI.SD_SEM      = SD_lowIBI.SD_SEM ;
        data.(L).lowIBI.SDP         = SD_lowIBI.SDPmean ;
        data.(L).lowIBI.SDPCL       = SD_lowIBI.SDPconf(1,:) ;
        data.(L).lowIBI.SDPCu       = SD_lowIBI.SDPconf(2,:) ;
        data.(L).lowIBI.sig_all     = SD_lowIBI.sig_all;
        data.(L).lowIBI.sig         = SD_lowIBI.sig;
        data.(L).lowIBI.sig_FR_diff = SD_lowIBI.sig_FR_diff;
        data.(L).lowIBI.sig_time    = SD_lowIBI.sig_time;
        data.(L).lowIBI.sig_n_bins  = SD_lowIBI.sig_n_bins;
        data.(L).lowIBI.sig_sign    = SD_lowIBI.sig_sign;
        data.(L).lowIBI.SDsubstractedSDP            = data.(L).lowIBI.SD - data.(L).lowIBI.SDP; % spikes/s, difference between mean and jittered data
        data.(L).lowIBI.SDsubstractedSDP_normalized = data.(L).lowIBI.SDsubstractedSDP ./ data.(L).lowIBI.SDP *100; % percent signal change
        data.(L).lowIBI.FR_ModIndex_SubtrSDP        = max(data.(L).lowIBI.SDsubstractedSDP) - min(data.(L).lowIBI.SDsubstractedSDP); % difference between max and min FR
        data.(L).lowIBI.FR_ModIndex_PcS             = max(data.(L).lowIBI.SDsubstractedSDP_normalized) - min(data.(L).lowIBI.SDsubstractedSDP_normalized); % difference between max and min % signal change
        
        real_data_highIBI = average_smooth_data(data.(L).highIBI_spike_phases_histogram2,cfg);
        SD_highIBI        = ecg_bna_do_statistics(real_data_highIBI,shuffled_SDF_highIBI,1:cfg.phase.N_phase_bins);
        
        data.(L).highIBI.SD          = SD_highIBI.SD_mean;
        data.(L).highIBI.SD_STD      = SD_highIBI.SD_STD;
        data.(L).highIBI.SD_SEM      = SD_highIBI.SD_SEM ;
        data.(L).highIBI.SDP         = SD_highIBI.SDPmean ;
        data.(L).highIBI.SDPCL       = SD_highIBI.SDPconf(1,:) ;
        data.(L).highIBI.SDPCu       = SD_highIBI.SDPconf(2,:) ;
        data.(L).highIBI.sig_all     = SD_highIBI.sig_all;
        data.(L).highIBI.sig         = SD_highIBI.sig;
        data.(L).highIBI.sig_FR_diff = SD_highIBI.sig_FR_diff;
        data.(L).highIBI.sig_time    = SD_highIBI.sig_time;
        data.(L).highIBI.sig_n_bins  = SD_highIBI.sig_n_bins;
        data.(L).highIBI.sig_sign    = SD_highIBI.sig_sign;
        data.(L).highIBI.SDsubstractedSDP            = data.(L).highIBI.SD - data.(L).highIBI.SDP; % spikes/s, difference between mean and jittered data
        data.(L).highIBI.SDsubstractedSDP_normalized = data.(L).highIBI.SDsubstractedSDP ./ data.(L).highIBI.SDP *100; % percent signal change
        data.(L).highIBI.FR_ModIndex_SubtrSDP        = max(data.(L).highIBI.SDsubstractedSDP) - min(data.(L).highIBI.SDsubstractedSDP); % difference between max and min FR
        data.(L).highIBI.FR_ModIndex_PcS             = max(data.(L).highIBI.SDsubstractedSDP_normalized) - min(data.(L).highIBI.SDsubstractedSDP_normalized); % difference between max and min % signal change
        
        % Mosher's cosine fit the smoothed phase histogram
        [modIndex_phase_hist,phase_hist] = ...
            fitCardiacModulation(cfg.phase.phase_bin_centers, ...
            data.(L).spike_phases_histogram, {'PSTH'}, 0, [221]);
        
        data.(L).spike_phases_histogram_smoothed = phase_hist;
        data.(L).histogram_MI                    = modIndex_phase_hist(1);
        data.(L).histogram_p                     = modIndex_phase_hist(2);
        data.(L).histogram_phase                 = modIndex_phase_hist(3);
        data.(L).rsquared                        = modIndex_phase_hist(4);
        
%         % lowIBI
%         [~,lowIBI_phase_hist] = ...
%             fitCardiacModulation(cfg.phase.phase_bin_centers, ...
%             data.(L).lowIBI_spike_phases_histogram, {'PSTH'}, 0, [221]);
%         data.(L).lowIBI_spike_phases_histogram_smoothed = lowIBI_phase_hist;
        
%         % highIBI
%         [~,highIBI_phase_hist] = ...
%             fitCardiacModulation(cfg.phase.phase_bin_centers, ...
%             data.(L).highIBI_spike_phases_histogram, {'PSTH'}, 0, [221]);
%         data.(L).highIBI_spike_phases_histogram_smoothed = highIBI_phase_hist;
        
        %% LINEAR FITS
        % linear fit with 'fitlm' for all the data
%         data.(L).linear                 = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'linear');
        
        % linear fit with 'fitlm' for smoothed data
%         data.(L).linear_smoothed        = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'linear');
        
        % linear fit for the lower IBI
%         data.(L).lowIBI_linear          = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'linear');
        
        % linear fit for smoothed data
%         data.(L).lowIBI_linear_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram, length(data.(L).lowIBI_timeRRstart), 'linear');
        
        % linear fit for the higher IBI
%         data.(L).highIBI_linear         = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'linear');
        
        % % linear fit for smoothed data
%         data.(L).highIBI_linear_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram, length(data.(L).highIBI_timeRRstart), 'linear');
        
        %% COSINE FITS
        % cosine fit with 'fit' on all data points instead of smoothed mean
%         data.(L).cosine              = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'cosine');
        
        % cosine fit with 'fit' for smoothed data
%         data.(L).cosine_smoothed     = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'cosine');
        
        % cosine fit for the lower IBI
%         data.(L).lowIBI_cosine       = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'cosine');
        
        % fit for smoothed data
%         data.(L).lowIBI_cosine_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram, length(data.(L).lowIBI_timeRRstart), 'cosine');
        
        % cosine fit for the higher IBI
%         data.(L).highIBI_cosine      = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'cosine');
        
        % fit for smoothed data
%         data.(L).highIBI_cosine_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram, length(data.(L).lowIBI_timeRRstart), 'cosine');
        
        %% POSITIVE VON MISES FITS
        % fit all the data
%         data.(L).vonMisesPos          = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', 1);
        
        % fit all the data - smoothed
%         data.(L).vonMisesPos_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'vonMises', 1);
        
        % pos. von Mises fits - low IBI
%         data.(L).lowIBI_vonMisesPos   = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'vonMises', 1);
        
        % fit for smoothed data
%         data.(L).lowIBI_vonMisesPos_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram, length(data.(L).lowIBI_timeRRstart), 'vonMises', 1);
        
        % pos. von Mises fits - high IBI
%         data.(L).highIBI_vonMisesPos  = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'vonMises', 1);
        
        % fit for smoothed data
%         data.(L).highIBI_vonMisesPos_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram, length(data.(L).lowIBI_timeRRstart), 'vonMises', 1);
        
        %% NEGATIVE VON MISES FITS
        % fit all the data
%         data.(L).vonMisesNeg          = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', -1);
        
        % fit all the data - smoothed
%         data.(L).vonMisesNeg_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'vonMises', -1);
        
        % neg. von Mises fits - low IBI
%         data.(L).lowIBI_vonMisesNeg   = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'vonMises', -1);
        
        % fit for smoothed data
%         data.(L).lowIBI_vonMisesNeg_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram, length(data.(L).lowIBI_timeRRstart), 'vonMises', -1);
        
        % neg. von Mises fits - high IBI
%         data.(L).highIBI_vonMisesNeg  = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'vonMises', -1);
        
        % fit for smoothed data
%         data.(L).highIBI_vonMisesNeg_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).highIBI_spike_phases_histogram, length(data.(L).lowIBI_timeRRstart), 'vonMises', -1);
        
        [~, ~, bin] = histcounts(data.(L).spike_phases_radians, cfg.phase.phase_bins);
        
        data.(L).waveforms_microvolts                 = 10^6 * WF_one_stream(eventsTaken,:);
        waveforms_upsampled                           = interpft(data.(L).waveforms_microvolts, 32*4, 2);
        data.(L).waveforms_upsampled_microvolts       = shift2peak(cfg.phase.wf_times_interp_ms, waveforms_upsampled, cfg);
        data.(L).waveforms_byBin_microvolts           = arrayfun(@(x) mean(data.(L).waveforms_upsampled_microvolts(bin == x,:),1), 1:cfg.phase.N_phase_bins, 'UniformOutput', false);
        data.(L).waveforms_byBin_microvolts           = cat(1,data.(L).waveforms_byBin_microvolts{:});
        % 7. Calculate spike features with Mosher's procedure
        tic
        sMetric = ...
            struct('extremAmp', cell(length(data.(L).spike_phases_radians),1), ...
            'widthHW', cell(length(data.(L).spike_phases_radians),1), ...
            'widthTP', cell(length(data.(L).spike_phases_radians),1), ...
            'repolTime', cell(length(data.(L).spike_phases_radians),1));
        parfor wfNum = 1:length(data.(L).spike_phases_radians)
            sMetric(wfNum)=spikeWaveMetrics(double(data.(L).waveforms_upsampled_microvolts(wfNum,:)), 37, cfg.phase.Fs*4, 0, [1 0 0 0]); % 37 - index of th peak for updsampled data
        end
        toc
        % 8. put the resulting data together
        data.(L).AMP_microV               = [sMetric.extremAmp];
        data.(L).HW_ms                    = 10^3 * [sMetric.widthHW];
        data.(L).TPW_ms                   = 10^3 * [sMetric.widthTP];
        data.(L).REP_ms                   = 10^3 * [sMetric.repolTime];
        
        data.(L).AMP_microV_byBin         = arrayfun(@(x) nanmean(data.(L).AMP_microV(bin == x)), 1:cfg.phase.N_phase_bins); % mean by phase
        data.(L).HW_ms_byBin              = arrayfun(@(x) nanmean(data.(L).HW_ms(bin == x)), 1:cfg.phase.N_phase_bins);
        data.(L).TPW_ms_byBin             = arrayfun(@(x) nanmean(data.(L).TPW_ms(bin == x)), 1:cfg.phase.N_phase_bins);
        data.(L).REP_ms_byBin             = arrayfun(@(x) nanmean(data.(L).REP_ms(bin == x)), 1:cfg.phase.N_phase_bins);
        
        % 9. estimate significance with bootstrapping
        [data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).AMP_microV, bin, cfg.phase);
        [data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).HW_ms, bin, cfg.phase);
        [data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).TPW_ms, bin, cfg.phase);
        [data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).REP_ms, bin, cfg.phase);
        
        % 10. compute cosine fits and put those into the resulting
        % structure
        featureMatrix = ...
            [data.(L).AMP_microV_byBin; data.(L).HW_ms_byBin; ...
            data.(L).TPW_ms_byBin; data.(L).REP_ms_byBin];
        if size(featureMatrix,1) > 4
            featureMatrix = featureMatrix';
            if size(featureMatrix,1) ~= 4
                error('Dimensions of feature matrix aren''t suitable for the analysis')
            end
        end
        
        if length(data.(L).AMP_microV) - sum(isnan(data.(L).AMP_microV)) < 3
            continue
        end
        
        if sum(isnan(data.(L).AMP_microV)) == length(data.(L).AMP_microV) % if they're all nan
            data.(L).AMP_MI                   = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).AMP_microV, eventPhases, cfg);
            data.(L).AMP_MI                   = coefs;
            clear coefs
        end
        
        if sum(isnan(data.(L).HW_ms)) == length(data.(L).HW_ms) % if they're all nan
            data.(L).HW_MI                    = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).HW_ms, eventPhases, cfg);
            data.(L).HW_MI                    = coefs;
            clear coefs
        end
        
        if sum(isnan(data.(L).TPW_ms)) == length(data.(L).TPW_ms) % if they're all nan
            data.(L).TPW_MI                   = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).TPW_ms, eventPhases, cfg);
            data.(L).TPW_MI                   = coefs;
            clear coefs
        end
        
        if sum(isnan(data.(L).REP_ms)) == length(data.(L).REP_ms) % if they're all nan
            data.(L).REP_MI                   = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).REP_ms, eventPhases, cfg);
            data.(L).REP_MI                   = coefs;
            clear coefs
        end
        
        % use Mosher's cosine fitting procedure with binned data
        modIndex        = fitCardiacModulation(cfg.phase.phase_bin_centers, featureMatrix, {'AMP', 'HW', 'TPW', 'REP'}, 0, 0);
        data.(L).AMP_MI_binned = modIndex(1,:);
        data.(L).HW_MI_binned  = modIndex(2,:);
        data.(L).TPW_MI_binned = modIndex(3,:);
        data.(L).REP_MI_binned = modIndex(4,:);
        
        % store smoothed data for each measure
        if sum(isnan(data.(L).AMP_microV_byBin)) == length(data.(L).AMP_microV_byBin)
            data.(L).AMP_microV_byBin_smoothed    = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).AMP_microV_byBin_smoothed    = DAG_circ_smooth(data.(L).AMP_microV_byBin);
        end
        if sum(isnan(data.(L).HW_ms_byBin)) == length(data.(L).HW_ms_byBin)
            data.(L).HW_ms_byBin_smoothed         = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).HW_ms_byBin_smoothed         = DAG_circ_smooth(data.(L).HW_ms_byBin);
        end
        if sum(isnan(data.(L).TPW_ms_byBin)) == length(data.(L).TPW_ms_byBin)
            data.(L).TPW_ms_byBin_smoothed        = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).TPW_ms_byBin_smoothed        = DAG_circ_smooth(data.(L).TPW_ms_byBin);
        end
        if sum(isnan(data.(L).REP_ms_byBin)) == length(data.(L).REP_ms_byBin)
            data.(L).REP_ms_byBin_smoothed        = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).REP_ms_byBin_smoothed        = DAG_circ_smooth(data.(L).REP_ms_byBin);
        end
        
        %         data.(L).allCorr                      = allCorr;
        %         data.(L).allLinMod                    = allLinMod;
        
        [data.(L).AMP_max_consec_bins, data.(L).AMP_modulation_index] = ...
            significant_bins(data.(L).AMP_microV_byBin_smoothed', data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg);
        [data.(L).HW_max_consec_bins, data.(L).HW_modulation_index] = ...
            significant_bins(data.(L).HW_ms_byBin_smoothed', data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg);
        [data.(L).TPW_max_consec_bins, data.(L).TPW_modulation_index] = ...
            significant_bins(data.(L).TPW_ms_byBin_smoothed', data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg);
        [data.(L).REP_max_consec_bins, data.(L).REP_modulation_index] = ...
            significant_bins(data.(L).REP_ms_byBin_smoothed', data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg);
        
        % 11. estimate how much spike amplitude is far from the
        % thresholds
        % bin from the closest threshold first
        dropnan = isnan(data.(L).AMP_microV);
        higher_threshold_detected = ...
            all(data.(L).AMP_microV(dropnan) > data.thresholds_microV(1));
        lower_threshold_detected = ...
            all(data.(L).AMP_microV(dropnan) > data.thresholds_microV(2));
        if higher_threshold_detected && lower_threshold_detected % both thresholds are lower than all the spike amplitudes
            % take the higher one
            data.(L).AMP_voltageBins   = data.thresholds_microV(1):500;
            data.(L).AMP_voltageBinned = histc(data.(L).AMP_microV, data.(L).AMP_voltageBins); % bin amplitudes
        elseif lower_threshold_detected && ~higher_threshold_detected
            % take the lower one
            data.(L).AMP_voltageBins   = data.thresholds_microV(2):500;
            data.(L).AMP_voltageBinned = histc(data.(L).AMP_microV, data.(L).AMP_voltageBins);
        elseif ~lower_threshold_detected && ~higher_threshold_detected
            error('It''s an unknown error: there are spikes below both thresholds')
        elseif ~lower_threshold_detected && higher_threshold_detected
            error('It''s an unknown error: there are spikes that are below only lower threshold')
        end
        
        first_nonzero_bin          = find(data.(L).AMP_voltageBinned,1,'first');
        
        data.(L).distance2thr = first_nonzero_bin; % basically, number of empty bins from the spike with the lowest amplitude and the closest threshold
        
        % 12. compute correlation between features phase dynamics and
        % spike dynamics
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).AMP_microV_byBin_smoothed);
        data.(L).cc_PSTH_feature(1)    = cc(2,1); % for spike AMP
        data.(L).pp_PSTH_feature(1)    = pp(2,1);
        data.(L).pperm_PSTH_feature(1) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).AMP_microV_byBin_smoothed, cfg.phase.n_shuffles, cfg.phase.tail, cfg.phase.alpha_level, cfg.phase.stat, cfg.phase.reports, cfg.phase.seed_state);
%         [c,p] = corrcoef(data.(L).spike_phases_radians, data.(L).AMP_microV);
        
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).HW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(2)    = cc(2,1); % for HW
        data.(L).pp_PSTH_feature(2)    = pp(2,1);
%         data.(L).pperm_PSTH_feature(2) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).HW_ms_byBin_smoothed, cfg.phase.correlation.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).TPW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(3)    = cc(2,1); % for spike TPW
        data.(L).pp_PSTH_feature(3)    = pp(2,1);
%         data.(L).pperm_PSTH_feature(3) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).TPW_ms_byBin_smoothed, cfg.phase.correlation.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).REP_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(4)    = cc(2,1); % for spike REP
        data.(L).pp_PSTH_feature(4)    = pp(2,1);
%         data.(L).pperm_PSTH_feature(4) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).REP_ms_byBin_smoothed, cfg.phase.correlation.n_permutations, 0, 0.05, 'linear', 0);
        
    end
    save([cfg.cardioballistic_folder filesep data.unitId '_' data.target '__spikes_ECGphase.mat'], 'data', '-v7.3')
    clear data
end
end

function real_data = average_smooth_data(input, cfg)
Kernel = normpdf(-5*cfg.time.gaussian_kernel:cfg.time.PSTH_binwidth:5*cfg.time.gaussian_kernel,0,cfg.time.gaussian_kernel); % build kernel
tmp = [ [zeros(1,80); input(:,1:end-1)'] ...
    input' ...
    [input(:,2:end)'; zeros(1,80)] ];
conv_data = conv2(tmp,Kernel,'same');
real_data = conv_data(:,cfg.phase.N_phase_bins+1:end-cfg.phase.N_phase_bins);
% average_data = mean(conv_data);
% std_data = std(conv_data);
% ste_data = sterr(conv_data);
end

function output = cardioballistic_fit(feature_data, eventPhases, cfg)
scaled_feature = feature_data / nanmean(feature_data);

% compute starting point for scaling factor
a1 = (max(scaled_feature) - min(scaled_feature))/2;

% compute starting point for phase
dropnan = isnan(scaled_feature);
y       = scaled_feature(~dropnan);
ph      = eventPhases(~dropnan);
b1      = mod(circ_mean(ph, y'), 2*pi); % add modulo by 2pi as circ_mean and circ_median can return negative output even having input within 0-2pi

% starting point for intercept
c1 = nanmean(scaled_feature);

startPoint = [a1 b1 c1];

[fittedmdl,gof,~,lin_mdl] = cosine_fit(eventPhases, scaled_feature', cfg, startPoint);

coefs = coeffvalues(fittedmdl);
coefs(2) = mod(coefs(2),2*pi);

% 5 coefficient related to cosine fitting
% - the modulation index, the slope of the cosine function
% - p-value of the modulation index
% - phase of modulation
% - R-squared
% - intercept from linear model
output = [coefs(1) lin_mdl.Coefficients.pValue(2) coefs(2) gof.rsquare lin_mdl.Coefficients.Estimate(1)];

end

function [max_consec_bins, feature_modulation_index] = significant_bins(average_real, lowerPercentile_2_5, upperPercentile_97_5, average_reshuffled)
% figure out significant differences of spike feature dynamics
sig_above = bsxfun(@gt, average_real, upperPercentile_97_5);
sig_below = bsxfun(@lt, average_real, lowerPercentile_2_5);

consec_above = diff([0 find(diff(sig_above)) numel(sig_above)]); % number of repetitions of the same element
consec_below =  diff([0 find(diff(sig_below)) numel(sig_below)]);

if sig_above(1) % if starts with zero
    clust_above = consec_above(1:2:end); % then take lengths of non-zero clusters
else
    clust_above = consec_above(2:2:end);
end

if sig_below(1) % if starts with zero
    clust_below = consec_below(1:2:end);
else
    clust_below = consec_below(2:2:end);
end

max_consec_bins = max([clust_above clust_below]); % max number of consecutive significant bins

if isempty(max_consec_bins)
    max_consec_bins = 0;
end

feature_modulation_index = ...
    (max(average_real) - min(average_real)) / nanmean(average_reshuffled);

end

function [fittedmdl,gof,yfit,lin_mdl] = cosine_fit(x, y, cfg, startPoint_cos)
% drop nans
dropnans = isnan(y);
y        = y(~dropnans);
x        = x(~dropnans);

% do the non-linear fit
[fittedmdl,gof] = fit(double(x),y,cfg.fit.cos_mod,'StartPoint', startPoint_cos, 'Lower', cfg.fit.cos_lower, 'Upper', cfg.fit.cos_upper);

coefs = coeffvalues(fittedmdl); % get model coefficients
coefs(2) = mod(coefs(2),2*pi);

yfit_all = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), x);
yfit     = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), cfg.phase.phase_bin_centers);

% employ a linear fit to get a p-value vs. fitting with a
% constant model
lin_mdl = fitlm(y(:), yfit_all);

end

function spikes_realigned_microV = shift2peak(wf_times_interp_ms, waveforms_upsampled,cfg)
peak_idx = ...
    wf_times_interp_ms > cfg.phase.peak_time_range(1) & wf_times_interp_ms < cfg.phase.peak_time_range(2);
[~, idx] = max(abs(waveforms_upsampled(:,peak_idx)), [], 2); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
spikes_realigned_microV = bsxfun(@(x,y) circshift(x,y,1), waveforms_upsampled', 37 - idx'); % shift to peak (idx 37 for the interpolated data), in microvolts
spikes_realigned_microV = spikes_realigned_microV';
end

function [lowerPercentile_2_5, upperPercentile_97_5, average_reshuffled] = compute_reshuffles(data, bin, cfg)
%% compute reshuffles
rng(0)
[~, reshuffled_spike_order] = sort(rand(cfg.n_permutations, length(data)), 2); % get random order of elements
data_reshuffled      = data(reshuffled_spike_order);
data_reshuffled      = arrayfun(@(x) nanmean(data_reshuffled(:, bin == x),2), 1:cfg.N_phase_bins, 'UniformOutput', false); % mean by phase
data_reshuffled      = cat(2, data_reshuffled{:});
average_reshuffled   = nanmean(data_reshuffled, 1);
lowerPercentile_2_5  = prctile(data_reshuffled, 2.5, 1);
upperPercentile_97_5 = prctile(data_reshuffled, 97.5, 1);
end

% function Y_microV = upsample_spikes(waveforms, wf_times, wf_times_interp)
% %% Luba's interpolation
% y = waveforms;
% y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
% peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
% [~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
% idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
% 
% Y_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts
% end

% function AMP = compute_AMP(waveforms, spike_phases, wf_times, wf_times_interp, phase_bins, nReshuffles)
% %% phase bins
% [AMP.counts,~,AMP.bin] = histcounts(spike_phases, phase_bins);
% AMP.bin = single(AMP.bin);
% 
% %% Luba's interpolation
% y = waveforms;
% y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
% peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
% [~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
% idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
% 
% AMP.waveforms_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts
% 
% AMP.real_microV                   = AMP.waveforms_microV(37,:); % real peak amplitude
% AMP.abs_microV                    = abs(AMP.real_microV); % absolute peak amplitude
% AMP.sign_peak                     = sign(mean(AMP.real_microV)); % signum at the peak
% AMP.mean_by_phase                 = arrayfun(@(x) mean(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % mean by phase
% AMP.mean_by_phase_smoothed        = circ_smooth(AMP.mean_by_phase, 16); % pi/128 * 16 --> pi/8
% AMP.std_by_phase                  = arrayfun(@(x) std(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % standard deviation by phase
% 
% %% compute reshuffles
% [~, reshuffled_spike_order] = sort(rand(nReshuffles, size(waveforms,1)), 2); % get random order of elements
% AMP_reshuffled = AMP.abs_microV(reshuffled_spike_order);
% AMP_reshuffled      = arrayfun(@(x) mean(AMP_reshuffled(:, AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false); % mean by phase
% AMP_reshuffled      = cat(2, AMP_reshuffled{:});
% AMP.lowerPrctile_2_5    = prctile(AMP_reshuffled, 2.5, 1);
% AMP.uppperPrctile_97_5  = prctile(AMP_reshuffled, 97.5, 1);
% 
% % average waveforms by phase and then figure out spike parameters
% AMP.WF_by_phase = arrayfun(@(x) mean(AMP.waveforms_microV(:,AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false);
% AMP.WF_by_phase = cat(2, AMP.WF_by_phase{:});
% AMP.mean_by_bin = AMP.WF_by_phase(37,:)';
% end
