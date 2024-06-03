function ecg_bna_compute_session_ECG_related_spikePhase(trials,population,Rpeaks,cfg)

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
    data.cc_lag_list       = cfg.spk.lag_list;
    data.criteria          = pop.criteria;
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        data.(L).meanHR_bpm                     = single(NaN);
        data.(L).medianHR_bpm                   = single(NaN);
        data.(L).stdHR_bpm                      = single(NaN);
        data.(L).SDNN_ms                        = single(NaN);
        data.(L).spike_phases_radians           = single(NaN);
        data.(L).spike_phases_histogram         = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).spike_phases_histogram2        = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).spike_phases_histogram_smoothed= single(nan(1,cfg.spk.N_phase_bins));
        data.(L).histogram_MI                   = single(NaN);
        data.(L).histogram_p                    = single(NaN);
        data.(L).histogram_phase                = single(NaN);
        data.(L).rsquared                       = single(NaN);
        
        data.(L).linear.yfit                    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).linear.coefs                   = single([NaN NaN]);
        data.(L).linear.rsquared                = single(NaN);
        data.(L).linear.adjrsquared             = single(NaN);
        data.(L).linear.sse                     = single(NaN);
        data.(L).linear.dfe                     = single(NaN);
        data.(L).linear.rmse                    = single(NaN);
        data.(L).linear.pvalue                  = single([NaN; NaN]);
        data.(L).linear.aic                     = single(NaN);
        data.(L).linear.bic                     = single(NaN);
        
        data.(L).linear_smoothed.yfit                    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).linear_smoothed.coefs                   = single([NaN NaN]);
        data.(L).linear_smoothed.rsquared                = single(NaN);
        data.(L).linear_smoothed.adjrsquared             = single(NaN);
        data.(L).linear_smoothed.sse                     = single(NaN);
        data.(L).linear_smoothed.dfe                     = single(NaN);
        data.(L).linear_smoothed.rmse                    = single(NaN);
        data.(L).linear_smoothed.pvalue                  = single([NaN; NaN]);
        data.(L).linear_smoothed.aic                     = single(NaN);
        data.(L).linear_smoothed.bic                     = single(NaN);
        
        data.(L).cosine.average                 = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).cosine.startPoint              = single([NaN NaN NaN]);
        data.(L).cosine.yfit                    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).cosine.coefs                   = single([NaN NaN NaN]);
        data.(L).cosine.rsquared                = single(NaN);
        data.(L).cosine.adjrsquared             = single(NaN);
        data.(L).cosine.sse                     = single(NaN);
        data.(L).cosine.dfe                     = single(NaN);
        data.(L).cosine.rmse                    = single(NaN);
        data.(L).cosine.pvalue                  = single(NaN);
        data.(L).cosine.aic                     = single(NaN);
        data.(L).cosine.bic                     = single(NaN);
        
        data.(L).cosine_smoothed.average                 = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).cosine_smoothed.startPoint              = single([NaN NaN NaN]);
        data.(L).cosine_smoothed.yfit                    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).cosine_smoothed.coefs                   = single([NaN NaN NaN]);
        data.(L).cosine_smoothed.rsquared                = single(NaN);
        data.(L).cosine_smoothed.adjrsquared             = single(NaN);
        data.(L).cosine_smoothed.sse                     = single(NaN);
        data.(L).cosine_smoothed.dfe                     = single(NaN);
        data.(L).cosine_smoothed.rmse                    = single(NaN);
        data.(L).cosine_smoothed.pvalue                  = single(NaN);
        data.(L).cosine_smoothed.aic                     = single(NaN);
        data.(L).cosine_smoothed.bic                     = single(NaN);
        
        data.(L).vonMisesPos.average            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesPos.startPoint         = single([NaN NaN NaN NaN]);
        data.(L).vonMisesPos.yfit               = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesPos.coefs              = single([NaN NaN NaN NaN]);
        data.(L).vonMisesPos.rsquared           = single(NaN);
        data.(L).vonMisesPos.adjrsquared        = single(NaN);
        data.(L).vonMisesPos.sse                = single(NaN);
        data.(L).vonMisesPos.dfe                = single(NaN);
        data.(L).vonMisesPos.rmse               = single(NaN);
        data.(L).vonMisesPos.pvalue             = single(NaN);
        data.(L).vonMisesPos.aic                = single(NaN);
        data.(L).vonMisesPos.bic                = single(NaN);
        
        data.(L).vonMisesPos_smoothed.average            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesPos_smoothed.startPoint         = single([NaN NaN NaN NaN]);
        data.(L).vonMisesPos_smoothed.yfit               = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesPos_smoothed.coefs              = single([NaN NaN NaN NaN]);
        data.(L).vonMisesPos_smoothed.rsquared           = single(NaN);
        data.(L).vonMisesPos_smoothed.adjrsquared        = single(NaN);
        data.(L).vonMisesPos_smoothed.sse                = single(NaN);
        data.(L).vonMisesPos_smoothed.dfe                = single(NaN);
        data.(L).vonMisesPos_smoothed.rmse               = single(NaN);
        data.(L).vonMisesPos_smoothed.pvalue             = single(NaN);
        data.(L).vonMisesPos_smoothed.aic                = single(NaN);
        data.(L).vonMisesPos_smoothed.bic                = single(NaN);
        
        data.(L).vonMisesNeg.average            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesNeg.startPoint         = single([NaN NaN NaN NaN]);
        data.(L).vonMisesNeg.yfit               = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesNeg.coefs              = single([NaN NaN NaN NaN]);
        data.(L).vonMisesNeg.rsquared           = single(NaN);
        data.(L).vonMisesNeg.adjrsquared        = single(NaN);
        data.(L).vonMisesNeg.sse                = single(NaN);
        data.(L).vonMisesNeg.dfe                = single(NaN);
        data.(L).vonMisesNeg.rmse               = single(NaN);
        data.(L).vonMisesNeg.pvalue             = single(NaN);
        data.(L).vonMisesNeg.aic                = single(NaN);
        data.(L).vonMisesNeg.bic                = single(NaN);
        
        data.(L).vonMisesNeg_smoothed.average            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesNeg_smoothed.startPoint         = single([NaN NaN NaN NaN]);
        data.(L).vonMisesNeg_smoothed.yfit               = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).vonMisesNeg_smoothed.coefs              = single([NaN NaN NaN NaN]);
        data.(L).vonMisesNeg_smoothed.rsquared           = single(NaN);
        data.(L).vonMisesNeg_smoothed.adjrsquared        = single(NaN);
        data.(L).vonMisesNeg_smoothed.sse                = single(NaN);
        data.(L).vonMisesNeg_smoothed.dfe                = single(NaN);
        data.(L).vonMisesNeg_smoothed.rmse               = single(NaN);
        data.(L).vonMisesNeg_smoothed.pvalue             = single(NaN);
        data.(L).vonMisesNeg_smoothed.aic                = single(NaN);
        data.(L).vonMisesNeg_smoothed.bic                = single(NaN);
        
        % median split - IBI low
        data.(L).IBI_median                     = single(NaN);
        
        data.(L).lowIBI_timeRRstart             = single(NaN);
        data.(L).lowIBI_timeRRend               = single(NaN);
        data.(L).lowIBI_starts                  = single(NaN);
        data.(L).lowIBI_meanHR_bpm              = single(NaN);
        data.(L).lowIBI_medianHR_bpm            = single(NaN);
        data.(L).lowIBI_stdHR_bpm               = single(NaN);
        data.(L).lowIBI_SDNN_ms                 = single(NaN);
        
        data.(L).lowIBI_linear.yfit             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_linear.coefs            = single([NaN NaN]);
        data.(L).lowIBI_linear.rsquared         = single(NaN);
        data.(L).lowIBI_linear.adjrsquared      = single(NaN);
        data.(L).lowIBI_linear.sse              = single(NaN);
        data.(L).lowIBI_linear.dfe              = single(NaN);
        data.(L).lowIBI_linear.rmse             = single(NaN);
        data.(L).lowIBI_linear.pvalue           = single([NaN; NaN]);
        data.(L).lowIBI_linear.aic              = single(NaN);
        data.(L).lowIBI_linear.bic              = single(NaN);
        
        data.(L).lowIBI_linear_smoothed.yfit             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_linear_smoothed.coefs            = single([NaN NaN]);
        data.(L).lowIBI_linear_smoothed.rsquared         = single(NaN);
        data.(L).lowIBI_linear_smoothed.adjrsquared      = single(NaN);
        data.(L).lowIBI_linear_smoothed.sse              = single(NaN);
        data.(L).lowIBI_linear_smoothed.dfe              = single(NaN);
        data.(L).lowIBI_linear_smoothed.rmse             = single(NaN);
        data.(L).lowIBI_linear_smoothed.pvalue           = single([NaN; NaN]);
        data.(L).lowIBI_linear_smoothed.aic              = single(NaN);
        data.(L).lowIBI_linear_smoothed.bic              = single(NaN);
        
        data.(L).lowIBI_cosine.average          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_cosine.startPoint       = single([NaN NaN NaN]);
        data.(L).lowIBI_cosine.yfit             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_cosine.coefs            = single([NaN NaN NaN]);
        data.(L).lowIBI_cosine.rsquared         = single(NaN);
        data.(L).lowIBI_cosine.adjrsquared      = single(NaN);
        data.(L).lowIBI_cosine.sse              = single(NaN);
        data.(L).lowIBI_cosine.dfe              = single(NaN);
        data.(L).lowIBI_cosine.rmse             = single(NaN);
        data.(L).lowIBI_cosine.pvalue           = single(NaN);
        data.(L).lowIBI_cosine.aic              = single(NaN);
        data.(L).lowIBI_cosine.bic              = single(NaN);
        
        data.(L).lowIBI_cosine_smoothed.average          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_cosine_smoothed.startPoint       = single([NaN NaN NaN]);
        data.(L).lowIBI_cosine_smoothed.yfit             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_cosine_smoothed.coefs            = single([NaN NaN NaN]);
        data.(L).lowIBI_cosine_smoothed.rsquared         = single(NaN);
        data.(L).lowIBI_cosine_smoothed.adjrsquared      = single(NaN);
        data.(L).lowIBI_cosine_smoothed.sse              = single(NaN);
        data.(L).lowIBI_cosine_smoothed.dfe              = single(NaN);
        data.(L).lowIBI_cosine_smoothed.rmse             = single(NaN);
        data.(L).lowIBI_cosine_smoothed.pvalue           = single(NaN);
        data.(L).lowIBI_cosine_smoothed.aic              = single(NaN);
        data.(L).lowIBI_cosine_smoothed.bic              = single(NaN);
        
        data.(L).lowIBI_vonMisesPos.average     = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesPos.startPoint  = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesPos.yfit        = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesPos.coefs       = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesPos.rsquared    = single(NaN);
        data.(L).lowIBI_vonMisesPos.adjrsquared = single(NaN);
        data.(L).lowIBI_vonMisesPos.sse         = single(NaN);
        data.(L).lowIBI_vonMisesPos.dfe         = single(NaN);
        data.(L).lowIBI_vonMisesPos.rmse        = single(NaN);
        data.(L).lowIBI_vonMisesPos.pvalue      = single(NaN);
        data.(L).lowIBI_vonMisesPos.aic         = single(NaN);
        data.(L).lowIBI_vonMisesPos.bic         = single(NaN);
        
        data.(L).lowIBI_vonMisesPos_smoothed.average     = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesPos_smoothed.startPoint  = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesPos_smoothed.yfit        = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesPos_smoothed.coefs       = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesPos_smoothed.rsquared    = single(NaN);
        data.(L).lowIBI_vonMisesPos_smoothed.adjrsquared = single(NaN);
        data.(L).lowIBI_vonMisesPos_smoothed.sse         = single(NaN);
        data.(L).lowIBI_vonMisesPos_smoothed.dfe         = single(NaN);
        data.(L).lowIBI_vonMisesPos_smoothed.rmse        = single(NaN);
        data.(L).lowIBI_vonMisesPos_smoothed.pvalue      = single(NaN);
        data.(L).lowIBI_vonMisesPos_smoothed.aic         = single(NaN);
        data.(L).lowIBI_vonMisesPos_smoothed.bic         = single(NaN);
        
        data.(L).lowIBI_vonMisesNeg.average     = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesNeg.startPoint  = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesNeg.yfit        = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesNeg.coefs       = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesNeg.rsquared    = single(NaN);
        data.(L).lowIBI_vonMisesNeg.adjrsquared = single(NaN);
        data.(L).lowIBI_vonMisesNeg.sse         = single(NaN);
        data.(L).lowIBI_vonMisesNeg.dfe         = single(NaN);
        data.(L).lowIBI_vonMisesNeg.rmse        = single(NaN);
        data.(L).lowIBI_vonMisesNeg.pvalue      = single(NaN);
        data.(L).lowIBI_vonMisesNeg.aic         = single(NaN);
        data.(L).lowIBI_vonMisesNeg.bic         = single(NaN);
        
        data.(L).lowIBI_vonMisesNeg_smoothed.average     = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesNeg_smoothed.startPoint  = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesNeg_smoothed.yfit        = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).lowIBI_vonMisesNeg_smoothed.coefs       = single([NaN NaN NaN NaN]);
        data.(L).lowIBI_vonMisesNeg_smoothed.rsquared    = single(NaN);
        data.(L).lowIBI_vonMisesNeg_smoothed.adjrsquared = single(NaN);
        data.(L).lowIBI_vonMisesNeg_smoothed.sse         = single(NaN);
        data.(L).lowIBI_vonMisesNeg_smoothed.dfe         = single(NaN);
        data.(L).lowIBI_vonMisesNeg_smoothed.rmse        = single(NaN);
        data.(L).lowIBI_vonMisesNeg_smoothed.pvalue      = single(NaN);
        data.(L).lowIBI_vonMisesNeg_smoothed.aic         = single(NaN);
        data.(L).lowIBI_vonMisesNeg_smoothed.bic         = single(NaN);
        
        % median split - IBI high
        data.(L).highIBI_timeRRstart            = single(NaN);
        data.(L).highIBI_timeRRend              = single(NaN);
        data.(L).highIBI_starts                 = single(NaN);
        data.(L).highIBI_meanHR_bpm             = single(NaN);
        data.(L).highIBI_medianHR_bpm           = single(NaN);
        data.(L).highIBI_stdHR_bpm              = single(NaN);
        data.(L).highIBI_SDNN_ms                = single(NaN);
        
        data.(L).highIBI_linear.yfit            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_linear.coefs           = single([NaN NaN]);
        data.(L).highIBI_linear.rsquared        = single(NaN);
        data.(L).highIBI_linear.adjrsquared     = single(NaN);
        data.(L).highIBI_linear.sse             = single(NaN);
        data.(L).highIBI_linear.dfe             = single(NaN);
        data.(L).highIBI_linear.rmse            = single(NaN);
        data.(L).highIBI_linear.pvalue          = single([NaN; NaN]);
        data.(L).highIBI_linear.aic             = single(NaN);
        data.(L).highIBI_linear.bic             = single(NaN);
        
        data.(L).highIBI_linear_smoothed.yfit            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_linear_smoothed.coefs           = single([NaN NaN]);
        data.(L).highIBI_linear_smoothed.rsquared        = single(NaN);
        data.(L).highIBI_linear_smoothed.adjrsquared     = single(NaN);
        data.(L).highIBI_linear_smoothed.sse             = single(NaN);
        data.(L).highIBI_linear_smoothed.dfe             = single(NaN);
        data.(L).highIBI_linear_smoothed.rmse            = single(NaN);
        data.(L).highIBI_linear_smoothed.pvalue          = single([NaN; NaN]);
        data.(L).highIBI_linear_smoothed.aic             = single(NaN);
        data.(L).highIBI_linear_smoothed.bic             = single(NaN);
        
        data.(L).highIBI_cosine.average         = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_cosine.startPoint      = single([NaN NaN NaN]);
        data.(L).highIBI_cosine.yfit            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_cosine.coefs           = single([NaN NaN NaN]);
        data.(L).highIBI_cosine.rsquared        = single(NaN);
        data.(L).highIBI_cosine.adjrsquared     = single(NaN);
        data.(L).highIBI_cosine.sse             = single(NaN);
        data.(L).highIBI_cosine.dfe             = single(NaN);
        data.(L).highIBI_cosine.rmse            = single(NaN);
        data.(L).highIBI_cosine.pvalue          = single(NaN);
        data.(L).highIBI_cosine.aic             = single(NaN);
        data.(L).highIBI_cosine.bic             = single(NaN);
        
        data.(L).highIBI_cosine_smoothed.average         = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_cosine_smoothed.startPoint      = single([NaN NaN NaN]);
        data.(L).highIBI_cosine_smoothed.yfit            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_cosine_smoothed.coefs           = single([NaN NaN NaN]);
        data.(L).highIBI_cosine_smoothed.rsquared        = single(NaN);
        data.(L).highIBI_cosine_smoothed.adjrsquared     = single(NaN);
        data.(L).highIBI_cosine_smoothed.sse             = single(NaN);
        data.(L).highIBI_cosine_smoothed.dfe             = single(NaN);
        data.(L).highIBI_cosine_smoothed.rmse            = single(NaN);
        data.(L).highIBI_cosine_smoothed.pvalue          = single(NaN);
        data.(L).highIBI_cosine_smoothed.aic             = single(NaN);
        data.(L).highIBI_cosine_smoothed.bic             = single(NaN);
        
        data.(L).highIBI_vonMisesPos.average    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesPos.startPoint = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesPos.yfit       = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesPos.coefs      = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesPos.rsquared   = single(NaN);
        data.(L).highIBI_vonMisesPos.adjrsquared= single(NaN);
        data.(L).highIBI_vonMisesPos.sse        = single(NaN);
        data.(L).highIBI_vonMisesPos.dfe        = single(NaN);
        data.(L).highIBI_vonMisesPos.rmse       = single(NaN);
        data.(L).highIBI_vonMisesPos.pvalue     = single(NaN);
        data.(L).highIBI_vonMisesPos.aic        = single(NaN);
        data.(L).highIBI_vonMisesPos.bic        = single(NaN);
        
        data.(L).highIBI_vonMisesPos_smoothed.average    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesPos_smoothed.startPoint = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesPos_smoothed.yfit       = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesPos_smoothed.coefs      = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesPos_smoothed.rsquared   = single(NaN);
        data.(L).highIBI_vonMisesPos_smoothed.adjrsquared= single(NaN);
        data.(L).highIBI_vonMisesPos_smoothed.sse        = single(NaN);
        data.(L).highIBI_vonMisesPos_smoothed.dfe        = single(NaN);
        data.(L).highIBI_vonMisesPos_smoothed.rmse       = single(NaN);
        data.(L).highIBI_vonMisesPos_smoothed.pvalue     = single(NaN);
        data.(L).highIBI_vonMisesPos_smoothed.aic        = single(NaN);
        data.(L).highIBI_vonMisesPos_smoothed.bic        = single(NaN);
        
        data.(L).highIBI_vonMisesNeg.average    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesNeg.startPoint = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesNeg.yfit       = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesNeg.coefs      = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesNeg.rsquared   = single(NaN);
        data.(L).highIBI_vonMisesNeg.adjrsquared= single(NaN);
        data.(L).highIBI_vonMisesNeg.sse        = single(NaN);
        data.(L).highIBI_vonMisesNeg.dfe        = single(NaN);
        data.(L).highIBI_vonMisesNeg.rmse       = single(NaN);
        data.(L).highIBI_vonMisesNeg.pvalue     = single(NaN);
        data.(L).highIBI_vonMisesNeg.aic        = single(NaN);
        data.(L).highIBI_vonMisesNeg.bic        = single(NaN);
        
        data.(L).highIBI_vonMisesNeg_smoothed.average    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesNeg_smoothed.startPoint = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesNeg_smoothed.yfit       = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).highIBI_vonMisesNeg_smoothed.coefs      = single([NaN NaN NaN NaN]);
        data.(L).highIBI_vonMisesNeg_smoothed.rsquared   = single(NaN);
        data.(L).highIBI_vonMisesNeg_smoothed.adjrsquared= single(NaN);
        data.(L).highIBI_vonMisesNeg_smoothed.sse        = single(NaN);
        data.(L).highIBI_vonMisesNeg_smoothed.dfe        = single(NaN);
        data.(L).highIBI_vonMisesNeg_smoothed.rmse       = single(NaN);
        data.(L).highIBI_vonMisesNeg_smoothed.pvalue     = single(NaN);
        data.(L).highIBI_vonMisesNeg_smoothed.aic        = single(NaN);
        data.(L).highIBI_vonMisesNeg_smoothed.bic        = single(NaN);
        
        data.(L).waveforms_microvolts           = single(nan(1,32)); % 32 - number of points per spike
        data.(L).waveforms_upsampled_microvolts = single(nan(1,128)); % 128 - number of points per upsampled spike 4*32
        data.(L).waveforms_byBin_microvolts     = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_microV                     = single(NaN);
        data.(L).AMP_voltageBinned              = single(nan(1,500));
        data.(L).AMP_voltageBins                = single(nan(1,500));
        data.(L).HW_ms                          = single(NaN);
        data.(L).TPW_ms                         = single(NaN);
        data.(L).REP_ms                         = single(NaN);
        data.(L).AMP_microV_byBin               = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_ms_byBin                    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_ms_byBin                   = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_ms_byBin                   = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_reshuffled_avg             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_lowerPrctile_2_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_upperPrctile_97_5          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_reshuffled_avg              = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_lowerPrctile_2_5            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_upperPrctile_97_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_reshuffled_avg             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_lowerPrctile_2_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_upperPrctile_97_5          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_reshuffled_avg             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_lowerPrctile_2_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_upperPrctile_97_5          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_MI                         = single(nan(1,5));
        data.(L).HW_MI                          = single(nan(1,5));
        data.(L).TPW_MI                         = single(nan(1,5));
        data.(L).REP_MI                         = single(nan(1,5));
        data.(L).AMP_MI_binned                  = single(nan(1,5));
        data.(L).HW_MI_binned                   = single(nan(1,5));
        data.(L).TPW_MI_binned                  = single(nan(1,5));
        data.(L).REP_MI_binned                  = single(nan(1,5));
        data.(L).AMP_microV_byBin_smoothed      = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_ms_byBin_smoothed           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_ms_byBin_smoothed          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_ms_byBin_smoothed          = single(nan(1,cfg.spk.N_phase_bins));
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
        data.(L).FRbyRR_Hz                      = single(nan(1,1));
        data.(L).cycleDurations_s               = single(nan(1,1));
        data.(L).pearson_r                      = single(nan(length(cfg.spk.lag_list), 1));
        data.(L).pearson_p                      = single(nan(length(cfg.spk.lag_list), 1));
        data.(L).permuted_p                     = single(nan(length(cfg.spk.lag_list), 1));
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
        % 0. figure out RR-intervals lying within trials
        trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
        trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
        
        RR_within_trial_idx = false(length(valid_RRinterval_starts),1);
        for RRnum = 1:length(RR_within_trial_idx)
            RR_within_trial_idx(RRnum) = ...
                any(valid_RRinterval_starts(RRnum)>trial_starts_one_stream & ...
                valid_RRinterval_ends(RRnum)<trial_ends_one_stream);
        end
        valid_RRinterval_ends      = valid_RRinterval_ends(RR_within_trial_idx); % get rid of RRs beyond the current set of trials
        valid_RRinterval_starts    = valid_RRinterval_starts(RR_within_trial_idx);
        
        % compute parameters of heart activity
        data.(L).meanHR_bpm   = mean(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).medianHR_bpm = median(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).stdHR_bpm    = std(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).SDNN_ms      = std(1000 * (valid_RRinterval_ends - valid_RRinterval_starts));
        
        % implement median split to heart-cycle durations
        RRs     = valid_RRinterval_ends - valid_RRinterval_starts;
        M_IBI   = median(valid_RRinterval_ends - valid_RRinterval_starts);
        
        data.(L).IBI_median           = M_IBI;
        
        data.(L).lowIBI_timeRRstart   = valid_RRinterval_starts(RRs <= M_IBI);
        data.(L).lowIBI_timeRRend     = valid_RRinterval_ends(RRs <= M_IBI);
        data.(L).lowIBI_starts        = RRs(RRs <= M_IBI);
        data.(L).lowIBI_meanHR_bpm    = mean(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_medianHR_bpm  = median(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_stdHR_bpm     = std(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_SDNN_ms       = std(1000 * data.(L).lowIBI_starts);
        
        data.(L).highIBI_timeRRstart  = valid_RRinterval_starts(RRs > M_IBI);
        data.(L).highIBI_timeRRend    = valid_RRinterval_ends(RRs > M_IBI);
        data.(L).highIBI_starts       = RRs(RRs > M_IBI);
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
        [eventPhases, eventsTaken, cycleNums_withSpikes] = DAG_eventPhase(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream); % all data
        [lowIBI_eventPhases, lowIBI_eventsTaken, lowIBI_cycleNums_withSpikes] = DAG_eventPhase(data.(L).lowIBI_timeRRstart, data.(L).lowIBI_timeRRend, AT_one_stream);
        [highIBI_eventPhases, highIBI_eventsTaken, highIBI_cycleNums_withSpikes] = DAG_eventPhase(data.(L).highIBI_timeRRstart, data.(L).highIBI_timeRRend, AT_one_stream);
        
        % check the number of spikes left after computing phases
        if length(eventPhases) < 3 || length(lowIBI_eventPhases) < 3 || length(highIBI_eventPhases) < 3 || ...
                length(valid_RRinterval_starts) < 3 || length(data.(L).lowIBI_timeRRstart) < 3 || length(data.(L).highIBI_timeRRstart) < 3
            continue
        end
        
        % 6. Put results for real data together
        data.(L).spike_phases_radians                 = eventPhases;
        data.(L).spike_phases_histogram               = hist(data.(L).spike_phases_radians, cfg.spk.phase_bin_centers); % compute overal phase histogram
        data.(L).spike_phases_histogram2              = hist3([eventPhases, cycleNums_withSpikes], 'ctrs', {cfg.spk.phase_bin_centers 1:length(valid_RRinterval_starts)});
        
        data.(L).lowIBI_spike_phases_radians          = lowIBI_eventPhases;
        data.(L).lowIBI_spike_phases_histogram        = hist(data.(L).lowIBI_spike_phases_radians, cfg.spk.phase_bin_centers);
        data.(L).lowIBI_spike_phases_histogram2       = hist3([lowIBI_eventPhases, lowIBI_cycleNums_withSpikes], 'ctrs', {cfg.spk.phase_bin_centers 1:length(data.(L).lowIBI_timeRRstart)});
        
        data.(L).highIBI_spike_phases_radians         = highIBI_eventPhases;
        data.(L).highIBI_spike_phases_histogram       = hist(data.(L).highIBI_spike_phases_radians, cfg.spk.phase_bin_centers);
        data.(L).highIBI_spike_phases_histogram2      = hist3([highIBI_eventPhases, highIBI_cycleNums_withSpikes], 'ctrs', {cfg.spk.phase_bin_centers 1:length(data.(L).highIBI_timeRRstart)});
        
        % Mosher's cosine fit the smoothed phase histogram
        [modIndex_phase_hist,phase_hist] = ...
            fitCardiacModulation(cfg.spk.phase_bin_centers, ...
            data.(L).spike_phases_histogram, {'PSTH'}, 0, [221]);
        
        data.(L).spike_phases_histogram_smoothed = phase_hist;
        data.(L).histogram_MI                    = modIndex_phase_hist(1);
        data.(L).histogram_p                     = modIndex_phase_hist(2);
        data.(L).histogram_phase                 = modIndex_phase_hist(3);
        data.(L).rsquared                        = modIndex_phase_hist(4);
        
        %% LINEAR FITS
        % linear fit with 'fitlm' for all the data
        data.(L).linear              = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'linear');
        
        % linear fit with 'fitlm' for smoothed data
        data.(L).linear_smoothed     = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'linear');
        
        % linear fit for the lower IBI
        data.(L).lowIBI_linear       = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'linear');
        
        % linear fit for the higher IBI
        data.(L).highIBI_linear      = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'linear');
        
        %% COSINE FITS
        % cosine fit with 'fit' on all data points instead of smoothed mean
        data.(L).cosine              = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'cosine');
        
        % cosine fit with 'fit' for smoothed data
        data.(L).cosine_smoothed     = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'cosine');
        
        % cosine fit for the lower IBI
        data.(L).lowIBI_cosine       = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'cosine');
        
        % cosine fit for the higher IBI
        data.(L).highIBI_cosine      = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'cosine');
        
        %% POSITIVE VON MISES FITS
        % fit all the data
        data.(L).vonMisesPos          = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', 1);
        
        % fit all the data - smoothed
        data.(L).vonMisesPos_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'vonMises', 1);
        
        % pos. von Mises fits - low IBI
        data.(L).lowIBI_vonMisesPos   = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'vonMises', 1);
        
        % pos. von Mises fits - high IBI
        data.(L).highIBI_vonMisesPos  = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'vonMises', 1);
        
        %% NEGATIVE VON MISES FITS
        % fit all the data
        data.(L).vonMisesNeg          = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', -1);
        
        % fit all the data - smoothed
        data.(L).vonMisesNeg_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'vonMises', -1);
        
        % neg. von Mises fits - low IBI
        data.(L).lowIBI_vonMisesNeg   = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'vonMises', -1);
        
        % neg. von Mises fits - high IBI
        data.(L).highIBI_vonMisesNeg  = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'vonMises', -1);
        
        
        [~, ~, bin] = histcounts(data.(L).spike_phases_radians, cfg.spk.phase_bins);
        
        data.(L).waveforms_microvolts                 = 10^6 * WF_one_stream(eventsTaken,:);
        waveforms_upsampled                           = interpft(data.(L).waveforms_microvolts, 32*4, 2);
        data.(L).waveforms_upsampled_microvolts       = shift2peak(cfg.spk.wf_times_interp_ms, waveforms_upsampled);
        data.(L).waveforms_byBin_microvolts           = arrayfun(@(x) mean(data.(L).waveforms_upsampled_microvolts(bin == x,:),1), 1:cfg.spk.N_phase_bins, 'UniformOutput', false);
        data.(L).waveforms_byBin_microvolts           = cat(1,data.(L).waveforms_byBin_microvolts{:});
        % 7. Calculate spike features with Mosher's procedure
        tic
        sMetric = ...
            struct('extremAmp', cell(length(data.(L).spike_phases_radians),1), ...
            'widthHW', cell(length(data.(L).spike_phases_radians),1), ...
            'widthTP', cell(length(data.(L).spike_phases_radians),1), ...
            'repolTime', cell(length(data.(L).spike_phases_radians),1));
        parfor wfNum = 1:length(data.(L).spike_phases_radians)
            sMetric(wfNum)=spikeWaveMetrics(double(data.(L).waveforms_upsampled_microvolts(wfNum,:)), 37, cfg.spk.Fs*4, 0, [1 0 0 0]); % 37 - index of th peak for updsampled data
        end
        toc
        % 8. put the resulting data together
        data.(L).AMP_microV               = [sMetric.extremAmp];
        data.(L).HW_ms                    = 10^3 * [sMetric.widthHW];
        data.(L).TPW_ms                   = 10^3 * [sMetric.widthTP];
        data.(L).REP_ms                   = 10^3 * [sMetric.repolTime];
        
        data.(L).AMP_microV_byBin         = arrayfun(@(x) nanmean(data.(L).AMP_microV(bin == x)), 1:cfg.spk.N_phase_bins); % mean by phase
        data.(L).HW_ms_byBin              = arrayfun(@(x) nanmean(data.(L).HW_ms(bin == x)), 1:cfg.spk.N_phase_bins);
        data.(L).TPW_ms_byBin             = arrayfun(@(x) nanmean(data.(L).TPW_ms(bin == x)), 1:cfg.spk.N_phase_bins);
        data.(L).REP_ms_byBin             = arrayfun(@(x) nanmean(data.(L).REP_ms(bin == x)), 1:cfg.spk.N_phase_bins);
        
        % 9. estimate significance with bootstrapping
        [data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).AMP_microV, bin, cfg);
        [data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).HW_ms, bin, cfg);
        [data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).TPW_ms, bin, cfg);
        [data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).REP_ms, bin, cfg);
        
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
        modIndex        = fitCardiacModulation(cfg.spk.phase_bin_centers, featureMatrix, {'AMP', 'HW', 'TPW', 'REP'}, 0, 0);
        data.(L).AMP_MI_binned = modIndex(1,:);
        data.(L).HW_MI_binned  = modIndex(2,:);
        data.(L).TPW_MI_binned = modIndex(3,:);
        data.(L).REP_MI_binned = modIndex(4,:);
        
        % store smoothed data for each measure
        if sum(isnan(data.(L).AMP_microV_byBin)) == length(data.(L).AMP_microV_byBin)
            data.(L).AMP_microV_byBin_smoothed    = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).AMP_microV_byBin_smoothed    = circ_smooth(data.(L).AMP_microV_byBin);
        end
        if sum(isnan(data.(L).HW_ms_byBin)) == length(data.(L).HW_ms_byBin)
            data.(L).HW_ms_byBin_smoothed         = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).HW_ms_byBin_smoothed         = circ_smooth(data.(L).HW_ms_byBin);
        end
        if sum(isnan(data.(L).TPW_ms_byBin)) == length(data.(L).TPW_ms_byBin)
            data.(L).TPW_ms_byBin_smoothed        = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).TPW_ms_byBin_smoothed        = circ_smooth(data.(L).TPW_ms_byBin);
        end
        if sum(isnan(data.(L).REP_ms_byBin)) == length(data.(L).REP_ms_byBin)
            data.(L).REP_ms_byBin_smoothed        = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).REP_ms_byBin_smoothed        = circ_smooth(data.(L).REP_ms_byBin);
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
        higher_threshold_detected = ...
            all(rmmissing(data.(L).AMP_microV) > data.thresholds_microV(1));
        lower_threshold_detected = ...
            all(rmmissing(data.(L).AMP_microV) > data.thresholds_microV(2));
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
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).AMP_microV_byBin_smoothed);
        data.(L).cc_PSTH_feature(1)    = cc(2,1); % for spike AMP
        data.(L).pp_PSTH_feature(1)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(1) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).AMP_microV_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).HW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(2)    = cc(2,1); % for HW
        data.(L).pp_PSTH_feature(2)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(2) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).HW_ms_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).TPW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(3)    = cc(2,1); % for spike TPW
        data.(L).pp_PSTH_feature(3)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(3) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).TPW_ms_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).REP_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(4)    = cc(2,1); % for spike REP
        data.(L).pp_PSTH_feature(4)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(4) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).REP_ms_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        % II. Compute unit firing rate per RR-interval
        data.(L).timeRRstart           = valid_RRinterval_starts;
        [data.(L).FRbyRR_Hz, ...
            data.(L).cycleDurations_s] = ...
            computeFRperCycle(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream);
        % compute correlation with different lag
        for lagNum = 1:length(cfg.spk.lag_list)
            % create data variables
            fr_hz = data.(L).FRbyRR_Hz;
            rr_s  = circshift(data.(L).cycleDurations_s, cfg.spk.lag_list(lagNum));
            
            % compute correlation coefficient
            [temp_r, temp_p] = corrcoef(fr_hz, rr_s);
            data.(L).pearson_r(lagNum) = temp_r(2,1);
            data.(L).pearson_p(lagNum) = temp_p(2,1);
            data.(L).permuted_p(lagNum) = mult_comp_perm_corr(fr_hz, rr_s, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
            clear temp_r temp_p
            
        end
    end
    save([cfg.cardioballistic_folder filesep data.unitId '_' data.target '__spikes_ECGphase.mat'], 'data', '-v7.3')
    clear data
end
end

function output = cardioballistic_fit(feature_data, eventPhases, cfg)
scaled_feature = feature_data / nanmean(feature_data);

% compute starting point for scaling factor
a1 = (max(scaled_feature) - min(scaled_feature))/2;

% compute starting point for phase
[y,dropnan]  = rmmissing(scaled_feature);
ph           = eventPhases(~dropnan);
b1 = mod(circ_mean(ph, y'), 2*pi); % add modulo by 2pi as circ_mean and circ_median can return negative output even having input within 0-2pi

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
sig_above = average_real > upperPercentile_97_5;
sig_below = average_real < lowerPercentile_2_5;

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
    (max(average_real) - min(average_real)) / mean(average_reshuffled, 'omitnan');

end

function [fittedmdl,gof,yfit,lin_mdl] = cosine_fit(x, y, cfg, startPoint_cos)
% drop nans
[y,dropnans] = rmmissing(y);
x            = x(~dropnans);

% do the non-linear fit
[fittedmdl,gof] = fit(double(x),y,cfg.fit.cos_mod,'StartPoint', startPoint_cos, 'Lower', cfg.fit.cos_lower, 'Upper', cfg.fit.cos_upper);

coefs = coeffvalues(fittedmdl); % get model coefficients
coefs(2) = mod(coefs(2),2*pi);

yfit_all = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), x);
yfit     = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), cfg.spk.phase_bin_centers);

% employ a linear fit to get a p-value vs. fitting with a
% constant model
lin_mdl = fitlm(y(:), yfit_all);

end

function spikes_realigned_microV = shift2peak(wf_times_interp_ms, waveforms_upsampled)
peak_idx = discretize(wf_times_interp_ms, [0.32 0.48]) == 1;
[~, idx] = max(abs(waveforms_upsampled(:,peak_idx)), [], 2); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
spikes_realigned_microV = bsxfun(@(x,y) circshift(x,y,1), waveforms_upsampled', 37 - idx'); % shift to peak (idx 37 for the interpolated data), in microvolts
spikes_realigned_microV = spikes_realigned_microV';
end

function [FRbyRR_Hz, cycleDurations] = computeFRperCycle(interval_starts, interval_ends, eventTimes)
% make inputs vertical vectors
interval_starts = interval_starts(:);
interval_ends = interval_ends(:);
eventTimes = eventTimes(:);

% Calculate ECG cycle durations
cycleDurations = interval_ends - interval_starts;

% cycleNums = arrayfun(@(x) sum(x > interval_starts & x < interval_ends), eventTimes, 'UniformOutput', false);
spikeCounts = arrayfun(@(x,y) sum(eventTimes > x & eventTimes < y), interval_starts, interval_ends);
FRbyRR_Hz = spikeCounts ./ cycleDurations;

% figure,
% scatter(cycleDurations, FRbyRR_Hz)
end

function [lowerPercentile_2_5, upperPercentile_97_5, average_reshuffled] = compute_reshuffles(data, bin, cfg)
%% compute reshuffles
rng(0)
[~, reshuffled_spike_order] = sort(rand(cfg.spk.n_permutations, length(data)), 2); % get random order of elements
data_reshuffled      = data(reshuffled_spike_order);
data_reshuffled      = arrayfun(@(x) nanmean(data_reshuffled(:, bin == x),2), 1:cfg.spk.N_phase_bins, 'UniformOutput', false); % mean by phase
data_reshuffled      = cat(2, data_reshuffled{:});
average_reshuffled = mean(data_reshuffled, 1, 'omitnan');
lowerPercentile_2_5    = prctile(data_reshuffled, 2.5, 1);
upperPercentile_97_5  = prctile(data_reshuffled, 97.5, 1);
end

function Y_microV = upsample_spikes(waveforms, wf_times, wf_times_interp)
%% Luba's interpolation
y = waveforms;
y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
[~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data

Y_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts
end

function AMP = compute_AMP(waveforms, spike_phases, wf_times, wf_times_interp, phase_bins, nReshuffles)
%% phase bins
[AMP.counts,~,AMP.bin] = histcounts(spike_phases, phase_bins);
AMP.bin = single(AMP.bin);

%% Luba's interpolation
y = waveforms;
y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
[~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data

AMP.waveforms_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts

AMP.real_microV                   = AMP.waveforms_microV(37,:); % real peak amplitude
AMP.abs_microV                    = abs(AMP.real_microV); % absolute peak amplitude
AMP.sign_peak                     = sign(mean(AMP.real_microV)); % signum at the peak
AMP.mean_by_phase                 = arrayfun(@(x) mean(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % mean by phase
AMP.mean_by_phase_smoothed        = circ_smooth(AMP.mean_by_phase, 16); % pi/128 * 16 --> pi/8
AMP.std_by_phase                  = arrayfun(@(x) std(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % standard deviation by phase

%% compute reshuffles
[~, reshuffled_spike_order] = sort(rand(nReshuffles, size(waveforms,1)), 2); % get random order of elements
AMP_reshuffled = AMP.abs_microV(reshuffled_spike_order);
AMP_reshuffled      = arrayfun(@(x) mean(AMP_reshuffled(:, AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false); % mean by phase
AMP_reshuffled      = cat(2, AMP_reshuffled{:});
AMP.lowerPrctile_2_5    = prctile(AMP_reshuffled, 2.5, 1);
AMP.uppperPrctile_97_5  = prctile(AMP_reshuffled, 97.5, 1);

% average waveforms by phase and then figure out spike parameters
AMP.WF_by_phase = arrayfun(@(x) mean(AMP.waveforms_microV(:,AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false);
AMP.WF_by_phase = cat(2, AMP.WF_by_phase{:});
AMP.mean_by_bin = AMP.WF_by_phase(37,:)';
end
