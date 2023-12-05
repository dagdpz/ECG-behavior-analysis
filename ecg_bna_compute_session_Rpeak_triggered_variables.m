function [ triggered_site_data ] = ecg_bna_compute_session_Rpeak_triggered_variables( site_proc , triggers_in, trials,cfg )
% ecg_bna_compute_session_Rpeak_triggered_variables  - compute evoked_LFP,
% Powspctrm, ITPC, and phaseBP variables for different conditions for each site of a session.
% A condition is a combination of
% possibly hand-space tuning (for consitency with LFP analysis),
% control/inactivation, choice/instructed, type-effector values.
%
% USAGE:
%	[ session_data ] = ecg_bna_compute_session_Rpeak_triggered_variables(
%	site_proc, events, cfg )
%
% INPUTS:
%		site_proc  	- 1xN struct containing raw LFP data for a
%		session,  output from ecg_bna_process_combined_LFP_ECG
%       events  - cell array containing states to be
%       analysed and corresponding time windows
%       cfg     - struct containing configuration settings
%           Required fields:
%               random_seed                 - random seed for
%               reproducibility of random shuffling of Rpeaks
%               session_results_fldr        - folder to which the
%               results of the session should be saved
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               diff_condition      - conditions to compare, the plot
%               for compared conditions would be shown one on top of the
%               other
%               ref_hemisphere      - reference hemisphere for ipsi- and
%               contra- hand and space labeling
%           Optional Fields:
%               diff_color          - color to be used for plotting the
%               compared conditions
%               diff_legend         - legend to be used while plotting the
%               compared conditions
%               random_permute_triggers     - flag indicating whether to
%               randomly permute the Rpeak triggers
%               n_shuffles                  - integer which specifies how
%               many times the Rpeaks need to be randomly shuffled to
%               compute statistics
%
%
% OUTPUTS:
%		session_data	- output structure which saves the average
%       evoked LFP in a time window around Rpeak for trials of given
%       conditions
%
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_get_condition_trials,
% ecg_bna_get_Rpeak_evoked_LFP, ecg_bna_get_shuffled_Rpeak_evoked_LFP,
% ecg_bna_plot_evoked_lfp, ecg_bna_compute_diff_condition_average
%
% See also ecg_bna_compute_session_evoked_ECG,
% ecg_bna_compute_session_Rpeak_evoked_TFS,
% ecg_bna_compute_session_evoked_ECG_R2Rt,
% ecg_bna_compute_session_Rpeak_evoked_state_onsets

% suppress warning for xticklabel
warning ('off', 'MATLAB:hg:willberemoved');

% folder to save sitewise results
site_results_folder = fullfile(cfg.sites_lfp_fldr);
if ~exist(site_results_folder, 'dir')
    mkdir(site_results_folder);
end
resampling_factor=site_proc.tfs.resampling_factor;
LFP_samples_end=round(cumsum(site_proc.LFP_samples)*resampling_factor);
LFP_samples_start=[0 LFP_samples_end(1:end-1)]+1;

%% get triggers for this site
trig_fn=fieldnames(triggers_in);
for f=1:numel(trig_fn)
    FN=trig_fn{f};
    current=triggers_in.(FN);
    triggers.([FN '_real'])=[current.real.event_sample];
    if isfield(current,'shuffled')
        triggers.([FN '_shuffled'])=[current.shuffled.event_sample];
    end
end

%% convert to ipsi/contra
sitetrials=ph_get_unit_trials(site_proc,trials);
positions=[sitetrials.tar_pos]-[sitetrials.fix_pos];
hemifields=num2cell(sign(real(positions)));
fixations=num2cell([sitetrials.fix_pos]);
positions=num2cell(positions);
[sitetrials.position]=deal(positions{:});
[sitetrials.hemifield]=deal(hemifields{:});
[sitetrials.fixation]=deal(fixations{:});
sitetrials=ph_LR_to_CI(cfg,site_proc,sitetrials);  %% convert...

%% take over IDs
trig.condition = struct();
trig.site_ID = site_proc.site_ID;
trig.session = site_proc.session;
trig.target  = site_proc.target;

% loop through conditions

for cn = 1:length(cfg.condition)
    % store details of analysed condition
    trig.condition(cn).label    = cfg.conditionname{cn};    
    
    % get trial indices for the given condition
    cond_trials = ecg_bna_get_condition_trials(sitetrials, cfg.condition(cn));
%     %% this needs to go because its WRONG!
%     for tr = 1:length(sitetrials)
%         state_onset_values = [sitetrials(tr).states_onset];
%         cond_trials_missing_timing(tr) = ~any(isnan(state_onset_values));
%     end
%     cond_trials=cond_trials&cond_trials_missing_timing;
    
    %         %% FIX THIS BS
    %         trig.condition(cn).ntrials(hs)     = sum(cond_trials);
    %         trig.condition(cn).noisytrials(hs) = sum(cond_trials & [site_proc.trials.noisy]);
    %         fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
    %         fprintf('Total number of trials %g\n', sum(cond_trials));
    %         fprintf('Number of noisy trials %g\n', sum(cond_trials & [site_proc.trials.noisy]));
    %         cond_trials = cond_trials & ~[site_proc.trials.noisy];
    
    if sum(cond_trials) == 0
        continue;
    end
    % fprintf('Condition Trials = %i , Num Trial_idx = %i, CIX = %i\n',numel(cond_trials),sum(trial_idx), sum(cix))
    
    % loop through time windows around the states to analyse
    for e = 1:size(cfg.analyse_states, 1)
        state=cfg.analyse_states(e,:);
        event_name=state{1};
        width_in_samples=[floor([state{3}]*site_proc.tfs.sr) ceil([state{4}]*site_proc.tfs.sr)]; %% time to samples!
        
        LFP_samples_end_con=LFP_samples_end(cond_trials)-width_in_samples(2);
        LFP_samples_start_con=LFP_samples_start(cond_trials)-width_in_samples(1);
        notenoughsamples=LFP_samples_end_con<LFP_samples_start_con;
        LFP_samples_start_con(notenoughsamples)=[];
        LFP_samples_end_con(notenoughsamples)=[];
        %         continous_ends=LFP_samples_start_con(2:end)-LFP_samples_end_con(1:end-1)==1;
        %         LFP_samples_end_con([continous_ends false])=[];
        %         LFP_samples_start_con([false continous_ends])=[];
        
        time=[width_in_samples(1):width_in_samples(2)]/site_proc.tfs.sr;
        
        
        trig.condition(cn).event(e).trials = find(cond_trials); %find(cond_trials);
        trig.condition(cn).event(e).ntrials = sum(cond_trials);
        trig.condition(cn).event(e).event_name = event_name;
        trig.condition(cn).event(e).time=time;
        trig.condition(cn).event(e).tfr_time=time;
        
        trig_con_s=triggers.([event_name '_real']);
        trig_con_s(trig_con_s<LFP_samples_start_con(1)-width_in_samples(1))=0;
        trig_con_s(trig_con_s>LFP_samples_end_con(end)-width_in_samples(2))=0;
        for t=1:numel(LFP_samples_end_con)-1
            inbetween=trig_con_s<LFP_samples_start_con(t+1) & trig_con_s>LFP_samples_end_con(t);
            trig_con_s(inbetween)=0;
        end
        realD = ecg_bna_get_triggered_parameters(site_proc,trig_con_s, width_in_samples);
        
        %% change to dependent on state
        if true % keep for now, but probably we are going to ALWAYS shuffle
            %% compute shuffled power spectra, ITPC spectra, lfp, and bandpassed ITPC:
            trig_con_s=triggers.([event_name '_shuffled']);
            trig_con_s(trig_con_s<LFP_samples_start_con(1)-width_in_samples(1))=0;
            trig_con_s(trig_con_s>LFP_samples_end_con(end)-width_in_samples(2))=0;
            for t=1:numel(LFP_samples_end_con)-1
                inbetween=trig_con_s<LFP_samples_start_con(t+1) & trig_con_s>LFP_samples_end_con(t);
                trig_con_s(inbetween)=0;
            end
            shuffledD = ecg_bna_get_triggered_parameters(site_proc,trig_con_s,width_in_samples);
        else % some sort of dummies
            shuffledD.pow.mean=zeros(size(realD.pow.mean));
            shuffledD.pow.std =zeros(size(realD.pow.std));
            shuffledD.itpc.mean=zeros(size(realD.itpc.mean));
            shuffledD.itpc.std =zeros(size(realD.itpc.std));
            shuffledD.lfp.mean=zeros(size(realD.lfp.mean));
            shuffledD.lfp.std =zeros(size(realD.lfp.std));
            shuffledD.itpcbp.mean=zeros(size(realD.itpcbp.mean));
            shuffledD.itpcbp.std =zeros(size(realD.itpcbp.std));
            shuffledD.powbp.mean=zeros(size(realD.powbp.mean));
            shuffledD.powbp.std =zeros(size(realD.powbp.std));
            shuffledD.ntriggers =0;
        end
        
        normalized = ecg_bna_compute_shufflePredictor_normalization_general(realD,shuffledD,cfg);
        significance = ecg_bna_compute_significance(realD,shuffledD,cfg);
        
        trig.condition(cn).event(e).real=realD;
        trig.condition(cn).event(e).shuffled=shuffledD;
        trig.condition(cn).event(e).normalized=normalized;
        trig.condition(cn).event(e).significance=significance;
        
    end
end

% plots - if we don't shuffle, there will be no shuffled!
methods= {'real','shuffled','normalized'};
for mt = 1: numel(methods)
    ecg_bna_plots_per_site( trig, cfg, methods{mt}) % per site!
    % Note: ===> last input could be 'real', 'shuffled', or 'normalized'
end

triggered_site_data = trig;
save(fullfile(site_results_folder, [trig.site_ID '.mat']), 'triggered_site_data');
close all;
end