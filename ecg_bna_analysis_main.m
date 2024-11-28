function ecg_bna_analysis_main(project,versions)
%ecg_bna_analysis_main('Pulv_bodysignal',{'ver_LS2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% whether the LFP should be processed (true) or not (false)
% if the LFP for the sessions to analyse has already been processed, and no
% settings need to be changed, this flag can be set to false to skip LFP
% being processed again
% If LFP was not previously processed and the flag is set to false,
% analysis won't happen
% TODO: check if LFP is processed, of not, process LFP even if flag is set
% to false

%% INITIALIZATION
% loop through settings file

cfg = [];
cfg.project = project;
cfg.results_folder = ['Y:\Projects\' cfg.project];
ecg_bna_location     =which('ecg_bna_define_folders');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_folders'])-1);

for v = 1:length(versions)
    cfg.version = versions{v};
    run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
    cfg = ecg_bna_define_folders(cfg);
    
    %% Get info about sessions to be analysed
    % Read the info about sessions to analyse
    sessions_info = cfg.session_info;
    
    %% tmp
    
    
    
    %% temporary - loading data for those and plotting
    %     ecg_bna_modulation_indices_vs_motion_indices(cfg)
    
    %% apply exclusion criteria and save lists of units - we do it once
    if cfg.spk.compute_unit_subsets
        ecg_bna_get_unit_list(cfg,1);
    end
    
    %% per session processing..
    if cfg.process_per_session
        %% %%%% Important for LFP IBI split! don't Change!
        if exist(fullfile(cfg.results_folder,filesep, 'all_sessions_IBIsplit_numRpeaks.mat'),'file')
            allsessions_IBIsplit = load(fullfile(cfg.results_folder,filesep, 'all_sessions_IBIsplit_numRpeaks.mat'));
            allsessions_IBIsplit = allsessions_IBIsplit.allsessions_IBIsplit;
        end
        %%
        for i = 1:length(sessions_info)
            java.lang.System.gc() % added by Luba to control the number of graphical device interface handles (hope to handle the problem with freezing plots while creating figures)
            
            session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
            load(sessions_info(i).Input_trials);
            
            % First make seed and ecg shuffles, then use those shuffles for all subfunctions
            seed_filename=[cfg.ECG_root_results_fldr filesep 'seed.mat']; %% not quite sure yet where to put seed
            if exist(seed_filename,'file')
                load(seed_filename, 'seed');
                rng(seed);
            else
                seed=rng;
                save(seed_filename,'seed');
            end
            
            
            if cfg.process_spikes
                
                %% do ECG spike analysis and computations related to cardioballistic effect
                if cfg.spk.compute_spike_histograms || cfg.spk.compute_spike_phase || cfg.spk.compute_correlation
                    Rpeaks = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),cfg.time);
                    %                     sandbox_bb_vs_waveform(cfg, Rpeaks)
                    %                     Rpeaks=ecg_bna_jitter(sessions_info(i),cfg.spk);
                    cfg.Input_WC=sessions_info(i).Input_WC;
                    load(sessions_info(i).Input_spikes);
                    
%                     ecg_bna_compute_session_phase_response_analysis(trials,population,Rpeaks,cfg)
                    
                    if cfg.spk.compute_spike_histograms
                        ecg_bna_compute_session_spike_histogram(trials,population,Rpeaks,sessions_info(i),cfg);
                    end
                    
                    if cfg.spk.compute_spike_phase
                        ecg_bna_compute_session_ECG_related_spikePhase(trials,population,Rpeaks,sessions_info(i),cfg)
                    end
                    
                    if cfg.spk.compute_correlation
                        ecg_bna_compute_session_correlation_analysis(trials,population,Rpeaks,cfg)
                    end
                end
                
                if cfg.spk.plot_spike_histograms
                    ecg_bna_plot_session_spike_histogram(sessions_info(i),cfg);
                end
                if cfg.spk.plot_spike_phase
                    ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),cfg)
                end
                if cfg.spk.plot_correlation
                    ecg_bna_plot_session_correlation(sessions_info(i),cfg)
                end
                
            end
            if cfg.process_LFP
                if isfield(cfg.lfp,'IBI') && ~isfield(cfg,'PSTH_binwidth') && cfg.lfp.IBI==1
                    [Rpeaks, IBIsplit_concat] = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),cfg.lfp);
                elseif isfield(cfg.lfp,'IBI') && ~isfield(cfg,'PSTH_binwidth') && cfg.lfp.IBI==0
                    Rpeaks= ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),cfg.lfp);
                end
                %Rpeaks=ecg_bna_jitter(sessions_info(i),cfg.lfp);
                
                % new condition to only produce IBI split table in the
                % Rpeak shuffles and save it, rather than compliting the
                % whole analysis:
                if isfield(cfg,'produce_List_only') && (cfg.produce_List_only == 1)
                    
                    allsessions_IBIsplit(i).session = sessions_info(i).Date;
                    
                    if cfg.lfp.IBI==1
                        if cfg.lfp.IBI_low == 1 || cfg.lfp.IBI_high == 0
                            allsessions_IBIsplit(i).Rest_IBI_low_total = sum(IBIsplit_concat.Rest_IBI_low);
                            allsessions_IBIsplit(i).Task_IBI_low_total = sum(IBIsplit_concat.Task_IBI_low);
                            switch length(fieldnames(IBIsplit_concat))
                                case 2
                                    allsessions_IBIsplit(i).Rest_IBI_low_total = sum(IBIsplit_concat.Rest_IBI_low);
                                    allsessions_IBIsplit(i).Task_IBI_low_total = sum(IBIsplit_concat.Task_IBI_low);
                                case 1 && isfield(IBIsplit_concat,'Rest_IBI_high')
                                    allsessions_IBIsplit(i).Rest_IBI_low_total = sum(IBIsplit_concat.Rest_IBI_low);
                                    allsessions_IBIsplit(i).Task_IBI_low_total = [];
                                case 1 && isfield(IBIsplit_concat,'Task_IBI_high')
                                    allsessions_IBIsplit(i).Rest_IBI_low_total = [];
                                    allsessions_IBIsplit(i).Task_IBI_low_total = sum(IBIsplit_concat.Task_IBI_low);
                                otherwise
                                    allsessions_IBIsplit(i).Rest_IBI_low_total = [];
                                    allsessions_IBIsplit(i).Task_IBI_low_total = [];
                            end
                            % Read LFP data
                            cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session_IBIlow');
                            cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site_IBIlow');
                            
                        elseif cfg.lfp.IBI_high == 1 || cfg.lfp.IBI_low == 0
                            switch length(fieldnames(IBIsplit_concat))
                                case 2
                                    allsessions_IBIsplit(i).Rest_IBI_high_total = sum(IBIsplit_concat.Rest_IBI_high);
                                    allsessions_IBIsplit(i).Task_IBI_high_total = sum(IBIsplit_concat.Task_IBI_high);
                                case 1 && isfield(IBIsplit_concat,'Rest_IBI_high')
                                    allsessions_IBIsplit(i).Rest_IBI_high_total = sum(IBIsplit_concat.Rest_IBI_high);
                                    allsessions_IBIsplit(i).Task_IBI_high_total = [];
                                case 1 && isfield(IBIsplit_concat,'Task_IBI_high')
                                    allsessions_IBIsplit(i).Rest_IBI_high_total = [];
                                    allsessions_IBIsplit(i).Task_IBI_high_total = sum(IBIsplit_concat.Task_IBI_high);
                                otherwise
                                    allsessions_IBIsplit(i).Rest_IBI_high_total = [];
                                    allsessions_IBIsplit(i).Task_IBI_high_total = [];
                            end
                            % Read LFP data
                            cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session_IBIhigh');
                            cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site_IBIhigh');
                        end
                    else
                        % Read LFP data
                        cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session');
                        cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site');
                    end
                    %% this is new
                    sitesdir=fileparts(sessions_info(i).Input_LFP{:});
                    [sitefiles]=dir(sessions_info(i).Input_LFP{:});
                    sr=unique([trials.TDT_LFPx_SR]);
                    ts_original=1/sr;
                    
                    continue;
                else
                    fprintf('Analysing for session %s\n', session_name);
                    if cfg.lfp.IBI==1
                        if cfg.lfp.IBI_low == 1 || cfg.lfp.IBI_high == 0
                            % Read LFP data
                            cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session_IBIlow');
                            cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site_IBIlow');
                            
                        elseif cfg.lfp.IBI_high == 1 || cfg.lfp.IBI_low == 0
                            % Read LFP data
                            cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session_IBIhigh');
                            cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site_IBIhigh');
                        end
                        
                    else
                        % Read LFP data
                        cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session');
                        cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site');
                    end
                    %% this is new
                    sitesdir=fileparts(sessions_info(i).Input_LFP{:});
                    [sitefiles]=dir(sessions_info(i).Input_LFP{:});
                    sr=unique([trials.TDT_LFPx_SR]);
                    ts_original=1/sr;
                    
                    for s = 1:length(sitefiles)
                        load([sitesdir filesep sitefiles(s).name], 'sites');
                        site_LFP= ecg_bna_process_LFP(sites, cfg, ts_original);
                        n_LFP_samples_per_block=site_LFP.tfs.n_samples_per_block;
                        
                        
                        %% CHECK BLOCKS OVERLAP - remove LFP for which there is no trigger block - we did the opposite already in resample_triggers!
                        blocks_not_present_in_triggers=n_LFP_samples_per_block(1, ~ismember(n_LFP_samples_per_block(1,:),[Rpeaks.block]));
                        if numel(blocks_not_present_in_triggers)==size(n_LFP_samples_per_block,2)
                            continue;
                        end
                        
                        tfr_samples_to_remove=[];
                        for b=1:numel(blocks_not_present_in_triggers)
                            B=blocks_not_present_in_triggers(b);
                            start_to_remove=sum(n_LFP_samples_per_block(2,n_LFP_samples_per_block(1,:)<B))+1;
                            end_to_remove=n_LFP_samples_per_block(2,n_LFP_samples_per_block(1,:)==B)+start_to_remove-1;
                            b_idx=[site_LFP.block]==B;
                            FN_tr={'LFP_samples','dataset','block','run','n','LFP_tStart','LFP_t0'};
                            for f=1:numel(FN_tr)
                                FN=FN_tr{f};
                                site_LFP.(FN)(b_idx)=[];
                            end
                            tfr_samples_to_remove=[tfr_samples_to_remove start_to_remove:end_to_remove];
                        end
                        FN_tfs={'phabp','powbp','pha','pow','lfp'};
                        for f=1:numel(FN_tfs)
                            FN=FN_tfs{f};
                            site_LFP.tfs.(FN)(:,tfr_samples_to_remove)=[];
                        end
                        
                        blockstarts=[1, find(diff([trials.block]))+1];
                        t_offset_per_block=[[trials(blockstarts).block];[trials(blockstarts).TDT_LFPx_tStart]];
                        triggers.ecg.shuffled = ecg_bna_resample_triggers(Rpeaks,'shuffled_ts',t_offset_per_block,n_LFP_samples_per_block,1/site_LFP.tfs.sr);
                        triggers.ecg.real     = ecg_bna_resample_triggers(Rpeaks,'RPEAK_ts',t_offset_per_block,n_LFP_samples_per_block,1/site_LFP.tfs.sr);
                        
                        site_data=ecg_bna_compute_session_Rpeak_triggered_variables( site_LFP,triggers,trials,cfg );
                        triggered_session_data.sites(s) = site_data;
                        triggered_session_data.session = site_data.session;
                    end
                    
                    % make a folder to save figures
                    session_result_folder = fullfile(cfg.session_lfp_fldr);
                    if ~exist(session_result_folder, 'dir')
                        mkdir(session_result_folder);
                    end
                    save(fullfile(session_result_folder, ['Triggered_session_' triggered_session_data.session '.mat']), 'triggered_session_data');
                    clear session_proc_lfp site_data triggered_session_data;
                    
                end
            end
        end
        allresult_folder = fullfile(cfg.results_folder);
        if ~exist(allresult_folder, 'dir')
            mkdir(allresult_folder);
        end
%         save(fullfile(allresult_folder,filesep, 'all_sessions_IBIsplit_numRpeaks.mat'), 'allsessions_IBIsplit');
    end
    
    if cfg.process_spikes
        %% additionaly exclude by R-peak number and (in the future) other heart-related criteria
        if cfg.spk.ecg_exclusion_criteria
            ecg_bna_get_unit_list_ecg_params(cfg)
        end
        
        %% copy selected units separately - we do it once
        if cfg.spk.move_files
            %             % move units for either task or rest
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_600', cfg.per_session_folder, cfg.per_session_selected, cfg)
            %             % move files for units with no cardioballistic effect
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_selected_600_noCB'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_noCB', cfg.per_session_folder, ...
            %                 output_folder, cfg)
            %
            %             % move files for units WITH the cardioballistic effect
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_selected_600_withCB'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_withCB', cfg.per_session_folder, ...
            %                 output_folder, cfg)
            %
            %             % move files for units with no significant cardioballitic
            %             % effect AND those that had huge amplitude
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_selected_600_noCB_corr'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_noCB_corr', cfg.per_session_folder, ...
            %                 output_folder, cfg)
            %
            %             % move units for both task and rest without cardioballistic
            %             % effect, high spike amplitude and non-significant cc between
            %             % phase PSTH and AMP phase dynamic
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_selected_600_noCB_corr_ccs'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_noCB_corr_ccs', cfg.per_session_folder, ...
            %                 output_folder, cfg)
            %
            %             % both task and rest - high amplitudes
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_selected_600_highAmp'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_highAmp', cfg.per_session_folder, ...
            %                 output_folder, cfg)
            
            
            % ><\\\'> ??? <'///><
            
            
            % move units for both task and rest
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_600', [cfg.per_session_folder '_-0.25-0.25s'], cfg.per_session_stable, cfg)
            
            % move units for both task and rest without cardioballistic
            % effect
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_noCB'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noCB', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % move units for both task and rest with cardioballistic effect
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_withCB'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_withCB', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % move units for both task and rest without cardioballistic
            % effect and with high spike amplitude
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_noCB_corr'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noCB_corr', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % move units for both task and rest without cardioballistic
            % effect, high spike amplitude and non-significant cc between
            % phase PSTH and AMP phase dynamic
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_noCB_corr_ccs'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noCB_corr_ccs', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % high amplitude
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_high_amplitude'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_high_amplitude', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % low amplitude
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_low_amplitude'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_low_amplitude', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % [any] low amp + ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_low_amplitude_ccs_any'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_low_amplitude_ccs_any', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % [both] low amp + ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_low_amplitude_ccs_both'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_low_amplitude_ccs_both', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % [any] no low amp + ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_noLow_amplitude_ccs_any'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_any', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            % [both] no low amp + ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'per_unit_stable_600_noLow_amplitude_ccs_both'];
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_both', [cfg.per_session_folder '_-0.25-0.25s'], ...
                output_folder, cfg)
            
            
            % ><\\\'> ??? <'///><
            
            %% move cardioballistic files
            %             % for either task or rest
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_600', cfg.cardioballistic_folder, cfg.cardioballistic_selected, cfg)
            %
            %             % for units with no cardioballistic effect
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_selected_600_noCB'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_noCB', cfg.cardioballistic_folder, output_folder, cfg)
            %
            %             % for units WITH cardioballistic effect
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_selected_600_withCB'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_withCB', cfg.cardioballistic_folder, output_folder, cfg)
            %
            %             % for units with no significant cardioballitic effect AND those that had huge amplitude
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_selected_600_noCB_corr'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_noCB_corr', cfg.cardioballistic_folder, output_folder, cfg)
            %
            %             % move units for both task and rest without cardioballistic
            %             % effect, high spike amplitude and non-significant cc between
            %             % phase PSTH and AMP phase dynamic
            %             output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_selected_600_noCB_corr_ccs'];
            %             ecg_bna_copy_selected_units('unitInfo_after_exclusion_noCB_corr_ccs', cfg.cardioballistic_folder, output_folder, cfg)
            
            % for both task and rest
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_600', cfg.cardioballistic_folder, cfg.cardioballistic_stable, cfg)
            
            % for units with no cardioballistic effect
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_noCB'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noCB', cfg.cardioballistic_folder, output_folder, cfg)
            
            % for units WITH cardioballistic effect
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_withCB'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_withCB', cfg.cardioballistic_folder, output_folder, cfg)
            
            % for units with no significant cardioballitic effect AND those that had huge amplitude
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_noCB_corr'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noCB_corr', cfg.cardioballistic_folder, output_folder, cfg)
            
            % move units for both task and rest without cardioballistic
            % effect, high spike amplitude and significant cc between phase
            % PSTH and AMP phase dynamic
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_noCB_corr_ccs'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noCB_corr_ccs', cfg.cardioballistic_folder, output_folder, cfg)
            
            % high amplitude
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_high_amplitude'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_high_amplitude', cfg.cardioballistic_folder, output_folder, cfg)
            
            % low amplitude
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_low_amplitude'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_low_amplitude', cfg.cardioballistic_folder, output_folder, cfg)
            
            % low amplitude + any ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_low_amplitude_ccs_any'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_low_amplitude_ccs_any', cfg.cardioballistic_folder, output_folder, cfg)
            
            % low amplitude + both ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_low_amplitude_ccs_both'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_low_amplitude_ccs_both', cfg.cardioballistic_folder, output_folder, cfg)
            
            
            % no low amplitude + any ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_noLow_amplitude_ccs_any'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_any', cfg.cardioballistic_folder, output_folder, cfg)
            
            % no low amplitude + both ccs
            output_folder = [cfg.SPK_root_results_fldr filesep 'cardioballistic_stable_600_noLow_amplitude_ccs_both'];
            if ~exist(output_folder, 'dir')
                mkdir(output_folder)
            end
            ecg_bna_copy_selected_units('unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_both', cfg.cardioballistic_folder, output_folder, cfg)
            
            
        end
    end
    %% average across sessions
    if cfg.process_population
        %
        %         % ==================================================================================================
        %         %            Temp, for checking similar selected units of Luba
        %         % ==================================================================================================
        %         monkeys = unique({cfg.session_info.Monkey});
        %         load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_',monkeys{1},...
        %             '_TaskRest\unit_lists\unitInfo_after_exclusion_stableTaskAndRest_noCB_corr.mat']);
        % %         load(['Y:\Projects\Pulv_bodysignal\ephys\ECG_TaskRest_',monkeys{1},...
        % %             '_merged\tuning_table_combined_CI.mat']);
        % %
        %         site_ID_idx  = contains(keys.tuning_table(2:end,1),unit_ids);
        %         cfg.site_IDS = keys.tuning_table(site_ID_idx,find_column_index(keys.tuning_table,'site_ID'));
        %         % ==================================================================================================
        %         % ==================================================================================================
        
        if cfg.process_LFP
            keys=ecg_bna_get_unit_list(cfg,0);
            cfg.site_IDS=keys.tuning_table(2:end,find_column_index(keys.tuning_table,'site_ID'));
            
            monkeys = unique({cfg.session_info.Monkey});
            cfg.monkey = [monkeys{:}];
            if isfield(cfg.lfp,'IBI') && cfg.lfp.IBI==1
                if cfg.lfp.IBI_low == 1 || cfg.lfp.IBI_high == 0
                    % Read LFP data
                    cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session_IBIlow');
                    cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site_IBIlow');
                    cfg.grand_avg_fldr   = fullfile(cfg.analyse_lfp_folder, 'grand_average_IBIlow');
                    
                elseif cfg.lfp.IBI_high == 1 || cfg.lfp.IBI_low == 0
                    % Read LFP data
                    cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session_IBIhigh');
                    cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site_IBIhigh');
                    cfg.grand_avg_fldr   = fullfile(cfg.analyse_lfp_folder, 'grand_average_IBIhigh');
                end
                
            else
                % Read LFP data
                cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session');
                cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site');
                cfg.grand_avg_fldr   = fullfile(cfg.analyse_lfp_folder, 'grand_average_all');
            end
            
%             grand_avg = ecg_bna_compute_grand_avg(cfg,'w_units');
%             grand_avg = ecg_bna_compute_grand_avg(cfg,'wo_units');
%             grand_avg = ecg_bna_compute_grand_avg(cfg,'all');
            withunits_list = {'all','w_units','wo_units'};
            shad_list = {'percentile', 'std'};
            for shad = 1:2
                for withunits = 1:length(withunits_list)
                    grand_avg = ecg_bna_compute_grand_avg_DiffNormMethods(cfg,withunits_list{withunits},'zscore','abs',shad_list{shad}); % shading options: 'std', 'percentile'
                    grand_avg = ecg_bna_compute_grand_avg_DiffNormMethods(cfg,withunits_list{withunits},'raw','raw',shad_list{shad});
                    grand_avg = ecg_bna_compute_grand_avg_DiffNormMethods(cfg,withunits_list{withunits},'subtraction','raw',shad_list{shad});
                    grand_avg = ecg_bna_compute_grand_avg_DiffNormMethods(cfg,withunits_list{withunits},'mean_division','raw',shad_list{shad});
                    grand_avg = ecg_bna_compute_grand_avg_DiffNormMethods(cfg,withunits_list{withunits},'std_division','raw',shad_list{shad});
                end
            end
            
            if isfield(cfg.lfp , 'TaskRest_diff') && cfg.lfp.TaskRest_diff==1
                Task_Rest_raw_diff = ecg_bna_compute_TaskRest_raw_diff(cfg, 'all');
                Task_Rest_raw_diff = ecg_bna_compute_TaskRest_raw_diff(cfg, 'w_units');
                Task_Rest_raw_diff = ecg_bna_compute_TaskRest_raw_diff(cfg, 'wo_units');
            end
%             
            IBI_high_grand_avg_results_file = fullfile([cfg.analyse_lfp_folder filesep 'grand_average_IBIhigh' filesep cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Condition_diff_grand_avg_sessions_sites','all','.mat']);
            IBI_low_grand_avg_results_file = fullfile([cfg.analyse_lfp_folder filesep 'grand_average_IBIlow' filesep cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Condition_diff_grand_avg_sessions_sites','all','.mat']);
            
            if isfield(cfg.lfp , 'IBI_diff') && exist(IBI_low_grand_avg_results_file,'file') && exist(IBI_high_grand_avg_results_file,'file')
               ecg_bna_compute_IBI_raw_diff(cfg, 'all', cfg.lfp.IBIdiff_type);
               ecg_bna_compute_IBI_raw_diff(cfg, 'w_units', cfg.lfp.IBIdiff_type);
               ecg_bna_compute_IBI_raw_diff(cfg, 'wo_units', cfg.lfp.IBIdiff_type);
            end
            
%             if isfield(cfg.lfp, 'cluster_perm') && cfg.lfp.cluster_perm==1 && isfield(cfg.lfp , 'IBI_diff') && exist(IBI_low_grand_avg_results_file,'file') && exist(IBI_high_grand_avg_results_file,'file')
            if isfield(cfg.lfp, 'cluster_perm') && cfg.lfp.cluster_perm==1 
               ecg_bna_cluster_perm_condition_diff(cfg,'all','itpc');
               ecg_bna_cluster_perm_condition_diff(cfg,'all','pow');
%                ecg_bna_cluster_perm_condition_diff(cfg,'all','lfp');
            end
        end
        
        if cfg.process_spikes
            
%             ecg_bna_avg_spike_histogram_clean(cfg, 'per_unit_-0.25-0.25s', 'unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_any')
            
%             ecg_bna_population_correlation_analysis(cfg, 'correlation_analysis', 'unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_any')
            
            % loop though all selection lists
            for listNum = 1:length(cfg.pop.unit_selection_lists)
                
                %% time-domain analysis
                ecg_bna_avg_spike_histogram_clean(cfg, 'per_unit_-0.25-0.25s', cfg.pop.unit_selection_lists{listNum})
                
                %% phase-domain analysis will be here
                
                
                %% correlation analysis
%                 ecg_bna_population_correlation_analysis(cfg, 'correlation_analysis', cfg.pop.unit_selection_lists{listNum})
                
                %% plot 
%                 ecg_bna_plot_averageHR(cfg, cfg.pop.unit_selection_lists{listNum})
                
            end
            
%             ecg_bna_population_cardioballistic(cfg, 'cardioballistic', 'unitInfo_after_SNR_exclusion_selected_noLow_amplitude_ccs_any', 'Population_noLowAmp_ccs_any_cardioballistic')
%             
%             ecg_bna_plot_circular_fits(cfg, 'per_unit_-0.25-0.25s', 'Circular_population_results_0-0.5s')
%             
%             ecg_bna_plot_venns(cfg)
%             
%             % stable 600
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_600');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg, 'SPK_root_results_fldr', '', 'cardioballistic_stable_600', 'data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 noCB
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_noCB');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_noCB','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 withCB
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_withCB');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_withCB','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 noCB corr
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_noCB_corr');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_noCB_corr','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 noCB corr ccs
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_noCB_corr_ccs');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_noCB_corr_ccs','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 + high amplitude
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_high_amplitude');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_high_amplitude','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 + low amplitude
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_low_amplitude');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_low_amplitude','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 + low amp + any ccs
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_low_amplitude_ccs_any');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_low_amplitude_ccs_any','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 + low amp + both ccs
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_low_amplitude_ccs_both');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_low_amplitude_ccs_both','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 + low amp + any ccs
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_noLow_amplitude_ccs_any');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_noLow_amplitude_ccs_any','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             % stable 600 + low amp + both ccs
%             output_folder = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic_stable_noLow_amplitude_ccs_both');
%             SPK_cardioballistic=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','', 'cardioballistic_stable_600_noLow_amplitude_ccs_both','data');
%             ecg_bna_population_cardioballistic(SPK_cardioballistic, output_folder, cfg)
%             
%             ecg_bna_plot_circular_fits(cfg, 'per_unit_0-0.5s', 'Circular_population_results_0-0.5s')
%             ecg_bna_plot_circular_fits(cfg, 'per_unit_-0.25-0.25s', 'Circular_population_results_-0.25-0.25s')
%             ecg_bna_plot_circular_fits(cfg, 'per_unit_-0.5-0s', 'Circular_population_results_-0.5-0s')
            
        end
        
    end
end
end

function Out = load_stuff(sessions_info,cfg,subfolder,namepart,per,varname)
out_idx=0;
for i = 1:length(sessions_info)
    monkey=sessions_info(i).Monkey(1:3);
    
    results_folder = fullfile(cfg.(subfolder),per);
    file=dir(fullfile(results_folder, [namepart monkey '_' sessions_info(i).Date '*.mat']));
    for f=1:numel(file)
        out_idx=out_idx+1;
        to_load = load([results_folder filesep file(f).name],varname);
        Out(out_idx)=to_load.(varname);
    end
end
end