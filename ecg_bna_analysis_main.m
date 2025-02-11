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

sanitycheck=0;
for v = 1:length(versions)
    cfg.version = versions{v};
    run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
    cfg = ecg_bna_define_folders(cfg);
    
    %% Get info about sessions to be analysed    
    sessions_info = cfg.session_info;
    
    %% per session processing..
    if cfg.process_per_session       
        for i = 1:length(sessions_info)
            % First make same seed for each session -> that way even if
            % the data changes in a subsequent re-running, shuffles in
            % unaffected sessions are unaffected...
            seed_filename=[cfg.ECG_root_results_fldr filesep 'seed.mat']; %% not quite sure yet where to put seed
            if exist(seed_filename,'file')
                load(seed_filename);
                rng(seed);
            else
                seed=rng;
                save(seed_filename,'seed');
            end
            
            java.lang.System.gc() % added by Luba to control the number of graphical device interface handles (hope to handle the problem with freezing plots while creating figures)
                        
            load(sessions_info(i).Input_trials);            
            
            % reading in actual TDT clock block starts (in seconds - inprecise, but that is irrelevant)
            blocks=unique([trials.block]);
            blockstart=ecg_bna_get_anchor_times(monkey,sessions_info(i).Date,blocks); 
                        
            cfg.event_types=cfg.analyse_states(:,2);
            cfg.events=cfg.analyse_states(:,1);
            for e=1:numel(cfg.events)
                current_event=cfg.analyse_states(e,:);
                Event=ecg_bna_get_events(sessions_info(i),trials,current_event,blockstart);
                Triggers.(cfg.events{e})=ecg_bna_jitter(Event,cfg.spk);
            end
                    
            
            if cfg.process_spikes
                cfg.Input_WC=sessions_info(i).Input_WC;
                
                %% apply exclusion criteria and save lists of units - we do it once
                if cfg.spk.compute_unit_subsets
                    ecg_bna_get_unit_list(cfg,1);
                end                
                %% copy selected units separately - we do it once
                if cfg.spk.move_files
                    ecg_bna_copy_selected_units(cfg)
                end
                
                %% do ECG spike analysis and computations related to cardioballistic effect
                if cfg.spk.compute_spike_histograms || cfg.spk.compute_spike_phase
                    if sanitycheck
                        load(strrep(sessions_info(i).Input_spikes,'population','by_block'));
                        ecg_bna_sanity_check(by_block,trials,Triggers,blockstart,cfg)
                    end
                    load(sessions_info(i).Input_spikes);
                    if cfg.spk.compute_spike_histograms
                        ecg_bna_PSTH(trials,population,Triggers,blockstart,cfg)
                    end
                    
%                     if cfg.spk.compute_spike_phase
%                         ecg_bna_compute_session_ECG_related_spikePhase(trials,population,Rpeaks,cfg)
%                     end
                end
                
                if cfg.spk.plot_spike_histograms
                    EPO=ecg_bna_get_plotoptions(cfg,'spk');
                    ecg_bna_plot_PSTH(sessions_info(i),EPO,cfg);
                    %ecg_bna_plot_session_spike_histogram(sessions_info(i),cfg);
                end
%                 if cfg.spk.plot_spike_phase
%                     ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),cfg)
%                 end
            end
            if cfg.process_LFP
                fprintf('Analysing for session %s\n', [sessions_info(i).Monkey '_' sessions_info(i).Date]);
                cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session');
                cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site');
                
                %% this is new
                sitesdir=fileparts(sessions_info(i).Input_LFP{:});
                [sitefiles]=dir(sessions_info(i).Input_LFP{:});
                sr=unique([trials.TDT_LFPx_SR]);
                ts_original=1/sr;
                
                %% load all sites
                tic
                allSitesData = ecg_bna_remove_rawLFP_outliers(sitesdir,sitefiles,cfg,ts_original);
                toc;
                ecg_bna_rawLFP_butterfly_plots(cfg,allSitesData,sr)
                %% exclude outliars (too many samples with too high/low voltage ? or other criteria)
                   %% threshold = 1V (fixed threshold), more than 1% of bins above that threshold
                %% exclude blocks with less than 3 sites
                %% compute common ground for remaining blocks
                %% identify valid blocks (because some BLOCKS will be excluded, not only sites)

                %% create corrected site LFP data for each valid site, and loop through that
                %% instead of looping through sitefiles and load each of them again
                
                for s = 1:length(allSitesData) %% loop only through valid sites
%                     load([sitesdir filesep sitefiles{s}.name], 'sites');
                    %% adress correct site
                    %% site = ...
                    sites = allSitesData(s).site;
                    site_LFP= ecg_bna_process_LFP(sites, cfg, ts_original);
                    n_LFP_samples_per_block=site_LFP.tfs.n_samples_per_block;
                    
                    site_triggers     = ecg_bna_resample_triggers2(Triggers,[blocks;blockstart(blocks)],n_LFP_samples_per_block,site_LFP.tfs.sr);                     
                    site_data=ecg_bna_compute_triggered_LFP_variables(site_LFP,site_triggers,trials,cfg );                    
                    
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
    
    %% average across sessions
    if cfg.process_population
        
        keys=ecg_bna_get_unit_list(cfg,0);
        cfg.site_IDS=keys.tuning_table(2:end,find_column_index(keys.tuning_table,'site_ID'));
        
        if cfg.process_LFP
            monkeys = unique({cfg.session_info.Monkey});
            cfg.monkey = [monkeys{:}];
            cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session');
            cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site');
            
            grand_avg = ecg_bna_compute_grand_avg(cfg,'w_units');
            grand_avg = ecg_bna_compute_grand_avg(cfg,'wo_units');
            grand_avg = ecg_bna_compute_grand_avg(cfg,'all');
        end
        
        if cfg.process_spikes
            SPK_PSTH=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','','per_unit','Output');
            %ecg_bna_avg_spike_histogram(SPK_PSTH,sessions_info, cfg);
            ecg_bna_avg_spike_histogram_clean(SPK_PSTH,cfg);
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

