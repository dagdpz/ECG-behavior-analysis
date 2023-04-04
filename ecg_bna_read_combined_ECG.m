function session_ecg = ecg_bna_read_combined_ECG( session_info, plottrials  )
% ecg_bna_read_combined_ECG - function to read the ECG raw data and
% timestamps and trial information for all trials in a session from
% combined files (For example, files stored inside
% Y:\Data\[Monkey]\[Date])
%
% USAGE:
%	session_ecg = ecg_bna_read_combined_ECG( session_info, plottrials  )
%
% INPUTS:
%   session_info        - struct containing information about all sessions
%   to be processed
%   plottrials          - whether to plot ECG data for individual trials
% OUTPUTS:
%   session_info        - a copy of input struct session_info
%
% REQUIRES: ecg_bna_get_block_Rpeak_times
%
% See also lfp_tfa_process_LFP, ecg_bna_read_preproc_ECG, ecg_bna_read_combined_ECG


close all;

if nargin < 2
    plottrials = 0;
end


% struct to save data for a site
session_ecg = struct();

if ~exist(session_info.Input_ECG_combined, 'dir')
    fprintf('No file with ECG data found in the specified directory \n%s\n', ...
        session_info.Input_ECG_combined);
    return;
end
block_files = dir(fullfile(session_info.Input_ECG_combined, '*.mat'));

% prepare results folder
results_fldr = fullfile(session_info.proc_ecg_fldr);
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

if isfield(session_info, 'Input_ECG')
    if ~exist(session_info.Input_ECG, 'file')
        fprintf('No file found \n%s\n', ...
            session_info.Input_ECG);
        return;
    end
    load(session_info.Input_ECG);
    if exist('out', 'var')
        block_Rpeaks = out;
        clear out;
    end
end

comp_trial = 0; % iterator for completed trials

% save data inside struct
% first loop through each block
for b = 1:length(block_files)
    
    % get info about site
    nrblock_combinedFile = str2num(block_files(b).name(end-5:end-4));
    run = str2num(block_files(b).name(end-8:end-7));
    
    % check if ECG peaks exists for this block
    if length(block_Rpeaks) >= b
        block_Rpeak = block_Rpeaks(b);
        if isempty(block_Rpeak) || isempty(block_Rpeak.Rpeak_t)
            continue;
        end
    else
        continue;
    end
    
    fprintf('=============================================================\n');
    fprintf('Reading ECG for block, %g\n', nrblock_combinedFile);
    
    %             % block ECG timestamps
    %             block_ecg_timestamps = [];
    %             trial_ecg_timestamps = [];
    
    session_ecg.session = session_info.session;
    
    load(fullfile(session_info.Input_ECG_combined, block_files(b).name));
    
    % add first trial INI to block timestamps
    ts = 1/(trial(1).TDT_ECG1_samplingrate);
    first_trial_INI_timestamps = linspace(0, ts*(length(First_trial_INI.ECG1)-1), length(First_trial_INI.ECG1));
    trial_ecg_timestamps = first_trial_INI_timestamps;
    
    %% get information common to all sites for a session
    
    % now loop through each trial for this site
    for t = 1:length(trial)
        completed = trial(t).completed;
        if true
            type = trial(t).type;
            effector = trial(t).effector;
            completed = trial(t).completed;
            choice_trial = trial(t).choice;
            perturbation = nan;
            
            % decide about perturbation ... this looks horrible!
            if isfield(session_info, 'Preinj_blocks') && ~isempty(session_info.Preinj_blocks) && ismember(nrblock_combinedFile, session_info.Preinj_blocks)
                perturbation = 0;
            elseif exist('ses', 'var') && isfield(ses, 'first_inj_block') && ~isempty(ses.first_inj_block) && nrblock_combinedFile < ses.first_inj_block
                perturbation = 0;
            end
            if isfield(session_info, 'Postinj_blocks') && ~isempty(session_info.Postinj_blocks) && ismember(nrblock_combinedFile, session_info.Postinj_blocks)
                perturbation = 1;
            elseif exist('ses', 'var') && isfield(ses, 'first_inj_block') && ~isempty(ses.first_inj_block) && nrblock_combinedFile >= ses.first_inj_block
                perturbation = 1;
            end
            if isnan(perturbation) && isfield(trial, 'perturbation') && isempty(trial(t).perturbation)
                perturbation = trial(t).perturbation;
            end
            
            
            start_time = 0; % trial start time
            fs = trial(t).TDT_ECG1_samplingrate; % sample rate
            ts = (1/fs); % sample time
            ECG = trial(t).TDT_ECG1; % ecg data
            nsamples = numel(ECG);
            end_time = start_time + (ts*(nsamples-1));
            timestamps = linspace(start_time, end_time, nsamples);
            trial_ecg_timestamps = timestamps + ts + trial_ecg_timestamps(end);
            % save retrieved data into struct
            comp_trial = comp_trial + 1;
            session_ecg.trials(comp_trial).completed = completed;
            session_ecg.trials(comp_trial).type = type;
            session_ecg.trials(comp_trial).effector = effector;
            session_ecg.trials(comp_trial).run = run;
            session_ecg.trials(comp_trial).block = nrblock_combinedFile;
            session_ecg.trials(comp_trial).dataset = [];
            session_ecg.trials(comp_trial).choice_trial = choice_trial;
            session_ecg.trials(comp_trial).reach_hand = 0;
            session_ecg.trials(comp_trial).reach_space = 0;
            session_ecg.trials(comp_trial).hndspc_lbl  = [];
            session_ecg.trials(comp_trial).time = timestamps;
            session_ecg.trials(comp_trial).ecg_data = ECG;
            session_ecg.trials(comp_trial).fsample  = fs;
            session_ecg.trials(comp_trial).tsample = ts;
            session_ecg.trials(comp_trial).tstart = start_time;
            session_ecg.trials(comp_trial).trialperiod = [trial_ecg_timestamps(1) trial_ecg_timestamps(end)];
            session_ecg.trials(comp_trial).perturbation  = perturbation;
            % flag to mark noisy trials, default False, filled in by lfp_tfa_reject_noisy_lfp.m
            session_ecg.trials(comp_trial).noisy = ~completed;
            
            % get state onset times and onset samples - test and delete
            session_ecg.trials(comp_trial).states = struct();
            st_idx = 0;
            for st = 1:length(trial(t).TDT_states)
                % get state ID
                state_id = trial(t).TDT_states(st);
                % get state onset time
                state_onset = trial(t).TDT_state_onsets(trial(t).TDT_states == state_id);
                % for combined file, first sample
                % correspond to timestamp = 0
                if state_onset < 0
                    continue;
                end
                % get sample number of state onset time
                state_onset_sample = find(abs(timestamps - state_onset(1)) == min(abs(timestamps - state_onset(1))), 1);
                st_idx = st_idx + 1;
                % save into struct
                session_ecg.trials(comp_trial).states(st_idx).id = state_id;
                session_ecg.trials(comp_trial).states(st_idx).onset_t  = state_onset(1);
                session_ecg.trials(comp_trial).states(st_idx).onset_s  = state_onset_sample;
            end
            %end
        end
    end
    session_ecg = ecg_bna_get_block_Rpeak_times( session_ecg, block_Rpeak, nrblock_combinedFile, plottrials, results_fldr );
    
    %%% Noise rejection - should this be included within processing check this? %%%
    %state_filt_lfp(i) = lfp_tfa_reject_noisy_lfp( state_lfp(i), lfp_tfa_cfg.noise );
end

% save allsites_lfp
results_mat = fullfile(results_fldr, ['session_ecg_' session_info.session '.mat']);
save(results_mat, 'session_ecg');
end

