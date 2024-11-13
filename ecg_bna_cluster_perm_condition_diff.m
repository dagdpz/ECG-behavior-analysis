function ecg_bna_cluster_perm_condition_diff(cfg,withunits,data_label)
if isfield(cfg.lfp,'IBI') && cfg.lfp.IBI==1
    if cfg.lfp.IBI_low == 1 || cfg.lfp.IBI_high == 0
        
        type = 'IBIlow';
        path_to_grand_avg   = fullfile(cfg.analyse_lfp_folder, 'grand_average_IBIlow');
        cd(path_to_grand_avg)
        load(fullfile(path_to_grand_avg,filesep,[cfg.monkey,...
            '_R peak_Triggered_target_wise_Grand_grand_avg_sessions_sites',withunits,'.mat']));
        
    elseif cfg.lfp.IBI_high == 1 || cfg.lfp.IBI_low == 0
        
        type = 'IBIhigh';
        path_to_grand_avg   = fullfile(cfg.analyse_lfp_folder, 'grand_average_IBIhigh');
        cd(path_to_grand_avg)
        load([cfg.monkey  ,'_R peak_Triggered_target_wise_Grand_grand_avg_sessions_sites',withunits,'.mat'])
    end
elseif isfield(cfg.lfp,'IBI') && cfg.lfp.IBI==0
    type = 'woIBIsplit';
    path_to_grand_avg   = fullfile(cfg.analyse_lfp_folder, 'grand_average_all');
    cd(path_to_grand_avg)
    load(fullfile(path_to_grand_avg,filesep,[cfg.monkey,...
        '_R peak_Triggered_target_wise_Grand_grand_avg_sessions_sites',withunits,'.mat']));
    
end
% Assuming the data is loaded in the following variables:
% data_task: data for the task condition (freq x time x trials)
% data_rest: data for the rest condition (freq x time x trials)
% tfr_time: A vector representing the time points (1 x time)
% tfr_freq: A vector representing the frequency points (1 x freq)

target = {grand_avg.target};
for tr = 1:length(target)
    nSites = grand_avg(tr).nSites;
    if (nSites==0)
        continue;
    else
        cond = {[grand_avg(tr).avg(1).cond_name,'-',grand_avg(tr).avg(2).cond_name]};
        
        if strcmp(data_label,'itpc') || strcmp(data_label,'pow')
            data_rest = grand_avg(tr).avg(1).(data_label); % grand adverage of target, Rest Condition , across grand_avg(2).nSites;
            data_task =  grand_avg(tr).avg(2).(data_label); % grand adverage of target, Task Condition ,across grand_avg(2).nSites;
            tfr_freq = cfg.lfp.foi; %logspace(log10(2),log10(120),60);
            
            % Step 1: Create FieldTrip-compliant structures for task and rest conditions
            freq_task = [];
            freq_task.label     = {data_label};  % Dummy label (no real channels here)
            freq_task.dimord    = 'rpt_freq_time';  % Trials, frequencies, and time
            freq_task.freq      = tfr_freq;  % Frequency vector
            freq_task.time      = tfr_time;  % Time vector
            freq_task.powspctrm = permute(data_task, [3, 1, 2]);  % [trials x frequencies x time]
            
            freq_rest = [];
            freq_rest.label     = {data_label};  % Dummy label (no real channels here)
            freq_rest.dimord    = 'rpt_freq_time';  % Trials, frequencies, and time
            freq_rest.freq      = tfr_freq;  % Frequency vector
            freq_rest.time      = tfr_time;  % Time vector
            freq_rest.powspctrm = permute(data_rest, [3, 1, 2]);  % [trials x frequencies x time]
            
            % Step 2: Prepare configuration for cluster-based permutation testing
            cfg_perm = [];
            cfg_perm.method           = 'montecarlo';             % Monte Carlo permutation method
            %         cfg_perm.statistic        = 'ft_statfun_indepsamplesT'; % Independent samples t-test
            cfg_perm.statistic        = 'ft_statfun_depsamplesT'; % Independent samples t-test
            cfg_perm.correctm         = 'cluster';                % Cluster-based correction
            cfg_perm.clusteralpha     = 0.01;                     % Cluster-forming threshold
            cfg_perm.clusterstatistic = 'maxsum';                 % Max-sum clustering or 'maxsize', 'wcm' (weighted cluster mass), 'wcs' (weighted cluster size)
            cfg_perm.minnbchan        = 0;                        % No neighboring channels (no spatial information)
            cfg_perm.tail             = 0;                        % Two-tailed test
            cfg_perm.clustertail      = 0;                        % Two-tailed for clustering
            cfg_perm.alpha            = 0.01;                     % Alpha level for the permutation test
            cfg.correcttail           = 'alpha';                  % FDR correction for two-tailed T-test
            cfg_perm.numrandomization = 5000;                     % Number of randomizations
            cfg_perm.neighbours            = [];
            
            % Define a dummy connectivity matrix for the frequencies
            num_freqs = length(tfr_freq);  % Number of frequency bins
            cfg_perm.connectivity = ones(num_freqs, num_freqs);  % Fully connected between frequencies
            
            
            % Step 3: Create the design matrix
            
            % %         n_trials_task = size(data_task, 3);  % Number of trials in task condition
            % %         n_trials_rest = size(data_rest, 3);  % Number of trials in rest condition
            % %
            % %         design = zeros(1, n_trials_task + n_trials_rest);
            % %         design(1, 1:n_trials_task) = 1;  % Label task trials as 1
            % %         design(1, (n_trials_task+1):(n_trials_task+n_trials_rest)) = 2;  % Label rest trials as 2
            % %
            % %         cfg_perm.design = design;        % Design: task vs. rest
            % %         cfg_perm.ivar   = 1;             % The independent variable (1: task, 2: rest)
            
            % Define the design matrix for dependent samples (Task vs. Rest)
            nSites = size(freq_task.powspctrm,1);
            
            design = zeros(2, 2 * nSites);
            design(1, 1:nSites)         = 1;          % Task condition
            design(1, nSites+1:end)     = 2;          % Rest condition
            design(2, :)                   = [1:nSites, 1:nSites];  % site indices
            
            cfg_perm.design = design;
            cfg_perm.uvar   = 2;                              % Unit variable (Site index)
            cfg_perm.ivar   = 1;                              % Independent variable (condition: Task vs. Rest)
            
            % Step 4: Perform the cluster-based permutation test
            [stat] = ft_freqstatistics(cfg_perm, freq_task, freq_rest);
            
        elseif strcmp(data_label,'lfp')
            data_rest = grand_avg(tr).avg(1).(data_label); % grand adverage of target, Rest Condition , across grand_avg(2).nSites;
            data_task =  grand_avg(tr).avg(2).(data_label); % grand adverage of target, Task Condition ,across grand_avg(2).nSites;
            
            
            % Step 1: Create FieldTrip-compliant structures for task and rest conditions
            LFP_task = [];
            LFP_task.label     = {data_label};  % Dummy label (no real channels here)
            LFP_task.dimord    = 'rpt_freq_time';  % Trials, frequencies, and time
            LFP_task.time      = time;  % Time vector
            LFP_task.powspctrm = permute(data_task, [3, 1, 2]);  % [trials x frequencies x time]
            
            freq_rest = [];
            freq_rest.label     = {data_label};  % Dummy label (no real channels here)
            freq_rest.dimord    = 'rpt_freq_time';  % Trials, frequencies, and time
            freq_rest.freq      = tfr_freq;  % Frequency vector
            freq_rest.time      = tfr_time;  % Time vector
            freq_rest.powspctrm = permute(data_rest, [3, 1, 2]);  % [trials x frequencies x time]
            
            % Step 2: Prepare configuration for cluster-based permutation testing
            cfg_perm = [];
            cfg_perm.method           = 'montecarlo';             % Monte Carlo permutation method
            cfg_perm.statistic        = 'ft_statfun_depsamplesT'; % Independent samples t-test
            cfg_perm.correctm         = 'cluster';                % Cluster-based correction
            cfg_perm.clusteralpha     = 0.01;                     % Cluster-forming threshold
            cfg_perm.clusterstatistic = 'maxsum';                 % Max-sum clustering or 'maxsize', 'wcm' (weighted cluster mass), 'wcs' (weighted cluster size)
            cfg_perm.minnbchan        = 0;                        % No neighboring channels (no spatial information)
            cfg_perm.tail             = 0;                        % Two-tailed test
            cfg_perm.clustertail      = 0;                        % Two-tailed for clustering
            cfg_perm.alpha            = 0.01;                     % Alpha level for the permutation test
            cfg.correcttail           = 'alpha';                  % FDR correction for two-tailed T-test
            cfg_perm.numrandomization = 5000;                     % Number of randomizations
            cfg_perm.neighbours            = [];
            
            % Define a dummy connectivity matrix for the frequencies
            num_freqs = length(tfr_freq);  % Number of frequency bins
            cfg_perm.connectivity = ones(num_freqs, num_freqs);  % Fully connected between frequencies
            
            
            % Step 3:  Define the design matrix for dependent samples (Task vs. Rest)
            nSites = size(LFP_task.powspctrm,1);
            
            design = zeros(2, 2 * nSites);
            design(1, 1:nSites)         = 1;          % Task condition
            design(1, nSites+1:end)     = 2;          % Rest condition
            design(2, :)                   = [1:nSites, 1:nSites];  % site indices
            
            cfg_perm.design = design;
            cfg_perm.uvar   = 2;                              % Unit variable (Site index)
            cfg_perm.ivar   = 1;                              % Independent variable (condition: Task vs. Rest)
            
            % Step 4: Perform the cluster-based permutation test
            [stat] = ft_freqstatistics(cfg_perm, LFP_task, freq_rest);
            
        end
        
        %     temp = cat(3, stat.stat, stat.stat );
        %     temp = permute(temp, [3, 1, 2]);
        %     temp = temp(1,:,:);
        %     stat_temp = stat;
        %     stat_temp.stat = temp;
        %     temp = [];
        %     temp = cat(3, stat.mask, stat.mask );
        %     temp = permute(temp, [3, 1, 2]);
        %     temp = temp(1,:,:);
        %     stat_temp.mask = temp;
        %     clear temp
        %
        %     % Step 5: Visualize the results (optional, using a time-frequency plot)
        %     cfg_perm_plot = [];
        %     cfg_perm_plot.channel = 1;
        %     cfg_perm_plot.parameter = 'stat';      % Use the t-statistics from the test
        %     cfg_perm_plot.maskparameter = 'mask';  % Show significant clusters
        %     cfg_perm_plot.maskstyle = 'outline';   % Outline significant clusters
        %     cfg_perm_plot.masknans = 'yes';         % Disable masking of NaNs if they exist
        %
        %     figure;
        %     ft_singleplotTFR(cfg_perm_plot, stat_temp); % Initial plot
        %     title('Significant Clusters between Task and Rest Conditions');
        %     xlabel('Time (s)');
        %
        h(tr) = figure;
        imagesc(stat.time, 1:numel(stat.freq),stat.stat);
        set(gca,'YDir','normal');
        % set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
        fbandstart = unique(cfg.lfp.frequency_bands(:))';
        set(gca,'Xlim',[-.25 .25]);
        line([0 0], ylim, 'color', 'k');
        % horizontal lines to separate frequency bands
        fbandstart_idx = zeros(size(fbandstart));
        for f = fbandstart
            f_idx = find(abs(stat.freq - f) == min(abs(stat.freq - f)), 1, 'first');
            yline(f_idx, 'color', 'k', 'linestyle', '--');
            fbandstart_idx(fbandstart == f) = f_idx;
        end
        set(gca,'TickDir','out')
        set(gca, 'ytick', fbandstart_idx);
        set(gca, 'yticklabel', fbandstart);
        set(gca, 'ylim', [0.5,numel(stat.freq) + 0.5]);
        hold on
        contour(stat.time, 1:numel(stat.freq), stat.mask, [0.5, 0.5], 'LineColor', 'k', 'LineWidth', 2); % Overlay significant clusters
        % clabel(C, h);
        colormap(jet);  % Change the colormap if desired
        colorbar;  % Show colorbar
        mtit([ cfg.monkey,'-',target{tr},'-avg of ',num2str(nSites),' sites ' ,...
            withunits, ' - ',cond{1},'-',data_label,'- Significant difference - numPerm =',...
            num2str(cfg_perm.numrandomization),'-',type],'xoff', 0.05,'yoff', 0.02,...
            'color', [0 0 0], 'fontsize', 8,'Interpreter', 'none');
        
        
        results_file{tr} = fullfile([path_to_grand_avg,filesep,cfg.monkey,'-',target{tr},' - ',data_label,'- Significant difference of ',cond{1},'_sites ', withunits]);
        export_fig(h(tr),[results_file{tr},'.pdf']);
        
    end
    close all,
end
end