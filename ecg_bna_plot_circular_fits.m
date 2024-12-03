function ecg_bna_plot_circular_fits(cfg, data_folder, destination_folder)

basepath_to_save = [cfg.SPK_root_results_fldr filesep destination_folder];
if ~isfolder(basepath_to_save)
    mkdir(basepath_to_save)
end

dataset_name = [basepath_to_save filesep 'data.mat'];

%% temporary - do plotting of fit parameters
phase_by_area_cond1_vMpos = [];
phase_by_area_cond2_vMpos = [];
area_id_vMpos             = [];

phase_by_area_cond1_vMneg = [];
phase_by_area_cond2_vMneg = [];
area_id_vMneg             = [];

C = colororder;
C(1,:) = [0 0 0];

var_list_ephys = {'NrEvents', 'sig_n_bins', ...
    ...
    'sig_FR_diff', 'sig_time', 'sig_n_bins', 'sig_sign', ...
    'FR_ModIndex_SubtrSDP', 'FR_ModIndex_PcS', ...
    ...
    'lowIBI_sig_FR_diff', 'lowIBI_sig_time', 'lowIBI_sig_n_bins', 'lowIBI_sig_sign', ...
    'lowIBI_FR_ModIndex_SubtrSDP', 'lowIBI_FR_ModIndex_PcS', ...
    ...
    'highIBI_sig_FR_diff', 'highIBI_sig_time', 'highIBI_sig_n_bins', 'highIBI_sig_sign', ...
    'highIBI_FR_ModIndex_SubtrSDP', 'highIBI_FR_ModIndex_PcS'};

var_list_cardioballistic = ...
    {'lowIBI_medianHR_bpm', 'highIBI_medianHR_bpm', ...
    ...
    'sig_FR_diff', 'sig_time', 'sig_n_bins', 'sig_sign', ...
    'FR_ModIndex_SubtrSDP', 'FR_ModIndex_PcS', ...
    ...
    'lowIBI_sig_FR_diff', 'lowIBI_sig_time', 'lowIBI_sig_n_bins', 'lowIBI_sig_sign', ...
    'lowIBI_FR_ModIndex_SubtrSDP', 'lowIBI_FR_ModIndex_PcS', ...
    ...
    'highIBI_sig_FR_diff', 'highIBI_sig_time', 'highIBI_sig_n_bins', 'highIBI_sig_sign', ...
    'highIBI_FR_ModIndex_SubtrSDP', 'highIBI_FR_ModIndex_PcS'}; % ,

load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG\unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_any.mat'], 'unit_ids', 'targets')

%     unqTargets = unique(targets);
unqTargets = {'VPL', 'dPul', 'MD'};

for groupNum = 1:length(cfg.spk.compare_conditions)
    cond1_num = cfg.spk.compare_conditions{groupNum}(1);
    cond2_num = cfg.spk.compare_conditions{groupNum}(2);
    
    if ~exist(dataset_name,'file')
        for targNum = 1:length(unqTargets)
            
            currTargIds = cellfun(@(x) strcmp(x, unqTargets{targNum}), targets);
            curr_unit_ids = unit_ids(currTargIds);
            dt.(unqTargets{targNum}) = ecg_bna_load_variables(cfg, curr_unit_ids, 'cardioballistic', 'data', var_list_cardioballistic);
            
            ephys_dt.(unqTargets{targNum}) = ecg_bna_load_variables(cfg, curr_unit_ids, data_folder, 'Output', var_list_ephys);
            
        end
        save(dataset_name,'dt','ephys_dt')
    else
        load(dataset_name,'dt','ephys_dt')
    end
end

for groupNum = 1:length(cfg.spk.compare_conditions)
    cond1_num = cfg.spk.compare_conditions{groupNum}(1);
    cond2_num = cfg.spk.compare_conditions{groupNum}(2);
    
    %% Peak Response Phases and Times BETWEEN LOW AND HIGH IBI
    % preallocate variabes
    meanPhase_lowIBI  = zeros(length(unqTargets),length(cfg.condition));
    meanPhase_highIBI = zeros(length(unqTargets),length(cfg.condition));
    circ_cc           = zeros(length(unqTargets),length(cfg.condition));
    circ_pp           = zeros(length(unqTargets),length(cfg.condition));
    obs_diff          = zeros(length(unqTargets),length(cfg.condition));
    p_value           = zeros(length(unqTargets),length(cfg.condition));
    p_wwtest          = zeros(length(unqTargets),length(cfg.condition));
    
    meanTime_lowIBI   = zeros(length(unqTargets),length(cfg.condition));
    meanTime_highIBI  = zeros(length(unqTargets),length(cfg.condition));
    cc                = zeros(length(unqTargets),length(cfg.condition));
    pp                = zeros(length(unqTargets),length(cfg.condition));
    p_wilcox          = zeros(length(unqTargets),length(cfg.condition));
    h_wilcox          = zeros(length(unqTargets),length(cfg.condition));
    
    for c = 1:length(cfg.condition)
        L = cfg.condition(c).name;
        
        for targNum = 1:length(unqTargets)
            T = unqTargets{targNum};
            
            % find units with significant responses in time domain in
            % either condition
            sig = ...
                ephys_dt.(T).(cfg.condition(cond1_num).name).sig_n_bins > cfg.time.n_sig_bins | ...
                ephys_dt.(T).(cfg.condition(cond2_num).name).sig_n_bins > cfg.time.n_sig_bins;
            % find units with consistent sign of phase responses in both
            % conditions
%             eq_sign_phase = eq(dt.(T).(cfg.condition(cond1_num).name).sig_sign, dt.(T).(cfg.condition(cond2_num).name).sig_sign);
            
            not_nan_sig = ~isnan(dt.(T).(L).lowIBI_sig_time) & ~isnan(dt.(T).(L).highIBI_sig_time) & sig; %  & eq_sign_phase
            lowIBI_peak_phases = cfg.phase.phase_bin_centers(dt.(T).(L).lowIBI_sig_time(not_nan_sig));
            highIBI_peak_phases = cfg.phase.phase_bin_centers(dt.(T).(L).highIBI_sig_time(not_nan_sig));
            
            if numel(lowIBI_peak_phases) < 2 || numel(highIBI_peak_phases) < 2
                continue
            end
            
            meanPhase_lowIBI(targNum,c)  = mod(circ_mean(lowIBI_peak_phases),2*pi);
            meanPhase_highIBI(targNum,c) = mod(circ_mean(highIBI_peak_phases),2*pi);
            [circ_cc(targNum,c), circ_pp(targNum,c)] = circ_corrcc(lowIBI_peak_phases, highIBI_peak_phases);
            [obs_diff(targNum,c), p_value(targNum,c)] = DAG_circ_perm_means(lowIBI_peak_phases, highIBI_peak_phases, 10000);
            % compute Watson-Williams test to compare circular means,
            % collect its warnings to know whether its assumptions were
            % violated by the dataset
            % Turn off all warnings and save the current warning state
            warningState = warning('off', 'all');
            warning('on', 'all');
            lastwarn('');
            [p_wwtest(targNum,c), tbl] = circ_wwtest(lowIBI_peak_phases, highIBI_peak_phases);
            lastWarnMsg = lastwarn();
            tbl{5,1} = lastWarnMsg;
            
            % Re-enable the original warning state
            warning(warningState);
            
            figure,
            set(gcf,'Position',[681   399   628   580])
            scatterhist(lowIBI_peak_phases, highIBI_peak_phases, ...
                'NBins', 10, ...
                'Location', 'SouthWest', ...
                'Direction', 'out', ...
                'Style', 'stairs', ...
                'Kernel', 'on', ...
                'Color', cfg.condition(c).color);
            hold on
            % plot means
            scatter(meanPhase_lowIBI(targNum,c), 0, 'Marker', 'v', 'MarkerFaceColor', cfg.condition(c).color, 'MarkerEdgeColor', 'none')
            scatter(0, meanPhase_highIBI(targNum,c), 'Marker', '<', 'MarkerFaceColor', cfg.condition(c).color, 'MarkerEdgeColor', 'none')
            hold off
            title({[T ': N = ' num2str(sum(not_nan_sig))], ...
                ['circ cc = ' num2str(circ_cc(targNum,c)) '; p = ' num2str(circ_pp(targNum,c))], ...
                ['Watson-Williams test: p = ' num2str(p_wwtest(targNum,c))]})
            xlim([0 2*pi])
            ylim([0 2*pi])
            xlabel(['lowIBI ' L ': Peak Phases, radians [0 2\pi]'])
            ylabel(['highIBI ' L ': Peak Phases, radians [0 2\pi]'])
            
            
            % save Watson-Williams test results
            filename = ['Watson-Williams_Test_PeakPhase_lowIBI_vs_highIBI_' T '_' L];
            writetable(cell2table(tbl), [basepath_to_save filesep filename '.xls'])
            
            % save scatter histogram
            filename = ['PeakPhase_lowIBI_vs_highIBI_' T '_' L];
            export_fig(gcf, [basepath_to_save, filesep ,filename], '-pdf');
            close(gcf);
            
            % the same for the time domain
            lowIBI_peak_times = ephys_dt.(T).(L).lowIBI_sig_time(not_nan_sig);
            highIBI_peak_times = ephys_dt.(T).(L).highIBI_sig_time(not_nan_sig);
            
            meanTime_lowIBI(targNum,c)  = mean(lowIBI_peak_times);
            meanTime_highIBI(targNum,c) = mean(highIBI_peak_times);
            
            [cc_tmp, pp_tmp] = corrcoef(lowIBI_peak_times, highIBI_peak_times);
            cc(targNum,c) = cc_tmp(2,1);
            pp(targNum,c) = pp_tmp(2,1);
            [p_wilcox(targNum,c), h_wilcox(targNum,c)] = signrank(lowIBI_peak_times, highIBI_peak_times);
            
            figure
            set(gcf,'Position',[681 381 560 598])
            s = scatterhist(lowIBI_peak_times, highIBI_peak_times, ...
                'NBins', 10, ...
                'Location', 'SouthWest', ...
                'Direction', 'out', ...
                'Style', 'stairs', ...
                'Kernel', 'on', ...
                'Color', cfg.condition(c).color);
            hold on
            % plot means
            scatter(meanTime_lowIBI(targNum,c), 0, 'Marker', 'v', 'MarkerFaceColor', cfg.condition(c).color, 'MarkerEdgeColor', 'none')
            scatter(0, meanTime_highIBI(targNum,c), 'Marker', '<', 'MarkerFaceColor', cfg.condition(c).color, 'MarkerEdgeColor', 'none')
            hold off
            title({[T ': N = ' num2str(sum(not_nan_sig))], ...
                ['cc = ' num2str(cc(targNum,c)) '; p = ' num2str(pp(targNum,c))], ...
                ['Paired Wilcoxon test: p = ' num2str(p_wilcox(targNum,c))]})
            xlim([0 400])
            ylim([0 400])
            xlabel(['lowIBI ' L ': Peak Times, ms after R-peak'])
            ylabel(['highIBI ' L ': Peak Times, ms after R-peak'])
            
            % save scatter histogram
            filename = ['PeakTimes_lowIBI_vs_highIBI_' T '_' L];
            export_fig(gcf, [basepath_to_save, filesep ,filename], '-pdf');
            close(gcf);
            
        end
        
    end
    
    % phase domain: construct resulting table
    h              = fdr_bky(circ_pp,0.05,'yes'); % multiple comparison correction
    ww_h_corrected = fdr_bky(p_wwtest,0.05,'yes');
    
    t = table(unqTargets', ...
        circ_cc(:,1), circ_cc(:,2), circ_pp(:,1), circ_pp(:,2), h(:,1), h(:,2), ...
        meanPhase_lowIBI(:,1), meanPhase_lowIBI(:,2), meanPhase_highIBI(:,1), meanPhase_highIBI(:,2), ww_h_corrected(:,1), ww_h_corrected(:,2), ...
        'VariableNames',{'Target', ...
        'Rest: Circ. Corr.',        'Task: Circ. Corr.', ...
        'Rest: Circ. p',            'Task: Circ. p', ...
        'Rest: After FDR',          'Task: After FDR', ...
        'Rest: Mean Phase lowIBI',  'Task: Mean Phase lowIBI', ...
        'Rest: Mean Phase highIBI', 'Task: Mean Phase highIBI', ...
        'Rest: WW After FDR',       'Task: WW After FDR'});
    filename = 'Table_PeakPhase_lowIBI_vs_highIBI';
    writetable(t, [basepath_to_save filesep filename '.xls'])
    clear t
    
    % permutation test results
    t = table(unqTargets', ...
        obs_diff(:,1), obs_diff(:,2), p_value(:,1), p_value(:,2), ...
        'VariableNames',{'Target', ...
        'Rest: Obs. Diff.', 'Task: Obs. Diff.', 'Rest: permuted p', 'Task: permuted p'});
    filename = 'Table_PeakPhasePermuted_lowIBI_vs_highIBI';
    writetable(t, [basepath_to_save filesep filename '.xls'])
    clear t
    
    % time domain: construct resulting table
    h              = fdr_bky(pp,0.05,'yes'); % multiple comparison correction
    wilcox_h_corrected = fdr_bky(p_wwtest,0.05,'yes');
    
    t = table(unqTargets', ...
        cc(:,1), cc(:,2), pp(:,1), pp(:,2), h(:,1), h(:,2), ...
        meanTime_lowIBI(:,1), meanTime_lowIBI(:,2), meanTime_highIBI(:,1), meanTime_highIBI(:,2), wilcox_h_corrected(:,1), wilcox_h_corrected(:,2), ...
        'VariableNames',{'Target', ...
        'Rest: Corr.',             'Task: Corr.', ...
        'Rest: Corr. p',           'Task: Corr. p', ...
        'Rest: After FDR',         'Task: After FDR', ...
        'Rest: Mean Time lowIBI',  'Task: Mean Time lowIBI', ...
        'Rest: Mean Time highIBI', 'Task: Mean Time highIBI', ...
        'Rest: Wilcox After FDR',  'Task: Wilcox After FDR'});
    filename = 'Table_PeakTime_lowIBI_vs_highIBI';
    writetable(t, [basepath_to_save filesep filename '.xls'])
    clear t
    
    %% RESPONSE PHASES BETWEEN REST AND TASK
%     figure,
%     set(gcf, 'Position', [331 423 1012 324])
    for targNum = 1:length(unqTargets)
        T = unqTargets{targNum};
        
        sig = ...
            ephys_dt.(T).(cfg.condition(cond1_num).name).sig_n_bins > cfg.time.n_sig_bins | ...
            ephys_dt.(T).(cfg.condition(cond2_num).name).sig_n_bins > cfg.time.n_sig_bins;
        not_nan_sig = ~isnan(dt.(T).Rest.sig_time) & ~isnan(dt.(T).Task.sig_time) & sig;
        cond1_peak_phases = cfg.phase.phase_bin_centers(dt.(T).Rest.sig_time(not_nan_sig));
        cond2_peak_phases = cfg.phase.phase_bin_centers(dt.(T).Task.sig_time(not_nan_sig));
        
%         subplot(1,3,targNum)
        figure
        set(gcf,'Position',[681 421 514 558])
        [circ_cc(targNum), circ_pp(targNum)] = circ_corrcc(cond1_peak_phases, cond2_peak_phases);
        [p(targNum), tbl] = circ_wwtest(cond1_peak_phases, cond2_peak_phases);
        
        scatterhistogram(cond1_peak_phases, cond2_peak_phases, ...
            'NumBins', 10, ...
            'HistogramDisplayStyle', 'smooth', ...
            'ScatterPlotLocation', 'NorthEast', ...
            'Title', {[T ': N = ' num2str(sum(not_nan_sig))], ...
            ['circ cc = ' num2str(circ_cc(targNum)) '; p = ' num2str(circ_pp(targNum))], ...
            ['Watson-Williams test: p = ' num2str(p(targNum))]})
        xlim([0 2*pi])
        ylim([0 2*pi])
        xlabel('Rest: Peak Phases, radians [0 2\pi]')
        ylabel('Task: Peak Phases, radians [0 2\pi]')
        
        % save Watson-Williams test results
        filename = ['Watson-Williams_Test_PeakPhase_' T];
        writetable(cell2table(tbl), [basepath_to_save filesep filename '.xls'])
        clear tbl 
        
        % save scatter histogram
        filename = ['PeakPhase_' T];
        export_fig(gcf, [basepath_to_save, filesep ,filename], '-pdf');
        close(gcf);
        
    end
    
    % construct results 
    h              = fdr_bky(circ_pp,0.05,'yes');
    ww_h_corrected = fdr_bky(p,0.05,'yes');
    
    t = table(unqTargets', circ_cc', circ_pp', h', ww_h_corrected', ...
        'VariableNames',{'Target', 'Circ. Corr.', 'Circ. p', 'After FDR', 'WW After FDR'});
    filename = 'Table_PeakPhase_Rest_vs_Task';
    writetable(t, [basepath_to_save filesep filename '.xls'])
    
    
    %% DIFFERENCE BETWEEN MEDIAN HR BETWEEN LOWIBI AND HIGHIBI
    figure
    set(gcf, 'Position', [557 545 350 451])
    
    for targNum = 1:length(unqTargets)
        T = unqTargets{targNum};
        
        for c = 1:length(cfg.condition)
            L = cfg.condition(c).name;
            
            s(c) = subplot(2,1,c);
            hist(1000*60 ./ dt.(T).(L).highIBI_medianHR_bpm - 1000*60 ./ dt.(T).(L).lowIBI_medianHR_bpm)
            xlabel('median(highIBI) - median(lowIBI), ms')
            ylabel('Unit Counts')
            title(L)
            
        end
        
        linkaxes(s, 'xy')
    
        sgtitle([cfg.version ' ' unqTargets{targNum}], 'interpreter', 'none')
        
        filename = [T '_MedianRRdurations'];
        export_fig(gcf, [basepath_to_save, filesep ,filename], '-pdf');
        close(gcf);
    
    end
    
    %% BUILD A MIXED EFFECTS MODEL TO FIGURE OUT R-SQUARED DIFFERENCES BETWEEN TIME AND PHASE DOMAINS
    
%     % preallocate data variables for mixed-effects model
%     data.Rsquared_log10 = [];
%     data.Domain         = {};
%     data.Condition      = {};
%     data.BrainNucleus   = {};
%     data.UnitID         = {};
%     
%     for targNum = 1:length(unqTargets)
%         T = unqTargets{targNum};
%         
%         % figure out significantly modulated units by spike analysis
%         sig = ...
%             ephys_dt.(T).(cfg.condition(cond1_num).name).sig_n_bins > cfg.time.n_sig_bins | ...
%             ephys_dt.(T).(cfg.condition(cond2_num).name).sig_n_bins > cfg.time.n_sig_bins;
%     
%         data.Rsquared_log10     = [data.Rsquared_log10; ...
%             dt.(T).(cfg.condition(cond1_num).name).cosine.rsquared(sig); ephys_dt.(T).(cfg.condition(cond1_num).name).cosine.rsquared(sig); ...
%             dt.(T).(cfg.condition(cond2_num).name).cosine.rsquared(sig); ephys_dt.(T).(cfg.condition(cond2_num).name).cosine.rsquared(sig)];
%     %     ephys_dt.(T).Rest.cosine.rsquared
%         data.Domain       = [data.Domain; ...
%             repmat({'Time'}, sum(sig), 1); repmat({'Phase'}, sum(sig), 1); ...
%             repmat({'Time'}, sum(sig), 1); repmat({'Phase'}, sum(sig), 1)];
%         
%         data.Condition    = [data.Condition; ...
%             repmat({'Rest'}, 2*sum(sig), 1); ...
%             repmat({'Task'}, 2*sum(sig), 1)];
%         
%         data.BrainNucleus = [data.BrainNucleus; repmat(unqTargets(targNum),4*length(dt.(T).unitId(sig)),1)];
%         data.UnitID       = [data.UnitID; ...
%             dt.(T).unitId(sig); ephys_dt.(T).unit_ID(sig); ...
%             dt.(T).unitId(sig); ephys_dt.(T).unit_ID(sig)];
%     
%         if targNum == length(unqTargets)
%             
%             data = struct2table(data);
%         
%             % Convert relevant columns to categorical
%             data.Domain       = categorical(data.Domain);
%             data.Condition    = categorical(data.Condition);
%             data.BrainNucleus = categorical(data.BrainNucleus);
%             data.UnitID       = categorical(data.UnitID);
%             
%             % Fit a mixed-effects model
%             lme = fitlme(data, 'Rsquared_log10 ~ 1 + BrainNucleus * Condition + (1|BrainNucleus) + (1|Condition) + (1|UnitID)');
%         
%             disp(lme);
%             
%             filename = 'Rsquared_Comparison_TimePhaseDomains';
%             save(filename,'lme')
%             clear lme
%             
%             rm = fitrm(data,'Rsquared_log10 ~ BrainNucleus * Domain * Condition + (1|UnitID)');
%             
%             r = ranova(data);
%         end
%         
%     end
    
    %% R-SQUARED FOR COSINE AND POS. VON MISES FITS IN TIME- AND PHASE DOMAINS
    dist_type = {'cosine', 'vonMisesPos'};
    bin_edges  = 10 .^ (-8:0.2:0);
    
    for dist_num = 1:length(dist_type)
        
        for targNum = 1:length(unqTargets)
            T = unqTargets{targNum};
            
            % figure out significantly modulated units by spike analysis
            sig = ...
                ephys_dt.(T).(cfg.condition(cond1_num).name).sig_n_bins > cfg.time.n_sig_bins | ...
                ephys_dt.(T).(cfg.condition(cond2_num).name).sig_n_bins > cfg.time.n_sig_bins;
            
            sig_cond1 = ephys_dt.(T).(cfg.condition(cond1_num).name).sig_n_bins > cfg.time.n_sig_bins;
            sig_cond2 = ephys_dt.(T).(cfg.condition(cond2_num).name).sig_n_bins > cfg.time.n_sig_bins;
            
            % extract data
            time_non_neg1 = ephys_dt.(T).(cfg.condition(cond1_num).name).(dist_type{dist_num}).rsquared > 0; % time domain
            time_non_neg2 = ephys_dt.(T).(cfg.condition(cond2_num).name).(dist_type{dist_num}).rsquared > 0; % time domain
            time_non_neg  = time_non_neg1 & time_non_neg2;
            
            phase_non_neg1 = dt.(T).(cfg.condition(cond1_num).name).(dist_type{dist_num}).rsquared > 0; % phase domain
            phase_non_neg2 = dt.(T).(cfg.condition(cond2_num).name).(dist_type{dist_num}).rsquared > 0; % phase domain
            phase_non_neg  = phase_non_neg1 & phase_non_neg2;
            
            vars(dist_num).(T).Rsq_time_cond1  = ephys_dt.(T).(cfg.condition(cond1_num).name).(dist_type{dist_num}).rsquared(sig_cond1 & time_non_neg & phase_non_neg); % time domain
            vars(dist_num).(T).Rsq_time_cond2  = ephys_dt.(T).(cfg.condition(cond2_num).name).(dist_type{dist_num}).rsquared(sig_cond2 & time_non_neg & phase_non_neg); % time domain
            
            vars(dist_num).(T).Rsq_phase_cond1 = dt.(T).(cfg.condition(cond1_num).name).(dist_type{dist_num}).rsquared(sig_cond1 & time_non_neg & phase_non_neg); % phase domain
            vars(dist_num).(T).Rsq_phase_cond2 = dt.(T).(cfg.condition(cond2_num).name).(dist_type{dist_num}).rsquared(sig_cond2 & time_non_neg & phase_non_neg); % phase domain
            
            % compute histograms
            vars(dist_num).(T).curr_Rsq_time_cond1  = histc(vars(dist_num).(T).Rsq_time_cond1, bin_edges);
            vars(dist_num).(T).curr_Rsq_time_cond2  = histc(vars(dist_num).(T).Rsq_time_cond2, bin_edges);
            
            vars(dist_num).(T).curr_Rsq_phase_cond1 = histc(vars(dist_num).(T).Rsq_phase_cond1, bin_edges);
            vars(dist_num).(T).curr_Rsq_phase_cond2 = histc(vars(dist_num).(T).Rsq_phase_cond2, bin_edges);
            
            % compute means
            vars(dist_num).(T).mean_Rsq_time_cond1  = nanmean(vars(dist_num).(T).Rsq_time_cond1);
            vars(dist_num).(T).mean_Rsq_time_cond2  = nanmean(vars(dist_num).(T).Rsq_time_cond2);
            vars(dist_num).(T).mean_Rsq_phase_cond1 = nanmean(vars(dist_num).(T).Rsq_phase_cond1);
            vars(dist_num).(T).mean_Rsq_phase_cond2 = nanmean(vars(dist_num).(T).Rsq_phase_cond2);
            
            % compute t-tests
            [~, vars(dist_num).p_time_vs_phase_cond1(targNum)]  = ttest (log10(vars(dist_num).(T).Rsq_time_cond1), log10(vars(dist_num).(T).Rsq_phase_cond1));
            [~, vars(dist_num).p_time_vs_phase_cond2(targNum)]  = ttest (log10(vars(dist_num).(T).Rsq_time_cond2), log10(vars(dist_num).(T).Rsq_phase_cond2));
            
            figure
            set(gcf, 'Position', [557 399 1150 597])
            
            for c = 1:length(cfg.condition)
                L = cfg.condition(c).name;
                
                % time domain
                s1 = subplot(2,3,c);
                stairs(bin_edges, vars(dist_num).(T).(['curr_Rsq_time_cond' num2str(c)]), 'Color', cfg.condition(c).color)
                hold on
                plot(vars(dist_num).(T).(['mean_Rsq_time_cond' num2str(c)]), 0, '.', 'MarkerSize', 15, 'Color', cfg.condition(c).color)
                set(gca,'XScale','log')
                title(['Time Domain: ' cfg.condition(c).name '; mean R^2 = ' num2str(vars(dist_num).(T).(['mean_Rsq_time_cond' num2str(c)]))])
                xlabel('R^2')
                ylabel('Unit Counts')
                legend('R^2 Histogram', 'Mean R^2', 'Location', 'best')
                clear curr_Rsq_hist
                
                % phase domain
                s2 = subplot(2,3,3+c);
                stairs(bin_edges, vars(dist_num).(T).(['curr_Rsq_phase_cond' num2str(c)]), 'Color', cfg.condition(c).color)
                hold on
                plot(vars(dist_num).(T).(['mean_Rsq_phase_cond' num2str(c)]), 0, '.', 'MarkerSize', 15, 'Color', cfg.condition(c).color)
                set(gca,'XScale','log')
                title(['Phase Domain: ' cfg.condition(c).name  '; mean R^2 = ' num2str(vars(dist_num).(T).(['mean_Rsq_phase_cond' num2str(c)]))])
                xlabel('R^2')
                ylabel('Unit Counts')
                clear curr_Rsq_hist
                
                linkaxes([s1 s2], 'xy')
                
            end
            
%             % scatter plot - time domain
%             subplot(2,3,3)
%             plot_Rsquared_scatter(log10(vars(dist_num).(T).Rsq_time_cond1), log10(vars(dist_num).(T).Rsq_time_cond2))
%             
%             % scatter plot - phase domain
%             subplot(2,3,6)
%             plot_Rsquared_scatter(log10(vars(dist_num).(T).Rsq_phase_cond1), log10(vars(dist_num).(T).Rsq_phase_cond2))
            
            sgtitle({[cfg.version ' ' T ' ' dist_type{dist_num}], ...
                ['Paired t-test p (Rest) = ' num2str(vars(dist_num).p_time_vs_phase_cond1(targNum)) '; ' ...
                'p (Task) = ' num2str(vars(dist_num).p_time_vs_phase_cond2(targNum))]}, ...
                'interpreter', 'none')
            
            filename = ['Rsquared_' dist_type{dist_num} '_time_vs_phase_' T '_' cfg.condition(c).name];
            export_fig(gcf, [basepath_to_save, filesep ,filename], '-pdf');
            close(gcf);
        
        end
        
        % create a table with p-values for t-test
        [h, crit_p] = fdr_bky([vars(dist_num).p_time_vs_phase_cond1 vars(dist_num).p_time_vs_phase_cond2], 0.05, 'yes');
        t = table(unqTargets', ...
            vars(dist_num).p_time_vs_phase_cond1', vars(dist_num).p_time_vs_phase_cond2', ...
            h(1:3)', h(4:6)', [crit_p; nan; nan], ...
            'VariableNames',{'Target', 'Rest: p', 'Task: p', 'Rest: h', 'Task: h', 'crit. p'});
        filename = ['Table_Rsquared_time_vs_phase_' dist_type{dist_num}];
        writetable(t, [basepath_to_save filesep filename '.xls'])
        
    end
    
    if targNum == length(unqTargets)
        save([basepath_to_save filesep 'Comparison_Rsq_time_vs_phase.mat'], 'vars')
        clear vars
    end
    
    %% PEAK PHASES FOR COSINE FITS - COMPARISON BETWEEN LOWIBI AND HIGHIBI
    phase_bin_edges = 0:0.2:2*pi;
    field_name = {'lowIBI_', 'highIBI_'};
    dist_type = {'cosine', 'vonMisesPos'};
    
    for dist_num = 1:length(dist_type)
    
        for targNum = 1:length(unqTargets)
            T = unqTargets{targNum};
            
            % figure out significantly modulated units by spike analysis
            sig = ...
                ephys_dt.(T).(cfg.condition(cond1_num).name).sig_n_bins > cfg.time.n_sig_bins | ...
                ephys_dt.(T).(cfg.condition(cond2_num).name).sig_n_bins > cfg.time.n_sig_bins;
            
            for lin_field_num = 1:length(field_name)
                
                if strcmp(dist_type{dist_num},'cosine')
                    coef_num = 2;
                elseif strcmp(dist_type{dist_num},'vonMisesPos')
                    coef_num = 4;
                end
                
                % extract data
                vars(targNum).peak_time_cond1{lin_field_num}  = mod(ephys_dt.(T).(cfg.condition(cond1_num).name).([field_name{lin_field_num} dist_type{dist_num}]).coefs(sig,coef_num), 2*pi);
                vars(targNum).peak_time_cond2{lin_field_num}  = mod(ephys_dt.(T).(cfg.condition(cond2_num).name).([field_name{lin_field_num} dist_type{dist_num}]).coefs(sig,coef_num), 2*pi);
                vars(targNum).peak_phase_cond1{lin_field_num} = dt.(T).(cfg.condition(cond1_num).name).([field_name{lin_field_num} dist_type{dist_num}]).coefs(sig,coef_num);
                vars(targNum).peak_phase_cond2{lin_field_num} = dt.(T).(cfg.condition(cond2_num).name).([field_name{lin_field_num} dist_type{dist_num}]).coefs(sig,coef_num);
                
                % mean peak phases
                vars(targNum).meanPeak_time_cond1{lin_field_num}  = circ_mean(rmmissing(vars(targNum).peak_time_cond1{lin_field_num}));
                vars(targNum).meanPeak_time_cond2{lin_field_num}  = circ_mean(rmmissing(vars(targNum).peak_time_cond2{lin_field_num}));
                vars(targNum).meanPeak_phase_cond1{lin_field_num} = circ_mean(rmmissing(vars(targNum).peak_phase_cond1{lin_field_num}));
                vars(targNum).meanPeak_phase_cond2{lin_field_num} = circ_mean(rmmissing(vars(targNum).peak_phase_cond2{lin_field_num}));
                
                % compute histograms
                vars(targNum).curr_peak_time_cond1{lin_field_num}  = histc(vars(targNum).peak_time_cond1{lin_field_num}, phase_bin_edges);
                vars(targNum).curr_peak_time_cond2{lin_field_num}  = histc(vars(targNum).peak_time_cond2{lin_field_num}, phase_bin_edges);
                
                vars(targNum).curr_peak_phase_cond1{lin_field_num} = histc(vars(targNum).peak_phase_cond1{lin_field_num}, phase_bin_edges);
                vars(targNum).curr_peak_phase_cond2{lin_field_num} = histc(vars(targNum).peak_phase_cond2{lin_field_num}, phase_bin_edges);
                
            end
            
            % compute circular distances
            vars(targNum).circ_dist_time_cond1  = circ_dist(vars(targNum).peak_time_cond1{1},vars(targNum).peak_time_cond1{2});
            vars(targNum).circ_dist_phase_cond1 = circ_dist(vars(targNum).peak_phase_cond1{1},vars(targNum).peak_phase_cond1{2});
            vars(targNum).circ_dist_time_cond2  = circ_dist(vars(targNum).peak_time_cond2{1},vars(targNum).peak_time_cond2{2});
            vars(targNum).circ_dist_phase_cond2 = circ_dist(vars(targNum).peak_phase_cond2{1},vars(targNum).peak_phase_cond2{2});
            
            % compute Watson-Williams Tests between lowIBI and highIBI
            % conditions
            non_nan_cond1 = ...
                ~isnan(vars(targNum).peak_time_cond1{1}) & ~isnan(vars(targNum).peak_time_cond1{2}) & ...
                ~isnan(vars(targNum).peak_phase_cond1{1}) & ~isnan(vars(targNum).peak_phase_cond1{1});
            
            non_nan_cond2 = ...
                ~isnan(vars(targNum).peak_time_cond2{1}) & ~isnan(vars(targNum).peak_time_cond2{2}) & ...
                ~isnan(vars(targNum).peak_phase_cond2{1}) & ~isnan(vars(targNum).peak_phase_cond2{1});
            
            vars(targNum).p_time_lowIBI_vs_highIBI(cond1_num)  = circ_wwtest(vars(targNum).peak_time_cond1{1}(non_nan_cond1), vars(targNum).peak_time_cond1{2}(non_nan_cond1));
            vars(targNum).p_phase_lowIBI_vs_highIBI(cond1_num) = circ_wwtest(vars(targNum).peak_phase_cond1{1}(non_nan_cond1), vars(targNum).peak_phase_cond1{2}(non_nan_cond1));
            vars(targNum).p_time_lowIBI_vs_highIBI(cond2_num)  = circ_wwtest(vars(targNum).peak_time_cond2{1}(non_nan_cond2), vars(targNum).peak_time_cond2{2}(non_nan_cond2));
            vars(targNum).p_phase_lowIBI_vs_highIBI(cond2_num) = circ_wwtest(vars(targNum).peak_phase_cond2{1}(non_nan_cond2), vars(targNum).peak_phase_cond2{2}(non_nan_cond2));
            
            % save data
            if targNum == length(unqTargets)
                % control for multiple comparisons - FDR
                A = ...
                    [vars(1).p_time_lowIBI_vs_highIBI;
                    vars(2).p_time_lowIBI_vs_highIBI;
                    vars(3).p_time_lowIBI_vs_highIBI;
                    ...
                    vars(1).p_phase_lowIBI_vs_highIBI;
                    vars(2).p_phase_lowIBI_vs_highIBI;
                    vars(3).p_phase_lowIBI_vs_highIBI];
                
                [h,p_crit] = fdr_bky(A);
                
                % create and save the table with p-values and significance
                t = table(unqTargets', ...
                    [vars(1).p_time_lowIBI_vs_highIBI(1) vars(2).p_time_lowIBI_vs_highIBI(1) vars(3).p_time_lowIBI_vs_highIBI(1)]', ...    % p-val in Rest: peak - time domain: lowIBI vs. highIBI
                    h(1:3,1), ...                                                                                                          % h in Rest: peak - time domain: lowIBI vs. highIBI
                    [vars(1).p_phase_lowIBI_vs_highIBI(1) vars(2).p_phase_lowIBI_vs_highIBI(1) vars(3).p_phase_lowIBI_vs_highIBI(1)]', ... % p-val in Rest: peak - phase domain: lowIBI vs. highIBI
                    h(4:end,1), ...                                                                                                        % h in Rest: peak - phase domain: lowIBI vs. highIBI
                    ...
                    [vars(1).p_time_lowIBI_vs_highIBI(2) vars(2).p_time_lowIBI_vs_highIBI(2) vars(3).p_time_lowIBI_vs_highIBI(2)]', ...    % p-val in Task: peak - time domain: lowIBI vs. highIBI
                    h(1:3,2), ...                                                                                                          % h in Task: peak - time domain: lowIBI vs. highIBI
                    [vars(1).p_phase_lowIBI_vs_highIBI(2) vars(2).p_phase_lowIBI_vs_highIBI(2) vars(3).p_phase_lowIBI_vs_highIBI(2)]', ... % p-val in Task: peak - phase domain: lowIBI vs. highIBI
                    h(4:end,2), ...                                                                                                        % h in Task: peak - phase domain: lowIBI vs. highIBI
                    [p_crit; nan; nan], ...
                    'VariableNames',...
                    {'Target', 'Rest/Time-Domain/p', 'Rest/Time-Domain/h', 'Rest/Phase-Domain/p', 'Rest/Phase-Domain/h', ...
                    'Task/Time-Domain/p', 'Task/Time_Domain/h', 'Task/Phase-Domain/p', 'Task/Phase-Domain/h', ...
                    'crit. p'});
                filename = ['Table_PeakPhase_time_vs_phase_' dist_type{dist_num}];
                writetable(t, [basepath_to_save filesep filename '.xls'])
                clear t A h p_crit
                
                % create and save table with mean angles
                t2 = ...
                    table(unqTargets', ...
                    [vars(1).meanPeak_time_cond1{1} vars(2).meanPeak_time_cond1{1} vars(3).meanPeak_time_cond1{1}]', ... % Time-Domain: Rest - lowIBI
                    [vars(1).meanPeak_time_cond1{2} vars(2).meanPeak_time_cond1{2} vars(3).meanPeak_time_cond1{2}]', ... % Time-Domain: Rest - highIBI
                    [vars(1).meanPeak_time_cond2{1} vars(2).meanPeak_time_cond2{1} vars(3).meanPeak_time_cond2{1}]', ... % Time-Domain: Task - lowIBI
                    [vars(1).meanPeak_time_cond2{2} vars(2).meanPeak_time_cond2{2} vars(3).meanPeak_time_cond2{2}]', ... % Time-Domain: Task - highIBI
                    ...
                    [vars(1).meanPeak_phase_cond1{1} vars(2).meanPeak_phase_cond1{1} vars(3).meanPeak_phase_cond1{1}]', ... % Phase-Domain: Rest - lowIBI
                    [vars(1).meanPeak_phase_cond1{2} vars(2).meanPeak_phase_cond1{2} vars(3).meanPeak_phase_cond1{2}]', ... % Phase-Domain: Rest - highIBI
                    [vars(1).meanPeak_phase_cond2{1} vars(2).meanPeak_phase_cond2{1} vars(3).meanPeak_phase_cond2{1}]', ... % Phase-Domain: Task - lowIBI
                    [vars(1).meanPeak_phase_cond2{2} vars(2).meanPeak_phase_cond2{2} vars(3).meanPeak_phase_cond2{2}]', ... % Phase-Domain: Task - highIBI
                    'VariableNames', ...
                    {'Target', ...
                    'Time-Domain/Rest/lowIBI', 'Time-Domain/Rest/highIBI', 'Time-Domain/Task/lowIBI', 'Time-Domain/Task/highIBI', ...
                    'Phase-Domain/Rest/lowIBI', 'Phase-Domain/Rest/highIBI', 'Phase-Domain/Task/lowIBI', 'Phase-Domain/Task/highIBI'});
                filename = ['Table_MeanPeakPhase_time_vs_phase_' dist_type{dist_num}];
                writetable(t2, [basepath_to_save filesep filename '.xls'])
                clear t2
                
                % save vars structure
                save([basepath_to_save filesep 'Comparison_phase_lowIBI_vs_highIBI_' dist_type{dist_num} '.mat'], 'vars')
                %             clear vars
            end
            
            
            f0 = figure;
            set(f0, 'Position', [557   620   700   400])
            
            for lin_field_num = 1:length(field_name)
                
                % time domain
                s1 = subplot(length(cfg.condition),length(field_name), 2*(lin_field_num-1)+1);
                stairs(phase_bin_edges, vars(targNum).(['curr_peak_time_cond' num2str(c)]){lin_field_num}, 'Color', cfg.condition(c).color)
                %                 set(gca,'XScale','log')
                if lin_field_num == 1
                    title({['Time Domain: ' cfg.condition(c).name ': ' field_name{lin_field_num}(1:end-1)], ...
                        ['WW-test p = ' num2str(vars(targNum).p_time_lowIBI_vs_highIBI(c))]})
                else
                    title(['Time Domain: ' cfg.condition(c).name ': ' field_name{lin_field_num}(1:end-1)])
                end
                xlabel('Peak Phase, radians')
                clear curr_Rsq_hist
                
                % phase domain
                s2 = subplot(length(cfg.condition),length(field_name),2*(lin_field_num-1)+2);
                stairs(phase_bin_edges, vars(targNum).(['curr_peak_phase_cond' num2str(c)]){lin_field_num}, 'Color', cfg.condition(c).color)
                %                 set(gca,'XScale','log')
                if lin_field_num == 1
                    title({['Phase Domain: ' cfg.condition(c).name ': ' field_name{lin_field_num}(1:end-1)], ...
                        ['WW-test p = ' num2str(vars(targNum).p_phase_lowIBI_vs_highIBI(c))]})
                else
                    title(['Phase Domain: ' cfg.condition(c).name ': ' field_name{lin_field_num}(1:end-1)])
                end
                xlabel('Peak Phase, radians')
                clear curr_Rsq_hist
                
                sgtitle([cfg.version ' ' T], 'interpreter', 'none')
                
                linkaxes([s1 s2], 'xy')
                
            end
            
            filename = ['Rsquared_time_vs_phase_' T '_' cfg.condition(c).name '_' dist_type{dist_num}];
            export_fig(f0, [basepath_to_save,filesep ,filename], '-pdf');
            close(f0);
            
            f1 = figure;
            set(f1, 'Position', [557   620   700   400])
            
            for c = 1:length(cfg.condition)
                L = cfg.condition(c).name;
                
                figure(f1);
                subplot(2,2,1+2*(c-1))
                hist([vars(targNum).(['circ_dist_time_cond' num2str(c)]) vars(targNum).(['circ_dist_phase_cond' num2str(c)])])
                C = colororder;
                hold on
                scatter(circ_mean(vars(targNum).(['circ_dist_time_cond' num2str(c)])), 0, [], 'MarkerFaceColor', 'b')
                scatter(circ_mean(vars(targNum).(['circ_dist_phase_cond' num2str(c)])), 0, [], 'MarkerFaceColor', 'y')
                xlabel('Phase Difference Related to R-peak')
                ylabel('Unit Counts')
                [~,p] = ttest(vars(targNum).(['circ_dist_time_cond' num2str(c)]), vars(targNum).(['circ_dist_phase_cond' num2str(c)]));
                title([T ': t-test p = ' num2str(p)])
                
                legend('Circular Distance (Time Domain)', 'Circular Distance (Phase Domain)', 'Mean Dist. (Time Domain)', 'Mean Dist. (Phase Domain)')
                
                
                subplot(2,2,2*c)
                scatter(vars(targNum).(['circ_dist_time_cond' num2str(c)]), vars(targNum).(['circ_dist_phase_cond' num2str(c)]))
                xlim([-pi pi])
                ylim([-pi pi])
                axis square
                box on
                hold on
                
                nan_ids1 = isnan(vars(targNum).(['circ_dist_time_cond' num2str(c)]));
                nan_ids2 = isnan(vars(targNum).(['circ_dist_phase_cond' num2str(c)]));
                not_nan  = ~nan_ids1 & ~nan_ids2;
                
                %         [cc1, pp1] = corrcoef(variables.(['circ_dist_time_cond' num2str(c)]), variables.(['circ_dist_phase_cond' num2str(c)]),'Row','pairwise');
                % compute circular correlation coefficients
                [cc, pp] = circ_corrcc(vars(targNum).(['circ_dist_time_cond' num2str(c)])(not_nan), vars(targNum).(['circ_dist_phase_cond' num2str(c)])(not_nan));
                
                p = polyfit(vars(targNum).(['circ_dist_time_cond' num2str(c)])(not_nan), vars(targNum).(['circ_dist_phase_cond' num2str(c)])(not_nan),1);
                
                x_lim = get(gca,'XLim');
                plot(x_lim, p(1)*x_lim + p(2))
                xlabel('Rest: Circular Differences')
                ylabel('Task: Circular Differences')
                title({['circ. cc = ' num2str(cc) '; circ. p = ' num2str(pp)], ...
                    ['y = ' num2str(p(1)) 'x + ' num2str(p(2))]})
                
                if c == 2
                    filename = ['PhaseDifferences_time_vs_phase_' T '_' dist_type{dist_num}];
                    export_fig(f1, [basepath_to_save,filesep ,filename], '-pdf');
                    close(f1);
                end
                
            end
            
        end
    
    end
    
    %% PARAMETERS BY CONDITION
    for targNum = 1:length(unqTargets)
            T = unqTargets{targNum};
            
            %% stability vs. linear gof
            plot_stability_vs_linear_gof(dt.(T), {'linear'}, cfg, cond1_num, cond2_num, T, 'Rsq_vs_Stability_', basepath_to_save)
            
            plot_stability_vs_linear_gof(dt.(T), {'lowIBI_linear', 'highIBI_linear'}, cfg, cond1_num, cond2_num, T, 'MedSplit_Rsq_vs_Stability_', basepath_to_save)
            
            %% R-interval numbers
            plot_RpeakNum_histograms(ephys_dt.(T), '', cfg, cond1_num, cond2_num, T, 'HeartCycleNumber_', basepath_to_save)
            
            %% Rsquared histograms
            plot_Rsquared_histograms(dt.(T), {''}, cfg, cond1_num, cond2_num, C, T, 'R-squared_', basepath_to_save)
            
            plot_Rsquared_histograms(dt.(T), {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, T, 'MedSplit_R-squared_', basepath_to_save)
            
            %% Phases of the peak from different fits
            plot_phase_histograms(dt.(T), {''}, cfg, cond1_num, cond2_num, C, T, 'PeakPhase_' , basepath_to_save)
            
            plot_phase_histograms(dt.(T), {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, T, 'MedSplit_PeakPhase_' , basepath_to_save)
            
            %% scaling factor for von Mises fits
            plot_scalingFactor_histograms(dt.(T), {''}, cfg, cond1_num, cond2_num, C, T, 'ScalingFactors_', basepath_to_save)
            
            plot_scalingFactor_histograms(dt.(T), {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, T, 'MedSplit_ScalingFactors_', basepath_to_save)
            
            %% kappa - concentration parameter for von Mises fits
            plot_kappa_histograms(dt.(T), {''}, cfg, cond1_num, cond2_num, C, T, 'Kappas_', basepath_to_save)
            
            plot_kappa_histograms(dt.(T), {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, T, 'MedSplit_Kappas_', basepath_to_save)
            
            %% kappas vs. phase for von Mises fits
            plot_kappas_vs_phase(dt.(T), {''}, cfg, cond1_num, cond2_num, C, T, 'Kappa_vs_Phase_', basepath_to_save)
    
    end
    %% figure out the better fit (von Mises Pos. vs. von Mises Neg.)
    % 2. sum of pvalues is lower for the better fit
    % 3. at least one of the conditions (Rest or Task) is fitted
    % significantly (p < 0.01)
    
    % for only one condition - significantly fitted with a linear
    % function
    % 1) pos linear fits - cond1 - rest
    lin_sig_pos_cond1 = ...
        ( dt.(cfg.condition(cond1_num).name).linear.pvalue(:,2) < 0.05 & ... % significant fit
        dt.(cfg.condition(cond1_num).name).linear.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond1_num).name).linear.coefs(:,2) > 0 & ...
        ...
        ( (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared) & ...
        (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared) );
    
    % 2) pos linear fits - cond2 - task
    lin_sig_pos_cond2 = ...
        ( dt.(cfg.condition(cond2_num).name).linear.pvalue(:,2) < 0.05 & ... % significant fit
        dt.(cfg.condition(cond2_num).name).linear.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond2_num).name).linear.coefs(:,2) > 0 & ...
        ...
        ( (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared) & ...
        (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared) );
    
    % 3) neg linear fits - cond1 - rest
    lin_sig_neg_cond1 = ...
        ( dt.(cfg.condition(cond1_num).name).linear.pvalue(:,2) < 0.05 & ... % significant fit
        dt.(cfg.condition(cond1_num).name).linear.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond1_num).name).linear.coefs(:,2) < 0 & ...
        ...
        ( (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared) & ...
        (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared) );
    
    % 4) neg linear fits - cond2 - task
    lin_sig_neg_cond2 = ...
        ( dt.(cfg.condition(cond2_num).name).linear.pvalue(:,2) < 0.05 & ... % significant fit
        dt.(cfg.condition(cond2_num).name).linear.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond2_num).name).linear.coefs(:,2) < 0 & ...
        ...
        ( (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared) & ...
        (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared) );
    
    % 5) von Mises pos - cond1 - rest
    % check that these units are fitted better with a circular fit in
    % both condition as I'm going to segregate units with at least one
    % linear fit
    vmpos_sig_cond1 = ...
        ( dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue < 0.05 & ...
        dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared;
    
    
    % 6) von Mises pos - cond2 - task
    % check that these units are fitted better with a circular fit in
    % both condition as I'm going to segregate units with at least one
    % linear fit
    vmpos_sig_cond2 = ...
        ( dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue < 0.05 & ...
        dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared;
    
    % 7) von Mises neg - cond1 - rest
    % check that these units are fitted better with a circular fit in
    % both condition as I'm going to segregate units with at least one
    % linear fit
    vmneg_sig_cond1 = ...
        ( dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue < 0.05 & ...
        dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared;
    
    % 8) von Mises neg - cond2 - task
    % check that these units are fitted better with a circular fit in
    % both condition as I'm going to segregate units with at least one
    % linear fit
    vmneg_sig_cond2 = ...
        ( dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue < 0.05 & ...
        dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > 0 ) & ...
        ...
        dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
        dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared;
    
    % all nonsig units
    selection.all_nonsig_ids = ...
        ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        ~vmpos_sig_cond1 & ~vmpos_sig_cond2 & ~vmneg_sig_cond1 & ~vmneg_sig_cond2;
    
    % pos von Mises
    selection.vmpos_sig_cond1_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        vmpos_sig_cond1 & ~vmpos_sig_cond2;
    selection.vmpos_sig_cond2_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        ~vmpos_sig_cond1 & vmpos_sig_cond2;
    selection.vmneg_sig_cond1_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        vmneg_sig_cond1 & ~vmneg_sig_cond2;
    selection.vmneg_sig_cond2_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        ~vmneg_sig_cond1 & vmneg_sig_cond2;
    
    selection.vmpos_sig_both_ids  = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        vmpos_sig_cond1 & vmpos_sig_cond2;
    selection.vmneg_sig_both_ids  = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        vmneg_sig_cond1 & vmneg_sig_cond2;
    
    selection.vmpos_nonsig_ids    = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        ~vmpos_sig_cond1 & ~vmpos_sig_cond2;
    selection.vmneg_nonsig_ids    = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
        ~vmneg_sig_cond1 & ~vmneg_sig_cond2;
    
    %% figure out the preferred type of fit
    vm_pos_ids = (dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared + dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared) > ...
        (dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared + dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared) & ...
        (min([dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue],[],2) < 0.05 | ...
        min([dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue],[],2) < 0.05);
    
    vm_neg_ids = (dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared + dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared) < ...
        (dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared + dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared) & ...
        (min([dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue],[],2) < 0.05 | ...
        min([dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue],[],2) < 0.05);
    
    % R-squared for significant units
    bin_edges       = 10 .^ (-6:0.2:0);
    % R-squared
    lin_Rsq_cond1   = dt.(cfg.condition(cond1_num).name).linear.rsquared(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2);
    lin_Rsq_cond2   = dt.(cfg.condition(cond1_num).name).linear.rsquared(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2);
    
    vMPos_Rsq_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids);
    vMNeg_Rsq_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids);
    
    vMPos_Rsq_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids);
    vMNeg_Rsq_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids);
    
    % bins
    lin_RsqBins_cond1   = histc(lin_Rsq_cond1, bin_edges);
    lin_RsqBins_cond2   = histc(lin_Rsq_cond2, bin_edges);
    
    vMPos_RsqBins_cond1 = histc(vMPos_Rsq_cond1, bin_edges);
    vMNeg_RsqBins_cond1 = histc(vMNeg_Rsq_cond1, bin_edges);
    
    vMPos_RsqBins_cond2 = histc(vMPos_Rsq_cond2, bin_edges);
    vMNeg_RsqBins_cond2 = histc(vMNeg_Rsq_cond2, bin_edges);
    
    % medians
    med_rest = [median(vMPos_Rsq_cond1) median(vMNeg_Rsq_cond1)];
    med_task = [median(dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids)) ...
        median(dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids))];
    
    vMPos_Rsq_cond1     = vMPos_Rsq_cond1(:);
    vMNeg_RsqBins_cond1 = vMNeg_RsqBins_cond1(:);
    vMNeg_RsqBins_cond2 = vMNeg_RsqBins_cond2(:);
    lin_RsqBins_cond1   = lin_RsqBins_cond1(:);
    lin_RsqBins_cond2   = lin_RsqBins_cond2(:);
    
    figure,
    set(gcf, 'Position', [557   620   754   206])
    colororder(C(2:end,:))
    s1 = subplot(1,2,1);
    stairs(bin_edges, [vMPos_RsqBins_cond1 vMNeg_RsqBins_cond1 lin_RsqBins_cond1], 'LineWidth', 2)
    hold on
    plot(med_rest(1), 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(med_rest(2), 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    set(gca,'XScale','log')
    xlabel('log10(R^2)')
    ylabel('Number of Units')
    title(cfg.condition(cond1_num).name)
    legend('von Mises Pos.', 'von Mises Neg.', 'Linear Fit', 'Location', 'Best')
    
    s2 = subplot(1,2,2);
    stairs(bin_edges, [vMPos_RsqBins_cond2 vMNeg_RsqBins_cond2 lin_RsqBins_cond2], 'LineWidth', 2)
    hold on
    plot(med_task(1), 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(med_task(2), 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    set(gca,'XScale','log')
    xlabel('log10(R^2)')
    title(cfg.condition(cond2_num).name)
    
    linkaxes([s1 s2], 'xy')
    
    sgtitle([cfg.version ' ' T], 'interpreter', 'none')
    
    filename = ['R-squaredSig_' T '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close all;
    
    if ~isempty(vMPos_Rsq_cond1) & ~isempty(vMPos_Rsq_cond2)
        p_wilcox(1) = signrank(vMPos_Rsq_cond1, vMPos_Rsq_cond2);
    else
        p_wilcox(1) = NaN;
    end
    
    if ~isempty(vMNeg_Rsq_cond1) & ~isempty(vMNeg_Rsq_cond2)
        p_wilcox(2) = signrank(vMNeg_Rsq_cond1, vMNeg_Rsq_cond2);
    else
        p_wilcox(2) = NaN;
    end
    
    T = table({'Pos. von Mises', 'Neg. von Mises'}', ...
        med_rest', med_task', p_wilcox', ...
        'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median R^2'], [cfg.condition(cond2_num).name ': median R^2'], 'Paired Wilcoxon''s p'});
    writetable(T, [basepath_to_save filesep filename '.xls'])
    clear p_wilcox
    
    %% scaling factors of significant fits
    plot_scalingFactor_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, T, 'SigScalingFactorDifference_', basepath_to_save, selection)
    
    plot_scalingFactor_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, T, 'MedSplit_SigScalingFactorDifference_', basepath_to_save, selection)
    
    % [Pos. von Mises] Phases of significant fits
    phase_bin_edges = [0:0.2:2*pi];
    
    cond1_pos.counts_sig_cond1       = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids,4), phase_bin_edges);
    cond1_pos.counts_sig_cond2       = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond2_ids,4), phase_bin_edges);
    cond1_pos.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_both_ids,4), phase_bin_edges);
    cond1_pos.counts_nonsig          = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_nonsig_ids,4), phase_bin_edges);
    
    cond2_pos.counts_sig_cond1       = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids,4), phase_bin_edges);
    cond2_pos.counts_sig_cond2       = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond2_ids,4), phase_bin_edges);
    cond2_pos.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_both_ids,4), phase_bin_edges);
    cond2_pos.counts_nonsig          = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_nonsig_ids,4), phase_bin_edges);
    
    cond1_pos.counts_nonsig    = cond1_pos.counts_nonsig(:);
    cond1_pos.counts_sig_cond1 = cond1_pos.counts_sig_cond1(:);
    cond1_pos.counts_sig_cond2 = cond1_pos.counts_sig_cond2(:);
    cond1_pos.counts_sig_cond1_cond2 = cond1_pos.counts_sig_cond1_cond2(:);
    
    cond2_pos.counts_sig_cond2 = cond2_pos.counts_sig_cond2(:);
    cond2_pos.counts_nonsig    = cond2_pos.counts_nonsig(:);
    cond2_pos.counts_sig_cond2 = cond2_pos.counts_sig_cond2(:);
    cond2_pos.counts_sig_cond1 = cond2_pos.counts_sig_cond1(:);
    cond2_pos.counts_sig_cond1_cond2 = cond2_pos.counts_sig_cond1_cond2(:);
    
    figure,
    set(gcf, 'Position', [557   620   754   206])
    s1 = subplot(1,2,1);
    colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
    bar([0:0.2:2*pi], [cond1_pos.counts_sig_cond1 cond1_pos.counts_sig_cond2 cond1_pos.counts_sig_cond1_cond2 cond1_pos.counts_nonsig], 'stacked')
    xlabel('Heart Cycle Phase [0-2\pi]')
    ylabel('# Units')
    legend({['Only ' cfg.condition(cond1_num).name], ['Only ' cfg.condition(cond2_num).name], ['Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name], 'None'}, 'Location', 'best')
    title(cfg.condition(cond1_num).name)
    
    s2 = subplot(1,2,2);
    colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
    bar([0:0.2:2*pi], [cond2_pos.counts_sig_cond1 cond2_pos.counts_sig_cond2 cond2_pos.counts_sig_cond1_cond2 cond2_pos.counts_nonsig], 'stacked')
    xlabel('Heart Cycle Phase [0-2\pi]')
    title(cfg.condition(cond2_num).name)
    
    linkaxes([s1 s2], 'xy')
    
    sgtitle([T ': Pos. von Mises'])
    
    filename = ['SigPhasesVMPos_' T '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close all;
    
    % [Neg. von Mises]
    cond1_neg.counts_sig_cond1       = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids,4), phase_bin_edges);
    cond1_neg.counts_sig_cond2       = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond2_ids,4), phase_bin_edges);
    cond1_neg.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_both_ids,4), phase_bin_edges);
    cond1_neg.counts_nonsig          = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_nonsig_ids,4), phase_bin_edges);
    
    cond2_neg.counts_sig_cond1       = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids,4), phase_bin_edges);
    cond2_neg.counts_sig_cond2       = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond2_ids,4), phase_bin_edges);
    cond2_neg.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_both_ids,4), phase_bin_edges);
    cond2_neg.counts_nonsig          = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_nonsig_ids,4), phase_bin_edges);
    
    cond1_neg.counts_nonsig          = cond1_neg.counts_nonsig(:);
    cond1_neg.counts_sig_cond1       = cond1_neg.counts_sig_cond1(:);
    cond1_neg.counts_sig_cond2       = cond1_neg.counts_sig_cond2(:);
    cond1_neg.counts_sig_cond1_cond2 = cond1_neg.counts_sig_cond1_cond2(:);
    cond2_neg.counts_sig_cond2       = cond2_neg.counts_sig_cond2(:);
    cond2_neg.counts_nonsig          = cond2_neg.counts_nonsig(:);
    cond2_neg.counts_sig_cond2       = cond2_neg.counts_sig_cond2(:);
    cond2_neg.counts_sig_cond1       = cond2_neg.counts_sig_cond1(:);
    cond2_neg.counts_sig_cond1_cond2 = cond2_neg.counts_sig_cond1_cond2(:);
    
    figure,
    set(gcf, 'Position', [557   620   754   206])
    s1 = subplot(1,2,1);
    colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
    bar(phase_bin_edges, [cond1_neg.counts_sig_cond1 cond1_neg.counts_sig_cond2 cond1_neg.counts_sig_cond1_cond2 cond1_neg.counts_nonsig], 'stacked')
    xlabel('Heart Cycle Phase [0-2\pi]')
    ylabel('# Units')
    legend({['Only ' cfg.condition(cond1_num).name], ['Only ' cfg.condition(cond2_num).name], ['Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name], 'None'}, 'Location', 'best')
    title(cfg.condition(cond1_num).name)
    
    s2 = subplot(1,2,2);
    colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
    bar(phase_bin_edges, [cond2_neg.counts_sig_cond1 cond2_neg.counts_sig_cond2 cond2_neg.counts_sig_cond1_cond2 cond2_neg.counts_nonsig], 'stacked')
    xlabel('Heart Cycle Phase [0-2\pi]')
    title(cfg.condition(cond2_num).name)
    
    linkaxes([s1 s2], 'xy')
    
    sgtitle([T ': Neg. von Mises'])
    
    filename = ['SigPhasesVMNeg_' T '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close all;
    
    
    % plot number of heart-responsive units in the rest and in the task
    pos_num_sig_cond1(targNum)          = 100 * sum(selection.vmpos_sig_cond1_ids) / length(curr_unit_ids);
    pos_num_sig_cond2(targNum)          = 100 * sum(selection.vmpos_sig_cond2_ids) / length(curr_unit_ids);
    pos_num_sig_cond1_cond2(targNum)    = 100 * sum(selection.vmpos_sig_both_ids) / length(curr_unit_ids);
    
    neg_num_sig_cond1(targNum)          = 100 * sum(selection.vmneg_sig_cond1_ids) / length(curr_unit_ids);
    neg_num_sig_cond2(targNum)          = 100 * sum(selection.vmneg_sig_cond2_ids) / length(curr_unit_ids);
    neg_num_sig_cond1_cond2(targNum)    = 100 * sum(selection.vmneg_sig_both_ids) / length(curr_unit_ids);
    
    all_num_nonsig(targNum)             = 100 * sum(selection.all_nonsig_ids) / length(curr_unit_ids);
    lin_num_sig_cond1_or_cond2(targNum) = 100 * sum(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) / length(curr_unit_ids);
    
    if targNum == length(unqTargets)
        figure,
        pos_colors = brighten([0 0 0.5; 0.5 0 0; 0.5 0 0.5], 0.5);
        neg_colors = brighten([0 0 0.5; 0.5 0 0; 0.5 0 0.5], -0.25);
        colororder([pos_colors; 1 1 1; 0 0 0; neg_colors]) % blue, red, magenta, white, black
        bar([pos_num_sig_cond1; pos_num_sig_cond2; pos_num_sig_cond1_cond2; ...
            all_num_nonsig; lin_num_sig_cond1_or_cond2; ...
            neg_num_sig_cond1; neg_num_sig_cond2; neg_num_sig_cond1_cond2]', 'stacked')
        set(gca, 'XTickLabel', unqTargets)
        legend({['Pos. von Mises: Only ' cfg.condition(cond1_num).name], ...
            ['Pos. von Mises: Only ' cfg.condition(cond2_num).name], ...
            ['Pos. von Mises: Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name], ...
            'None', ...
            'Linear Slopes', ...
            ['Neg. von Mises: Only ' cfg.condition(cond1_num).name], ...
            ['Neg. von Mises: Only ' cfg.condition(cond2_num).name], ...
            ['Neg. von Mises: Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name]
            }, 'Location', 'southoutside')
        ylabel('Percentage of Units')
        
        filename = ['HistogramPos_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
        
    end
    
    %% plot significant kappas
    plot_kappa_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, T, 'SigKappas_', basepath_to_save, selection)
    
    plot_kappa_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, T, 'MedSplit_SigKappas_', basepath_to_save, selection)
    
    %% plot significant kappas vs phases
    plot_kappas_vs_phase(dt, {''}, cfg, cond1_num, cond2_num, C, T, 'SigKappa_vs_SigPhase_', basepath_to_save)
    
    
    %% compare circular distribution parameters
    disp(['Compare Means: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Pos. von Mises'])
    circ_wwtest(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4), dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(vm_pos_ids,4))
    disp(['Compare Concentrations: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Pos. von Mises'])
    circ_ktest(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4)', dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(vm_pos_ids,4)')
    
    disp(['Compare Means: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Neg. von Mises'])
    circ_wwtest(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4), dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(vm_neg_ids,4))
    disp(['Compare Concentrations: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Neg. von Mises'])
    circ_ktest(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4), dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(vm_neg_ids,4))
    
    % compare mean phases by area
    phase_by_area_cond1_vMpos = [phase_by_area_cond1_vMpos; dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4)];
    phase_by_area_cond2_vMpos = [phase_by_area_cond2_vMpos; dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(vm_pos_ids,4)];
    
    phase_by_area_cond1_vMneg = [phase_by_area_cond1_vMneg; dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4)];
    phase_by_area_cond2_vMneg = [phase_by_area_cond2_vMneg; dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(vm_neg_ids,4)];
    
    area_id_vMpos            = [area_id_vMpos targNum*ones(1, length(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4)))];
    area_id_vMneg            = [area_id_vMneg targNum*ones(1, length(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4)))];
    
    clear cond1 cond2
end

%% 2-way ANOVA for pos. von Mises distributions
condition_id_vMpos = [zeros(1, length(phase_by_area_cond1_vMpos)) ones(1, length(phase_by_area_cond2_vMpos))];
[pval_vMpos, stats_vMpos] = circ_hktest([phase_by_area_cond1_vMpos; phase_by_area_cond2_vMpos], [area_id_vMpos area_id_vMpos], condition_id_vMpos, 1, {'Area: VPL, dPul, MD', 'Condition: Rest, Task', 'Area * Condition'});

T = table(stats_vMpos);
writetable(T, [basepath_to_save filesep 'PosVM_Harrison-Kanji.xls'])

figure, boxplot([phase_by_area_cond1_vMpos; phase_by_area_cond2_vMpos],[area_id_vMpos area_id_vMpos])
set(gca, 'XTickLabel', {'VPL', 'dPul', 'MD'})
ylim([0 2*pi])
xlabel('Area')
ylabel('Heart-Cycle Phase [0-2\pi]')

filename = ['WhiskerBoxPlot_byAreaPos_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf);

figure, boxplot([phase_by_area_cond1_vMpos; phase_by_area_cond2_vMpos],condition_id_vMpos)
set(gca, 'XTickLabel', {'Rest', 'Task'})
ylim([0 2*pi])
xlabel('Condition')
ylabel('Heart-Cycle Phase [0-2\pi]')

filename = ['WhiskerBoxPlot_byConditionPos_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf);

%% 2-way ANOVA for neg. von Mises distributions
condition_id_vMneg = [zeros(1, length(phase_by_area_cond1_vMneg)) ones(1, length(phase_by_area_cond2_vMneg))];
[pval_vMneg, stats_vMneg] = circ_hktest([phase_by_area_cond1_vMneg; phase_by_area_cond2_vMneg], [area_id_vMneg area_id_vMneg], condition_id_vMneg, 1, {'Area: VPL, dPul, MD', 'Condition: Rest, Task', 'Area * Condition'});

T = table(stats_vMneg);
writetable(T, [basepath_to_save filesep 'NegVM_Harrison-Kanji.xls'])

figure, boxplot([phase_by_area_cond1_vMneg; phase_by_area_cond2_vMneg],[area_id_vMneg area_id_vMneg])
set(gca, 'XTickLabel', {'VPL', 'dPul', 'MD'})
ylim([0 2*pi])
xlabel('Area')
ylabel('Heart-Cycle Phase [0-2\pi]')

filename = ['WhiskerBoxPlot_byAreaNeg_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf);

figure, boxplot([phase_by_area_cond1_vMneg; phase_by_area_cond2_vMneg],condition_id_vMneg)
set(gca, 'XTickLabel', {'Rest', 'Task'})
ylim([0 2*pi])
xlabel('Condition')
ylabel('Heart-Cycle Phase [0-2\pi]')

filename = ['WhiskerBoxPlot_byConditionNeg_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf);

end

function plot_stability_vs_linear_gof(data, field_name, cfg, cond1_num, cond2_num, targName, file_prefix, basepath_to_save)

figure,
set(gcf, 'Position', [557   620   754   200*length(field_name)])

for lin_field_num = 1:length(field_name)
    
    Rsq_cond1{lin_field_num} = data.(cfg.condition(cond1_num).name).(field_name{lin_field_num}).rsquared;
    Rsq_cond2{lin_field_num} = data.(cfg.condition(cond2_num).name).(field_name{lin_field_num}).rsquared;
    
    [cc_cond1_tmp, p_cond1_tmp] = corrcoef(Rsq_cond1{lin_field_num}, data.criteria.stability_F);
    cc_cond1{lin_field_num}     = cc_cond1_tmp(1,2);
    p_cond1{lin_field_num}      = p_cond1_tmp(1,2);
    [cc_cond2_tmp, p_cond2_tmp] = corrcoef(Rsq_cond2{lin_field_num}, data.criteria.stability_V);
    cc_cond2{lin_field_num}     = cc_cond2_tmp(1,2);
    p_cond2{lin_field_num}      = p_cond2_tmp(1,2);
    clear cc_cond1_tmp p_cond1_tmp cc_cond2_tmp p_cond2_tmp
    
    s1(lin_field_num) = subplot(length(field_name),length(cfg.condition),1+2*(lin_field_num-1));
    scatter(log10(Rsq_cond1{lin_field_num}), data.criteria.stability_F, '.', 'MarkerEdgeColor', cfg.condition(cond1_num).color)
    xlabel('log10(R^2 of Linear Fit)')
    ylabel('Stability')
    title({[cfg.condition(cond1_num).name ': ' field_name{lin_field_num}],
        ['cc = ' num2str(cc_cond1{lin_field_num}) '; p = ' num2str(p_cond1{lin_field_num})]}, 'interpreter', 'none')
    box on
    
    s2(lin_field_num) = subplot(length(field_name),length(cfg.condition),2*lin_field_num);
    scatter(log10(Rsq_cond2{lin_field_num}), data.criteria.stability_V, '.', 'MarkerEdgeColor', cfg.condition(cond2_num).color)
    title({[cfg.condition(cond2_num).name ': ' field_name{lin_field_num}], ...
        ['cc = ' num2str(cc_cond2{lin_field_num}) '; p = ' num2str(p_cond2{lin_field_num})]}, 'interpreter', 'none')
    box on
    
end

linkaxes([s1 s2], 'xy')

sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf);

T = table({[cfg.condition(cond1_num).name ': ' field_name{1}], [cfg.condition(cond2_num).name]}', ...
    [cc_cond1; cc_cond2], [p_cond1; p_cond2], ...
    'VariableNames', {'Condition', 'Pearson''s Rho', 'Pearson''s p-value'});
writetable(T, [basepath_to_save filesep filename '.xls'])

end

function plot_RpeakNum_histograms(ephys_dt, var_prefix, cfg, cond1_num, cond2_num, targName, file_prefix, basepath_to_save)
bin_centers = 250:500:11750;

figure,
set(gcf, 'Position', [557   620   754   206])

M1 = nanmedian(ephys_dt.(cfg.condition(cond1_num).name).NrEvents);
M2 = nanmedian(ephys_dt.(cfg.condition(cond2_num).name).NrEvents);

s1 = subplot(1,2,1);
hist(ephys_dt.(cfg.condition(cond1_num).name).NrEvents, bin_centers);
title([cfg.condition(cond1_num).name '-' var_prefix ': median = ' num2str(M1) ' heart cycles'])
xlim([0 12000])
xlabel('Number of Heart Cycles')
ylabel('Number of Units')

s2 = subplot(1,2,2);
hist(ephys_dt.(cfg.condition(cond2_num).name).NrEvents, bin_centers);
xlim([0 12000])
title([cfg.condition(cond2_num).name '-' var_prefix ': median = ' num2str(M2) ' heart cycles'])

linkaxes([s1 s2], 'y')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf)

if ~isempty(rmmissing(ephys_dt.(cfg.condition(cond1_num).name).NrEvents)) && ~isempty(rmmissing(ephys_dt.(cfg.condition(cond2_num).name).NrEvents))
    p_wilcox = signrank(ephys_dt.(cfg.condition(cond1_num).name).NrEvents, ephys_dt.(cfg.condition(cond2_num).name).NrEvents);
else
    p_wilcox = NaN;
end

T = table(M1, M2, p_wilcox','VariableNames', {[cfg.condition(cond1_num).name ': median NrEvents'], [cfg.condition(cond2_num).name ': median NrEvents'], 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])

end

function plot_Rsquared_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save)

bin_edges              = 10 .^ (-8:0.2:0);

figure,
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C)

for field_num = 1:length(var_prefix)
    
    % put R-squared into separate cells
    cos_Rsq_cond1{field_num}   = dt.(cfg.condition(cond1_num).name).([var_prefix{field_num} 'cosine']).rsquared;
    vMPos_Rsq_cond1{field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{field_num} 'vonMisesPos']).rsquared;
    vMNeg_Rsq_cond1{field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{field_num} 'vonMisesNeg']).rsquared;
    
    cos_Rsq_cond2{field_num}   = dt.(cfg.condition(cond2_num).name).([var_prefix{field_num} 'cosine']).rsquared;
    vMPos_Rsq_cond2{field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{field_num} 'vonMisesPos']).rsquared;
    vMNeg_Rsq_cond2{field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{field_num} 'vonMisesNeg']).rsquared;
    
    % exclude negative R-squared
    cos_Rsq_cond1{field_num}   = cos_Rsq_cond1{field_num}(cos_Rsq_cond1{field_num}>0);
    vMPos_Rsq_cond1{field_num} = vMPos_Rsq_cond1{field_num}(vMPos_Rsq_cond1{field_num}>0);
    vMNeg_Rsq_cond1{field_num} = vMNeg_Rsq_cond1{field_num}(vMNeg_Rsq_cond1{field_num}>0);
    
    cos_Rsq_cond2{field_num}   = cos_Rsq_cond2{field_num}(cos_Rsq_cond2{field_num}>0);
    vMPos_Rsq_cond2{field_num} = vMPos_Rsq_cond2{field_num}(vMPos_Rsq_cond2{field_num}>0);
    vMNeg_Rsq_cond2{field_num} = vMNeg_Rsq_cond2{field_num}(vMNeg_Rsq_cond2{field_num}>0);
    
    % bin data
    cos_Rsq_cond1_binned{field_num}   = histc(cos_Rsq_cond1{field_num}, bin_edges);
    vMPos_Rsq_cond1_binned{field_num} = histc(vMPos_Rsq_cond1{field_num}, bin_edges);
    vMNeg_Rsq_cond1_binned{field_num} = histc(vMNeg_Rsq_cond1{field_num}, bin_edges);
    
    cos_Rsq_cond2_binned{field_num}   = histc(cos_Rsq_cond2{field_num}, bin_edges);
    vMPos_Rsq_cond2_binned{field_num} = histc(vMPos_Rsq_cond2{field_num}, bin_edges);
    vMNeg_Rsq_cond2_binned{field_num} = histc(vMNeg_Rsq_cond2{field_num}, bin_edges);
    
    % compute medians
    med_rest(:, field_num) = {median(cos_Rsq_cond1{field_num}()); median(vMPos_Rsq_cond1{field_num}); median(vMNeg_Rsq_cond1{field_num})};
    med_task(:, field_num) = {median(cos_Rsq_cond2{field_num}); median(vMPos_Rsq_cond2{field_num}); median(vMNeg_Rsq_cond2{field_num})};
    
    % plot histgrams
    s1(field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(field_num-1));
    stairs(bin_edges, [cos_Rsq_cond1_binned{field_num} vMPos_Rsq_cond1_binned{field_num} vMNeg_Rsq_cond1_binned{field_num}], 'LineWidth', 2)
    hold on
    plot(med_rest{1,field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(med_rest{2,field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(med_rest{3,field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    set(gca,'XScale','log')
    xlabel('log10(R^2)')
    ylabel('Number of Units')
    title(cfg.condition(cond1_num).name)
    legend('Cosine', 'von Mises Pos.', 'von Mises Neg.', 'Location', 'Best')
    
    s2(field_num) = subplot(length(var_prefix),length(cfg.condition),2*field_num);
    stairs(bin_edges, [cos_Rsq_cond2_binned{field_num} vMPos_Rsq_cond2_binned{field_num} vMNeg_Rsq_cond2_binned{field_num}], 'LineWidth', 2)
    hold on
    plot(med_task{1,field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(med_task{2,field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(med_task{3,field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    set(gca,'XScale','log')
    xlabel('log10(R^2)')
    title(cfg.condition(cond2_num).name)
    
    p_wilcox{1,field_num} = signrank(dt.(cfg.condition(cond1_num).name).([var_prefix{1} 'cosine']).rsquared, dt.(cfg.condition(cond2_num).name).([var_prefix{field_num} 'cosine']).rsquared);
    p_wilcox{2,field_num} = signrank(dt.(cfg.condition(cond1_num).name).([var_prefix{1} 'vonMisesPos']).rsquared, dt.(cfg.condition(cond2_num).name).([var_prefix{field_num} 'vonMisesPos']).rsquared);
    p_wilcox{3,field_num} = signrank(dt.(cfg.condition(cond1_num).name).([var_prefix{1} 'vonMisesNeg']).rsquared, dt.(cfg.condition(cond2_num).name).([var_prefix{field_num} 'vonMisesNeg']).rsquared);
    
end

linkaxes([s1 s2], 'xy')

sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close(gcf)

T = table({'Cosine', 'Pos. von Mises', 'Neg. von Mises'}', ...
    med_rest, med_task, p_wilcox, ...
    'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median R^2'], [cfg.condition(cond2_num).name ': median R^2'], 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])

end

function plot_phase_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save, selection)

phase_bin_edges = [0:0.2:2*pi];

f1 = figure;
set(f1, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C)

for lin_field_num = 1:length(var_prefix)
    
    cos_Phase_cond1{lin_field_num}   = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'cosine']).coefs(:,2);
    cos_Phase_cond2{lin_field_num}   = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'cosine']).coefs(:,2);
    
    cos_Phase_cond1_binned{lin_field_num}   = histc(cos_Phase_cond1{lin_field_num}, phase_bin_edges);
    cos_Phase_cond2_binned{lin_field_num}   = histc(cos_Phase_cond2{lin_field_num}, phase_bin_edges);
    
    if exist('selection', 'var')
        % select only significant units
        vMPos_Phase_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
        vMNeg_Phase_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
        
        vMPos_Phase_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
        vMNeg_Phase_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
    else
        vMPos_Phase_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
        vMNeg_Phase_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);
        
        vMPos_Phase_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
        vMNeg_Phase_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);
    end
    
    vMPos_Phase_cond1_binned{lin_field_num} = histc(vMPos_Phase_cond1{lin_field_num}, phase_bin_edges);
    vMNeg_Phase_cond1_binned{lin_field_num} = histc(vMNeg_Phase_cond1{lin_field_num}, phase_bin_edges);
    
    vMPos_Phase_cond2_binned{lin_field_num} = histc(vMPos_Phase_cond2{lin_field_num}, phase_bin_edges);
    vMNeg_Phase_cond2_binned{lin_field_num} = histc(vMNeg_Phase_cond2{lin_field_num}, phase_bin_edges);
    
    circ_med_rest(:,lin_field_num) = ...
        {mod(circ_median(cos_Phase_cond1{lin_field_num}), 2*pi);
        mod(circ_median(vMPos_Phase_cond1{lin_field_num}), 2*pi);
        mod(circ_median(vMNeg_Phase_cond1{lin_field_num}), 2*pi)};
    circ_med_task(:,lin_field_num) = ...
        {mod(circ_median(cos_Phase_cond2{lin_field_num}), 2*pi);
        mod(circ_median(vMPos_Phase_cond2{lin_field_num}), 2*pi);
        mod(circ_median(vMNeg_Phase_cond2{lin_field_num}), 2*pi)};
    
    cos_Phase_cond1_binned{lin_field_num} = cos_Phase_cond1_binned{lin_field_num}(:);
    
    figure(f1)
    
    s1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    stairs(phase_bin_edges, [cos_Phase_cond1_binned{lin_field_num} vMPos_Phase_cond1_binned{lin_field_num} vMNeg_Phase_cond1_binned{lin_field_num}], 'LineWidth', 2)
    hold on
    plot(circ_med_rest{1,lin_field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(circ_med_rest{2,lin_field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(circ_med_rest{3,lin_field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    hold on
    
    xlabel('Heart-Cycle Phase of Peak')
    title(cfg.condition(cond1_num).name)
    legend('Cosine', 'von Mises Pos.', 'von Mises Neg.')
    xlim([0 2*pi])
    
    s2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    stairs(phase_bin_edges, [cos_Phase_cond2_binned{lin_field_num} vMPos_Phase_cond2_binned{lin_field_num} vMNeg_Phase_cond2_binned{lin_field_num}], 'LineWidth', 2)
    hold on
    plot(circ_med_task{1,lin_field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(circ_med_task{2,lin_field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(circ_med_task{3,lin_field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    xlabel('Heart-Cycle Phase of Peak')
    title(cfg.condition(cond2_num).name)
    xlim([0 2*pi])
    
    p_wwtest{lin_field_num,1} = circ_wwtest(cos_Phase_cond1{lin_field_num}, cos_Phase_cond2{lin_field_num});
    p_wwtest{lin_field_num,2} = circ_wwtest(vMPos_Phase_cond1{lin_field_num}, vMPos_Phase_cond2{lin_field_num});
    p_wwtest{lin_field_num,3} = circ_wwtest(vMNeg_Phase_cond1{lin_field_num}, vMNeg_Phase_cond2{lin_field_num});
    
    if length(var_prefix) <= 1
        
        f2 = figure;
        set(f2, 'Position', [557   516   903   300*length(var_prefix)])
        colororder(C)
        
        subplot(1,3,1);
        scatter(cos_Phase_cond1{lin_field_num}, cos_Phase_cond2{lin_field_num}, [], 'MarkerFaceColor', C(1,:))
        box on
        xlabel([cfg.condition(cond1_num).name ': Cardiac Phase, radians'])
        ylabel([cfg.condition(cond2_num).name ': Cardiac Phase, radians'])
        sgtitle([targName ': Peak Phases'])
        axis square
        [circ_cc, circ_pp] = circ_corrcc(cos_Phase_cond1{lin_field_num}, cos_Phase_cond2{lin_field_num});
        title({'Cosine Fit: ', ['circ cc = ' num2str(circ_cc) '; circ pp = ' num2str(circ_pp)]})
        
        subplot(1,3,2)
        scatter(vMPos_Phase_cond1{lin_field_num}, vMPos_Phase_cond2{lin_field_num}, [], 'MarkerFaceColor', C(2,:))
        box on
        axis square
        [circ_cc, circ_pp] = circ_corrcc(vMPos_Phase_cond1{lin_field_num}, vMPos_Phase_cond2{lin_field_num});
        title({'Pos. von Mises Fit: ', ['circ cc = ' num2str(circ_cc) '; circ pp = ' num2str(circ_pp)]})
        
        subplot(1,3,3)
        scatter(vMNeg_Phase_cond1{lin_field_num}, vMNeg_Phase_cond2{lin_field_num}, [], 'MarkerFaceColor', C(3,:))
        box on
        axis square
        [circ_cc, circ_pp] = circ_corrcc(vMNeg_Phase_cond1{lin_field_num}, vMNeg_Phase_cond2{lin_field_num});
        title({'Neg. von Mises Fit: ', ['circ cc = ' num2str(circ_cc) '; circ pp = ' num2str(circ_pp)]})
        
        filename = [file_prefix targName '_CircScatter_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(f2, [basepath_to_save filesep filename], '-pdf');
        close(f2);
        
    end
    
end

linkaxes([s1 s2], 'xy')

sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(f1, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close(f1);

T = table({'Cosine', 'Pos. von Mises', 'Neg. von Mises'}', circ_med_rest, circ_med_task, p_wwtest', ...
    'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median phase'], [cfg.condition(cond2_num).name ': median phase'], 'Watson-Williams''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])
clear T

% compare median splits
if length(var_prefix) == 2
    
    % compare for rest
    pp_medSplit_rest(1) = circ_wwtest(cos_Phase_cond1{1}, cos_Phase_cond1{2});
    pp_medSplit_rest(2) = circ_wwtest(vMPos_Phase_cond1{1}, vMPos_Phase_cond1{2});
    pp_medSplit_rest(3) = circ_wwtest(vMNeg_Phase_cond1{1}, vMNeg_Phase_cond1{2});
    
    %compare for task
    pp_medSplit_task(1) = circ_wwtest(cos_Phase_cond2{1}, cos_Phase_cond2{2});
    pp_medSplit_task(2) = circ_wwtest(vMPos_Phase_cond2{1}, vMPos_Phase_cond2{2});
    pp_medSplit_task(3) = circ_wwtest(vMNeg_Phase_cond2{1}, vMNeg_Phase_cond2{2});
    
    T = table({'Cosine', 'Pos. von Mises', 'Neg. von Mises'}', pp_medSplit_rest', pp_medSplit_task', ...
        'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': WW-Test'], [cfg.condition(cond2_num).name ': WW-Test']});
    
    writetable(T, [basepath_to_save filesep 'CompareSplits_' filename '.xls'])
    
end

end

function plot_scalingFactor_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save, selection)

bin_centers_a    = -0.195:0.01:0.195;

bin_centers_diff = -0.145:0.01:0.145;

figure
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C(2:end,:))

for lin_field_num = 1:length(var_prefix)
    
    cos_a_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'cosine']).coefs(:,1);
    cos_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'cosine']).coefs(:,1);
    
    if exist('selection', 'var')
        % select only significant units
        vMPos_a_cond1{lin_field_num} = ...
            dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,1);
        vMNeg_a_cond1{lin_field_num} = ...
            dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,1);
        
        vMPos_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,1);
        vMNeg_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,1);
    else
        vMPos_a_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,1);
        vMNeg_a_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,1);
        
        vMPos_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,1);
        vMNeg_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,1);
    end
    
    % bin scaling factors
    vMPos_a_cond1_binned{lin_field_num} = histc(vMPos_a_cond1{lin_field_num}, bin_centers_a);
    vMNeg_a_cond1_binned{lin_field_num} = histc(vMNeg_a_cond1{lin_field_num}, bin_centers_a);
    
    vMPos_a_cond2_binned{lin_field_num} = histc(vMPos_a_cond2{lin_field_num}, bin_centers_a);
    vMNeg_a_cond2_binned{lin_field_num} = histc(vMNeg_a_cond2{lin_field_num}, bin_centers_a);
    
    % plot scaling factor histograms
    vMNeg_a_cond1_binned{lin_field_num} = vMNeg_a_cond1_binned{lin_field_num}(:);
    vMNeg_a_cond2_binned{lin_field_num} = vMNeg_a_cond2_binned{lin_field_num}(:);
    
    s1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    stairs(bin_centers_a, [vMPos_a_cond1_binned{lin_field_num} vMNeg_a_cond1_binned{lin_field_num}], 'LineWidth', 2)
    legend('von Mises Pos.', 'von Mises Neg.')
    title(cfg.condition(cond1_num).name)
    xlabel('Scaling Factor')
    ylabel('Unit Counts')
    
    s2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    stairs(bin_centers_a, [vMPos_a_cond2_binned{lin_field_num} vMNeg_a_cond2_binned{lin_field_num}], 'LineWidth', 2)
    title(cfg.condition(cond2_num).name)
    
    sgtitle([targName ': Scaling Factors'])
    
    if length(var_prefix) <= 1
        
        f2 = figure;
        set(f2, 'Position', [557   516   903   300*length(var_prefix)])
        colororder(C)
        
        subplot(1,3,1);
        scatter(cos_a_cond1{lin_field_num}, cos_a_cond2{lin_field_num}, [], 'MarkerFaceColor', C(1,:))
        box on
        xlabel([cfg.condition(cond1_num).name ': Scaling Factor, spk/bin/RR'])
        ylabel([cfg.condition(cond2_num).name ': Scaling Factor, spk/bin/RR'])
        sgtitle([targName ': Scaling Factors'])
        [cc, pp] = corrcoef(cos_a_cond1{lin_field_num}, cos_a_cond2{lin_field_num});
        [p,~] = ... % p(1) is linear slope, p(2) is b-member
            polyfit(cos_a_cond1{lin_field_num}, cos_a_cond2{lin_field_num},1);
        title({'Cosine Fit: ', ...
            ['circ cc = ' num2str(cc(2,1)) '; circ pp = ' num2str(pp(2,1))], ...
            ['lin.fit: slope = ' num2str(p(1))]})
        
        subplot(1,3,2)
        scatter(vMPos_a_cond1{lin_field_num}, vMPos_a_cond2{lin_field_num}, [], 'MarkerFaceColor', C(2,:))
        box on
        [cc, pp] = corrcoef(vMPos_a_cond1{lin_field_num}, vMPos_a_cond2{lin_field_num});
        [p,~] = ... % p(1) is linear slope, p(2) is b-member
            polyfit(vMPos_a_cond1{lin_field_num}, vMPos_a_cond2{lin_field_num},1);
        title({'Pos. von Mises Fit: ', ...
            ['circ cc = ' num2str(cc(2,1)) '; circ pp = ' num2str(pp(2,1))], ...
            ['lin.fit: slope = ' num2str(p(1))]})
        
        subplot(1,3,3)
        scatter(vMNeg_a_cond1{lin_field_num}, vMNeg_a_cond2{lin_field_num}, [], 'MarkerFaceColor', C(3,:))
        box on
        [cc, pp] = corrcoef(vMNeg_a_cond1{lin_field_num}, vMNeg_a_cond2{lin_field_num});
        [p,~] = ... % p(1) is linear slope, p(2) is b-member
            polyfit(vMNeg_a_cond1{lin_field_num}, vMNeg_a_cond2{lin_field_num},1);
        if size(cc,1) == 2
            title({'Neg. von Mises Fit: ', ...
                ['circ cc = ' num2str(cc(2,1)) '; circ pp = ' num2str(pp(2,1))], ...
                ['lin.fit: slope = ' num2str(p(1))]})
        end
        
        filename = [file_prefix targName '_CircScatter_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(f2, [basepath_to_save filesep filename], '-pdf');
        close(f2);
        
    end
    
end

linkaxes([s1 s2], 'y')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close(gcf);

figure
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])

for lin_field_num = 1:length(var_prefix)
    
    % compute differences of absolute scaling factors
    vMPos_a_absdiff = abs(vMPos_a_cond1{lin_field_num}) - abs(vMPos_a_cond2{lin_field_num});
    
    vMNeg_a_absdiff = abs(vMNeg_a_cond1{lin_field_num}) - abs(vMNeg_a_cond2{lin_field_num});
    
    % median difference in scaling factors
    M_pos{lin_field_num} = median(vMPos_a_absdiff);
    M_neg{lin_field_num} = median(vMNeg_a_absdiff);
    
    % binned scaling factor differences
    vMPos_a_absdiff_binned = histc(vMPos_a_absdiff, bin_centers_diff);
    vMNeg_a_absdiff_binned = histc(vMNeg_a_absdiff, bin_centers_diff);
    
    S1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    stairs(bin_centers_diff, vMPos_a_absdiff_binned, 'LineWidth', 2)
    hold on
    plot(M_pos{lin_field_num}, 0, '.', 'MarkerSize', 15)
    title('Pos. von Mises: Scaling Factor Difference')
    xlabel('Scaling Factor (abs(Rest) - abs(Task))')
    
    S2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    stairs(bin_centers_diff, vMNeg_a_absdiff_binned, 'LineWidth', 2)
    hold on
    plot(M_neg{lin_field_num}, 0, '.', 'MarkerSize', 15)
    title('Neg. von Mises: Scaling Factor Difference')
    
    sgtitle([targName ': Scaling Factor Differences'])
    
    % check if means are zeros
    if ~isempty(vMPos_a_cond1{lin_field_num}) && ~isempty(vMPos_a_cond2{lin_field_num})
        p_wilcoxon_paired{1,lin_field_num} = signrank(abs(vMPos_a_cond1{lin_field_num}),abs(vMPos_a_cond2{lin_field_num}));
    else
        p_wilcoxon_paired{1,lin_field_num} = NaN;
    end
    
    if ~isempty(vMNeg_a_cond1{lin_field_num}) && ~isempty(vMNeg_a_cond2{lin_field_num})
        p_wilcoxon_paired{2,lin_field_num} = signrank(abs(vMNeg_a_cond1{lin_field_num}),abs(vMNeg_a_cond2{lin_field_num}));
    else
        p_wilcoxon_paired{2,lin_field_num} = NaN;
    end
    
end

linkaxes([S1 S2], 'xy')

if strfind(file_prefix, 'MedSplit_')
    filename = [file_prefix(1:9) 'Diff' file_prefix(10:end) targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
else
    filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
end
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf);

T = table(M_pos', M_neg', p_wilcoxon_paired', ...
    'VariableNames', {'vonMisesPos Scaling Factor Median Diff.', 'vonMisesNeg Scaling Factor Median Diff.', 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])

% compare median splits
if length(var_prefix) == 2
    
    % compare for rest
    pp_medSplit_rest(1) = signrank(cos_a_cond1{1}, cos_a_cond1{2});
    pp_medSplit_rest(2) = signrank(vMPos_a_cond1{1}, vMPos_a_cond1{2});
    pp_medSplit_rest(3) = signrank(vMNeg_a_cond1{1}, vMNeg_a_cond1{2});
    
    %compare for task
    pp_medSplit_task(1) = signrank(cos_a_cond2{1}, cos_a_cond2{2});
    pp_medSplit_task(2) = signrank(vMPos_a_cond2{1}, vMPos_a_cond2{2});
    pp_medSplit_task(3) = signrank(vMNeg_a_cond2{1}, vMNeg_a_cond2{2});
    
    T = table({'Cosine', 'Pos. von Mises', 'Neg. von Mises'}', pp_medSplit_rest', pp_medSplit_task', ...
        'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': Wilcox-Test'], [cfg.condition(cond2_num).name ': Wilcox-Test']});
    
    writetable(T, [basepath_to_save filesep 'CompareSplits_' filename '.xls'])
    
end

end

function plot_kappa_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save, selection)

kappa_bin_edges = exp(-1:2);

figure
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C(2:end,:))

for lin_field_num = 1:length(var_prefix)
    
    if exist('selection', 'var')
        % select only significant units
        vMPos_kappa_cond1{lin_field_num} = ...
            log(dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_kappa_cond1{lin_field_num} = ...
            log(dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
        vMPos_kappa_cond2{lin_field_num} = ...
            log(dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_kappa_cond2{lin_field_num} = ...
            log(dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
    else
        vMPos_kappa_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3);
        vMNeg_kappa_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3);
        
        vMPos_kappa_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3);
        vMNeg_kappa_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3);
    end
    
    vMPos_kappa_cond1_binned{lin_field_num} = histc(vMPos_kappa_cond1{lin_field_num}, kappa_bin_edges);
    vMNeg_kappa_cond1_binned{lin_field_num} = histc(vMNeg_kappa_cond1{lin_field_num}, kappa_bin_edges);
    
    vMPos_kappa_cond2_binned{lin_field_num} = histc(vMPos_kappa_cond2{lin_field_num}, kappa_bin_edges);
    vMNeg_kappa_cond2_binned{lin_field_num} = histc(vMNeg_kappa_cond2{lin_field_num}, kappa_bin_edges);
    
    M_cond1(lin_field_num,:) = ...
        {median(vMPos_kappa_cond1{lin_field_num}); median(vMNeg_kappa_cond1{lin_field_num})};
    M_cond2(lin_field_num,:) = ...
        {median(vMPos_kappa_cond2{lin_field_num}); median(vMNeg_kappa_cond2{lin_field_num})};
    
    vMNeg_kappa_cond1_binned{lin_field_num} = ...
        vMNeg_kappa_cond1_binned{lin_field_num}(:);
    vMNeg_kappa_cond2_binned{lin_field_num} = ...
        vMNeg_kappa_cond2_binned{lin_field_num}(:);
    
    s1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    bar(kappa_bin_edges,[vMPos_kappa_cond1_binned{lin_field_num} vMNeg_kappa_cond1_binned{lin_field_num}])
    xlabel('log(\kappa)')
    title(cfg.condition(cond1_num).name)
    legend({'von Mises Pos.', 'von Mises Neg.'})
    xlim([-2 4])
    
    s2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    bar(kappa_bin_edges,[vMPos_kappa_cond2_binned{lin_field_num} vMNeg_kappa_cond2_binned{lin_field_num}])
    xlabel('log(\kappa)')
    title(cfg.condition(cond2_num).name)
    xlim([-2 4])
    
    p_wilcoxon_paired{lin_field_num,1} = signrank(vMPos_kappa_cond1{lin_field_num}, vMPos_kappa_cond2{lin_field_num});
    
    if ~isempty(vMNeg_kappa_cond1{lin_field_num}) && ~isempty(vMNeg_kappa_cond2{lin_field_num})
        p_wilcoxon_paired{2,lin_field_num} = signrank(abs(vMNeg_kappa_cond1{lin_field_num}),abs(vMNeg_kappa_cond2{lin_field_num}));
    else
        p_wilcoxon_paired{2,lin_field_num} = NaN;
    end
    
end

linkaxes([s1 s2], 'xy')

sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close all;

T = table({'Pos. von Mises', 'Neg. von Mises'}', M_cond1', M_cond2', p_wilcoxon_paired, ...
    'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median \kappa'], [cfg.condition(cond2_num).name ': median \kappa'], 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])

% compare median splits
if length(var_prefix) == 2
    
    % compare for rest
    pp_medSplit_rest(1) = signrank(vMPos_kappa_cond1{1}, vMPos_kappa_cond1{2});
    pp_medSplit_rest(2) = signrank(vMNeg_kappa_cond1{1}, vMNeg_kappa_cond1{2});
    
    %compare for task
    pp_medSplit_task(1) = signrank(vMPos_kappa_cond2{1}, vMPos_kappa_cond2{2});
    pp_medSplit_task(2) = signrank(vMNeg_kappa_cond2{1}, vMNeg_kappa_cond2{2});
    
    T = table({'Pos. von Mises', 'Neg. von Mises'}', pp_medSplit_rest', pp_medSplit_task', ...
        'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': Wilcox-Test'], [cfg.condition(cond2_num).name ': Wilcox-Test']});
    
    writetable(T, [basepath_to_save filesep 'CompareSplits_' filename '.xls'])
    
end

end

function plot_kappas_vs_phase(data, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save, selection)

for lin_field_num = 1:length(var_prefix)
    
    if exist('selection', 'var')
        % take kappas
        vMPos_kappa_cond1 = log(dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_kappa_cond1 = log(dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
        vMPos_kappa_cond2 = log(dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_kappa_cond2 = log(dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
        % take phases
        vMPos_Phase_cond1 = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
        vMNeg_Phase_cond1 = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
        
        vMPos_Phase_cond2 = dt.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
        vMNeg_Phase_cond2 = dt.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
    else
        % take kappas
        vMPos_kappa_cond1 = log(data.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3));
        vMNeg_kappa_cond1 = log(data.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3));
        
        vMPos_kappa_cond2 = log(data.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3));
        vMNeg_kappa_cond2 = log(data.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3));
        
        % take phases
        vMPos_Phase_cond1 = data.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
        vMNeg_Phase_cond1 = data.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);
        
        vMPos_Phase_cond2 = data.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
        vMNeg_Phase_cond2 = data.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);
    end
    
    figure,
    set(gcf, 'Position', [557   620   754   206])
    colororder(C(2:end,:))
    
    subplot(1,2,1)
    scatter(vMPos_kappa_cond1,vMPos_Phase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
    hold on
    scatter(vMPos_kappa_cond2,vMPos_Phase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
    box on
    ylim([0 2*pi])
    title('Pos. von Mises')
    xlabel('log(\kappa)')
    ylabel('Peak Phase [0-2\pi]')
    legend({cfg.condition.name}, 'Location', 'Best')
    
    subplot(1,2,2)
    scatter(vMNeg_kappa_cond1,vMNeg_Phase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
    hold on
    scatter(vMNeg_kappa_cond2,vMNeg_Phase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
    box on
    ylim([0 2*pi])
    title('Neg. von Mises')
    
    filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(gcf)
    
end

end

function plot_Rsquared_scatter(a,b)

scatter(a, b)

xlim([-7 -2])
ylim([-7 -2])
axis square
box on
hold on

nan_ids1 = isnan(a);
nan_ids2 = isnan(b);
not_nan  = ~nan_ids1 & ~nan_ids2;

[cc, pp] = corrcoef(a, b,'Row','pairwise');
p = polyfit(a(not_nan),b(not_nan),1);

mdl = fitlm(a(not_nan),b(not_nan))

x_lim = get(gca,'XLim');
plot(x_lim, p(1)*x_lim + p(2))
xlabel('Rest: log_1_0(R^2)')
ylabel('Task: log_1_0(R^2)')
title({['cc = ' num2str(cc(2,1)) '; p = ' num2str(pp(2,1))], ...
    ['y = ' num2str(p(1)) 'x + ' num2str(p(2))]})

end

