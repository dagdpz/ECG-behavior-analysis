function ecg_bna_population_cardioballistic(cfg, data_folder, unitList, basepath_to_save)

if ~exist(basepath_to_save, 'dir')
    mkdir(basepath_to_save)
end

%% load data
output_folder = [cfg.SPK_root_results_fldr filesep basepath_to_save];
load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG\' unitList '.mat'], 'unit_ids', 'targets', 'ids_both')

var_list = ...
    {'histogram_MI', 'histogram_p', 'histogram_phase', ...
    'AMP_MI', 'HW_MI', 'TPW_MI', 'REP_MI', ...
    'AMP_max_consec_bins', 'HW_max_consec_bins', 'TPW_max_consec_bins', 'REP_max_consec_bins', ...
    'AMP_modulation_index', 'HW_modulation_index', 'TPW_modulation_index', 'REP_modulation_index', ...
    'cc_PSTH_feature', 'pp_PSTH_feature', 'distance2thr'};

% TargetBrainArea = {SPK_cardioballistic.target};
% 
% Ana_TargetBrainArea = cfg.targets;
% if cfg.combine_hemispheres
%     TargetBrainArea=cellfun(@(x) x(1:end-2),TargetBrainArea,'UniformOutput',false);
%     Ana_TargetBrainArea=cellfun(@(x) x(1:end-2),Ana_TargetBrainArea,'UniformOutput',false);
% end
Ana_TargetBrainArea = {'VPL', 'dPul', 'MD'}; % this is temporary

N_Areas=numel(Ana_TargetBrainArea);
N_conditions=numel(cfg.condition);

% %% load data
% 
% for a = 1: N_Areas
%     
%     T=unqTargets{a};
%     currTargIds = cellfun(@(x) strcmp(x, unqTargets{a}), targets);
%     curr_unit_ids  = unit_ids(currTargIds);
%     curr_ids_both = ids_both(currTargIds);
%     Out.(T) = ecg_bna_load_variables(cfg, curr_unit_ids, data_folder, 'Output', var_list, curr_ids_both);
%     
% end
% 
% 
% for a = 1: N_Areas
%     T=Ana_TargetBrainArea{a};
%     for c=1:N_conditions
%         L=cfg.condition(c).name;
% %         dat=[SPK_cardioballistic(ismember(TargetBrainArea,T)).(L)];
% %         dat_fieldnames=fieldnames(dat);
%         dat_fieldnames = ...
%             {'histogram_MI', 'histogram_p', 'histogram_phase', ...
%             'AMP_MI', 'HW_MI', 'TPW_MI', 'REP_MI', ...
%             'AMP_max_consec_bins', 'HW_max_consec_bins', 'TPW_max_consec_bins', 'REP_max_consec_bins', ...
%             'AMP_modulation_index', 'HW_modulation_index', 'TPW_modulation_index', 'REP_modulation_index', ...
%             'cc_PSTH_feature', 'pp_PSTH_feature', 'distance2thr'};
%         
%         for fn=1:numel(dat_fieldnames)
%             N=dat_fieldnames{fn};
%             tmp = arrayfun(@(x) x.(L).(N), SPK_cardioballistic(ismember(TargetBrainArea,T)), 'UniformOutput', false);
%             if strcmp(N, 'AMP_MI') || strcmp(N, 'HW_MI') || strcmp(N, 'TPW_MI') || strcmp(N, 'REP_MI') || ...
%                 strcmp(N, 'cc_PSTH_feature') || strcmp(N, 'pp_PSTH_feature')
%                 Out.(T).(L).(N) = vertcat(tmp{:});
% %             elseif strcmp(N, 'pearson_r') || strcmp(N, 'permuted_p')
%             elseif strcmp(N, 'pearson_r') || strcmp(N, 'pearson_p')
%                 Out.(T).(L).(N) = [tmp{:}];
%             else
%                 Out.(T).(L).(N) = horzcat(tmp{:});
%             end
%         end
%         
%         dat_fieldnames={'unitId','target','channel','unit','quantSNR','Single_rating','stability_rating','FR'};
%         for fn=1:numel(dat_fieldnames)
%             N=dat_fieldnames{fn};
%             if strcmp(N, 'unitId')
%                 Out.(T).(L).(N)=vertcat({SPK_cardioballistic(ismember(TargetBrainArea,T)).(N)}');
%             else
%                 Out.(T).(L).(N)=vertcat(SPK_cardioballistic(ismember(TargetBrainArea,T)).(N));
%             end
%         end
%     end
% end

for a = 1: N_Areas
    
    T=Ana_TargetBrainArea{a};
    currTargIds = cellfun(@(x) strcmp(x, Ana_TargetBrainArea{a}), targets);
    curr_unit_ids  = unit_ids(currTargIds);
    curr_ids_both = ids_both(currTargIds);
    Out.(T) = ecg_bna_load_variables(cfg, curr_unit_ids, data_folder, 'data', var_list, curr_ids_both);
    
end

%% number of 1 microVolt bins between the threshold and the first bin with at least 1 spike
figure,
set(gcf, 'Position', [1 471 1917 420])
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        counts = histc(Out.(T).(L).distance2thr, 0:20);
        subplot(1,N_Areas,a)
        stairs(0:20, counts, 'Color', cfg.condition(c).color, 'LineWidth', 2)
        hold on
        
    end
    title(T)
    if a == 1
        xlabel('Number of 1-\muV Bins Between the Threshold and the 1st Bin with Spikes')
        ylabel('Unit Count')
        legend({cfg.condition.name})
    end
end

save_figure_as(gcf,'Histograms_Number_of_Free_Bins',basepath_to_save,1)

%% number of consecutive bins that exceed 95% threshold
figure,
set(gcf, 'Position', [849 42 930 954])
field_names_mosher = ...
    {'AMP_max_consec_bins', 'HW_max_consec_bins', 'TPW_max_consec_bins', 'REP_max_consec_bins'};
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for feature_nums = 1:length(field_names_mosher)
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            counts = histc(Out.(T).(L).(field_names_mosher{feature_nums}), 0:6);
            subplot(length(field_names_mosher),N_Areas, (feature_nums-1)*N_Areas + a)
            stairs(0:6, counts, 'Color', cfg.condition(c).color, 'LineWidth', 2)
            
            if c == 1
                hold on
            else
                hold off
            end
            
            if feature_nums == 1
                title(T)
            end
            
            if a == 1
                ylabel(field_names_mosher{feature_nums}, 'Interpreter', 'none')
            end
            
            if a == 1 && feature_nums == 1
                xlabel('Number of Consecutive Bins Exceeding 95%')
            end
            
            if feature_nums == 1 && a == 1 && c == 2
                legend({cfg.condition.name})
            end
            
        end
    end
end

save_figure_as(gcf,'Histograms_Number_of_Consecutive_Bins',basepath_to_save,1)

%% correlation coefficients between feature phase dynamics and phase PSTH
figure,
set(gcf, 'Position', [849 42 930 954])
feature_names = {'AMP', 'HW', 'TPW', 'REP'};
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for feature_nums = 1:4
        subplot(length(feature_names), N_Areas, (feature_nums-1)*N_Areas + a)
        for c=1:N_conditions
            
            L=cfg.condition(c).name;
            
            if isempty(Out.(T).(L).cc_PSTH_feature)
                continue
            end
            counts = histc(Out.(T).(L).cc_PSTH_feature(feature_nums,:), -1:0.05:1);
            stairs(-1:0.05:1, counts, 'Color', cfg.condition(c).color, 'LineWidth', 2)
            
            if c == 1
                hold on
            else
                hold off
            end
            
            if feature_nums == 1
                title(T)
            end
            
            if a == 1
                ylabel(feature_names{feature_nums}, 'Interpreter', 'none')
            end
            
            if feature_nums == 1 && a == 1 && c == 2
                legend({cfg.condition.name})
                xlabel('CC between phase PSTH and feature phase dymanic')
            end
            
        end
    end
end

save_figure_as(gcf,'Histograms_CC_Phase_PSTH_&_Feature_Dynamic',basepath_to_save,1)

%% correlation coefficients between feature phase dynamic and phase PSTH - distinguished by significance
bar_colors = [0.8500 0.3250 0.0980; 0 0.4470 0.7410; 1 1 1];
feature_names = {'AMP', 'HW', 'TPW', 'REP'};
for feature_nums = 1:4
    figure,
    set(gcf, 'Position', [667 519 930 477])
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            if isempty(Out.(T).(L).cc_PSTH_feature)
                continue
            end
            
            curr_cc = Out.(T).(L).cc_PSTH_feature(feature_nums,:);
            curr_pp = Out.(T).(L).pp_PSTH_feature(feature_nums,:);
            
            if length(curr_cc) < 3
                continue
            end
            
            [cor_p, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
            pos_idx = curr_cc > 0;
            neg_idx = curr_cc < 0;
            
            % prepare data for plotting
            counts_pos    = histc(curr_cc(sig_idx & pos_idx)', -1:0.05:1); % significant and positive
            counts_neg    = histc(curr_cc(sig_idx & neg_idx)', -1:0.05:1); % significant and negative
            counts_nonsig = histc(curr_cc(~sig_idx)', -1:0.05:1);
            
            % plot
            plot_data = [counts_pos(:) counts_neg(:) counts_nonsig(:)];
            subplot(N_conditions, N_Areas, (c-1)*N_Areas + a)
            b = bar(-1:0.05:1, plot_data, 'stacked');
            set(b, 'FaceColor', 'Flat')
            set(b, {'CData'}, {bar_colors(1,:); bar_colors(2,:); bar_colors(3,:)})
            title([T ': ' L])
            
            if a == 1 && c == 1
                xlabel('CC between phase PSTH and feature phase dymanic')
                legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'})
            end
        end
            
    end
    sgtitle(feature_names{feature_nums})
    save_figure_as(gcf,['Histograms_CC_Phase_PSTH_&_' feature_names{feature_nums} '_Dynamic'], basepath_to_save, 1)
end

%% Mosher's features
figure,
set(gcf, 'Position', [849 42 930 954])
field_names_mosher = ...
    {'AMP_MI', 'HW_MI', 'TPW_MI', 'REP_MI'};
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for feature_nums = 1:length(field_names_mosher)
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            if isempty(Out.(T).(L).(field_names_mosher{feature_nums}))
                continue
            end
            
            counts = histc(100*Out.(T).(L).(field_names_mosher{feature_nums})(:, 1), 0:100); % modulation indices are stored in the 1st line; in % now
            subplot(length(field_names_mosher),N_Areas, (feature_nums-1)*N_Areas + a)
            stairs(0:100, counts, 'Color', cfg.condition(c).color, 'LineWidth', 2)
            
            xlim([0 20])
            
            if c == 1
                hold on
            else
                hold off
            end
            
            if feature_nums == 1
                title(T)
            end
            
            if a == 1
                ylabel(field_names_mosher{feature_nums}, 'Interpreter', 'none')
            end
            
            if feature_nums == 1 && a == 1 && c == 2
                xlabel('Mosher''s Modulation Index, %')
                legend({cfg.condition.name})
            end
            
        end
    end
end

save_figure_as(gcf,'Histograms_Mosher_Modulation_Indices', basepath_to_save, 1)

%% Mosher's features - significance
field_names_mosher = ...
    {'AMP_MI', 'HW_MI', 'TPW_MI', 'REP_MI'};
for feature_nums = 1:length(field_names_mosher)
    figure,
    set(gcf, 'Position', [667 519 930 477])
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            curr_ff = 100*Out.(T).(L).(field_names_mosher{feature_nums})(:, 1); % in % now
            curr_pp = Out.(T).(L).(field_names_mosher{feature_nums})(:, 2);
            
            if length(curr_ff) < 3
                continue
            end
            
            [cor_p, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
            pos_idx = curr_ff > 0;
            
            % prepare data for plotting
            counts_pos    = histc(curr_ff(sig_idx & pos_idx)', 0:100); % significant and positive
            counts_nonsig = histc(curr_ff(~sig_idx)', 0:100);
            
            % plot
            plot_data = [counts_pos; counts_nonsig];
            subplot(N_conditions, N_Areas, (c-1)*N_Areas + a)
            b = bar(0:100, plot_data, 'stacked');
            set(b, 'FaceColor', 'Flat')
            set(b, {'CData'}, {[cfg.condition(c).color]; [1 1 1]})
            
            xlim([0 20])
            
            title([T ': ' L])
            
            if a == 1
                legend({'Sig.', 'Non-Sig.'})
            end
            
            if a == 1 && c == 1
                xlabel('Mosher''s Motion Index, %')
                ylabel('Unit Counts')
            end
        end
    end
    sgtitle(feature_names{feature_nums})
    save_figure_as(gcf,['Histograms_' feature_names{feature_nums} '_Distribution'], basepath_to_save, 1)
end

%% Mosher's modulation phase
figure,
set(gcf, 'Position', [849 42 930 954])
field_names_mosher = ...
    {'AMP_MI', 'HW_MI', 'TPW_MI', 'REP_MI'};
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for feature_nums = 1:length(field_names_mosher)
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            curr_ph = Out.(T).(L).(field_names_mosher{feature_nums})(:, 3); % modulation phases are stored in the 3rd line
            curr_pp = Out.(T).(L).(field_names_mosher{feature_nums})(:, 2);
            
            if length(curr_ph) < 3
                continue
            end
            
            [cor_p, sig_idx] = bonf_holm(curr_pp, 0.05);
            
            curr_ph = curr_ph(sig_idx); % only significant
            counts = histc(curr_ph, 0:0.2:2*pi); 
            subplot(length(field_names_mosher),N_Areas, (feature_nums-1)*N_Areas + a)
            stairs(0:0.2:2*pi, counts, 'Color', cfg.condition(c).color);
            hold on
            l(c) = line(0:0.2:2*pi, smooth(counts, 'rlowess', 1)', 'LineWidth', 2, 'Color', cfg.condition(c).color);
            
            if c == 1
                hold on
            else
                hold off
            end
            
            if feature_nums == 1
                title(T)
            end
            
            if a == 1
                ylabel(field_names_mosher{feature_nums}, 'Interpreter', 'none')
            end
            
            if feature_nums == 1 && a == 1 && c == 2
                xlabel('Mosher''s Modulation Phase (0-2\pi)')
                legend(l, {cfg.condition.name})
            end
            
        end
    end
end
save_figure_as(gcf,'Histograms_Mosher_Modulation_Phases_Significant', basepath_to_save, 1)

%% Mosher's phases - significance
field_names_mosher = ...
    {'AMP_MI', 'HW_MI', 'TPW_MI', 'REP_MI'};
for feature_nums = 1:length(field_names_mosher)
    figure,
    set(gcf, 'Position', [667 519 930 477])
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            curr_ff = Out.(T).(L).(field_names_mosher{feature_nums})(:, 3);
            curr_pp = Out.(T).(L).(field_names_mosher{feature_nums})(:, 2);
            
            if length(curr_ff) < 3
                continue
            end
            
            [cor_p, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
            pos_idx = curr_ff > 0;
            
            % prepare data for plotting
            counts_pos    = histc(curr_ff(sig_idx & pos_idx)', 0:0.2:2*pi); % significant and positive
            counts_nonsig = histc(curr_ff(~sig_idx)', 0:0.2:2*pi);
            
            % plot
            plot_data = [counts_pos; counts_nonsig];
            subplot(N_conditions, N_Areas, (c-1)*N_Areas + a)
            b = bar(0:0.2:2*pi, plot_data, 'stacked');
            set(b, 'FaceColor', 'Flat')
            set(b, {'CData'}, {[cfg.condition(c).color]; [1 1 1]})
            
            title([T ': ' L])
            
            if a == 1
                legend({'Sig.', 'Non-Sig.'})
            end
            
            if a == 1 && c == 1
                xlabel('Modulation Index')
                ylabel('Unit Counts')
            end
        end
    end
    sgtitle(feature_names{feature_nums})
    save_figure_as(gcf,['Histograms_' feature_names{feature_nums} '_Phase_Distribution'], basepath_to_save, 1)
end

%% modulation index for features: (max - min)/mean_reshuffled
figure,
set(gcf, 'Position', [849 42 930 954])
field_names_mosher = ...
    {'AMP_modulation_index', 'HW_modulation_index', 'TPW_modulation_index', 'REP_modulation_index'};
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for feature_nums = 1:length(field_names_mosher)
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            if isempty(Out.(T).(L).(field_names_mosher{feature_nums}))
                continue
            end
            
            counts = histc(Out.(T).(L).(field_names_mosher{feature_nums}), 0:0.05:1);
            subplot(length(field_names_mosher),N_Areas, (feature_nums-1)*N_Areas + a)
            stairs(0:0.05:1, counts, 'Color', cfg.condition(c).color, 'LineWidth', 2)
            
            if c == 1
                hold on
            else
                hold off
            end
            
            if feature_nums == 1
                title(T)
            end
            
            if a == 1
                ylabel(field_names_mosher{feature_nums}, 'Interpreter', 'none')
            end
            
            if feature_nums == 1 && a == 1 && c == 2
                xlabel('Modulation Index: (max - min)/mean reshuffled')
                legend({cfg.condition.name})
            end
            
        end
    end
end
save_figure_as(gcf,'Histograms_DAG_Modulation_Indices', basepath_to_save, 1)

%% Mosher's modulation index vs. Igor's modulation index
figure,
set(gcf, 'Position', [849 42 930 954])
field_names_mosher = ...
    {'AMP_modulation_index', 'HW_modulation_index', 'TPW_modulation_index', 'REP_modulation_index'};
field_names_dag    = ...
    {'AMP_MI', 'HW_MI', 'TPW_MI', 'REP_MI'};
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for feature_nums = 1:length(field_names_mosher)
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            if isempty(Out.(T).(L).(field_names_mosher{feature_nums}))
                continue
            end
            
            curr_mosher = 100*Out.(T).(L).(field_names_mosher{feature_nums}); % in % now
            curr_dag    = Out.(T).(L).(field_names_dag{feature_nums})(:,1)';
            subplot(length(field_names_mosher),N_Areas, (feature_nums-1)*N_Areas + a)
            scatter(curr_mosher, curr_dag, [], cfg.condition(c).color)
            
            box on
            
            if c == 1
                hold on
            else
                hold off
            end
            
            if feature_nums == 1
                title(T)
            end
            
            if a == 1
                ylabel({'Modulation Index:', '(max - min)/mean reshuffled'})
            end
            
            if feature_nums == 1 && a == 1 && c == 2
                xlabel('Mosher''s Modulation Index, %')
                legend({cfg.condition.name})
            end
        end
    end
end
save_figure_as(gcf,'Scatters_Modulation_Indices', basepath_to_save, 1)

%% plot PSTH modulation phases
figure,
set(gcf, 'Position', [849 42 930 954])
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        curr_ff = Out.(T).(L).histogram_phase;
        curr_pp = Out.(T).(L).histogram_p;
        
        if numel(curr_ff) < 3
            continue
        end
        
        [cor_p, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
        
        % save data variables to use in the later analysis
        save([basepath_to_save filesep 'Phase_' T '_' L '.mat'], 'curr_ff', 'curr_pp')
        pos_idx = curr_ff > 0;
        
        % prepare data for plotting
        curr_ff_sig   = curr_ff(sig_idx & pos_idx)';
        counts_pos    = histc(curr_ff_sig, 0:0.2:2*pi); % significant and positive
        counts_nonsig = histc(curr_ff(~sig_idx)', 0:0.2:2*pi);
        
        counts_nonsig = counts_nonsig(:);
        counts_pos    = counts_pos(:);
        
        % plot
        plot_data = [counts_pos counts_nonsig];
        subplot(N_conditions, N_Areas, (c-1)*N_Areas + a)
        b = bar(0:0.2:2*pi, plot_data, 'stacked');
        set(b, 'FaceColor', 'Flat')
        set(b, {'CData'}, {[cfg.condition(c).color]; [1 1 1]})
        xlim([0 2*pi])
        
        if a == 1 && c == 1
            xlabel('Phase of the Heart Cycle, radians')
            ylabel('Unit Counts')
        end
        
        title([L ':' T])
    end
end
save_figure_as(gcf,'Histograms_PSTH_Phases', basepath_to_save, 1)

%% scatters for PSTH phases
figure,
set(gcf, 'Position', [63 525 1785 421])
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    
    task_ff = Out.(T).Task.histogram_phase;
    rest_ff = Out.(T).Rest.histogram_phase;
    
    task_p = Out.(T).Task.histogram_p;
    rest_p = Out.(T).Rest.histogram_p;
    
    [~, task_sig_idx] = bonf_holm(task_p, 0.05); % correct for multiple comparisons
    [~, rest_sig_idx] = bonf_holm(rest_p, 0.05); % correct for multiple comparisons
    
    [rho(a), pval(a)] = circ_corrcc(rest_ff(rest_sig_idx & task_sig_idx), task_ff(rest_sig_idx & task_sig_idx));
    
    subplot(1, N_Areas, a)
    scatter(rest_ff(rest_sig_idx & task_sig_idx), task_ff(rest_sig_idx & task_sig_idx), [], [0 0 0])
    xlim([0 2*pi])
    ylim([0 2*pi])
    lsline
    box on
    title([T ': circ_cc = ' num2str(rho(a)) ', p = ' num2str(pval(a))], 'interpreter', 'none')
    
    if a == 1
        xlabel('Rest: Average Modulation Phase, radians')
        ylabel('Task: Average Modulation Phase, radians')
    end
end
[~, sig_idx] = bonf_holm(pval, 0.05);
save_figure_as(gcf,'Scatters_Modulation_Phases', basepath_to_save, 1)

%% scatters for PSTH modulation strength
figure,
set(gcf, 'Position', [63 525 1785 421])
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    
    task_mi = Out.(T).Task.histogram_MI;
    rest_mi = Out.(T).Rest.histogram_MI;
    
    task_p = Out.(T).Task.histogram_p;
    rest_p = Out.(T).Rest.histogram_p;
    
    [~, task_sig_idx] = bonf_holm(task_p, 0.05); % correct for multiple comparisons
    [~, rest_sig_idx] = bonf_holm(rest_p, 0.05); % correct for multiple comparisons
    
    if numel(task_mi) > 3
        [rho_tmp, pval_tmp] = corrcoef(rest_mi(rest_sig_idx & task_sig_idx), task_mi(rest_sig_idx & task_sig_idx));
    else
        [rho_tmp, pval_tmp] = deal(nan(2));
    end
    rho(a) = rho_tmp(2,1);
    pval(a) = pval_tmp(2,1);
    
    subplot(1, N_Areas, a)
    scatter(100*rest_mi(rest_sig_idx & task_sig_idx), 100*task_mi(rest_sig_idx & task_sig_idx), [], [0 0 0])
    xlim([0 60])
    ylim([0 60])
    lsline
    box on
    title([T ': cc = ' num2str(rho(a)) ', p = ' num2str(pval(a))])
    
    if a == 1
        xlabel('Rest: Average Modulation Phase, radians')
        ylabel('Task: Average Modulation Phase, radians')
    end
end
save_figure_as(gcf,'Scatters_Modulation_Strengths', basepath_to_save, 1)

%% stem plots for average modulation phases vs. modulation strength
figure,
set(gcf, 'Position', [190 106 1533 884])
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        curr_mi = Out.(T).(L).histogram_MI;
        curr_ff = Out.(T).(L).histogram_phase;
        curr_pp = Out.(T).(L).histogram_p;
        [cor_p, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
        
        subplot(N_conditions, N_Areas, (c-1)*N_Areas + a)
        stem(curr_ff(sig_idx), 100*curr_mi(sig_idx), 'Color', cfg.condition(c).color) % 
        xlim([0 2*pi])
        
        title([L ':' T])
        
        if a == 1 && c == 1
            xlabel('Phase of the Heart Cycle, radians')
            ylabel('% Modulation')
        end
        
    end
end

save_figure_as(gcf,'Stem_ModulationPhases_vs_ModulationStrengths', basepath_to_save, 1)

function save_figure_as(fig_id,filename,basepath_to_save,savePlot)
if savePlot;
    export_fig(fig_id,[basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close all;
end
end

end