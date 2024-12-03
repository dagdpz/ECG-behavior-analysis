function ecg_bna_population_correlation_analysis(cfg, data_folder, unitList)

for targGrNum = 1:length(cfg.targets_spike_data)
    
    unqTargets = cfg.targets_spike_data{targGrNum};
    N_Areas = length(unqTargets);
    N_conditions = length(cfg.condition);
    
    var_list = {'pearson_r', 'pearson_p', 'n_cycles', ...
        'FR_pearson_r', 'FR_pearson_p', 'FR_n_cycles', ...
        'RR_pearson_r', 'RR_pearson_p', 'RR_n_cycles'};
    
    var_prefix = {'', 'FR_', 'RR_'};
    
    bar_colors = [244 149 173; ... % pink  - pos
        126 221 95; ...    % green - neg
        127 127 127]/255;  % grey  - non-significant
    
    pos_neg_lag_colors = [170 60 0;         % pos
        255 153 85;         % neg
        255 255 255]/255;   % zero
    
    % output folder
    targSuff = [];
    for targNum = 1:N_Areas
        targSuff = [targSuff '_' unqTargets{targNum}];
    end
    basepath_to_save = [cfg.SPK_root_results_fldr filesep 'Population_correlation_analysis' unitList(9:end) targSuff];
    if ~exist(basepath_to_save, 'dir')
        mkdir(basepath_to_save)
    end
    
    %% load data
    load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG' filesep unitList], 'unit_ids', 'targets', 'ids_both')
    
    dataset_name       = [basepath_to_save filesep 'data.mat'];
    ephys_dataset_name = [cfg.SPK_root_results_fldr filesep 'Population_time_domain' unitList(9:end) targSuff filesep 'Output.mat'];
    
    if ~exist(dataset_name,'file')
        for a = 1: N_Areas
            T=unqTargets{a};
            currTargIds = cellfun(@(x) strcmp(x, unqTargets{a}), targets);
            curr_unit_ids = unit_ids(currTargIds);
            curr_ids_both = ids_both(currTargIds);
            dt.(T) = ecg_bna_load_variables(cfg, curr_unit_ids, data_folder, 'data', var_list, curr_ids_both);
        end
        save(dataset_name,'dt')
    else
        load(dataset_name,'dt')
    end
    load(ephys_dataset_name,'Out') % load ephys dataset anyway
    
    % Sankey plot options
    default_colormap  = [0.7 0.7 0.7; 0.8500 0.3250 0.0980; 0 0.4470 0.7410; ...
        0.7 0.7 0.7; bar_colors(1,:); bar_colors(2,:)];
    options.color_map         = default_colormap;
    options.flow_transparency = 0.8;    % opacity of the flow paths
    options.bar_width         = 100;    % width of the category blocks
    options.show_perc         = true;   % show percentage over the blocks
    options.text_color        = [0 0 0];% text color for the percentages
    options.show_layer_labels = true;   % show layer names under the chart
    options.show_cat_labels   = true;   % show categories over the blocks.
    options.show_legend       = false;  % show legend with the category names.
    % if the data is not a table, then the
    % categories are labeled as catX-layerY
    
    options2.color_map         = [0.7 0.7 0.7; ([113 113 255])/255; 0 0 0; ...
        0.7 0.7 0.7; bar_colors(1,:); bar_colors(2,:)];
    options2.flow_transparency = 0.8;    % opacity of the flow paths
    options2.bar_width         = 100;    % width of the category blocks
    options2.show_perc         = true;   % show percentage over the blocks
    options2.text_color        = [0 0 0];% text color for the percentages
    options2.show_layer_labels = true;   % show layer names under the chart
    options2.show_cat_labels   = true;   % show categories over the blocks.
    options2.show_legend       = false;  % show legend with the category names.
    % if the data is not a table, then the
    % categories are labeled as catX-layerY
    
    %% autocorrelation & correlation coefs between FR and RR
    for prefixNum = 1:length(var_prefix)
        
        f0 = figure;
        set(f0, 'Position', [22 441 1857 491])
        
        if isempty(var_prefix{prefixNum})
            sgttl = 'Correlation: FR vs. RR';
        elseif contains(var_prefix{prefixNum}, 'FR_')
            sgttl = 'Autocorrelation FR';
        elseif contains(var_prefix{prefixNum}, 'RR_')
            sgttl = 'Autocorrelation RR';
        end
        
        for a = 1: N_Areas
            
            T=unqTargets{a};
            
            for c=1:N_conditions
                
                L=cfg.condition(c).name;
                
                % take data
                R=dt.(T).(L).pearson_r;
                P=dt.(T).(L).pearson_p;
                
                % figure out significance at lag zero
                [sig_idx,~] = fdr_bky(P(ceil(end/2),:), 0.05);
                
                pos_sig = R(ceil(end/2),:) > 0 & sig_idx;
                neg_sig = R(ceil(end/2),:) < 0 & sig_idx;
                non_sig = ~(pos_sig | neg_sig);
                
                splot_num = (c-1)*N_Areas + a;
                subplot(N_conditions, N_Areas, splot_num)
                hold on
                box on
                if sum(non_sig)>0
                    p1 = plot(dt.(T).cc_lag_list(1,:), dt.(T).(L).([var_prefix{prefixNum} 'pearson_r'])(:, non_sig), 'Color', bar_colors(3,:));
                else
                    p1 = plot(NaN,NaN);
                end
                if sum(pos_sig)>0
                    p2 = plot(dt.(T).cc_lag_list(1,:), dt.(T).(L).([var_prefix{prefixNum} 'pearson_r'])(:, pos_sig), 'Color', bar_colors(1,:));
                else
                    p2 = plot(NaN,NaN);
                end
                if sum(neg_sig)>0
                    p3 = plot(dt.(T).cc_lag_list(1,:), dt.(T).(L).([var_prefix{prefixNum} 'pearson_r'])(:, neg_sig), 'Color', bar_colors(2,:));
                else
                    p3 = plot(NaN,NaN);
                end
                
                xlim([-12 12])
                ylim([-1 1])
                
                if splot_num == 1
                    legend([p1(1) p2(1) p3(1)],{'Lag 0: Non-sig.', 'Lag 0: Sig.Pos.', 'Lag 0: Sig.Neg.'}, 'Location', 'best')
                end
                if splot_num == 2
                    title([sgttl ': ' T ' - ' L])
                else
                    title([T ':' L])
                end
                
            end
            
        end
        
        subplot(N_conditions, N_Areas, 1)
        xlabel('Shift of RR sequence relative to FR, # cardiac cycles')
        ylabel('Correlation Coefficients')
        
        save_figure_as(f0, [var_prefix{prefixNum} 'CC_Lags_AllAreas'], basepath_to_save, 1)
        
    end
    
    %% (1) CC histograms for one lag, for areas and conditions
    %% (2) pie charts
    %% (3) histograms of unit fractions
    %% (4) sankey plots between R-peak responsive and FR-RR correlated units
    bin_resolution=0.05;
    
    f1 = figure;
    set(f1, 'Position', [667 519 930 477])
    lag_to_plot=0;
    currLag=cfg.correlation.lag_list==lag_to_plot;
    
    f2 = figure;
    set(f2, 'Position', [301 229 1009 655])
    
    f3 = figure;
    set(f3, 'Position', [702 411 1026 449])
    
    % preallocate variables
    p_sign = nan(N_Areas,N_conditions);
    
    p_fisher_sig_nonsig_cc   = nan(N_Areas,1);
    odds_ratio_sig_nonsig_cc = nan(N_Areas,1);
    
    for a = 1: N_Areas
        
        T=unqTargets{a};
        colormap(bar_colors)
        
        % clean up this variable
        responsive_rest_task = [];
        
        %     y_lim=1;
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            curr_cc = dt.(T).(L).pearson_r(currLag, :);
            curr_pp = dt.(T).(L).pearson_p(currLag, :);
            %         [~, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
            [sig_idx,~] = fdr_bky(curr_pp, 0.05);
            pos_idx = curr_cc > 0;
            neg_idx = curr_cc < 0;
            nan_idx = isnan(curr_cc);
            
            % prepare data for plotting
            counts_pos    = histc(curr_cc(sig_idx & pos_idx), -1:bin_resolution:1); % significant and positive
            counts_neg    = histc(curr_cc(sig_idx & neg_idx), -1:bin_resolution:1); % significant and negative
            counts_nonsig = histc(curr_cc(~sig_idx & ~nan_idx), -1:bin_resolution:1);
            
            counts_pos    = counts_pos(:);
            counts_neg    = counts_neg(:);
            counts_nonsig = counts_nonsig(:);
            
            % proportions for histograms
            proportions_pos    = counts_pos / length(curr_cc);
            proportions_neg    = counts_neg / length(curr_cc);
            proportions_nonsig = counts_nonsig / length(curr_cc);
            
            % compute medians and IQRs
            M_median(a,c) = nanmedian(curr_cc);
            I(a,c) = iqr(curr_cc);
            
            % compute means and stds
            M_mean(a,c) = nanmean(curr_cc);
            S_STD(a,c)  = std(curr_cc);
            
            % compute sign test to know whether medians differ from 0
            p_sign(a,c) = signtest(curr_cc);
            
            % before computing t-test, make sure those distributions are normal
            % with Kolmogorov-Smirnov test
            [~,p_kstest(a,c)] = kstest(zscore(rmmissing(curr_cc)));
            
            % compute one-sample t-test to know whether averages differ from 0
            [~,p_ttest(a,c)] = ttest(rmmissing(curr_cc));
            
            % (1) plot histograms for the current lag
            figure(f1);
            
            splot_num = (c-1)*N_Areas + a;
            plot_data = [proportions_pos proportions_neg proportions_nonsig];
            s(a,c) = subplot(N_conditions, N_Areas, splot_num);
            b = bar(-1:bin_resolution:1, plot_data, 'stacked');
            set(b, 'FaceColor', 'Flat')
            for ii = 1:size(bar_colors,1)
                b(ii).FaceColor = bar_colors(ii,:);
            end
            xlim([-0.5 0.5])
            title([T ': ' L ', ' num2str(lag_to_plot) ' lag'])
            
            if a == 1 && c == 1
                legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'})
                xlabel('CC between FR and RR duration')
                ylabel('Unit Counts')
            end
            
            % (2) pie charts
            figure(f2);
            
            splot_num = (c-1)*N_Areas + a;
            subplot(N_conditions, N_Areas, splot_num)
            
            unit_counts_sig_nonsig_by_condition(c,:,a) = ...
                [sum(counts_pos)+sum(counts_neg) sum(counts_nonsig)]; % consider significant alltogether
            
            % prepare responsiveness table to compute Cohen's kappa further
            responsive_rest_task(:,c) = sig_idx';
            
            unit_counts_by_condition(c,:) = [sum(counts_pos) sum(counts_neg) sum(counts_nonsig)];
            
            unit_counts_by_area(c,a,:) = permute([sum(counts_pos) sum(counts_neg) sum(counts_nonsig)], [1 3 2]);
            
            pie_id = pie(unit_counts_by_condition(c,:));
            for ii = 1:2:length(pie_id)
                set(pie_id(ii),'FaceColor', bar_colors(mod(ii,2)+floor(ii/2),:))
            end
            
            title([T ':' L])
            
            % (4) Sankey plot
            % prepare the categorical data for plotting
            t = [];
            t(:,1) = Out.(T).(L).HeartResponseType_byUnit; % time-domain data
            t(t(:,1) == -1,1) = 2; % replace -1 with 2 to sort things properly
            t(:,2) = NaN;
            t(pos_idx & sig_idx,2)   = 1;
            t(neg_idx & sig_idx,2)   = 2;
            t(~sig_idx & ~nan_idx,2) = 0;
            
            nan_ids      = any(isnan(t),2);
            t(nan_ids,:) = [];
            
            t = sortrows(t,[1 2],'ascend');
            
            % set up proper colormap - if there are little data, category order
            % can be different from the sorted
            unqItems_2nd_col_1cat = unique(t(t(:,1) == 0,2)); % figure out how many unique categories from the 2nd column can be found in the 1st category of the 1st column
            if length(unqItems_2nd_col_1cat) == 1
                
                if unqItems_2nd_col_1cat == 2
                    color_id = 3;
                elseif unqItems_2nd_col_1cat == 1
                    color_id = 2;
                elseif unqItems_2nd_col_1cat == 0
                    color_id = 1;
                end
                
                other_color_ids = find(~ismember([1 2 3],color_id));
                new_colormap      = [options.color_map(1:3,:); options.color_map(3+color_id,:); options.color_map(3+other_color_ids,:)];
                options.color_map = new_colormap;
                
            elseif length(unqItems_2nd_col_1cat) == 2
                
            end
            %         options.color_map = [0.7 0.7 0.7; 0.8500 0.3250 0.0980; 0 0.4470 0.7410; ...
            %                     0.7 0.7 0.7; bar_colors(1,:); bar_colors(2,:)];
            
            tmp = cell(size(t));
            tmp(t(:,1) == 0,1) = deal({'no response'});
            tmp(t(:,1) == 1,1) = deal({'increase'});
            tmp(t(:,1) == 2,1) = deal({'decrease'});
            
            tmp(t(:,2) == 0,2) = deal({'no corr.'});
            tmp(t(:,2) == 1,2) = deal({'pos.corr.'});
            tmp(t(:,2) == 2,2) = deal({'neg.corr.'});
            
            tmp = cell2table(tmp,'VariableNames',{'R-peak_Resp.', 'FR-RR_Corr.'});
            
            tmp.('R-peak_Resp.') = categorical(tmp.('R-peak_Resp.'));
            tmp.('FR-RR_Corr.') = categorical(tmp.('FR-RR_Corr.'));
            
            %         tmp = sortrows(tmp,[1,2],'descend');
            
            f4 = figure;
            plotSankeyFlowChart(tmp,options)
            title([L ': ' T])
            save_figure_as(f4,['Rpeak_vs_FR-RR-Corr_SankeyPlot_' L '_' T],basepath_to_save,1)
            clear tmp
            options.color_map = default_colormap;
            
            % compute Fisher's test to figure out contingency between R-peak
            % responsiveness and FR-RR correlation
            noRpeak_noCorr = sum(t(:,1) == 0 & t(:,2) == 0);
            noRpeak_Corr   = sum(t(:,1) == 0 & t(:,2) ~= 0);
            Rpeak_noCorr   = sum(t(:,1) ~= 0 & t(:,2) == 0);
            Rpeak_Corr     = sum(t(:,1) ~= 0 & t(:,2) ~= 0);
            con_tbl = [noRpeak_noCorr noRpeak_Corr;
                Rpeak_noCorr Rpeak_Corr];
            [~,p_fisher_Rpeak_Corr(a,c),stats] = fishertest(con_tbl);
            oddsratio_fisher_Rpeak_Corr(a,c) = stats.OddsRatio;
            clear stats
            
            % create and save contingency table
            filename = ['Rpeak_FR-RR_Contingency_' L '_' T];
            TBL = table({'';'noRpeak';'Rpeak'}, [NaN; con_tbl(:,1)], [NaN; con_tbl(:,2)], ...
                'VariableNames', {'-', 'noCorr', 'Corr'});
            writetable(TBL, [basepath_to_save filesep filename '.xls'])
            clear TBL con_tbl
            
            % create Sankey plot for R-peak related increase and decrease
            % together
            f5 = figure;
            t = [];
            t(:,1) = Out.(T).(L).HeartResponseType_byUnit; % time-domain data
            t(t(:,1) ~= 0,1) = 1; % replace sig modulated with 1
            t(:,2) = NaN;
            t(pos_idx & sig_idx,2)   = 1;
            t(neg_idx & sig_idx,2)   = 2;
            t(~sig_idx & ~nan_idx,2) = 0;
            
            nan_ids      = any(isnan(t),2);
            t(nan_ids,:) = [];
            
            t = sortrows(t,[1 2],'ascend');
            
            tmp = cell(size(t));
            tmp(t(:,1) == 0,1) = deal({'no response'});
            tmp(t(:,1) == 1,1) = deal({'response'});
            
            tmp(t(:,2) == 0,2) = deal({'no corr.'});
            tmp(t(:,2) == 1,2) = deal({'pos.corr.'});
            tmp(t(:,2) == 2,2) = deal({'neg.corr.'});
            
            tmp = cell2table(tmp,'VariableNames',{'R-peak_Resp.', 'FR-RR_Corr.'});
            
            tmp.('R-peak_Resp.') = categorical(tmp.('R-peak_Resp.'));
            tmp.('FR-RR_Corr.') = categorical(tmp.('FR-RR_Corr.'));
            
            plotSankeyFlowChart(tmp,options2)
            title([L ': ' T])
            save_figure_as(f5,['Rpeak(Merged)_vs_FR-RR-Corr_SankeyPlot_' L '_' T],basepath_to_save,1)
            clear tmp
            
            % compute chi2-test and create cross-tabulation table
            [tbl,chi2_Rpeak_Corr(a,c),p_chi2_Rpeak_Corr(a,c)] = crosstab(t(:,1), t(:,2));
            con_tbl = table({'No R-peak Resp.', 'R-peak Resp.'}', tbl(:,1), tbl(:,2), tbl(:,3), ...
                'VariableNames', {'-', 'No FR-RR Corr.', 'Pos. FR-RR Corr.', 'Neg. FR-RR Corr.'});
            
        end
        
        % compute exact Fisher's test to test for unit count differences
        % between conditions
        p_fisher(:,a) = ecg_bna_fisher_test(unit_counts_by_condition(1,:), unit_counts_by_condition(2,:));
        p_fisher = round(p_fisher,10);
        
        % create table with unit counts for pos., neg. and non-correlated units
        filename = ['PosNegCorr_UnitCounts_' T];
        TBL = table({'Pos. Corr.'; 'Neg. Corr.'; 'Non-Sig.'}, unit_counts_by_condition(1,:)', unit_counts_by_condition(2,:)', ...
            'VariableNames', {'FR-RR Corr.', ...
            'Rest: Unit Counts', 'Task: Unit Counts'});
        writetable(TBL, [basepath_to_save filesep filename '.xls'])
        clear TBL
        
        % compute exact Fisher's test to figure out association between
        % significantly modulated units in the rest and in the task
        [~,p_fisher_sig_nonsig_cc(a),stats] = fishertest(unit_counts_sig_nonsig_by_condition(:,:,a));
        odds_ratio_sig_nonsig_cc(a)         = stats.OddsRatio;
        
        % save contingency tables for each area
        filename = ['SigCC_NonSigCC_Rest_Task_' T];
        TBL = table({'Rest';'Task'}, [unit_counts_sig_nonsig_by_condition(:,1,a)], [unit_counts_sig_nonsig_by_condition(:,2,a)], ...
            'VariableNames', {'-', 'SigCC', 'NonSigCC'});
        writetable(TBL, [basepath_to_save filesep filename '.xls'])
        clear TBL
        
        % compute chi-squared for cross-tabulated results of responsive vs
        % non-responsive
        [~,chi2(a),p_chi2(a)] = crosstab(responsive_rest_task(:,1), responsive_rest_task(:,2));
        
    end
    
    % create table with Fisher's test results checking for contingency
    % between R-peak responsiveness and FR-RR correlation
    p_fisher_Rpeak_Corr = round(p_fisher_Rpeak_Corr,10);
    [h_fisher, p_crit_fisher] = fdr_bky(p_fisher_Rpeak_Corr); % correct for multiple comparisons
    filename = 'Fisher_Test_Rpeak_FR-RR_Contingency';
    TBL = table(unqTargets', ...
        h_fisher(:,1), p_fisher_Rpeak_Corr(:,1), oddsratio_fisher_Rpeak_Corr(:,1), ...
        h_fisher(:,2), p_fisher_Rpeak_Corr(:,2), oddsratio_fisher_Rpeak_Corr(:,2), [p_crit_fisher; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: Fisher''s h', 'Rest: Fisher''s p', 'Rest: Fisher''s Odds Ratio', ...
        'Task: Fisher''s h', 'Task: Fisher''s p', 'Task: Fisher''s Odds Ratio', 'FDR Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL h_fisher p_crit_fisher
    
    % create table with chi2 results checking for contingency
    % between R-peak responsiveness and FR-RR correlation
    p_chi2_Rpeak_Corr = round(p_chi2_Rpeak_Corr,10);
    [h_chi2, p_crit_chi2] = fdr_bky(p_fisher_Rpeak_Corr); % correct for multiple comparisons
    filename = 'Chi2_Test_Rpeak_FR-RR_Contingency';
    TBL = table(unqTargets', ...
        h_chi2(:,1), p_chi2_Rpeak_Corr(:,1), chi2_Rpeak_Corr(:,1), ...
        h_chi2(:,2), p_chi2_Rpeak_Corr(:,2), chi2_Rpeak_Corr(:,2), [p_crit_chi2; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: Chi2 h', 'Rest: Chi2 p', 'Rest: Chi2', ...
        'Task: Chi2 h', 'Task: Chi2 p', 'Task: Chi2 Odds Ratio', 'FDR Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL h_fisher p_crit_fisher
    
    % (3) plot unit fraction histograms
    figure(f3);
    for c = 1:N_conditions
        
        % proportion histograms
        subplot(1,N_conditions,c)
        tmp = permute(unit_counts_by_area(c,:,:), [2 3 1]);
        unit_percentages = 100 * tmp ./ sum(tmp,2);
        b = bar(unit_percentages,'stacked');
        for ii = 1:size(bar_colors,1)
            b(ii).FaceColor = bar_colors(ii,:);
        end
        title(cfg.condition(c).name)
        set(gca,'XTickLabel',unqTargets)
        legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'},'Location','southoutside')
        
    end
    
    subplot(1,N_conditions,1)
    ylabel('Unit fractions, %')
    
    % link histogram axes
    linkaxes(s)
    
    save_figure_as(f1, ['Histograms_CC_FR_&_RR_at_lag_' num2str(lag_to_plot)], basepath_to_save, 1)
    save_figure_as(f2, 'Percentages_Pie_Charts', basepath_to_save, 1)
    save_figure_as(f3, 'Unit_Pcs', basepath_to_save, 1)
    
    % correct sig test results for multiple comparisons
    [h_sign, p_crit_sign] = fdr_bky(p_sign);
    
    % construct table with medians and their signtest results
    filename = 'Median_CCs';
    TBL = table(unqTargets', M_median(:,1), p_sign(:,1), h_sign(:,1), I(:,1), M_median(:,2), p_sign(:,2), h_sign(:,2), I(:,2), [p_crit_sign; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: Median cc', 'Rest: Sign Test p', 'Rest: Sign Test h', 'Rest: IQR', ...
        'Task: Median cc', 'Task: Sign Test p', 'Task: Sign Test h', 'Task: IQR', 'FDR: Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % correct Kolmogorov-Smirnov test results for multiple comparisons
    [h_kstest, p_crit_kstest] = fdr_bky(p_kstest);
    
    % table with Kolmogorov-Smirnov test results
    filename = 'Kolmogorov-Smirnov_test_Mean_ccs';
    TBL = table(unqTargets', p_kstest(:,1), h_kstest(:,1), p_kstest(:,2), h_kstest(:,2), [p_crit_kstest; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: KS p-value', 'Rest: KS h', ...
        'Task: KS p-value', 'Task: KS h', 'FDR: Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % correct ttest results for multiple comparisons
    [h_ttest, p_crit_ttest] = fdr_bky(p_ttest);
    
    % table with one-sample t-test results
    filename = 'One-Sample_t-test_results_Mean_ccs';
    TBL = table(unqTargets', M_mean(:,1), S_STD(:,1), p_ttest(:,1), h_ttest(:,1), M_mean(:,2), S_STD(:,2), p_ttest(:,2), h_ttest(:,2), [p_crit_ttest; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: cc mean', 'Rest: cc std', 'Rest: t-test p-value', 'Rest: t-test h', ...
        'Task: cc mean', 'Task: cc std', 'Task: t-test p-value', 'Task: t-test h', 'FDR: Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    
    % correct Fisher's test results for multiple comparisons
    [h_fisher, p_crit_fisher] = fdr_bky(p_fisher);
    
    % create table with exact Fisher's test results: rest vs. task
    filename = 'Exact_Fisher_Test_Rest_vs_Task';
    TBL = table(unqTargets', p_fisher(1,:)', h_fisher(1,:)', p_fisher(2,:)', h_fisher(2,:)', [p_crit_fisher; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Pos.Corr.: Fisher''s p', 'Pos.Corr.: Fisher''s h', 'Neg.Corr.: Fisher''s p', 'Neg.Corr.: Fisher'' h', 'FDR: Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % compute Fisher's test between areas
    id_VPL  = 1;
    id_dPul = 2;
    id_MD   = 3;
    for c = 1: N_conditions
        
        % rest - VPL  vs. dPul
        p_fisher_VPL_vs_dPul(:,c) = ecg_bna_fisher_test(unit_counts_by_area(c,id_VPL,:), unit_counts_by_area(c,id_dPul,:));
        
        % VPL  vs. MD
        p_fisher_VPL_vs_MD(:,c) = ecg_bna_fisher_test(unit_counts_by_area(c,id_VPL,:), unit_counts_by_area(c,id_MD,:));
        
        % dPul vs. MD
        p_fisher_dPul_vs_MD(:,c) = ecg_bna_fisher_test(unit_counts_by_area(c,id_dPul,:), unit_counts_by_area(c,id_MD,:));
        
    end
    
    % group Fisher's results by condition
    p_fisher_rest = [p_fisher_VPL_vs_dPul(:,1) p_fisher_VPL_vs_MD(:,1) p_fisher_dPul_vs_MD(:,1)];
    p_fisher_rest = round(p_fisher_rest,10);
    p_fisher_task = [p_fisher_VPL_vs_dPul(:,2) p_fisher_VPL_vs_MD(:,2) p_fisher_dPul_vs_MD(:,2)];
    p_fisher_task = round(p_fisher_task,10);
    
    % apply multiple comparison correction
    [h_fisher_rest, p_crit_fisher_rest] = fdr_bky(p_fisher_rest);
    [h_fisher_task, p_crit_fisher_task] = fdr_bky(p_fisher_task);
    
    % [rest] create tables with Fisher's test results between areas
    filename = 'Exact_Fisher_Test_Rest';
    TBL = table({'VPL-dPul'; 'VPL-MD'; 'dPul-MD'}, p_fisher_rest(1,:)', h_fisher_rest(1,:)', p_fisher_rest(2,:)', h_fisher_rest(2,:)', [p_crit_fisher_rest; NaN; NaN], ...
        'VariableNames', {'Area Combination', ...
        'Pos.Corr.: Fisher''s p', 'Pos.Corr.: Fisher''s h', 'Neg.Corr.: Fisher''s p', 'Neg.Corr.: Fisher'' h', 'FDR: Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % [rest] create tables with Fisher's test results between areas
    filename = 'Exact_Fisher_Test_Task';
    TBL = table({'VPL-dPul'; 'VPL-MD'; 'dPul-MD'}, p_fisher_task(1,:)', h_fisher_task(1,:)', p_fisher_task(2,:)', h_fisher_task(2,:)', [p_crit_fisher_task; NaN; NaN], ...
        'VariableNames', {'Area Combination', ...
        'Pos.Corr.: Fisher''s p', 'Pos.Corr.: Fisher''s h', 'Neg.Corr.: Fisher''s p', 'Neg.Corr.: Fisher'' h', 'FDR: Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % compare association between units counts of sig / non-sig CC in rest /
    % task
    p_fisher_sig_nonsig_cc = round(p_fisher_sig_nonsig_cc,10);
    [h_fisher_sig_nonsig_cc, p_crit_fisher_sig_nonsig_cc] = fdr_bky(p_fisher_sig_nonsig_cc);
    
    % create a table with Fisher's test results for association between
    % significance and condition
    filename = 'Fisher_Test_corrSign_condition';
    TBL = table(unqTargets', ...
        h_fisher_sig_nonsig_cc, p_fisher_sig_nonsig_cc, odds_ratio_sig_nonsig_cc, [p_crit_fisher_sig_nonsig_cc; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: Fisher''s h', 'Rest: Fisher''s p', 'Rest: Fisher''s Odds Ratio', 'FDR Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % create a table with chi-squared results
    [h_chi2, p_crit_chi2] = fdr_bky(p_chi2);
    
    chi2(a),p_chi2(a)
    filename = 'Chi2_FR-RR_Correlated_Rest_Task';
    TBL = table(unqTargets', ...
        chi2', p_chi2', h_chi2', [p_crit_chi2; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Chi2', 'p for chi2', 'h', 'FDR Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % %% create Sankey plots between FR-RR correlation sign unit proportions between rest and task
    %
    % % options
    % options.color_map = [0.7 0.7 0.7; 0.8500 0.3250 0.0980; 0 0.4470 0.7410; ...
    %                     0.7 0.7 0.7; 0.8500 0.3250 0.0980; 0 0.4470 0.7410];
    % options.flow_transparency = 0.8;   % opacity of the flow paths
    % options.bar_width = 100;            % width of the category blocks
    % options.show_perc = true;          % show percentage over the blocks
    % options.text_color = [0 0 0];      % text color for the percentages
    % options.show_layer_labels = true;  % show layer names under the chart
    % options.show_cat_labels = true;    % show categories over the blocks.
    % options.show_legend = false;       % show legend with the category names.
    %                                    % if the data is not a table, then the
    %                                    % categories are labeled as catX-layerY
    %
    % %                                    L1 = cfg.con
    % %                                    L2 = [];
    %
    % for a = 1: N_Areas
    % 	T=Ana_TargetBrainArea{a};
    %
    %
    %
    %     curr_cc = dt.(T).(L).pearson_r(currLag, :);
    %     curr_pp = dt.(T).(L).pearson_p(currLag, :);
    %     %         [~, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
    %     [sig_idx,~] = fdr_bky(curr_pp, 0.05);
    %     pos_idx = curr_cc > 0;
    %     neg_idx = curr_cc < 0;
    %     nan_idx = isnan(curr_cc);
    %
    %
    %
    %     t(:,2) = NaN;
    %     t(pos_idx & sig_idx,2)   = 1;
    %     t(neg_idx & sig_idx,2)   = 2;
    %     t(~sig_idx & ~nan_idx,2) = 0;
    %
    %     nan_ids      = any(isnan(t),2);
    %     t(nan_ids,:) = [];
    %
    %     t = sortrows(t,[1 2],'ascend');
    %
    %
    %
    %
    %
    %
    %     % prepare the categorical data for plotting
    %     t = [Out.(T).Rest.HeartResponseType_byUnit, Out.(T).Task.HeartResponseType_byUnit];
    %
    %     nan_ids    = any(isnan(t),2);
    %     t(nan_ids,:) = [];
    %
    %     tmp = cell(size(t));
    %     tmp(t == 0)  = deal({'no response'});
    %     tmp(t == 1)  = deal({'increase'});
    %     tmp(t == -1) = deal({'decrease'});
    %
    %     tmp = cell2table(tmp,'VariableNames',{'Rest','Task'});
    %
    %     tmp.Rest = categorical(tmp.Rest);
    %     tmp.Task = categorical(tmp.Task);
    %
    %     tmp = sortrows(tmp,[1,2],'descend');
    %
    %     figure,
    %     plotSankeyFlowChart(tmp,options)
    %     title(T)
    %     save_figure_as(['SankeyPlot_' T],output_folder,savePlot)
    % end
    
    %% explained variance
    f0 = figure;
    set(f0, 'Position', [301 229 1009 655])
    
    % preallocate variables
    [p_mann_whitney_pos, p_mann_whitney_neg] = deal(nan(N_Areas,1));
    
    for a = 1: N_Areas
        
        T=unqTargets{a};
        
        for c=1:N_conditions
            
            L=cfg.condition(c).name;
            
            % take data
            R=dt.(T).(L).pearson_r(ceil(end/2),:);
            P=dt.(T).(L).pearson_p(ceil(end/2),:);
            
            % figure out significance at lag zero
            [sig_idx,crit_p] = fdr_bky(P, 0.05);
            pos_sig = R > 0 & sig_idx;
            neg_sig = R < 0 & sig_idx;
            non_sig = ~(pos_sig | neg_sig);
            
            % compute fraction of explained variance
            var_pos_sig = R(pos_sig) .^ 2 * 100;
            var_neg_sig = R(neg_sig) .^ 2 * 100;
            var_non_sig = R(non_sig) .^ 2 * 100;
            
            % compute median variances
            expl_var_pos{a,c} = var_pos_sig;
            expl_var_neg{a,c} = var_neg_sig;
            M_all(a,c)    = nanmedian([var_pos_sig var_neg_sig var_non_sig]);
            M_sig(a,c)    = nanmedian([var_pos_sig var_neg_sig]);
            
            % compute histograms
            h_pos_sig = histc(var_pos_sig, 0:0.5:15);
            h_neg_sig = histc(var_neg_sig, 0:0.5:15);
            h_non_sig = histc(var_non_sig, 0:0.5:15);
            
            splot_num = (c-1)*N_Areas + a;
            subplot(N_conditions, N_Areas, splot_num)
            hold on
            box on
            b = bar(0:0.5:15, [h_pos_sig(:) h_neg_sig(:) h_non_sig(:)],'stacked');
            for ii = 1:size(bar_colors,1)
                b(ii).FaceColor = bar_colors(ii,:);
            end
            plot(M_all(a,c), 0, '^', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
            plot(M_sig(a,c), 0, '^', 'MarkerFaceColor', [153 101 21]/256, 'MarkerEdgeColor', [0 0 0])
            
            if splot_num == 1
                xlabel('Explained Variance, %')
                ylabel('Unit Counts')
                legend({'Sig.Pos.', 'Sig.Neg.', 'Non-sig.', 'Median: all', 'Median: sig.'},'Location','best')
            end
            
            title([T ':' L])
            
        end
        
        if ~isempty(expl_var_pos{a,1}) & ~isempty(expl_var_pos{a,2})
            p_mann_whitney_pos(a) = ranksum(expl_var_pos{a,1},expl_var_pos{a,2});
        end
        if ~isempty(expl_var_neg{a,1}) & ~isempty(expl_var_neg{a,2})
            p_mann_whitney_neg(a) = ranksum(expl_var_neg{a,1},expl_var_neg{a,2});
        end
        
    end
    save_figure_as(f0, 'Explained_Variance_Bars', basepath_to_save, 1)
    
    % create table with median variances
    filename = 'Median_Explained_Variances';
    TBL = table(unqTargets', M_all(:,1), M_all(:,2), M_sig(:,1), M_sig(:,2), ...
        'VariableNames', {'Target Area','Rest: Median Expl. Var., % (all)', 'Task: Median Expl. Var., % (all)', ...
        'Rest: Median Expl. Var., % (sig)', 'Task: Median Expl. Var., % (sig)'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % correct Mann-Whitney test results for multiple comparisons
    [h_mann_whitney, p_mann_whitney] = fdr_bky([p_mann_whitney_pos p_mann_whitney_neg]);
    filename = 'Explained_Variances_Mann-Whitney_U-Test';
    TBL = table(unqTargets', ...
        p_mann_whitney_pos, h_mann_whitney(:,1), ...
        p_mann_whitney_neg, h_mann_whitney(:,2), ...
        [p_mann_whitney; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Mann-Whitney p (pos.)', 'Mann-Whitney h (pos.)', ...
        'Mann-Whitney p (neg.)', 'Mann-Whitney h (neg.)', 'FDR crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    %% normalized correlation between FR and RR
    for prefixNum = 1:length(var_prefix)
        
        f0 = figure;
        set(f0, 'Position', [22 441 1857 491])
        
        %     sgttl = 'Correlation: FR vs. RR';
        
        for a = 1: N_Areas
            
            T=unqTargets{a};
            
            for c=1:N_conditions
                
                L=cfg.condition(c).name;
                
                figure(f0);
                %             compute normalized ccs
                R_FR_RR = dt.(T).(L).pearson_r(ceil(end/2),:);
                R       = dt.(T).(L).([var_prefix{prefixNum} 'pearson_r']);
                P       = dt.(T).(L).pearson_p(ceil(end/2),:);
                
                norm_r  = R ./ repmat(R(ceil(end/2),:),size(R,1),1);
                [sig_idx,crit_p] = fdr_bky(P, 0.05);
                pos_sig = R_FR_RR > 0 & sig_idx;
                neg_sig = R_FR_RR < 0 & sig_idx;
                non_sig = ~(pos_sig|neg_sig);
                
                splot_num = (c-1)*N_Areas + a;
                subplot(N_conditions, N_Areas, splot_num)
                hold on
                if sum(non_sig)>0
                    p1 = plot(dt.(T).cc_lag_list(1,:), norm_r(:, non_sig), 'Color', bar_colors(3,:), 'LineWidth', 0.5);
                else
                    p1 = plot(NaN,NaN);
                end
                
                
                if sum(non_sig)>0
                    p1 = plot(dt.(T).cc_lag_list(1,:), dt.(T).(L).([var_prefix{prefixNum} 'pearson_r'])(:, non_sig), 'Color', bar_colors(3,:));
                else
                    p1 = plot(NaN,NaN);
                end
                
                
                
                if sum(pos_sig)>0
                    plot(dt.(T).cc_lag_list(1,:), norm_r(:, pos_sig), 'Color', [0 0 0], 'LineWidth', 2)
                    p2 = plot(dt.(T).cc_lag_list(1,:), norm_r(:, pos_sig), 'Color', bar_colors(1,:), 'LineWidth', 0.5);
                else
                    p2 = plot(NaN,NaN);
                end
                if sum(neg_sig)>0
                    plot(dt.(T).cc_lag_list(1,:), norm_r(:, neg_sig), 'Color', [0 0 0], 'LineWidth', 2)
                    p3 = plot(dt.(T).cc_lag_list(1,:), norm_r(:, neg_sig), 'Color', bar_colors(2,:), 'LineWidth', 0.5);
                else
                    p3 = plot(NaN,NaN);
                end
                hline(-1)
                hline(0)
                hline(1)
                
                ylim([-3 3])
                if splot_num == 1 & prefixNum == 1
                    %                 legend([p1(1) p2(1) p3(1)],{'Lag 0: Non-sig.', 'Lag 0: Sig.Pos.', 'Lag 0: Sig.Neg.'}, 'Location', 'best')
                end
                if splot_num == 2
                    %                 title([sgttl ': ' T ' - ' L])
                else
                    title([T ':' L])
                end
                box on
            end
        end
        
        subplot(N_conditions, N_Areas, 1)
        xlabel('Shift of RR sequence relative to FR, # cardiac cycles')
        ylabel('cc normalized by cc(lag 0)')
        
        save_figure_as(f0, [var_prefix{prefixNum} 'CC_Lags_normalized_AllAreas'], basepath_to_save, 1)
        clear R P p1 p2 p3
    end
    
    %% lag scatters
    f1 = figure;
    set(f1, 'Position', [22 441 1857 491])
    
    sgttl = 'Lag of Most Effective Correlation';
    
    for a = 1: N_Areas
        T=unqTargets{a};
        
        for c=1:N_conditions
            
            L=cfg.condition(c).name;
            
            % plot optimal correlation coefficients
            % extract data
            R=dt.(T).(L).('pearson_r');
            P=dt.(T).(L).('pearson_p');
            [~, ids]  = max(abs(R), [], 1); % find max abs ccs; I will need 'ids' to figure out significance of corresponding ccs
            ccs=R(sub2ind(size(R),ids,1:numel(ids)));
            invalid=isnan(ccs);
            P=P(:,~invalid);
            R=R(:,~invalid);
            ccs=ccs(~invalid);
            ids=ids(~invalid);
            P=P(ceil(end/2),:);
            R=R(ceil(end/2),:);
            sig_idx = fdr_bky(P, 0.05); % figure out significant ones
            sig=sign(R).*sig_idx;
            
            splot_num = (c-1)*N_Areas + a;
            subplot(N_conditions, N_Areas, splot_num)
            hold on
            
            if sum(sig==0)
                p1 = plot(dt.(T).cc_lag_list(1,ids(sig==0)),  ccs(sig==0),  'o', 'Color', bar_colors(3, :)); % non-sig.
            else
                p1 = plot(NaN,NaN);
            end
            if sum(sig==-1)
                p2 = plot(dt.(T).cc_lag_list(1,ids(sig==-1)), ccs(sig==-1), 'o', 'Color', bar_colors(2, :)); % neg. sig.
            else
                p2 = plot(NaN,NaN);
            end
            if sum(sig==1)
                p3 = plot(dt.(T).cc_lag_list(1,ids(sig==1)),  ccs(sig==1),  'o', 'Color', bar_colors(1, :)); % pos. sig.
            else
                p3 = plot(NaN,NaN);
            end
            xlim([cfg.correlation.lag_list(1) cfg.correlation.lag_list(end)])
            ylim([-0.5 0.5])
            box on
            
            if splot_num == 1
                legend([p1(1) p2(1) p3(1)],{'Lag 0: Non-sig.', 'Lag 0: Sig.Pos.', 'Lag 0: Sig.Neg.'}, 'Location', 'best')
            end
            if splot_num == 2
                title([sgttl ': ' T ' - ' L])
            else
                title([T ':' L])
            end
        end
        
    end
    
    subplot(N_conditions, N_Areas, 1)
    xlabel('Shift of RR sequence relative to FR, # cardiac cycles')
    ylabel('Correlation Coefficient')
    
    save_figure_as(f1, 'CC_MaxAbsLags_AllAreas', basepath_to_save, 1)
    
    
    %% histograms of lag of max cc
    f9 = figure;
    set(f9, 'Position', [22 441 1857 491])
    
    f10 = figure;
    set(f10, 'Position', [22 441 1857 491])
    
    % preallocate variables
    M_median = nan(N_Areas,N_conditions);
    
    f11 = figure;
    % set(f11,'Position',[])
    
    for a = 1: N_Areas
        T=unqTargets{a};
        
        for c=1:N_conditions
            
            L=cfg.condition(c).name;
            
            % plot optimal correlation coefficients
            % extract data
            R=dt.(T).(L).('pearson_r');
            P=dt.(T).(L).('pearson_p');
            [~, ids]  = max(abs(R), [], 1); % find max abs ccs; I will need 'ids' to figure out significance of corresponding ccs
            ccs=R(sub2ind(size(R),ids,1:numel(ids)));
            invalid=isnan(ccs);
            P=P(:,~invalid);R=R(:,~invalid);ccs=ccs(~invalid);ids=ids(~invalid);
            P=P(ceil(end/2),:);
            R=R(ceil(end/2),:);
            sig_idx = fdr_bky(P, 0.05); % figure out significant ones
            sig=sign(R).*sig_idx;
            
            laglist=dt.(T).cc_lag_list(1,:);
            %         N=numel(ids);
            
            % compute lag medians only sig.
            M_median(a,c) = median(laglist(ids(sig_idx)));
            
            % compute sign test over sig lags
            [p_signtest_lag(a,c), h_signtest_lag(a,c)] = signtest(laglist(ids(sig_idx)));
            
            % comptute histogram counts
            h0=hist(laglist(ids(sig==0)),laglist);
            h1=hist(laglist(ids(sig==-1)),laglist);
            h2=hist(laglist(ids(sig==1)),laglist);
            
            % compute medians by subgroup
            m0(a,c)=median(laglist(ids(sig==0)));
            m1(a,c)=median(laglist(ids(sig==-1)));
            m2(a,c)=median(laglist(ids(sig==1)));
            
            % compute sign tests for pos and neg lags
            p_sign_pos(a,c) = signtest(laglist(ids(sig==1)));
            p_sign_neg(a,c) = signtest(laglist(ids(sig==-1)));
            
            % compute unit counts and fractions for pos and neg corr coef
            nonsig_counts(a,c) = sum(sig == 0);
            pos_counts(a,c)    = sum(sig == 1);
            neg_counts(a,c)    = sum(sig == -1);
            
            nonsig_prc(a,c) = nonsig_counts(a,c) / (nonsig_counts(a,c) + pos_counts(a,c) + neg_counts(a,c));
            pos_prc(a,c)    = pos_counts(a,c) / (nonsig_counts(a,c) + pos_counts(a,c) + neg_counts(a,c));
            neg_prc(a,c)    = neg_counts(a,c) /  (nonsig_counts(a,c) + pos_counts(a,c) + neg_counts(a,c));
            
            % compute unit counts and fractions with lag <> 0
            lag_zero_counts(a,c) = sum(laglist(ids) == 0 | ~sig_idx);
            pos_lag_counts(a,c)  = sum(laglist(ids) > 0 & sig_idx);
            neg_lag_counts(a,c)  = sum(laglist(ids) < 0 & sig_idx);
            
            lag_zero_prc(a,c) = lag_zero_counts(a,c) / (lag_zero_counts(a,c) + pos_lag_counts(a,c) + neg_lag_counts(a,c));
            pos_lag_prc(a,c)  = pos_lag_counts(a,c) / (lag_zero_counts(a,c) + pos_lag_counts(a,c) + neg_lag_counts(a,c));
            neg_lag_prc(a,c)  = neg_lag_counts(a,c) / (lag_zero_counts(a,c) + pos_lag_counts(a,c) + neg_lag_counts(a,c));
            
            % compute Fisher's test to figure out contingency between pos/neg
            % correlated units and pos/neg lag units
            posCorr_posLag = sum(sig == 1  & (laglist(ids) > 0 & sig_idx));
            posCorr_negLag = sum(sig == 1  & (laglist(ids) < 0 & sig_idx));
            negCorr_posLag = sum(sig == -1 & (laglist(ids) > 0 & sig_idx));
            negCorr_negLag = sum(sig == -1 & (laglist(ids) < 0 & sig_idx));
            
            con_tbl = [posCorr_posLag posCorr_negLag;
                negCorr_posLag negCorr_negLag];
            [~,p_fisher_corrSign_lagSign(a,c),stats] = fishertest(con_tbl);
            oddsratio_fisher_corrSign_lagSign(a,c) = stats.OddsRatio;
            clear stats
            
            % create and save contingency table
            filename = ['corrSign_vs_lagSign_' L '_' T];
            TBL = table({'';'posCorr';'negCorr'}, [NaN; con_tbl(:,1)], [NaN; con_tbl(:,2)], ...
                'VariableNames', {'-', 'posLag', 'negLag'});
            writetable(TBL, [basepath_to_save filesep filename '.xls'])
            clear TBL con_tbl
            
            % prepare ccs at peak lag for bars plots
            cc_neg_counts    = histc(ccs(ccs < 0 & sig_idx),-1:0.05:1);
            cc_pos_counts    = histc(ccs(ccs > 0 & sig_idx),-1:0.05:1);
            cc_nonsig_counts = histc(ccs(~sig_idx),-1:0.05:1);
            
            % compute medians for absolute ccs for pos and neg ccs
            M_median_ccs_neg(a,c) = median(ccs(ccs < 0 & sig_idx));
            M_median_ccs_pos(a,c) = median(ccs(ccs > 0 & sig_idx));
            
            % compute Mann-Whitney test to figure out differences between
            % absolute medians of pos and neg
            sig_pos_ids = ccs > 0 & sig_idx;
            sig_neg_ids = ccs < 0 & sig_idx;
            if sum(sig_pos_ids) & sum(sig_neg_ids)
                p_ccs_mann_whitney(a,c) = ranksum(abs(ccs(sig_pos_ids)),abs(ccs(sig_neg_ids)));
            end
            
            figure(f9),
            splot_num = (c-1)*N_Areas + a;
            subplot(N_conditions, N_Areas, splot_num)
            b = bar(-1:0.05:1,[cc_pos_counts' cc_neg_counts' cc_nonsig_counts'],'stacked');
            for ii = 1:size(bar_colors,1)
                b(ii).FaceColor = bar_colors(ii,:);
            end
            title([T ': ' L])
            clear b
            
            figure(f10);
            splot_num = (c-1)*N_Areas + a;
            subplot(N_conditions, N_Areas, splot_num)
            hold on
            box on
            
            b = bar(laglist, [h1(:) h2(:) h0(:)], 'stacked');
            colors_here = [bar_colors(2,:); bar_colors(1,:); bar_colors(3,:)];
            for ii = 1:size(bar_colors,1)
                b(ii).FaceColor = colors_here(ii,:);
            end
            clear b
            
            plot(m0(a,c), 0, '^', 'MarkerFaceColor', bar_colors(3,:), 'MarkerEdgeColor', [0 0 0]) % median non-sig
            plot(m1(a,c), 0, '^', 'MarkerFaceColor', bar_colors(2,:), 'MarkerEdgeColor', [0 0 0]) % median sig pos
            plot(m2(a,c), 0, '^', 'MarkerFaceColor', bar_colors(1,:), 'MarkerEdgeColor', [0 0 0]) % median sig neg
            plot(M_median(a,c),  0, '^', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])         % median overall
            
            xlabel('Lag for Most Effective Correlation, # cardiac cycles')
            ylabel('N units')
            title([T ': ' L])
            xlim([cfg.correlation.lag_list(1) cfg.correlation.lag_list(end)])
            
            if splot_num == 1
                legend({'Sig.Neg.', 'Sig.Pos.', 'Non-Sig.', 'Median non-sig.', 'Median sig.neg.', 'Median sig.pos.', 'Median only sig.'}, 'Location', 'best')
            end
        end
    end
    save_figure_as(f9, 'Histograms_MaxAbsCCs_AllAreas', basepath_to_save, 1)
    save_figure_as(f10, 'Histograms_MaxAbsLags_AllAreas', basepath_to_save, 1)
    
    % create table with median lags overall and sign test results
    [h_signtest_lag, p_crit_signtest_lag] = fdr_bky(p_signtest_lag); % correct for multiple comparisons with FDR-BKY
    filename = 'Median_Lags_Only_Sig';
    T = table(unqTargets', ...
        M_median(:,1), M_median(:,2), ...
        p_signtest_lag(:,1), p_signtest_lag(:,2), ...
        h_signtest_lag(:,1), h_signtest_lag(:,2), ...
        [p_crit_signtest_lag; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: Median Lag', 'Task: Median Lag', ...
        'Rest: Sign Test p-value', 'Task: Sign Test p-value', ...
        'Rest: Sign Test h', 'Task: Sign Test h', ...
        'FDR crit. p'});
    writetable(T, [basepath_to_save filesep filename '.xls'])
    clear T
    
    % create table with median lags by subgroup
    [h_sign_pos, sign_pos_crit_p] = fdr_bky(round(p_sign_pos,10));
    [h_sign_neg, sign_neg_crit_p] = fdr_bky(round(p_sign_neg,10));
    
    filename = 'Median_Lags';
    T = table(unqTargets', ...
        m1(:,1), m1(:,2), h_sign_neg(:,1), h_sign_neg(:,2), p_sign_neg(:,1), p_sign_neg(:,2), [sign_neg_crit_p; NaN; NaN], ...
        m2(:,1), m2(:,2), h_sign_pos(:,1), h_sign_pos(:,2), p_sign_pos(:,1), p_sign_pos(:,2), [sign_pos_crit_p; NaN; NaN], ...
        'VariableNames', ...
        {'Target Area',...
        'Rest: Median Lag, Neg. Corr.', 'Task: Median Lag, Neg. Corr.', 'Rest: h neg', 'Task: h neg', 'Rest: p neg', 'Task: p neg', 'Crit p neg', ...
        'Rest: Median Lag, Pos. Corr.', 'Task: Median Lag, Pos. Corr.', 'Rest: h pos', 'Task: h pos', 'Rest: p pos', 'Task: p pos', 'Crit p pos'});
    writetable(T, [basepath_to_save filesep filename '.xls'])
    clear T
    
    % create table with median peak ccs Mann-Whitney test results
    [h_ccs_mann_whitney, p_crit_ccs_mann_whitney] = fdr_bky(p_ccs_mann_whitney); % correct for multiple comparisons with FDR-BKY
    filename = 'Median_CCs_Pos_Neg';
    T = table(unqTargets', ...
        M_median_ccs_neg(:,1), M_median_ccs_neg(:,2), ...
        M_median_ccs_pos(:,1), M_median_ccs_pos(:,2), ...
        p_ccs_mann_whitney(:,1), p_ccs_mann_whitney(:,2), ...
        h_ccs_mann_whitney(:,1), h_ccs_mann_whitney(:,2), ...
        [p_crit_ccs_mann_whitney; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: Median Neg CC', 'Task: Median Neg CC', ...
        'Rest: Median Pos CC', 'Task: Median Pos CC', ...
        'Rest: MW Test p-value', 'Task: MW Test p-value', ...
        'Rest: MW Test h', 'Task: MW Test h', ...
        'FDR crit. p'});
    writetable(T, [basepath_to_save filesep filename '.xls'])
    clear T
    
    % plot histograms with unit percentages
    for c = 1:N_conditions
        
        figure(f11);
        subplot(1, N_conditions, c)
        hold on
        box on
        b = bar(100 * [pos_lag_prc(:,c) neg_lag_prc(:,c) lag_zero_prc(:,c)], 'stacked');
        for ii = 1:size(bar_colors,1)
            b(ii).FaceColor = pos_neg_lag_colors(ii,:);
        end
        clear b
        set(gca, 'XTick', [1 2 3], 'XTickLabel', unqTargets)
        legend({'Lag > 0', 'Lag < 0', 'Lag = 0 or non-sig.'},'Location','southoutside')
        title(cfg.condition(c).name)
        ylim([0 100])
        
        if c == 1
            ylabel('Proportion of Units, %')
        end
        
    end
    save_figure_as(f11, 'Histograms_Pos_Neg_Lags', basepath_to_save, 1)
    
    % compute Fisher's test to know significant differences in unit counts
    % between rest and task
    for a = 1: N_Areas
        T=unqTargets{a};
        
        disp(T)
        
        p_fishertest_lag(a,:) = ...
            ecg_bna_fisher_test([pos_lag_counts(a,1) neg_lag_counts(a,1) lag_zero_counts(a,1)], ...
            [pos_lag_counts(a,2) neg_lag_counts(a,2) lag_zero_counts(a,2)]);
        
    end
    
    % find corrected p-values
    [h_fishertest_lag, p_crit_fishertest_lag] = fdr_bky(round(p_fishertest_lag,10), 0.05);
    
    filename = 'Fisher_Test_posLag_negLag';
    TBL = table(unqTargets', ...
        p_fishertest_lag(:,1), p_fishertest_lag(:,2), ...
        h_fishertest_lag(:,1), h_fishertest_lag(:,2), [p_crit_fishertest_lag; NaN; NaN], ...
        'VariableNames', {'Brain Area', 'Fisher''s p, pos lag', 'Fisher''s p, neg lag', ...
        'Fisher''s h, pos lag', 'Fisher''s h, neg lag', 'FDR-BKY p crit.'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL
    
    % create table with Fisher's test results checking for contingency
    % between FR-RR correlation sign and lag sign
    p_fisher_corrSign_lagSign = round(p_fisher_corrSign_lagSign,10);
    [h_fisher, p_crit_fisher] = fdr_bky(p_fisher_corrSign_lagSign); % correct for multiple comparisons
    filename = 'Fisher_Test_corrSign_lagSign';
    TBL = table(unqTargets', ...
        h_fisher(:,1), p_fisher_corrSign_lagSign(:,1), oddsratio_fisher_corrSign_lagSign(:,1), ...
        h_fisher(:,2), p_fisher_corrSign_lagSign(:,2), oddsratio_fisher_corrSign_lagSign(:,2), [p_crit_fisher; NaN; NaN], ...
        'VariableNames', {'Target Area', ...
        'Rest: Fisher''s h', 'Rest: Fisher''s p', 'Rest: Fisher''s Odds Ratio', ...
        'Task: Fisher''s h', 'Fisher''s p', 'Task: Fisher''s Odds Ratio', 'FDR Crit. p'});
    writetable(TBL, [basepath_to_save filesep filename '.xls'])
    clear TBL h_fisher p_crit_fisher
    
    %% compare median lags between conditions
    % preallocate variables
    [h_ab_neg, p_ab_neg, h_ab_pos, p_ab_pos, p_mw_neg, h_mw_neg, p_mw_pos, h_mw_pos] = deal(nan(1,3));
    
    for a = 1: N_Areas
        T=unqTargets{a};
        
        R_rest=dt.(T).Rest.('pearson_r');
        R_task=dt.(T).Task.('pearson_r');
        
        P_rest=dt.(T).Rest.('pearson_p');
        P_task=dt.(T).Task.('pearson_p');
        
        [~, ids_rest]  = max(abs(R_rest), [], 1); % find max abs ccs; I will need 'ids' to figure out significance of corresponding ccs
        [~, ids_task]  = max(abs(R_task), [], 1);
        
        ccs_rest=R_rest(sub2ind(size(R_rest),ids_rest,1:numel(ids_rest)));
        ccs_task=R_task(sub2ind(size(R_task),ids_task,1:numel(ids_task)));
        
        invalid_rest=isnan(ccs_rest);
        invalid_task=isnan(ccs_task);
        
        P_rest=P_rest(:,~invalid_rest);
        R_rest=R_rest(:,~invalid_rest);
        %     ccs_rest=ccs_rest(~invalid_rest);
        ids_rest=ids_rest(~invalid_rest);
        
        P_task=P_task(:,~invalid_task);
        R_task=R_task(:,~invalid_task);
        %     ccs_task=ccs_task(~invalid_task);
        ids_task=ids_task(~invalid_task);
        
        P_rest       = P_rest(ceil(end/2),:);
        R_rest       = R_rest(ceil(end/2),:);
        sig_idx_rest = fdr_bky(P_rest, 0.05); % figure out significant ones
        sig_rest     = sign(R_rest).*sig_idx_rest;
        
        P_task       = P_task(ceil(end/2),:);
        R_task       = R_task(ceil(end/2),:);
        sig_idx_task = fdr_bky(P_task, 0.05); % figure out significant ones
        sig_task     = sign(R_task).*sig_idx_task;
        
        laglist=dt.(T).cc_lag_list(1,:);
        %     N=numel(ids);
        
        h1_rest=hist(laglist(ids_rest(sig_rest==-1)),laglist);
        h2_rest=hist(laglist(ids_rest(sig_rest==1)),laglist);
        
        h1_task=hist(laglist(ids_task(sig_task==-1)),laglist);
        h2_task=hist(laglist(ids_task(sig_task==1)),laglist);
        
        % compute Mann-Whitney Test to know whether there are differences in
        % median lags
        [p_mw_sig(a),h_mw_sig(a)] = ranksum(laglist(ids_rest(sig_rest~=0)), laglist(ids_task(sig_task~=0)));
        if sum(sig_rest==-1) & sum(sig_task==-1)
            [p_mw_neg(a),h_mw_neg(a)] = ranksum(laglist(ids_rest(sig_rest==-1)), laglist(ids_task(sig_task==-1)));
        end
        if sum(sig_rest==1) & sum(sig_task==1)
            [p_mw_pos(a),h_mw_pos(a)] = ranksum(laglist(ids_rest(sig_rest==1)), laglist(ids_task(sig_task==1)));
        end
        
        % compute Ansari-Brandley Test to knwo whether there are differences in
        % lag variances
        [h_ab_sig(a),p_ab_sig(a)] = ansaribradley(laglist(ids_rest(sig_rest~=0)), laglist(ids_task(sig_task~=0)));
        if sum(sig_rest==-1) > 1 & sum(sig_task==-1) > 1
            [h_ab_neg(a),p_ab_neg(a)] = ansaribradley(laglist(ids_rest(sig_rest==-1)), laglist(ids_task(sig_task==-1)));
        end
        if sum(sig_rest==1) > 1
            [h_ab_pos(a),p_ab_pos(a)] = ansaribradley(laglist(ids_rest(sig_rest==1)), laglist(ids_task(sig_task==1)));
        end
        
        
    end
    
    [h_mw_sig_cor,p_mw_sig_crit] = fdr_bky(p_mw_sig,0.05);
    [h_ab_sig_cor,p_ab_sig_crit] = fdr_bky(p_ab_sig,0.05);
    
    % create table with results of median differences and significant variances
    % for only significant lags
    filename = 'Differences_Between_Medians_Mann-Whitney_Variances_Ansari-Bradley_Only_Sig';
    T = table(unqTargets', ...
        p_mw_sig', h_mw_sig_cor', [p_mw_sig_crit; NaN; NaN], ...
        p_ab_sig', h_ab_sig_cor', [p_ab_sig_crit; NaN; NaN], ...
        'VariableNames', ...
        {'Target Area', ...
        'Mann-Whitney P-value', 'Mann-Whitney h', 'MW: FDR Crit. p', ...
        'Ansari-Bradley P-value', 'Ansari-Bradley h', 'AB: FDR Crit. p'});
    writetable(T, [basepath_to_save filesep filename '.xls'])
    clear T
    
    [h_mw_cor,p_mw_crit] = fdr_bky([p_mw_neg; p_mw_pos],0.05);
    [h_ab_cor,p_ab_crit] = fdr_bky([p_ab_neg; p_ab_pos],0.05);
    
    % create table with results of median differences
    filename = 'Differences_Between_Medians_Mann-Whitney';
    T = table(unqTargets', p_mw_neg', h_mw_cor(1,:)', p_mw_pos', h_mw_cor(2,:)', [p_mw_crit; NaN; NaN], ...
        'VariableNames', ...
        {'Target Area','Sig.Neg.: P-value', 'Sig.Neg.: h', ...
        'Sig.Pos.: P-value', 'Sig.Pos.: h', 'Crit. p'});
    writetable(T, [basepath_to_save filesep filename '.xls'])
    clear T
    
    % create table with results of variance differences
    filename = 'Differences_Between_Variances_Ansari-Bradley';
    T = table(unqTargets', p_ab_neg', h_ab_cor(1,:)', p_ab_pos', h_ab_cor(2,:)', [p_ab_crit; NaN; NaN], ...
        'VariableNames', ...
        {'Target Area','Sig.Neg.: P-value', 'Sig.Neg.: h', ...
        'Sig.Pos.: P-value', 'Sig.Pos.: h', 'Crit. p'});
    writetable(T, [basepath_to_save filesep filename '.xls'])
    clear T
    
    %% lag scatters task versus rest
    f2 = figure;
    set(gcf,'Position',[281 138 1347 843])
    
    % preallocate
    [h_ab,p_ab,p_wilcox,h_mann_whitney] = deal(nan(1,3));
    
    for a = 1: N_Areas
        
        subplot(3, N_Areas, a)
        
        T=unqTargets{a};
        L1=cfg.condition(1).name;
        L2=cfg.condition(2).name;
        
        [rc, lag_rest] = max(abs(dt.(T).(L1).pearson_r));
        [tc, lag_task] = max(abs(dt.(T).(L2).pearson_r));
        invalid=isnan(rc)|isnan(tc);
        lag_rest(invalid)=[];
        lag_task(invalid)=[];
        
        % compute correlation
        tmp=[dt.(T).cc_lag_list(1,lag_rest);dt.(T).cc_lag_list(1,lag_task)]';
        [rtaskrest,ptaskrest] = corrcoef(tmp(:,1),tmp(:,2));
        
        % compute (['SNR_' cfg.condition(c).saccadeTask]) paired test
        [p_wilcox(a),h_mann_whitney(a)] = signrank(tmp(:,1), tmp(:,2));
        
        % compute Ansari-Brandley Test to know whether there are differences in
        % lag variances
        [h_ab(a),p_ab(a)] = ansaribradley(tmp(:,1), tmp(:,2));
        
        % plot lag scatter
        scatter(dt.(T).cc_lag_list(1,lag_rest),dt.(T).cc_lag_list(1,lag_task), 20, cfg.area_colors{a}, 'filled', 'o');
        xlim([cfg.correlation.lag_list(1) cfg.correlation.lag_list(end)])
        ylim([cfg.correlation.lag_list(1) cfg.correlation.lag_list(end)])
        
        xlabel('Rest: Max. CC lag')
        ylabel('Task: Max. CC lag')
        title([T ', r=' num2str(rtaskrest(2,1)) ', p=' num2str(ptaskrest(2,1))]);
        box on
        axis square
        
        subplot(3, N_Areas, 3+a)
        hist(tmp(:,1));
        h = findobj(gca,'Type','patch');
        h.FaceColor = cfg.condition(1).color;
        hold on
        plot(mean(tmp(:,1)),0,'^','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',cfg.condition(1).color)
        title({'Rest: Marginal Histogram',['Wilcox p = ' num2str(p_wilcox(a))]})
        xlim([cfg.correlation.lag_list(1) cfg.correlation.lag_list(end)])
        if a == 1
            xlabel('Lag for Most Efficient Correlation')
            ylabel('Unit Counts')
        end
        
        subplot(3, N_Areas, 6+a)
        hist(tmp(:,2))
        h = findobj(gca,'Type','patch');
        h.FaceColor = cfg.condition(2).color;
        hold on
        plot(mean(tmp(:,2)),0,'^','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',cfg.condition(2).color)
        title('Task: Marginal Histogram')
        xlim([cfg.correlation.lag_list(1) cfg.correlation.lag_list(end)])
        if a == 1
            xlabel('Lag for Most Efficient Correlation')
            ylabel('Unit Counts')
        end
        
    end
    save_figure_as(f2, ['Scatterhist_Lags_Rest_vs_Task'], basepath_to_save, 1)
    
    %% CC histograms per area, for all lags and conditions
    for a = 1: N_Areas
        
        figure,
        set(gcf, 'Position', [22 441 1857 491])
        colormap(bar_colors)
        bin_resolution=0.1;
        T=unqTargets{a};
        
        for c=1:N_conditions
            N_lags = size(dt.(T).cc_lag_list,2);
            L=cfg.condition(c).name;
            y_lim=1;
            
            for lag_num = 1:N_lags
                curr_cc = dt.(T).(L).pearson_r(lag_num, :);
                curr_pp = dt.(T).(L).pearson_p(lag_num, :);
                [~, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
                pos_idx = curr_cc > 0;
                neg_idx = curr_cc < 0;
                
                % prepare data for plotting
                counts_pos    = histc(curr_cc(sig_idx & pos_idx)', -1:bin_resolution:1); % significant and positive
                counts_neg    = histc(curr_cc(sig_idx & neg_idx)', -1:bin_resolution:1); % significant and negative
                counts_nonsig = histc(curr_cc(~sig_idx)', -1:bin_resolution:1);
                counts_pos    = counts_pos(:);
                counts_neg    = counts_neg(:);
                counts_nonsig = counts_nonsig(:);
                
                % plot
                splot_num = (c-1)*N_lags + lag_num;
                plot_data = [counts_pos counts_neg counts_nonsig];
                subplot(N_conditions, N_lags, splot_num)
                b = bar(-1:bin_resolution:1, plot_data, 'stacked');
                for ii = 1:size(bar_colors,1)
                    b(ii).FaceColor = bar_colors(ii,:);
                end
                %             set(b, 'LineWidth',0.000000000000000000000000000000000001);
                xlim([-0.5 0.5])
                title(num2str(dt.(T).cc_lag_list(lag_num)))
                if lag_num > 1
                    set(gca,'yticklabel',[]);
                else
                    ylabel([T ': ' L ': N neurons']);
                end
                set(gca,'xtick',0,'xticklabel',0);
                xlabel('CC');
                if c == 1 && lag_num == 1
                    legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'})
                end
                
                y_lim=max(y_lim,max(sum(plot_data,2)));
            end
            
            for lag_num = 1:N_lags
                splot_num = (c-1)*N_lags + lag_num;
                subplot(N_conditions, N_lags, splot_num)
                set(gca,'ylim',[0,y_lim]);
            end
        end
        save_figure_as(gcf, ['Histograms_CC_FR_&_RR_' T], basepath_to_save, 1)
    end
    
    %% plot number of cardiac cycles used at lag 0
    figure,
    set(gcf, 'Position', [667 519 930 477])
    for a = 1: N_Areas
        T=unqTargets{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            N_cycles=dt.(T).(L).n_cycles;
            laglist=dt.(T).cc_lag_list;
            zero_lag_idx=find(laglist==0);
            N_cycles=N_cycles(zero_lag_idx);
            N_cycles(isnan(N_cycles))=[];
            splot_num = (c-1)*N_Areas + a;
            subplot(N_conditions,N_Areas,splot_num)
            hist(N_cycles, 20,'o-', 'Color', [0.5 0.5 0.5])
            title([T ': ' L ', min: ' num2str(min(N_cycles))])
        end
    end
    save_figure_as(gcf, 'RR_Counts_lag_0', basepath_to_save, 1)
    
    
    %% plot number of cardiac cycles used
    for prefixNum = 1:length(var_prefix)
        figure,
        set(gcf, 'Position', [667 519 930 477])
        for a = 1: N_Areas
            T=unqTargets{a};
            for c=1:N_conditions
                L=cfg.condition(c).name;
                N_cycles=dt.(T).(L).([var_prefix{prefixNum} 'n_cycles']);
                N_cycles(:,all(isnan(N_cycles),1))=[];
                splot_num = (c-1)*N_Areas + a;
                subplot(N_conditions,N_Areas,splot_num)
                plot(cfg.correlation.lag_list, 100 * N_cycles ./ repmat(N_cycles(ceil(end/2),:),size(N_cycles,1),1), 'o-', 'Color', [0.5 0.5 0.5])
                title([T ': ' L])
            end
        end
        
        subplot(N_conditions,N_Areas,1)
        xlabel('Lag, # Cardiac Cycles')
        ylabel('% of RRs Taken into Analysis')
        
        save_figure_as(gcf, [var_prefix{prefixNum} 'Percentages_RR_Counts_' T], basepath_to_save, 1)
    end
    clear Out dt
end
end

function save_figure_as(fig_id, filename,basepath_to_save,savePlot)
if savePlot;
    export_fig(fig_id,[basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(fig_id)
end
end