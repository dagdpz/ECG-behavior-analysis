function ecg_bna_avg_spike_histogram_clean(cfg, data_folder, unitList)
% Here comes some sort of across population plot i assume?

%% load data
load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG\' unitList '.mat'], 'unit_ids', 'targets', 'ids_both')

var_list = ...
    {'SD', 'SD_STD', 'SD_SEM', 'SDP', 'SDPCL', 'SDPCu', ...
    'sig_all', 'sig', 'sig_FR_diff', 'sig_time', 'sig_n_bins', 'sig_sign', ...
    'NrTrials', 'NrEvents', 'FR', 'SDsubstractedSDP', 'SDsubstractedSDP_normalized', 'FR_ModIndex_SubtrSDP', 'FR_ModIndex_PcS'};
% , ...
%     'lowIBI.SD', 'lowIBI.SD_STD', 'lowIBI.SD_SEM', 'lowIBI.SDP', 'lowIBI.SDPCL', 'lowIBI.SDPCu', ...
%     'lowIBI.sig_all', 'lowIBI.sig', 'lowIBI.sig_FR_diff', 'lowIBI.sig_time', 'lowIBI.sig_n_bins', 'lowIBI.sig_sign', ...
%     'lowIBI.NrEvents', 'lowIBI.SDsubstractedSDP', 'lowIBI.SDsubstractedSDP_normalized', 'lowIBI.FR_ModIndex_SubtrSDP', 'lowIBI.FR_ModIndex_PcS', ...
%     'highIBI.SD', 'highIBI.SD_STD', 'highIBI.SD_SEM', 'highIBI.SDP', 'highIBI.SDPCL', 'highIBI.SDPCu', ...
%     'highIBI.sig_all', 'highIBI.sig', 'highIBI.sig_FR_diff', 'highIBI.sig_time', 'highIBI.sig_n_bins', 'highIBI.sig_sign', ...
%     'highIBI.NrEvents', 'highIBI.SDsubstractedSDP', 'highIBI.SDsubstractedSDP_normalized', 'highIBI.FR_ModIndex_SubtrSDP', 'highIBI.FR_ModIndex_PcS'

BINS = linspace(cfg.analyse_states{1}{3},cfg.analyse_states{1}{4}, 101);

for targetGrNum = 1:length(cfg.targets_spike_data)
    
    unqTargets = cfg.targets_spike_data{targetGrNum};
    N_Areas    = length(unqTargets);
    N_conditions=numel(cfg.condition);
    n_sig_bins = cfg.time.n_sig_bins;
    
    % output folder
    targSuff = [];
    for targNum = 1:N_Areas
        targSuff = [targSuff '_' unqTargets{targNum}];
    end
    output_folder = [cfg.SPK_root_results_fldr filesep 'Population_time_domain' unitList(9:end) targSuff];
    if ~exist(output_folder,'dir')
        mkdir(output_folder)
    end
    
    dataset_name = [output_folder filesep 'Output.mat'];
    
    if ~exist(dataset_name,'file')
        for a = 1: N_Areas
            
            T=unqTargets{a};
            currTargIds   = cellfun(@(x) strcmp(x, unqTargets{a}), targets);
            curr_unit_ids = unit_ids(currTargIds);
            curr_ids_both = ids_both(currTargIds);
            Out.(T)       = ecg_bna_load_variables(cfg, curr_unit_ids, data_folder, 'Output', var_list, curr_ids_both);
            
            for c=1:N_conditions
                L=cfg.condition(c).name;
                Out.(T).(L).FR_perECGTriggeredAverage = nanmean(Out.(T).(L).SD,1);
                
                %% compute histograms - unit fractions and counts
                out = [Out.(T).(L)];
                
                Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(end,:))';
                Idx_Units_NaN =  sum(~Idx_Units_NonNaN);
                sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
                
                % decrease, increase, non-sign
                Out.(T).(L).Pc_SignFR = ([sum(out.sig_sign(sig) == -1), sum(out.sig_sign(sig) == 1) ,(sum(~sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
                Out.(T).(L).Nb_SignFR = ([sum(out.sig_sign(sig) == -1), sum(out.sig_sign(sig) == 1) ,(sum(~sig) -Idx_Units_NaN),] ) ;
                
                % [by unit] decrease, increase, non-sign
                %             HeartResponseType_byUnit = out.sig_sign;
                %             HeartResponseType_byUnit(~sig | ~Idx_Units_NonNaN) = 0;
                Out.(T).(L).HeartResponseType_byUnit = out.sig_sign; % increase, decrease, no response
                
                %% compute modulation indices with sign
                curr_sign = Out.(T).(L).sig_sign;
                curr_sign(curr_sign == 0) = 1; % to define real sign but not to lose non-significant ones
                Out.(T).(L).FR_ModIndex_SubtrSDP_signed = Out.(T).(L).FR_ModIndex_SubtrSDP .* curr_sign;
                Out.(T).(L).FR_ModIndex_PcS_signed      = Out.(T).(L).FR_ModIndex_PcS .* curr_sign;
                
                %% compute modulation directionality index
                Out.(T).(L).Mod_Directionality_Index = ...
                    (abs(max(Out.(T).(L).SDsubstractedSDP)) - abs(min(Out.(T).(L).SDsubstractedSDP))) ./ (abs(max(Out.(T).(L).SDsubstractedSDP)) + abs(min(Out.(T).(L).SDsubstractedSDP)));
                [tmp_max, max_id] = max(Out.(T).(L).SDsubstractedSDP);
                [tmp_min, min_id] = min(Out.(T).(L).SDsubstractedSDP);
                Out.(T).(L).RespMag_at_Max  = abs(tmp_max);
                Out.(T).(L).RespMag_at_Min  = abs(tmp_min);
                Out.(T).(L).RespTime_at_Max = BINS(max_id);
                Out.(T).(L).RespTime_at_Min = BINS(min_id);
                
            end
            
            % find shift of the max and shift of the min between conditions
            Out.(T).Shift_of_Max = Out.(T).Rest.RespTime_at_Max - Out.(T).Task.RespTime_at_Max;
            Out.(T).Shift_of_Min = Out.(T).Rest.RespTime_at_Min - Out.(T).Task.RespTime_at_Min;
            
            % preallocate
            [Out.(T).pp_rest_vs_task, Out.(T).cc_rest_vs_task] = deal(nan(length(Out.(T).unit_ID),1));
            for u = 1:length(Out.(T).unit_ID)
                
                if ~isnan(Out.(T).Rest.sig_FR_diff(u)) & ~isnan(Out.(T).Task.sig_FR_diff(u))
                    
                    [Out.(T).pp_rest_vs_task(u), Out.(T).cc_rest_vs_task(u)] = ...
                        mult_comp_perm_corr(Out.(T).Rest.SD(:,u), Out.(T).Task.SD(:,u), cfg.time.n_shuffles, cfg.time.tail, cfg.time.alpha_level, cfg.time.stat, cfg.time.reports, cfg.time.seed_state);
                    
                end
                
            end
            
        end
        save(dataset_name,'Out')
    else
        load(dataset_name,'Out')
    end
    
    bar_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 1 1 1];
    bar_colors_merged = [0 0 0.5; 0.7 0.7 0.7];
    
    savePlot = 1;
    saveTable= 0;
    OnlyUnits_withRestANDTask = 1;
    Graph_SelectionCriterion = 1;
    % colors = distinguishable_colors(25);
    curr_analyse_states = cfg.analyse_states{1};
    PSTH_bins = (curr_analyse_states{1,3}:cfg.time.PSTH_binwidth:curr_analyse_states{1,4})*1000;
    
    if ~exist(output_folder,'dir')
        mkdir(output_folder);
    end
    
    ConNames = {cfg.condition(:).name};
    
    lineProps = cell(length(cfg.condition), 1);
    for conNum = 1:length(cfg.condition)
        lineProps{conNum} = {'color',cfg.condition(conNum).color,'linewidth',4};
    end
    
    TargetBrainArea = fieldnames(Out);
    Ana_TargetBrainArea = TargetBrainArea;
    
    % N_Areas=numel(Ana_TargetBrainArea);
    
    Color_BrainArea = cfg.area_colors;
    
    % % create unit subset
    % for a = 1: N_Areas
    %     T=Ana_TargetBrainArea{a};
    %     for c=1:N_conditions
    %         L=cfg.condition(c).name;
    %
    %         Idx_Units_NonNaN_rest = ~isnan(Out.(T).Rest.SDsubstractedSDP(:,end)); % figure out units with nan responses
    %         Idx_Units_NonNaN_task = ~isnan(Out.(T).Task.SDsubstractedSDP(:,end));
    %
    %         sig_rest = ~isnan(Out.(T).Rest.sig_FR_diff) & (Out.(T).Rest.sig_n_bins > n_sig_bins);
    %         sig_task = ~isnan(Out.(T).Task.sig_FR_diff) & (Out.(T).Task.sig_n_bins > n_sig_bins);
    %
    %         unit_ids = ...
    %             Out.(T).Rest.unit_ID(Idx_Units_NonNaN_rest & Idx_Units_NonNaN_task & ...
    %             (sig_rest | sig_task));
    %         targets  = [];
    %
    %         if a == N_Areas
    %             save()
    %         end
    %     end
    % end
    
    %% plot histograms of modulation directionality indices
    f1 = figure;
    set(f1,'Position',[68 286 1422 695])
    
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            
            out = Out.(T).(L);
            
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            %         dec = (out.sig_sign == -1);
            %         inc = (out.sig_sign == 1);
            
            counts_sig    = histc(out.Mod_Directionality_Index(sig),-1:0.05:1);
            counts_nonsig = histc(out.Mod_Directionality_Index(~sig),-1:0.05:1);
            
            subplot(N_conditions,N_Areas,N_Areas*(c-1)+a)
            %         hist(out.Mod_Directionality_Index)
            bar(-1:0.05:1,[counts_sig; counts_nonsig]','stacked')
            title([T ': ' L])
            xlim([-1 1])
            
            if c == 1 && a == 1
                legend({'Sig.','Non-Sig.'},'Location','best')
                xlabel('Modulation Sign Index')
                ylabel('Unit Counts')
            end
            
        end
    end
    save_figure_as('Sign_Index_Histograms',output_folder,savePlot)
    
    %% stability of modulation directionality indices
    
    f1 = figure;
    set(f1,'Position',[3 268 1918 700])
    
    L1=cfg.condition(1).name;
    L2=cfg.condition(2).name;
    
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        
        out1 = Out.(T).(L1);
        out2 = Out.(T).(L2);
        
        sig1 = ~isnan(out1.sig_FR_diff) & (out1.sig_n_bins > n_sig_bins) ;
        sig2 = ~isnan(out2.sig_FR_diff) & (out2.sig_n_bins > n_sig_bins) ;
        
        nonsig1 = ~isnan(out1.sig_FR_diff) & (out1.sig_n_bins < n_sig_bins) ;
        nonsig2 = ~isnan(out2.sig_FR_diff) & (out2.sig_n_bins < n_sig_bins) ;
        
        sig_both = ~isnan(out1.sig_FR_diff) & ~isnan(out2.sig_FR_diff) & (out1.sig_n_bins > n_sig_bins | out2.sig_n_bins > n_sig_bins);
        
        counts_sig1 = histc(out1.Mod_Directionality_Index(sig1),-1:0.05:1);
        counts_sig2 = histc(out2.Mod_Directionality_Index(sig2),-1:0.05:1);
        
        subplot(1,N_Areas,a)
        hold on
        box on
        % plot all not nans
        p_nonsig1 = plot(1, out1.Mod_Directionality_Index(nonsig1),'o','MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor',cfg.condition(1).color);
        if sum(nonsig2)
            p_nonsig2 = plot(2, out2.Mod_Directionality_Index(nonsig2),'o','MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor',cfg.condition(2).color);
        end
        % plot all significant
        if sum(sig1)
            p_sig1 = plot(1, out1.Mod_Directionality_Index(sig1),'o','MarkerSize',7,'MarkerFaceColor',cfg.condition(1).color,'MarkerEdgeColor',cfg.condition(1).color);
        end
        if sum(sig2)
            p_sig2 = plot(2, out2.Mod_Directionality_Index(sig2),'o','MarkerSize',7,'MarkerFaceColor',cfg.condition(2).color,'MarkerEdgeColor',cfg.condition(2).color);
        end
        % plot significant at least in one condition
        if sum(sig_both)
            p_both = plot([1 2],[out1.Mod_Directionality_Index(sig_both); out2.Mod_Directionality_Index(sig_both)],'o-','Color',[1 0 1]/2);
        end
        
        xlim([0 3])
        
        set(gca,'XTick',[1 2],'XTickLabel',{'Rest','Task'})
        
        % marginal histogram for task
        fill1 = fill([2.1*ones(1,41) 2.1+counts_sig2/sum(counts_sig2)],[fliplr(-1:0.05:1) -1:0.05:1],'r','EdgeColor','none','FaceAlpha',0.3);
        % marginal histogram for rest
        fill2 = fill([0.9*ones(1,41) 0.9-counts_sig1/sum(counts_sig1)],[fliplr(-1:0.05:1) -1:0.05:1],'b','EdgeColor','none','FaceAlpha',0.3);
        
        title(T)
        
        if a == 1
            
            xlabel('Condition')
            ylabel('Modulation Sign Index')
            legend([p_nonsig1(1) p_nonsig2(1) p_sig1(1) p_sig2(1) p_both(1) fill1(1) fill2(1)], ...
                {'Rest: any', 'Task: any', 'Rest: Sig.', 'Task: Sig.', 'Sig. Any Condition', 'Rest: Histogram', 'Task: Histogram'}, ...
                'Location','best')
            
        end
        
    end
    
    save_figure_as('Sign_Index_Marginal_Histograms',output_folder,savePlot)
    
    %% Decrease and Increase grouped for brain region
    
    maxmin={'Max','Min'};
    symbolstoplot={'o','s','v'};
    for groupNum = 1:length(cfg.spk.compare_conditions)
        
        cond1_num = cfg.spk.compare_conditions{groupNum}(1);
        cond2_num = cfg.spk.compare_conditions{groupNum}(2);
        L1=cfg.condition(cond1_num).name;
        L2=cfg.condition(cond2_num).name;
        
        %m=colormap(bluewhitered);
        
        logscale_127=(0:126)'/127;
        colneg=[0 0 255];
        colpos=[255 0 0];
        
        m=[logscale_127*([255 255 255]-colneg)+repmat(colneg,size(logscale_127)); 255 255 255;...
            flipud(logscale_127*([255 255 255]-colpos)+repmat(colpos,size(logscale_127)))]/255;
        
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            
            f1 = figure('Name',sprintf(['ConditionsPeakComparison_' T]),'Position',[317 266 1400 585],'PaperPositionMode', 'auto');
            
            for maxmini=1:2
                subplot(1,2,maxmini);
                hold on
                
                toplot1=Out.(T).(L1).(['RespTime_at_' maxmin{maxmini}]);
                toplot2=Out.(T).(L2).(['RespTime_at_' maxmin{maxmini}]);
                cols=mean([Out.(T).(L1).Mod_Directionality_Index; Out.(T).(L1).Mod_Directionality_Index],1);
                colidx=round((cols+1)*(size(m,1)-1)/2)+1;
                sigidx=Out.(T).(L1).sig_n_bins>cfg.time.n_sig_bins | Out.(T).(L2).sig_n_bins>cfg.time.n_sig_bins;
                ix=~isnan(cols) & ~sigidx';
                colstoplot=m(colidx(ix),:);
                
                if any(ix)
                    sc=scatter(toplot1(ix),toplot2(ix),40,colstoplot,'marker',symbolstoplot{a});
                end
                ix=~isnan(cols) & sigidx';
                colstoplot=m(colidx(ix),:);
                if any(ix)
                    sc=scatter(toplot1(ix),toplot2(ix),40,colstoplot,'filled');
                    set(sc,'marker',symbolstoplot{a})
                end
                xlabel(['peak in ' L1 ', ms']);
                ylabel(['peak in ' L2 ', ms']);
                title(maxmin{maxmini})
                box on
                axis square
            end
            save_figure_as(['ConditionsPeakComparison_' T],output_folder,savePlot)
        end
    end
    
    %% percentages of units with different responsivenes - pos. / neg. / non-sig.
    figure('Name',sprintf('BarPlot_Pc'),'Position',[846 552 600 239],'PaperPositionMode', 'auto');
    for c=1:N_conditions
        L=cfg.condition(c).name;
        % prepare data matrix for plotting
        pc_mat = [];
        nb_mat = [];
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            pc_mat = vertcat(pc_mat, Out.(T).(L).Pc_SignFR);
            nb_mat = vertcat(nb_mat, Out.(T).(L).Nb_SignFR);
        end
        subplot(1,N_conditions,c);
        b = bar(pc_mat,'stacked');
        for ii = 1:length(b)
            b(ii).FaceColor = bar_colors(ii,:);
        end
        title(L,'interpreter','none');
        set(gca,'XTickLabel',Ana_TargetBrainArea,'fontsize',10);
        ylim([0 100])
        ylabel('Percentage of Units, %')
        % legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')
        ax = gca;
        ax.FontSize = 12;
        
        Tab = table(Ana_TargetBrainArea, nb_mat(:,1), nb_mat(:,2), nb_mat(:,3), 'VariableNames', {'Brain Area', 'Dec', 'Inc', 'No Resp'});
        savename = [output_folder filesep 'UnitCounts_' L '.xlsx'];
        writetable(Tab, savename)
        clear Tab
        
    end
    save_figure_as('Pc_CardiacRelatedUnits',output_folder,savePlot)
    
    %% Fisher's exact test to compare unit fractions between task and rest
    for groupNum = 1:length(cfg.spk.compare_conditions)
        
        cond1_num = cfg.spk.compare_conditions{groupNum}(1);
        cond2_num = cfg.spk.compare_conditions{groupNum}(2);
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            
            disp(T)
            
            p(a,:) = ecg_bna_fisher_test(Out.(T).(cfg.condition(cond2_num).name).Nb_SignFR([2 1 3]), Out.(T).(cfg.condition(cond1_num).name).Nb_SignFR([2 1 3]));
            
        end
        
        % find corrected p-values
        %     [p_corr, h] = bonf_holm(p, 0.05);
        h = fdr_bky(round(p,10), 0.05);
        
        Tab = table(Ana_TargetBrainArea, p(:,1), p(:,2), h(:,1), h(:,2), 'VariableNames', {'Brain Area', 'p_corr inc', 'p_corr dec', 'h inc', 'h dec'});
        savename = [output_folder filesep (cfg.condition(cond2_num).name) '_vs_' (cfg.condition(cond1_num).name) 'Table_Prevalences_pvalues_corrected.xlsx'];
        writetable(Tab, savename)
        clear Tab
    end
    clear p
    
    %% plot Sankey plots
    
    % options
    options.color_map = [0.7 0.7 0.7; 0.8500 0.3250 0.0980; 0 0.4470 0.7410; ...
        0.7 0.7 0.7; 0.8500 0.3250 0.0980; 0 0.4470 0.7410];
    options.flow_transparency = 0.8;   % opacity of the flow paths
    options.bar_width = 100;            % width of the category blocks
    options.show_perc = true;          % show percentage over the blocks
    options.text_color = [0 0 0];      % text color for the percentages
    options.show_layer_labels = true;  % show layer names under the chart
    options.show_cat_labels = true;    % show categories over the blocks.
    options.show_legend = false;       % show legend with the category names.
    % if the data is not a table, then the
    % categories are labeled as catX-layerY
    
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        
        % prepare the categorical data for plotting
        t      = cat(2,~isnan(Out.(T).Rest.sig_FR_diff) & Out.(T).Rest.sig_n_bins > n_sig_bins, ...
            ~isnan(Out.(T).Task.sig_FR_diff) & Out.(T).Task.sig_n_bins > n_sig_bins);
        t_sign = [Out.(T).Rest.HeartResponseType_byUnit, Out.(T).Task.HeartResponseType_byUnit];
        
        t      = t .* t_sign;
        
        nan_ids    = any(isnan(t),2);
        t(nan_ids,:) = [];
        disp(['Number of NaN units: ' num2str(sum(nan_ids))])
        
        tmp = cell(size(t));
        tmp(t == 0)  = deal({'no response'});
        tmp(t == 1)  = deal({'increase'});
        tmp(t == -1) = deal({'decrease'});
        
        tmp = cell2table(tmp,'VariableNames',{'Rest','Task'});
        
        tmp.Rest = categorical(tmp.Rest);
        tmp.Task = categorical(tmp.Task);
        
        tmp = sortrows(tmp,[1,2],'descend');
        
        figure,
        plotSankeyFlowChart(tmp,options)
        title(T)
        save_figure_as(['SankeyPlot_' T],output_folder,savePlot)
    end
    
    %% percentages of R-peak responsive / non-responsive
    figure('Name',sprintf('BarPlot_Pc'),'Position',[846 552 600 239],'PaperPositionMode', 'auto');
    for c=1:N_conditions
        L=cfg.condition(c).name;
        % prepare data matrix for plotting
        pc_mat = [];
        nb_mat = [];
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            pc_mat = vertcat(pc_mat, [Out.(T).(L).Pc_SignFR(1)+Out.(T).(L).Pc_SignFR(2) Out.(T).(L).Pc_SignFR(3)]);
            nb_mat = vertcat(nb_mat, [Out.(T).(L).Nb_SignFR(1)+Out.(T).(L).Nb_SignFR(2) Out.(T).(L).Nb_SignFR(3)]);
        end
        subplot(1,N_conditions,c);
        b = bar(pc_mat,'stacked');
        for ii = 1:length(b)
            b(ii).FaceColor = bar_colors_merged(ii,:);
        end
        title(L,'interpreter','none');
        set(gca,'XTickLabel',Ana_TargetBrainArea,'fontsize',10);
        ylim([0 100])
        ylabel('Percentage of Units, %')
        % legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')
        ax = gca;
        ax.FontSize = 12;
        
        Tab = table(Ana_TargetBrainArea, nb_mat(:,1), nb_mat(:,2), 'VariableNames', {'Brain Area', 'Resp', 'No Resp'});
        savename = [output_folder filesep 'UnitCounts_' L '.xlsx'];
        writetable(Tab, savename)
        clear Tab
        
    end
    save_figure_as('Pc_CardiacRelatedUnits_Merged',output_folder,savePlot)
    
    %% Fisher's exact test to compare unit fractions between task and rest (inc and dec together)
    if contains(unitList,'stable')
        L1=cfg.condition(1).name;
        L2=cfg.condition(2).name;
        
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            
            % create contingency table for Fisher's test
            con_tbl = [Out.(T).(L1).Nb_SignFR(1)+Out.(T).(L1).Nb_SignFR(2) Out.(T).(L2).Nb_SignFR(1)+Out.(T).(L2).Nb_SignFR(2); ...
                Out.(T).(L1).Nb_SignFR(3) Out.(T).(L2).Nb_SignFR(3)];
            
            % compute Fisher's test per area
            [fisher_h(a),fisher_p(a),stats] = fishertest(con_tbl);
            odds_ratio(a) = stats.OddsRatio;
            
            % create and save contingency table
            filename = ['R-peakResponsiveness_Condition_Contingency_' T];
            TBL = table({'R-peak Resp.';'R-peak No Resp.'}, con_tbl(:,1), con_tbl(:,2), ...
                'VariableNames', {'A', 'Rest', 'Task'});
            writetable(TBL, [output_folder filesep filename '.xls'])
            clear TBL con_tbl
            
        end
        
        % correct Fisher's test for multiple comparisons
        fisher_p = round(fisher_p,10);
        [fisher_h, fisher_crit_p] = fdr_bky(fisher_p);
        
        Tab = table(Ana_TargetBrainArea, odds_ratio', fisher_p', fisher_h', [fisher_crit_p; NaN; NaN], 'VariableNames', {'Brain Area', 'Fisher''s Odds Ratio', 'Fisher p', 'Fisher h', 'Crit. p'});
        savename = [output_folder filesep 'R-peak_responsiveness_Fisher_test.xlsx'];
        writetable(Tab, savename)
        clear Tab
    end
    
    %% compute chi2 test - only for stable dataset
    if contains(unitList,'stable')
        
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            
            out1 = [Out.(T).(L1)];
            out2 = [Out.(T).(L2)];
            
            % index for both conditions not nan
            %         idx_both_nonnan = ~isnan(out1.sig_FR_diff) & ~isnan(out2.sig_FR_diff);
            
            % significant in condition1 or condition2
            sig1    = ~isnan(out1.sig_FR_diff) & (out1.sig_n_bins > n_sig_bins) ;
            sig2    = ~isnan(out2.sig_FR_diff) & (out2.sig_n_bins > n_sig_bins) ;
            
            %         nonsig1 = ~isnan(out1.sig_FR_diff) & (out1.sig_n_bins <= n_sig_bins) ;
            %         nonsig2 = ~isnan(out2.sig_FR_diff) & (out2.sig_n_bins <= n_sig_bins) ;
            
            if sum(sig1) || sum(sig2)
                [con_tbl, chi2(a), chi2_p(a)] = crosstab(sig1,sig2);
            
                % create and save cross-tabulation
                filename = ['R-peakResponsiveness_Condition_Cross-Tabulation_' T];
                TBL = table({'Task: R-peak No Resp.';'Task: R-peak Resp.'}, con_tbl(:,1), con_tbl(:,2), ...
                    'VariableNames', {'A', 'Rest: R-peak No Resp.', 'Rest: R-peak Resp.'});
                writetable(TBL, [output_folder filesep filename '.xls'])
                clear TBL con_tbl
            end
            
        end
        
        % correct chi-squared p
        chi2_p = round(chi2_p,10);
        [chi2_h, chi2_crit_p] = fdr_bky(chi2_p);
        
        Tab = table(Ana_TargetBrainArea, chi2', chi2_p', chi2_h', [chi2_crit_p; NaN; NaN], 'VariableNames', {'Brain Area', 'Chi^2', 'Chi^2 p', 'Chi^2 h', 'Crit. p'});
        savename = [output_folder filesep 'Chi-Squared_Results.xlsx'];
        writetable(Tab, savename)
        clear Tab
        
    end
    
    %% plot Sankey plots (for inc and dec together)
    
    % options
    options.color_map = [0.4 0.4 0.4; 0 0 0.5;  ...
        0.4 0.4 0.4; 0 0 0.5];
    options.flow_transparency = 0.8;   % opacity of the flow paths
    options.bar_width = 100;            % width of the category blocks
    options.show_perc = true;          % show percentage over the blocks
    options.text_color = [1 1 1];      % text color for the percentages
    options.show_layer_labels = true;  % show layer names under the chart
    options.show_cat_labels = true;    % show categories over the blocks.
    options.show_legend = false;       % show legend with the category names.
    % if the data is not a table, then the
    % categories are labeled as catX-layerY
    
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        
        % prepare the categorical data for plotting
        %     t = [Out.(T).Rest.HeartResponseType_byUnit, Out.(T).Task.HeartResponseType_byUnit];
        
        t      = cat(2,~isnan(Out.(T).Rest.sig_FR_diff) & Out.(T).Rest.sig_n_bins > n_sig_bins, ...
            ~isnan(Out.(T).Task.sig_FR_diff) & Out.(T).Task.sig_n_bins > n_sig_bins);
        t_sign = [Out.(T).Rest.HeartResponseType_byUnit, Out.(T).Task.HeartResponseType_byUnit];
        
        t      = t .* t_sign;
        
        nan_ids    = any(isnan(t),2);
        t(nan_ids,:) = [];
        
        tmp = cell(size(t));
        tmp(t == 0)  = deal({'no response'});
        tmp(t ~= 0)  = deal({'response'});
        
        tmp = cell2table(tmp,'VariableNames',{'Rest','Task'});
        
        tmp.Rest = categorical(tmp.Rest);
        tmp.Task = categorical(tmp.Task);
        
        tmp = sortrows(tmp,[1,2],'ascend');
        
        figure,
        plotSankeyFlowChart(tmp,options)
        title(T)
        save_figure_as(['SankeyPlot_MergedIncDec_' T],output_folder,savePlot)
    end
    
    %% How much does the surrogate and mean firing diverge
    figure('Name',sprintf('Surrogate_MeanFr'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            out = [Out.(T).(L)];
            
            subplot(2,N_Areas,a);
            hold on;
            if numel(out.FR) > 1 && sum(~isnan(out.FR))
                scatter(out.FR, nanmean(out.SDP,1), 'filled', 'MarkerFaceColor',ccol)  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
                M=max([xlim,ylim]);
                line([0 M],[0 M], 'color','k');
                xlabel('average Firing rate','fontsize',14 );
                ylabel('mean FR of jittered data','fontsize',14 );
                axis square; box on;
                title(T);
            end
            
            subplot(2,N_Areas,a + N_Areas);
            hold on;
            if size(out.SD, 1) > 1  && sum(~isnan(nanmean(out.SD,2)))
                scatter(nanmean(out.SD,1) , nanmean(out.SDP,1), 'filled', 'MarkerFaceColor',ccol)  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
                M=max([xlim,ylim]);
                line([0 M],[0 M], 'color','k');
                xlabel('average Firing rate within analysis window','fontsize',14 );
                ylabel('mean FR of jittered data','fontsize',14 );
                axis square;box on;
                title(T);
            end
        end
    end
    save_figure_as('MeanSurrogate_MeanFr',output_folder,savePlot)
    
    %% Example for
    for a = 1:N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            figure('Name',sprintf(T),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
            
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(2,4,1);
            title([T ' ' L]);
            xlabel('Signal-to-Noise','fontsize',14 );
            ylabel('Modulation index (%)','fontsize',14 );
            
            box on;hold on;
            if ~isempty(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]))
                scatter(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask])(sig), out.FR_ModIndex_PcS(sig),  'filled', 'MarkerFaceColor', ccol);
                scatter(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask])(~sig) , out.FR_ModIndex_PcS(~sig), 'filled', 'MarkerFaceColor', ccol/2);
                axis square;
                
                %xf = [min(out.(['SNR_' cfg.condition(c).saccadeTask])), max(out.(['SNR_' cfg.condition(c).saccadeTask]))];
                [p,S] = polyfit(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask])(~isnan(out.FR_ModIndex_PcS)),out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_PcS)),1); %
                [y_fit,delta] = polyval(p,Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask])(~isnan(out.FR_ModIndex_PcS)),S);
                [coef, pval]  = corr(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_ModIndex_PcS, 'rows','complete') ;
                plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask])(~isnan(out.FR_ModIndex_PcS)), y_fit,'LineWidth', 2, 'Color', ccol);
                title({[T ' ' L], ['coef, p ', num2str([roundto(coef,2), roundto(pval,4)])]})
                %legend({'Significant Units', 'Non-Significant Units', 'Linear Fit (Overall)'}, 'Location', 'Best')
            end
            
            subplot(2,4,2);
            box on;hold on;
            scatter(out.FR(sig)  , out.FR_ModIndex_SubtrSDP(sig),  'filled', 'MarkerFaceColor', ccol)
            scatter(out.FR(~sig) , out.FR_ModIndex_SubtrSDP(~sig), 'filled', 'MarkerFaceColor', ccol/2)
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('mean FR','fontsize',14 );
            axis square;
            %legend({'Significant Units', 'Non-Significant Units'}, 'Location', 'Best')
            
            subplot(2,4,3);
            box on;hold on;
            scatter(out.sig_n_bins(sig)  , out.FR_ModIndex_SubtrSDP(sig), 'filled', 'MarkerFaceColor', ccol)
            scatter(out.sig_n_bins(~sig) , out.FR_ModIndex_SubtrSDP(~sig), 'filled', 'MarkerFaceColor', ccol/2)
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('sig. Nr. bins','fontsize',14 );
            axis square;
            %legend({'Significant Units', 'Non-Significant Units'}, 'Location', 'Best')
            
            subplot(2,4,5:6);
            box on;hold on;
            UnitSig_Rest = out.FR_ModIndex_SubtrSDP >30 & sig;
            UnitNotSign_Rest = out.FR_ModIndex_PcS > 30 & ~sig;
            text(-400,-15, Out.(T).unit_ID(UnitSig_Rest),'Color', ccol);
            if sum(UnitSig_Rest)
                line(PSTH_bins, out.SDsubstractedSDP_normalized(:,UnitSig_Rest), 'color', ccol, 'LineWidth', 1);
            end
            text(300,20, Out.(T).unit_ID(UnitNotSign_Rest),'Color','k');
            xlabel('Time from R-peak, ms')
            ylabel('% signal change')
            title('Units with Modulation strength > 30%');
            
            subplot(2,4,7:8);
            box on; hold on;
            text(300,20, Out.(T).unit_ID(UnitNotSign_Rest),'Color','k');
            if sum(UnitNotSign_Rest)
                line(PSTH_bins, out.SDsubstractedSDP_normalized(:,UnitNotSign_Rest), 'color', ccol,'LineWidth', 4);
            end
            xlabel('Time from R-peak, ms')
            ylabel('% signal change')
            title('Non-Sign units with Modulation strength > 30%');
            
            save_figure_as(['Check_ModulationIndex_',T,'_',(L) ],output_folder,savePlot)
        end
    end
    
    if Graph_SelectionCriterion
        %%  normalizing the Firing rate calcuations
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            figure('Name',sprintf(T),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
            for c=1:N_conditions
                L=cfg.condition(c).name;
                out = [Out.(T).(L)];
                
                sig =  ~isnan(out.sig_FR_diff);
                
                subplot(2,2,1:2);
                hold on
                box on
                % bar plot how many are significant & positive and negative FR?
                SDmean_SEM = nanstd(out.SDsubstractedSDP(:,sig),[],2)/ sqrt(sum(sig)) ;
                shadedErrorBar(PSTH_bins, nanmean(out.SDsubstractedSDP(:,sig),2), SDmean_SEM, lineProps{c}, 1);
                ylabel('Signal Change, spikes/s')
                title(['Pop:  (all significant) Cal Per Unit:SD-SDP ' T ' units'],'interpreter','none');
                legend(ConNames)
                
                subplot(2,2,3:4);
                hold on
                box on
                SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:,sig),[],2)/ sqrt(sum(sig)) ;
                shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:,sig),2), SDmean_SEM, lineProps{c},1);
                title(['Population:  (all significant)' (T) ' units'],'interpreter','none');
                ylabel('% signal change','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );
                title(['Pop:  (all significant) Cal Per Unit:(SD-SDP)/SDP ' T ' units'],'interpreter','none');
                legend(ConNames)
            end
            
            save_figure_as(['Normalization_FR_ECG_triggered_spike_' ,T],output_folder,savePlot)
        end
        
        %% Plot the dataset - check up
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            for c=1:N_conditions
                L=cfg.condition(c).name;
                ccol=cfg.condition(c).color;
                out = [Out.(T).(L)];
                
                figure('Name',sprintf([T '  ' L]),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
                rows_plot = 3;
                colums_plot = 3;
                
                for I = 1: 5
                    % How many different groups
                    if I == 1
                        idx1 = ([out.sig_n_bins] > 0 );
                        idx2 = ([out.FR_ModIndex_SubtrSDP] > 5 );
                        if sum(idx2) ~= sum(idx1&idx2)
                            disp('incorrect Selection') %% ???
                        end
                        subplot(rows_plot,colums_plot,4);
                        title('FR_ModIndex_SubtrSDP > 5 ','interpreter','none');
                        
                    elseif I == 2
                        idx1 = ([out.sig_n_bins] > 0 );
                        idx2 = ([out.FR_ModIndex_SubtrSDP] < 5 );
                        subplot(rows_plot,colums_plot,5);
                        title('FR_ModIndex_SubtrSDP < 5 ','interpreter','none');
                        
                    elseif I == 3
                        idx1 = ([out.sig_n_bins] > 0 );
                        idx2 = ([out.sig_n_bins] <= 4 );
                        subplot(rows_plot,colums_plot,7);
                        title('Nr. of bin: 1 - 4 bins for sign. Intervals','interpreter','none');
                    elseif I == 4
                        idx1 = ([out.sig_n_bins] >= 5);
                        idx2 = ([out.sig_n_bins] < 7 );
                        subplot(rows_plot,colums_plot,8);
                        title('Nr. of bin: 5 - 7 bins for sign. Intervals','interpreter','none');
                        
                    elseif I == 5
                        idx1 = ([out.sig_n_bins] > 6 );
                        idx2 = ([out.sig_n_bins] > 6 );
                        subplot(rows_plot,colums_plot,9);
                        title('Nr. of bin: above 6 bins for sign. Intervals','interpreter','none');
                        
                    end
                    Y1 = out.SD(:, idx1 & idx2)  ;
                    Y2 = out.SDP(:, idx1 & idx2)  ;
                    
                    if ~isempty(Y1)
                        cmap=jet(size(Y1,1));
                        colormap(cmap)
                        line(PSTH_bins, Y1 - Y2, 'linewidth', 3)
                        set(gca,'ylim',[-9, 9]);
                        ylabel('Firing rate (%)','fontsize',12 );
                        xlabel('Time relative to ECG peak (ms)','fontsize',12);
                    end
                end
                
                idx_ex = ([out.sig_n_bins] <= 4 );
                
                subplot(rows_plot,colums_plot,1);
                scatter(out.sig_n_bins , Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]), 'filled', 'MarkerFaceColor', ccol); hold on;
                scatter(out.sig_n_bins(idx_ex) , Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask])(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
                ylabel('Signal-to-Noise Ratio','fontsize',14 );
                xlabel('Nr. of bins of sig. Interval','fontsize',14 );
                axis square;
                
                subplot(rows_plot,colums_plot,2);
                scatter(out.sig_n_bins , out.FR_ModIndex_SubtrSDP, 'filled', 'MarkerFaceColor', ccol); hold on;
                scatter(out.sig_n_bins(idx_ex) , out.FR_ModIndex_SubtrSDP(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
                ylabel('FR Modulation (spike/s)','fontsize',14 );
                xlabel('Nr. of bins of sig. Interval','fontsize',14 );
                axis square;
                
                subplot(rows_plot,colums_plot,3);
                scatter(out.sig_time , out.FR_ModIndex_SubtrSDP, 'filled', 'MarkerFaceColor', ccol); hold on;
                scatter(out.sig_time(idx_ex) , out.FR_ModIndex_SubtrSDP(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
                ylabel('FR Modulation (spike/s)','fontsize',10);
                xlabel('TimePoint of sig. highest diff in FR','fontsize',10 );
                axis square;
                
                save_figure_as(['ModulationInRelationToBins_ECG_triggered_spike_' (L) , '_' ,(T)],output_folder,savePlot)
            end
        end
    end
    
    %% Plot the averages
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        figure('Name',sprintf(T),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
        
        %Bar graph how many
        subplot(2,4,3);
        toplot=arrayfun(@(x) Out.(T).(x.name).Pc_SignFR, cfg.condition, 'uniformoutput', false);
        barpairs =  vertcat(toplot{:});
        b = bar(barpairs,'stacked', 'Facecolor','flat' );
        b(3).FaceColor = [1 1 1];
        legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')
        set(gca,'XTickLabel',ConNames,'fontsize',10);
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            
            out = [Out.(T).(L)];
            % from all not NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(end,:));
            Idx_Units_NaN    = sum(isnan(out.SDsubstractedSDP(end,:)));
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            dec = (out.sig_sign == -1);
            inc = (out.sig_sign == 1);
            
            % bar plot how many are significant & positive and negative FR?
            subplot(2,4,1:2);
            hold on;
            if sum(sig)
                line(PSTH_bins, out.SDsubstractedSDP_normalized(:,sig), 'color', ccol, 'LineWidth', 1);
            end
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:, sig), [], 2)/ sqrt(sum(sig)) ;
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:, sig), 2) ,SDmean_SEM,lineProps{c},1);
            ylabel('normalized Firing rate (%)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            box on
            if c == 1
                ttl_sub12 = {['Population (all significant): ' (T) ' units']};
                ttl_sub12 = [ttl_sub12 ['\color[rgb]{' num2str(ccol) '}' L ', units = ' ,num2str(sum(sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN))]];
            elseif c == N_conditions
                ttl_sub12 = [ttl_sub12 ['\color[rgb]{' num2str(ccol) '}' L ', units = ' ,num2str(sum(sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN))]];
                title(ttl_sub12)
            end
            
            % display only significant units showing a increase in FR
            subplot(2,4,5); %
            hold on;
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:, inc & sig),0,2)/ sqrt(sum(inc & sig)) ;
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:, inc & sig),2), SDmean_SEM ,lineProps{c},1);
            set(gca,'ylim',[-10, 10]);
            ylabel('normalized Firing rate (%)','fontsize',14);
            xlabel('Time relative to R-peak (ms)','fontsize',14);
            box on
            if c == 1
                ttl_sub5 = {'units showing a sig. INCREASE in FR'};
                ttl_sub5 = [ttl_sub5 ['\color[rgb]{' num2str(ccol) '}' L ': units = ' num2str(sum(inc & sig))]];
            elseif c == N_conditions
                ttl_sub5 = [ttl_sub5 ['\color[rgb]{' num2str(ccol) '}' L ': units = ' num2str(sum(inc & sig))]];
                title(ttl_sub5)
            end
            
            % display only significant units showing a decrease in FR
            subplot(2,4,6);
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:, dec & sig),0,2)/ sqrt(sum(dec & sig)) ;
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:, dec & sig),2), SDmean_SEM, lineProps{c},1);
            set(gca,'ylim',[-10, 10]);
            ylabel('normalized Firing rate (spike/s)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            box on
            if c == 1
                ttl_sub6 = {'units showing a sig. DECREASE in FR','interpreter'};
                ttl_sub6 = [ttl_sub6 ['\color[rgb]{' num2str(ccol) '}' L ': units = ' num2str(sum(dec & sig))]];
            elseif c == N_conditions
                ttl_sub6 = [ttl_sub6 ['\color[rgb]{' num2str(ccol) '}' L ': units = ' num2str(sum(dec & sig))]];
                title(ttl_sub6)
            end
            
            
            %% MODULATION INDEX
            subplot(2,4,4);
            hold on;
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('Signal-to-Noise','fontsize',14 );
            axis square;
            %         SNR_fieldname = (['SNR_' cfg.condition(c).saccadeTask]);
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_ModIndex_PcS,sig,0,ccol,'raw') %% or FR_ModIndex_SubtrSDP ??
            
            subplot(2,4,7);
            hold on;
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('Signal-to-Noise','fontsize',14 );
            axis square;
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_ModIndex_SubtrSDP,sig,0,ccol,'-Surrogate')
            
            subplot(2,4,8);
            hold on
            ylabel('Modulation index(%)','fontsize',14 );
            xlabel('mean firing rate','fontsize',14);
            axis square;
            make_correlation_plot(out.FR, out.FR_ModIndex_SubtrSDP,sig,1,ccol,'-Surrogate, sig correlated')
            
        end
        save_figure_as(['Average_ECG_triggered_spike_' (T)],output_folder,savePlot)
    end
    
    %% MODULATION STRENGTH - for the percentage signal change
    figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            %% Modulation Index
            subplot(2,N_Areas,a);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('SNR','fontsize',14 );
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_ModIndex_PcS,sig,0,ccol,[T 'all units'])
            
            subplot(2,N_Areas,a +N_Areas);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('SNR','fontsize',14 );
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_ModIndex_PcS,sig,1,ccol,[T ' signficant only'])
        end
    end
    save_figure_as('ModulationIndex_Npsc',output_folder,savePlot)
    
    %% MODULATION STRENGTH - for SDP subtraction
    figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(2,N_Areas,a);
            ylabel('Modulation index (spike/s)','fontsize',14 );
            xlabel('Signal-to-Noise ratio (mV)','fontsize',14 );
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_ModIndex_SubtrSDP,sig,0,ccol,[T 'all units'])
            
            subplot(2,N_Areas,a +N_Areas);
            ylabel('Modulation index (spike/s)','fontsize',14 );
            xlabel('Signal-to-Noise ratio (mV)','fontsize',14 );
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_ModIndex_SubtrSDP,sig,1,ccol,[T ' signficicant only'])
        end
    end
    save_figure_as('ModulationIndex_NSubtr',output_folder,savePlot)
    
    %% MODULATION STRENGTH - compare both Normalizations ... SDP subtraction
    figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(2,N_Areas,a);
            ylabel('Modulation index (spike/s)','fontsize',14 );
            xlabel('Modulation index(%pSc)','fontsize',14 );
            make_correlation_plot(out.FR_ModIndex_PcS, out.FR_ModIndex_SubtrSDP,sig,0,ccol,[T 'all units']) %% was FR_ModIndex_PcS...
            
            subplot(2,N_Areas,a +N_Areas);
            ylabel('Modulation index (spike/s)','fontsize',14 );
            xlabel('Modulation index(%pSc)','fontsize',14 );
            make_correlation_plot(out.FR_ModIndex_PcS, out.FR_ModIndex_SubtrSDP,sig,1,ccol,[T ' signficicant only'])
        end
    end
    save_figure_as('ModulationIndex_Cmp_Npsc_NSubtr',output_folder,savePlot)
    
    %% MODULATION STRENGTH - FR as explanation - compare both Normalizations ... SDP subtraction
    figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(N_conditions,N_Areas,a + (c-1)*N_Areas);
            hold on;
            %% Modulation Index
            % all not significant & NaN units
            if any(sig)
                % scatter(out.FR_ModIndex_PcS(sig) , out.FR_ModIndex_SubtrSDP(sig), 'filled', 'MarkerFaceColor',Color)
                scatter(out.FR_ModIndex_PcS(sig) , out.FR_ModIndex_SubtrSDP(sig),30, out.FR_perECGTriggeredAverage(sig)/max(out.FR_perECGTriggeredAverage(sig)) , 'filled')
                %scatter(out.FR_ModIndex_PcS(~sig) , out.FR_ModIndex_SubtrSDP(~sig),30, out.FR(~sig)/max(out.FR) , 'filled')
                
                colorbar();
                colormap jet
                ylabel('Modulation index (spike/s)','fontsize',14 );
                xlabel('Modulation index(%pSc)','fontsize',14 );
                axis square; box on;
                title([L ' ' T]);
                if strcmp(T, 'VPL_R') || strcmp(T, 'dPul_R') || strcmp(T, 'MD')
                    set(gca,'xlim',[0, 60]);
                end
                set(gca,'ylim',[0, 15]);
            end
        end
    end
    save_figure_as('ModulationIndex_Compare_NpSc_Subtr_WithFR_ONLYSignificantUnits',output_folder,savePlot)
    
    %% MODULATION STRENGTH - FR as explanation - compare both Normalizations ... SDP subtraction
    figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            %% Modulation Index
            subplot(N_conditions,N_Areas,a + (c-1)*N_Areas);
            hold on;
            if any(sig)
                % all not significant & NaN units
                % scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
                scatter(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask])(sig),out.FR_ModIndex_PcS(sig) , 30, out.FR_perECGTriggeredAverage(sig)/max(out.FR_perECGTriggeredAverage(sig)) , 'filled')
                %scatter(out.FR_ModIndex_PcS(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig),30, out.FR(~idx_sig)/max(out.FR) , 'filled')
                
                colorbar();
                colormap jet
                ylabel('Modulation index (%sc)','fontsize',14 );
                xlabel('signal-to-noise ratio','fontsize',14 );
                axis square; box on;
                title([L ' ' T]);
                set(gca,'ylim',[0, 80]);
                set(gca,'xlim',[0, 25]);
            end
        end
    end
    save_figure_as('ModulationIndex_NpSc_SNR_FR_ONLYSignificantUnits',output_folder,savePlot)
    
    %% MODULATION STRENGTH pSC- FR
    figure('Name',sprintf('ModulationIndex_FiringRate'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(2,N_Areas,a);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );
            make_correlation_plot(out.FR,out.FR_ModIndex_PcS,sig,0,ccol,[T ' all units'])
            
            subplot(2,N_Areas,a +N_Areas);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );
            make_correlation_plot(out.FR,out.FR_ModIndex_PcS,sig,1,ccol,[T ' only significant'])
        end
    end
    save_figure_as('ModulationIndex_Npsc_FiringRate',output_folder,savePlot)
    
    %% Correlation between SNR - FR
    figure('Name',sprintf('Correlation_SNR_FiringRate'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            ccol=cfg.condition(c).color;
            
            out = [Out.(T).(L)];
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(2,N_Areas,a);
            xlabel('Signal-to-Noise ratio','fontsize',14 );
            ylabel('average Firing rate','fontsize',14 );
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_perECGTriggeredAverage,sig,0,ccol,[T ' all units'])
            
            subplot(2,N_Areas,a +N_Areas);
            xlabel('Signal-to-Noise ratio','fontsize',14 );
            ylabel('average Firing rate','fontsize',14 );
            make_correlation_plot(Out.(T).criteria.(['SNR_' cfg.condition(c).saccadeTask]),out.FR_perECGTriggeredAverage,sig,1,ccol,[T ' only significant'])
        end
    end
    save_figure_as('SNR_FiringRate',output_folder,savePlot)
    
    %% MODULATION STRENGTH - substractiveNormalization - FR
    % Shows that higher average firing rate is related to a higher modulation index
    figure('Name',sprintf('ModulationIndex_FiringRate'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            out = [Out.(T).(L)];
            ccol=cfg.condition(c).color;
            
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(2,N_Areas,a);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );
            make_correlation_plot(out.FR,out.FR_ModIndex_SubtrSDP,sig,0,ccol,[T ' all units'])
            
            subplot(2,N_Areas,a +N_Areas);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );
            make_correlation_plot(out.FR,out.FR_ModIndex_SubtrSDP,sig,1,ccol,[T ' only significant'])
        end
    end
    save_figure_as('ModulationIndex_NSubt_FR',output_folder,savePlot)
    
    %% MODULATION STRENGTH -  Bin size
    figure('Name',sprintf('ModulationIndex_BinSize'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        for c=1:N_conditions
            L=cfg.condition(c).name;
            out = [Out.(T).(L)];
            ccol=cfg.condition(c).color;
            
            sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
            
            subplot(2,N_Areas,a);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('Bin size of sig. Interval','fontsize',14 );
            make_correlation_plot(out.sig_n_bins,out.FR_ModIndex_PcS,sig,0,ccol,[T ' all units'])
            
            subplot(2,N_Areas,a +N_Areas);
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('Bin size of sig. Interval','fontsize',14 );
            make_correlation_plot(out.sig_n_bins,out.FR_ModIndex_PcS,sig,1,ccol,[T ' only significant'])
        end
    end
    save_figure_as('ModulationIndex_Nspc_Binsize',output_folder,savePlot)
    
    %% Decrease and Increase grouped for brain region
    figure('Name',sprintf('CardiacRelated_Change_FR'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for c=1:N_conditions
        L=cfg.condition(c).name;
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            acol=Color_BrainArea{a};
            out = [Out.(T).(L)];
            
            sig      =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins);
            dec      = (out.sig_sign == -1) & sig;
            inc      = (out.sig_sign == 1) & sig;
            %         idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
            %         idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
            %         idx_SigTime_After50       = (out.sig_time > 50);
            
            subplot(N_conditions,4,1+(c-1)*4);
            if a == 1
                title({[L ' Increase']})
            end
            ttl_set1 = get(gca, 'Title');
            ttl_str1 = ttl_set1.String;
            ttl_str1 = [ttl_str1; {['\color[rgb]{' num2str(acol) '}' T ' ' L ': units = ' num2str(sum(inc))]}];
            title(ttl_str1)
            hold on; axis square; box on;
            lineProps={'color',acol,'linewidth',3};
            %         text(-400,-1* a, [(T), ' ' L ': units = ' ,num2str(sum(inc))],'Color',acol)
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:, inc),0,2)/ sqrt(sum(inc));
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:, inc),2), SDmean_SEM ,lineProps,1);
            %set(gca,'ylim',[-10, 10]);
            ylabel('normalized Firing rate (%)','fontsize',12 );
            xlabel('Time relative to R-peak (ms)','fontsize',12 );
            % set(gca, 'XTick', PSTH_bins)
            
            % display only significant units showing a decrease in FR
            subplot(N_conditions,4,3+(c-1)*4);
            if a == 1
                title({[L ' Decrease' ]})
            end
            ttl_set2 = get(gca, 'Title');
            ttl_str2 = ttl_set2.String;
            ttl_str2 = [ttl_str2; {['\color[rgb]{' num2str(acol) '}' T ' ' L ': units = ' num2str(sum(dec))]}];
            title(ttl_str2)
            hold on; axis square; box on;
            lineProps={'color',acol,'linewidth',3};
            text(-400,1* a, [T ' ' L ': units = ' ,num2str(sum(dec)) ],'Color',acol)
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:,dec),0,2)/ sqrt(sum(dec)) ;
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:,dec),2), SDmean_SEM ,lineProps,1);
            %set(gca,'ylim',[-10, 10]);
            %         title(,'interpreter','none')
            ylabel('normalized Firing rate (spike/s)','fontsize',12 );
            xlabel('Time relative to R-peak (ms)','fontsize',12 );
            % set(gca, 'XTick', PSTH_bins)
            
            subplot(N_conditions,4,2+(c-1)*4);
            hold on; axis square; box on;
            scatter(out.sig_time(inc) , out.FR_ModIndex_SubtrSDP(inc), 'filled', 'MarkerFaceColor',acol);
            ylabel('Modulation strength (% pSc)','fontsize',12);
            xlabel('TimePoint of sig. highest diff in FR','fontsize',12 );
            title([L ' Increase' ],'interpreter','none')
            
            subplot(N_conditions,4,4+(c-1)*4);
            hold on; axis square; box on;
            scatter(out.sig_time(dec) , out.FR_ModIndex_SubtrSDP(dec), 'filled', 'MarkerFaceColor',acol);
            ylabel('Modulation strength (% pSc)','fontsize',12);
            xlabel('TimePoint of sig. highest diff in FR','fontsize',12 );
            title([L ' Decrease' ],'interpreter','none')
        end
    end
    save_figure_as('Suppression_Enhancement',output_folder,savePlot)
    
    %% Decrease and Increase grouped for brain region
    for groupNum = 1:length(cfg.spk.compare_conditions)
        
        cond1_num = cfg.spk.compare_conditions{groupNum}(1);
        cond2_num = cfg.spk.compare_conditions{groupNum}(2);
        
        figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
        ThreeTiming = {'T<0',  'T>0'} ; %{'T<-50', '-50>T<50', 'T>50'} ;
        for c=1:length(cfg.spk.compare_conditions{groupNum})
            L=cfg.condition(c).name;
            for i_Time = 1: length(ThreeTiming)
                for a = 1: N_Areas
                    T=Ana_TargetBrainArea{a};
                    acol=Color_BrainArea{a};
                    lineProps={'color',acol,'linewidth',3};
                    out = [Out.(T).(L)];
                    
                    sig      =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
                    dec      = (out.sig_sign == -1) & sig;
                    inc      = (out.sig_sign == 1)  & sig;
                    
                    %                 idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
                    %                 idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
                    %                 idx_SigTime_After50       = (out.sig_time > 50);
                    
                    switch i_Time
                        case 1 %{ 'BeforeMinus50'}
                            tim = out.sig_time < 0 ;
                            WindowIdx = find(PSTH_bins < 0);
                        case 2 %'After50'
                            tim = out.sig_time > 0;
                            WindowIdx = find(PSTH_bins > 0);
                    end
                    
                    %% INCREASE
                    subplot(length(cfg.spk.compare_conditions{groupNum}),6,2*(i_Time-1)+1 +(c-1)*6);
                    hold on; axis square; box on;
                    
                    if i_Time == 1 && c == 1
                        xlabel('Time from the R-peak, ms','fontsize',11 );
                        ylabel('% FR change','fontsize',11 );
                    end
                    
                    if a == 1
                        if i_Time == 1 % diastole
                            fill([0 -225 -225 0 0], [-30 -30 30 30 -30], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2)
                        elseif i_Time == 2 % systole
                            fill([0 225 225 0 0], [-30 -30 30 30 -30], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2)
                        end
                    end
                    
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:, inc & tim,:),0,2)/ sqrt(sum(inc & tim)) ;
                    shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:,inc & tim),2), SDmean_SEM, lineProps,1);
                    MI_groups =  max(nanmean(out.SDsubstractedSDP_normalized(WindowIdx,inc & tim),2)) - min(nanmean(out.SDsubstractedSDP_normalized(:,inc & tim),2)) ;
                    
                    if a == 1
                        title({[L ' INCREASE ' ThreeTiming{i_Time}]})
                    end
                    ttl_set1 = get(gca, 'Title');
                    ttl_str1 = ttl_set1.String;
                    ttl_str1 = [ttl_str1; {['\color[rgb]{' num2str(acol) '}' T ' : units = ' num2str(sum(inc & tim))]}];
                    title(ttl_str1)
                    
                    ylim([-30 30])
                    vline(0, 'k')
                    
                    %% Table stuff
                    %             Table_Units = table(T,L,{'sigINCREASE'}, ThreeTiming(i_Time),  sum(inc & sig & tim), roundto(MI_groups,2) );
                    %             Table_NrUnits = [Table_NrUnits; Table_Units];
                    
                    %% DECREASE
                    subplot(length(cfg.spk.compare_conditions{groupNum}),6,2*(i_Time-1)+2 +(c-1)*6);
                    hold on; axis square; box on;
                    
                    if a == 1
                        if i_Time == 1 % diastole
                            fill([0 -225 -225 0 0], [-30 -30 30 30 -30], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2)
                        elseif i_Time == 2 % systole
                            fill([0 225 225 0 0], [-30 -30 30 30 -30], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2)
                        end
                    end
                    
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:,dec & tim),0,2)/ sqrt(sum(dec & tim)) ;
                    shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:,dec & tim),2), SDmean_SEM ,lineProps,1);
                    MI_groups =   max(nanmean(out.SDsubstractedSDP_normalized(:,dec & tim),2)) - min(nanmean(out.SDsubstractedSDP_normalized(WindowIdx,dec & tim),2));
                    
                    if a == 1
                        title({[L ' DECREASE ' ThreeTiming{i_Time}]})
                    end
                    ttl_set1 = get(gca, 'Title');
                    ttl_str1 = ttl_set1.String;
                    ttl_str1 = [ttl_str1; {['\color[rgb]{' num2str(acol) '}' T ' : units = ' num2str(sum(dec & tim))]}];
                    title(ttl_str1)
                    
                    ylim([-30 30])
                    vline(0, 'k')
                    
                    %% Table stuff
                    %             Table_Units = table(T,L,{'sigDECREASE'}, ThreeTiming(i_Time),  sum(idx_SigDec & idx_sig & idx_Time), roundto(MI_groups,2)  );
                    %             Table_NrUnits = [Table_NrUnits; Table_Units];
                    
                    subplot(length(cfg.spk.compare_conditions{groupNum}),6,5 +(c-1)*6);
                    hold on; axis square; box on;
                    scatter(out.sig_time(dec) , out.FR_ModIndex_SubtrSDP(dec), 'filled', 'MarkerFaceColor',acol);
                    title([L ' Decrease' ])
                    if i_Time == 1 && c == 1
                        xlabel('Response Latency, ms','fontsize',11 );
                        ylabel('% FR Change','fontsize',11);
                    end
                    
                    subplot(length(cfg.spk.compare_conditions{groupNum}),6,6 +(c-1)*6);
                    hold on; axis square; box on;
                    scatter(out.sig_time(inc) , out.FR_ModIndex_SubtrSDP(inc), 'filled', 'MarkerFaceColor',acol);
                    title([L ' Increase' ])
                    
                end
            end
        end
        save_figure_as('Suppression_Enhancement_SeparatedForTime_NoLINES',output_folder,savePlot)
        
        if saveTable;
            filename = [output_folder , 'Table_NbUnits_IncreaseDecreaseAtSpecificTime.xlsx' ] ;
            writetable(Table_NrUnits,filename,'Sheet',1,  'Range' ,'A1' )
        end
    end
    
    %% check significance of heart response times
    pp = nan(12,1);
    hh = nan(12,1);
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        sig_targ1  = ~isnan(Out.(unqTargets{1}).(L).sig_FR_diff) & (Out.(unqTargets{1}).(L).sig_n_bins > n_sig_bins);
        sig_targ2 = ~isnan(Out.(unqTargets{2}).(L).sig_FR_diff) & (Out.(unqTargets{2}).(L).sig_n_bins > n_sig_bins);
        sig_targ3   = ~isnan(Out.(unqTargets{3}).(L).sig_FR_diff) & (Out.(unqTargets{3}).(L).sig_n_bins > n_sig_bins);
        
        inc_targ1  = Out.(unqTargets{1}).(L).sig_time((Out.(unqTargets{1}).(L).sig_sign == 1) & sig_targ1);
        mm(1+6*(c-1)) = median(inc_targ1);
        ci(1+6*(c-1)) = iqr(inc_targ1);
        inc_targ2 = Out.(unqTargets{2}).(L).sig_time((Out.(unqTargets{2}).(L).sig_sign == 1) & sig_targ2);
        mm(2+6*(c-1)) = median(inc_targ2);
        ci(2+6*(c-1)) = iqr(inc_targ2);
        inc_targ3   = Out.(unqTargets{3}).(L).sig_time((Out.(unqTargets{3}).(L).sig_sign == 1) & sig_targ3);
        mm(3+6*(c-1)) = median(sig_targ3);
        ci(3+6*(c-1)) = iqr(inc_targ3);
        
        dec_targ1  = Out.(unqTargets{1}).(L).sig_time((Out.(unqTargets{1}).(L).sig_sign == -1) & sig_targ1);
        mm(4+6*(c-1)) = median(dec_targ1);
        ci(4+6*(c-1)) = iqr(dec_targ1);
        dec_targ2 = Out.(unqTargets{2}).(L).sig_time((Out.(unqTargets{2}).(L).sig_sign == 1) & sig_targ2);
        mm(5+6*(c-1)) = median(dec_targ2);
        ci(5+6*(c-1)) = iqr(dec_targ2);
        dec_targ3   = Out.(unqTargets{3}).(L).sig_time((Out.(unqTargets{3}).(L).sig_sign == 1) & sig_targ3);
        mm(6+6*(c-1)) = median(dec_targ3);
        ci(6+6*(c-1)) = iqr(dec_targ3);
        
        if ~isempty(inc_targ1) && ~isempty(inc_targ2)
            disp(['Increase: ' unqTargets{1} ' vs. ' unqTargets{2}])
            [pp(1+6*(c-1)), hh(1+6*(c-1))] = ranksum(inc_targ1, inc_targ2)
        end
        if ~isempty(inc_targ1) && ~isempty(inc_targ3)
                disp(['Increase: ' unqTargets{1} ' vs. ' unqTargets{3}])
                [pp(2+6*(c-1)), hh(2+6*(c-1))] = ranksum(inc_targ1, inc_targ3)
            end
        if ~isempty(inc_targ2) && ~isempty(inc_targ3)
            disp(['Increase: ' unqTargets{2} ' vs. ' unqTargets{3}])
            [pp(3+6*(c-1)), hh(3+6*(c-1))] = ranksum(inc_targ2, inc_targ3)
        end
        
        if ~isempty(dec_targ1) && ~isempty(dec_targ2)
            disp('Decrease: unqTargets{1} vs. unqTargets{2}')
            [pp(4+6*(c-1)), hh(4+6*(c-1))] = ranksum(dec_targ1, dec_targ2)
        end
        if ~isempty(dec_targ1) & ~isempty(dec_targ3)
            disp(['Decrease: ' unqTargets{1} ' vs. ' unqTargets{3}])
            [pp(5+6*(c-1)), hh(5+6*(c-1))] = ranksum(dec_targ1, dec_targ3)
        end
        if ~isempty(dec_targ2) & ~isempty(dec_targ3)
            disp(['Decrease: ' unqTargets{2} ' vs. ' unqTargets{3}])
            [pp(6+6*(c-1)), hh(6+6*(c-1))] = ranksum(dec_targ1, dec_targ3)
        end
    end
    
    %% GENERATE A TABLE FOR ALL THE NUMBER OF UNITS PER CATEGORY
    figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    ThreeTiming = {'inc T<-50 or dec T>50', '-50>T<50', 'inc T>50 or dec T<-50'} ;
    
    for c=1:N_conditions
        L=cfg.condition(c).name;
        for i_Time = 1: length(ThreeTiming)
            for a = 1: N_Areas
                T=Ana_TargetBrainArea{a};
                acol=Color_BrainArea{a};
                out = [Out.(T).(L)];
                lineProps={'color',acol,'linewidth',3};
                
                sig      =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
                dec      = (out.sig_sign == -1) & sig;
                inc      = (out.sig_sign == 1) & sig;
                idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
                idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
                idx_SigTime_After50       = (out.sig_time > 50);
                
                switch i_Time
                    case 1 %{ 'BeforeMinus50'}
                        tim = (idx_SigTime_BeforeMinus50 & inc) |  (idx_SigTime_After50 & dec) ;
                    case 2 %'Around0'
                        tim = idx_SigTime_Around0 ;
                    case 3 %'After50'
                        tim = (idx_SigTime_BeforeMinus50 & dec) | (idx_SigTime_After50 & inc) ;
                end
                
                subplot(N_conditions,3,i_Time +(c-1)*3);
                hold on; axis square; box on;
                SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:,tim),0,2)/ sqrt(sum(tim)) ;
                shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:, tim),2), SDmean_SEM ,lineProps,1);
                
                if a == 1
                    title({[L ' ' ThreeTiming{i_Time}]})
                end
                ttl_set1 = get(gca, 'Title');
                ttl_str1 = ttl_set1.String;
                ttl_str1 = [ttl_str1; {['\color[rgb]{' num2str(acol) '}' T ' : units = ' num2str(sum(tim))]}];
                title(ttl_str1)
                
                if a == N_Areas
                    ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                    xlabel('Time relative to R-peak (ms)','fontsize',14 );
                    xlim(1000 * [curr_analyse_states{3} curr_analyse_states{4}])
                    ylim([-20 20])
                    vline(0, 'k')
                end
                ax = gca;
                ax.FontSize = 14;
            end
        end
    end
    save_figure_as('Suppression_Enhancement_SeparatedForTime_GroupedAccordingTo',output_folder,savePlot)
    
    %% split units into 2 groups: 1 - inc before the R-peak or dec after R-peak; 2 - inc or dec after the R-peak
    figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 196 824 800],'PaperPositionMode', 'auto');
    ThreeTiming = {'inc T<0 or dec T>0', 'inc T>0 or dec T<0'};
    ThreeTimingTitles = {'Inc. before or dec. after', 'Dec. before or inc. after'};
    N_Timing = length(ThreeTiming);
    
    for c=1:N_conditions
        L=cfg.condition(c).name;
        for i_Time = 1:N_Timing
            for a = 1: N_Areas
                T=Ana_TargetBrainArea{a};
                acol=Color_BrainArea{a};
                out = [Out.(T).(L)];
                lineProps={'color',acol,'linewidth',3};
                
                sig      =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
                dec      = (out.sig_sign == -1) & sig;
                inc      = (out.sig_sign == 1) & sig;
                idx_SigTime_BeforeMinus50 = (out.sig_time < 0);
                idx_SigTime_After50       = (out.sig_time > 0);
                
                switch i_Time
                    case 1 %{ 'Before R'}
                        tim = (idx_SigTime_BeforeMinus50 & inc) |  (idx_SigTime_After50 & dec) ;
                    case 2 %'After R'
                        tim = (idx_SigTime_BeforeMinus50 & dec) | (idx_SigTime_After50 & inc) ;
                end
                
                subplot(N_conditions,N_Timing,i_Time +(c-1)*N_Timing);
                hold on; axis square; box on;
                SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:,tim),0,2)/ sqrt(sum(tim)) ;
                shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:,tim),2), SDmean_SEM ,lineProps,1);
                
                if a == 1
                    title({[L ': ' ThreeTimingTitles{i_Time}]})
                end
%                 ttl_set1 = get(gca, 'Title');
%                 ttl_str1 = ttl_set1.String;
%                 ttl_str1 = [ttl_str1; {['\color[rgb]{' num2str(acol) '}' T ' : units = ' num2str(sum(tim))]}];
%                 title(ttl_str1)
                
                sum_tim(a) = sum(tim);
                
                if a == N_Areas && i_Time == 1 && c == 1
                    xlabel('Time relative to R-peak (ms)','fontsize',14);
                    ylabel('% Signal Change','fontsize',14);
                end
                
                if a == N_Areas
                    xlim(1000 * [curr_analyse_states{3} curr_analyse_states{4}])
                    ylim(cfg.time.y_lims)
                    vline(0, 'k')
                    legend({['VPL: ' num2str(sum_tim(1))], ['dPul: ' num2str(sum_tim(2))], ['MD: ' num2str(sum_tim(3))]}, 'Location', 'southwest')
                end
                ax = gca;
                ax.FontSize = 14;
            end
        end
    end
    save_figure_as('Suppression_Enhancement_SeparatedForTime2_GroupedAccordingTo',output_folder,savePlot)
    
    %% same plot type as above for R-peak non-responsive units
    figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 196 824 800],'PaperPositionMode', 'auto');
    ThreeTiming = {'inc T<0 or dec T>0', 'inc T>0 or dec T<0'};
    ThreeTimingTitles = {'Inc. before or dec. after', 'Inc. after or dec. before'};
    N_Timing = length(ThreeTiming);
    
    for c=1:N_conditions
        L=cfg.condition(c).name;
        for i_Time = 1:N_Timing
            for a = 1: N_Areas
                T=Ana_TargetBrainArea{a};
                acol=Color_BrainArea{a};
                out = [Out.(T).(L)];
                lineProps={'color',acol,'linewidth',3};
                
                sig      =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins <= n_sig_bins) ;
                dec      = (out.sig_sign == -1) & sig;
                inc      = (out.sig_sign == 1) & sig;
                idx_SigTime_BeforeMinus50 = (out.sig_time < 0);
                idx_SigTime_After50       = (out.sig_time > 0);
                
                switch i_Time
                    case 1 %{ 'Before R'}
                        tim = (idx_SigTime_BeforeMinus50 & inc) |  (idx_SigTime_After50 & dec) ;
                    case 2 %'After R'
                        tim = (idx_SigTime_BeforeMinus50 & dec) | (idx_SigTime_After50 & inc) ;
                end
                
                subplot(N_conditions,N_Timing,i_Time +(c-1)*N_Timing);
                hold on; axis square; box on;
                SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(:,tim),0,2)/ sqrt(sum(tim)) ;
                shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(:,tim),2), SDmean_SEM ,lineProps,1);
                
                if a == 1
                    title({[L ': ' ThreeTimingTitles{i_Time}]})
                end
                
                sum_tim(a) = sum(tim);
                
                if a == N_Areas && i_Time == 1 && c == 1
                    xlabel('Time relative to R-peak (ms)','fontsize',14);
                    ylabel('% Signal Change','fontsize',14);
                end
                
                if a == N_Areas
                    xlim(1000 * [curr_analyse_states{3} curr_analyse_states{4}])
                    ylim(cfg.time.y_lims)
                    vline(0, 'k')
                    legend({['VPL: ' num2str(sum_tim(1))], ['dPul: ' num2str(sum_tim(2))], ['MD: ' num2str(sum_tim(3))]}, 'Location', 'southwest')
                end
                ax = gca;
                ax.FontSize = 14;
            end
        end
    end
    save_figure_as('Suppression_Enhancement_SeparatedForTime2_NonResponsive',output_folder,savePlot)
    
    %% below are plots for units having BOTH task and rest - they actually will work out anyway
    if OnlyUnits_withRestANDTask
        %% modulation index - percent change
        for groupNum = 1:length(cfg.spk.compare_conditions)
            
            cond1_num = cfg.spk.compare_conditions{groupNum}(1);
            cond2_num = cfg.spk.compare_conditions{groupNum}(2);
            
            cond1_cond2_color = (cfg.condition(1).color + cfg.condition(2).color) / 2;
            
            figure;
            set(gcf, 'Position', [2 381 1914 553])
            for a = 1: N_Areas
                T = Ana_TargetBrainArea{a};
                
                subplot(1,3,a)
                
                unit_count = ~isnan([Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS_signed, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS_signed]');
                unit_count = sum(any(unit_count));
                sig_cond2_cond1 = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins > n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins > n_sig_bins;
                sig_cond2      = Out.(T).(cfg.condition(cond2_num).name).sig_n_bins >= n_sig_bins & Out.(T).(cfg.condition(cond1_num).name).sig_n_bins < n_sig_bins; % in condition 2 >= sig_bins, in condition 1 < sig bins
                sig_cond1      = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins >= n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins < n_sig_bins;
                %     [cc, pp] = corrcoef(Out.(T).Rest.FR_ModIndex_PcS, Out.(T).Task.FR_ModIndex_PcS,'Rows','complete');
                hold on
                scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS_signed, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS_signed, [], [0.7 0.7 0.7], 'filled')
                scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS_signed(sig_cond2_cond1), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS_signed(sig_cond2_cond1), [], cond1_cond2_color, 'filled')
                scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS_signed(sig_cond2), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS_signed(sig_cond2), [], cfg.condition(cond1_num).color, 'filled')
                scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS_signed(sig_cond1), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS_signed(sig_cond1), [], cfg.condition(cond2_num).color, 'filled')
                %     ls_line = lsline(gca);
                %     ls_line.Color = 'b';
                box on
                grid on
                axis square
                xlim([-60 60])
                ylim([-60 60])
                %     title({[T ': N = ' num2str(unit_count)], ...
                %         ['cc = ' num2str(cc(2)) '; p = ' num2str(pp(2))]})
                title([T ': N = ' num2str(unit_count)])
                if a == 1
                    xlabel(['% signal change in ' cfg.condition(cond1_num).name])
                    ylabel(['% signal change in ' cfg.condition(cond2_num).name])
                    legend({'non-significant', ['sig. ' cfg.condition(cond2_num).name ' & sig. ' cfg.condition(cond1_num).name], ['sig. only ' cfg.condition(cond2_num).name], ['sig. only ' cfg.condition(cond1_num).name]})
                end
                ax = gca;
                ax.FontSize = 16;
            end
            save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name '_Scatter_Pc_signal_change_signed'],output_folder,savePlot)
        end
        
        %% modulation index - absolute percent change
        for groupNum = 1:length(cfg.spk.compare_conditions)
            sig = struct;
            cond1_num = cfg.spk.compare_conditions{groupNum}(1);
            cond2_num = cfg.spk.compare_conditions{groupNum}(2);
            cond1_cond2_color = (cfg.condition(1).color + cfg.condition(2).color) / 2;
            
            figure
            set(gcf, 'Position', [2 332 1914 602])
            for a = 1: N_Areas
                T = Ana_TargetBrainArea{a};
                
                subplot(1,3,a)
                unit_count = ~isnan([Out.(T).((cfg.condition(cond1_num).name)).FR_ModIndex_PcS, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS]');
                unit_count = sum(any(unit_count));
                sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins > n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins > n_sig_bins; % both conditions
                %         sig_task_rest = Out.(T).(c_names{1}).sig_n_bins > n_sig_bins & Out.(T).(c_names{2}).sig_n_bins > n_sig_bins;
                sig.(cfg.condition(cond1_num).name)                  = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins <= n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins > n_sig_bins;
                %         sig_task      = Out.(T).(c_names{1}).sig_n_bins <= n_sig_bins & Out.(T).(c_names{2}).sig_n_bins > n_sig_bins;
                sig.(cfg.condition(cond2_num).name)                  = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins > n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins <= n_sig_bins;
                %         sig_rest      = Out.(T).(c_names{1}).sig_n_bins > n_sig_bins & Out.(T).(c_names{2}).sig_n_bins <= n_sig_bins;
                within_range  = Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS < 60 & Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS < 60;
                hold on
                % plot only responsive for either rest, task, or both and lsline
                % this scatter
                scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), [], [0.7 0.7 0.7], 'filled')
                ls_line = lsline(gca);
                ls_line.Color = [0.3 0.3 0.3];
                if sum(sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name])) + ...
                        sum(sig.(cfg.condition(cond1_num).name)) + ...
                        sum(sig.(cfg.condition(cond2_num).name)) < 2
                    cc = nan(2);
                    pp = nan(2);
                    
                    p = nan(1,2);
                else
                    [cc, pp] = ...
                        corrcoef(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                        Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                        'Rows','complete');
                    [p,~] = ... % p(1) is linear slope, p(2) is b-member
                        polyfit(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                        Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range),1);
                end
                
                fitlm(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range))
                % plot all in grey and then by group
                s1 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS, [], [0.7 0.7 0.7], 'filled');
                s2 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS(sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name])), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS(sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name])), [], cond1_cond2_color, 'filled');
                s3 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond1_num).name)), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond1_num).name)), [], cfg.condition(cond2_num).color, 'filled');
                s4 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond2_num).name)), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond2_num).name)), [], cfg.condition(cond1_num).color, 'filled');
                box on
                axis square
                xlim([0 60])
                ylim([0 60])
                if numel(cc) > 1
                    title({[T ': N = ' num2str(unit_count)], ...
                        ['cc = ' num2str(cc(1,2)) '; p = ' num2str(pp(1,2))], ...
                        ['y = ' num2str(p(1)) '*x + ' num2str(p(2))]})
                end
                if a == 1
                    xlabel(['% signal change in ' cfg.condition(cond1_num).name])
                    ylabel(['% signal change in ' cfg.condition(cond2_num).name])
                    legend([ls_line s1 s2 s3 s4], ...
                        {'linear fit', 'non-significant', ['sig. ' cfg.condition(cond2_num).name ' & sig. ' cfg.condition(cond1_num).name], ['sig. only ' cfg.condition(cond2_num).name], ['sig. only ' cfg.condition(cond1_num).name]})
                end
                ax = gca;
                ax.FontSize = 16;
                set(gca,'XTick',0:10:60,'XTickLabel',0:10:60)
                set(gca,'YTick',0:10:60,'YTickLabel',0:10:60)
            end
            save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name '_Scatter_Pc_signal_change'],output_folder,savePlot)
        end
        clear cc pp p
        
        %% modulation index - absolute percent change - histograms
        for groupNum = 1:length(cfg.spk.compare_conditions)
            sig = struct;
            cond1_num = cfg.spk.compare_conditions{groupNum}(1);
            cond2_num = cfg.spk.compare_conditions{groupNum}(2);
            cond1_cond2_color = (cfg.condition(1).color + cfg.condition(2).color) / 2;
            
            figure
            set(gcf, 'Position', [2 332 1914 602])
            for a = 1: N_Areas
                T = Ana_TargetBrainArea{a};
                
                subplot(2,3,a)
                
                sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins > n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins > n_sig_bins; % both conditions
                %         sig_task_rest = Out.(T).(c_names{1}).sig_n_bins > n_sig_bins & Out.(T).(c_names{2}).sig_n_bins > n_sig_bins;
                sig.(cfg.condition(cond1_num).name)                  = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins <= n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins > n_sig_bins;
                %         sig_task      = Out.(T).(c_names{1}).sig_n_bins <= n_sig_bins & Out.(T).(c_names{2}).sig_n_bins > n_sig_bins;
                sig.(cfg.condition(cond2_num).name)                  = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins > n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins <= n_sig_bins;
                %         sig_rest      = Out.(T).(c_names{1}).sig_n_bins > n_sig_bins & Out.(T).(c_names{2}).sig_n_bins <= n_sig_bins;
                within_range  = Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS < 60 & Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS < 60;
                hold on
                % plot only responsive for either rest, task, or both and lsline
                % this scatter
                scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), [], [0.7 0.7 0.7], 'filled')
                ls_line = lsline(gca);
                ls_line.Color = [0.3 0.3 0.3];
                if sum(sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name])) + ...
                        sum(sig.(cfg.condition(cond1_num).name)) + ...
                        sum(sig.(cfg.condition(cond2_num).name)) < 2
                    cc = nan(2);
                    pp = nan(2);
                    
                    p = nan(1,2);
                else
                    [cc, pp] = ...
                        corrcoef(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                        Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                        'Rows','complete');
                    [p,~] = ... % p(1) is linear slope, p(2) is b-member
                        polyfit(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                        Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range),1);
                end
                
                fitlm(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range))
                % plot all in grey and then by group
                s1 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS, [], [0.7 0.7 0.7], 'filled');
                s2 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS(sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name])), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS(sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name])), [], cond1_cond2_color, 'filled');
                s3 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond1_num).name)), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond1_num).name)), [], cfg.condition(cond2_num).color, 'filled');
                s4 = scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond2_num).name)), Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS(sig.(cfg.condition(cond2_num).name)), [], cfg.condition(cond1_num).color, 'filled');
                box on
                axis square
                xlim([0 60])
                ylim([0 60])
                if numel(cc) > 1
                    title({[T ': N = ' num2str(unit_count)], ...
                        ['cc = ' num2str(cc(1,2)) '; p = ' num2str(pp(1,2))], ...
                        ['y = ' num2str(p(1)) '*x + ' num2str(p(2))]})
                end
                if a == 1
                    xlabel(['% signal change in ' cfg.condition(cond1_num).name])
                    ylabel(['% signal change in ' cfg.condition(cond2_num).name])
                    legend([ls_line s1 s2 s3 s4], ...
                        {'linear fit', 'non-significant', ['sig. ' cfg.condition(cond2_num).name ' & sig. ' cfg.condition(cond1_num).name], ['sig. only ' cfg.condition(cond2_num).name], ['sig. only ' cfg.condition(cond1_num).name]})
                end
                ax = gca;
                ax.FontSize = 16;
                set(gca,'XTick',0:10:60,'XTickLabel',0:10:60)
                set(gca,'YTick',0:10:60,'YTickLabel',0:10:60)
            end
            
            
        end
        
        %% modulation index - FR
        for groupNum = 1:length(cfg.spk.compare_conditions)
            cond1_num = cfg.spk.compare_conditions{groupNum}(1);
            cond2_num = cfg.spk.compare_conditions{groupNum}(2);
            figure
            set(gcf, 'Position', [2 381 1914 553])
            for a = 1: N_Areas
                T = Ana_TargetBrainArea{a};
                
                subplot(1,3,a)
                unit_count = ~isnan([Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_SubtrSDP_signed; Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_SubtrSDP_signed]');
                unit_count = sum(any(unit_count));
                %     [cc, pp] = corrcoef(Out.(T).Rest.FR_ModIndex_SubtrSDP_signed, Out.(T).Task.FR_ModIndex_SubtrSDP_signed,'Rows','complete');
                scatter(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_SubtrSDP_signed, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_SubtrSDP_signed)
                %     ls_line = lsline(gca);
                %     ls_line.Color = 'b';
                box on
                axis square
                xlim([-7 7])
                ylim([-7 7])
                %     title({[T ': N = ' num2str(unit_count)], ...
                %         ['cc = ' num2str(cc(2)) '; p = ' num2str(pp(2))]})
                title([T ': N = ' num2str(unit_count)])
                if a == 1
                    xlabel(['\Delta Firing Rate in ' cfg.condition(cond1_num).name ', Hz'])
                    ylabel(['\Delta Firing Rate in ' cfg.condition(cond2_num).name ', Hz'])
                    legend({'', '', '', ''})
                end
                
            end
            save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name '_Scatter_Modulation_Magnitude_signed'],output_folder,savePlot)
        end
        
        %% time of max
        for groupNum = 1:length(cfg.spk.compare_conditions)
            cond1_num = cfg.spk.compare_conditions{groupNum}(1);
            cond2_num = cfg.spk.compare_conditions{groupNum}(2);
            
            colors = [0 0 0; 0 0 1; 1 0 0; 1 0 1]; % black, blue, red, magenta
            
            for a = 1: N_Areas
                T = Ana_TargetBrainArea{a};
                
                unit_count = ~isnan([Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_SubtrSDP_signed, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_SubtrSDP_signed, ...
                    Out.(T).(cfg.condition(cond1_num).name).sig_time, Out.(T).(cfg.condition(cond2_num).name).sig_time]');
                valid_unit_ids = all(unit_count,1);
                unit_count = sum(valid_unit_ids);
                
                valid_cond1_times      = Out.(T).(cfg.condition(cond1_num).name).sig_time(valid_unit_ids);
                valid_cond2_times      = Out.(T).(cfg.condition(cond2_num).name).sig_time(valid_unit_ids);
                
                valid_cond1_bins       = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins(valid_unit_ids);
                valid_cond2_bins       = Out.(T).(cfg.condition(cond2_num).name).sig_n_bins(valid_unit_ids);
                
                non_sig       = valid_cond1_bins <= n_sig_bins & valid_cond2_bins <= n_sig_bins & ...
                    ~isnan(valid_cond1_times) & ~isnan(valid_cond2_times); % I have to search for non-significant explicitly as there are some nan units
                sig_cond2_cond1 = valid_cond1_bins > n_sig_bins & valid_cond2_bins > n_sig_bins;
                sig_cond2      = valid_cond1_bins <= n_sig_bins & valid_cond2_bins > n_sig_bins;
                sig_cond1      = valid_cond1_bins > n_sig_bins & valid_cond2_bins <= n_sig_bins;
                groups        = zeros(unit_count,1);
                
                groups(sig_cond1)      = 1;
                groups(sig_cond2)      = 2;
                groups(sig_cond2_cond1) = 3;
                
                % units counts
                non_sig_count       = sum(non_sig);
                sig_cond2_cond1_count = sum(sig_cond2_cond1);
                sig_cond2_count      = sum(sig_cond2);
                sig_cond1_count      = sum(sig_cond1);
                
                % compute histograms - rest
                [dt.(cfg.condition(cond1_num).name).task_rest_bins, centers] = hist(valid_cond1_times(sig_cond2_cond1), -250:50:250);%/sum(sig_task_rest);
                dt.(cfg.condition(cond1_num).name).rest_bins      = hist(valid_cond1_times(sig_cond1), -250:50:250);%/sum(sig_rest);
                dt.(cfg.condition(cond1_num).name).task_bins      = hist(valid_cond1_times(sig_cond2), -250:50:250);%/sum(sig_task);
                dt.(cfg.condition(cond1_num).name).nonsig_bins    = hist(valid_cond1_times(non_sig), -250:50:250);%/sum(sig_task);
                
                % compute histograms - task
                [dt.(cfg.condition(cond2_num).name).task_rest_bins, centers] = hist(valid_cond2_times(sig_cond2_cond1), -250:50:250);%/sum(sig_task_rest);
                dt.(cfg.condition(cond2_num).name).rest_bins      = hist(valid_cond2_times(sig_cond1), -250:50:250);%/sum(sig_rest);
                dt.(cfg.condition(cond2_num).name).task_bins      = hist(valid_cond2_times(sig_cond2), -250:50:250);%/sum(sig_task);
                dt.(cfg.condition(cond2_num).name).nonsig_bins    = hist(valid_cond2_times(non_sig), -250:50:250);%/sum(sig_task);
                
                if sum(valid_cond1_times) & sum(valid_cond2_times)
                
                    figure,
                    set(gcf, 'Position', [723 152 799 832])
                    
                    % plot the main scatter
                    s = scatterhist(valid_cond1_times, valid_cond2_times, ...
                        'Group', groups, ...
                        ...%'Style', 'bar', ...
                        'Marker', '.o+x', ...
                        'MarkerSize', 6, ...
                        'Location', 'SouthWest', ...
                        'Direction', 'out', ...
                        'Color', colors);
                    legend('Non Significant', [cfg.condition(cond1_num).name ' Sig.'], [cfg.condition(cond2_num).name ' Sig'], [cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'])
                    
                    title({[T ': N = ' num2str(unit_count)], ...
                        ['Non-sig / {\color{blue}' cfg.condition(cond1_num).name ' sig} / {\color{red}' cfg.condition(cond2_num).name ' sig} / {\color{magenta}' cfg.condition(cond2_num).name '&' cfg.condition(cond1_num).name ' sig}'], ...
                        [num2str(non_sig_count) ' / {\color{blue}' num2str(sig_cond1_count) '} / {\color{red}' num2str(sig_cond2_count) '} / {\color{magenta}' num2str(sig_cond2_cond1_count) '}']})
                    xlabel([cfg.condition(cond1_num).name ': Time of Max \DeltaFR, ms from R-peak'])
                    ylabel([cfg.condition(cond2_num).name ': Time of Max \DeltaFR, ms from R-peak'])
                    
                    xlim([-250 250])
                    ylim([-250 250])
                    axis square
                    box on
                    
                    % plot the bottom marginal histogram
                    bar_data = [dt.(cfg.condition(cond1_num).name).nonsig_bins; dt.(cfg.condition(cond1_num).name).rest_bins; dt.(cfg.condition(cond1_num).name).task_bins; dt.(cfg.condition(cond1_num).name).task_rest_bins]';
                    st_bar = bar(s(2), centers, bar_data);
                    ylabel(s(2), 'Unit Count')
                    xlim(s(2), [-250 250])
                    for ii = 1:4
                        st_bar(ii).FaceColor = colors(ii,:);
                    end
                    clear st_bar
                    
                    % plot the side marginal histogram
                    bar_data = [dt.(cfg.condition(cond2_num).name).nonsig_bins; dt.(cfg.condition(cond2_num).name).rest_bins; dt.(cfg.condition(cond2_num).name).task_bins; dt.(cfg.condition(cond2_num).name).task_rest_bins]';
                    st_bar = barh(s(3), centers, bar_data);
                    xlabel(s(3), 'Unit Count')
                    %     xlim(s(3), [0 max(bar_data, [], 'all')])
                    ylim(s(3), [-250 250])
                    for ii = 1:4
                        st_bar(ii).FaceColor = colors(ii,:);
                    end
                    clear st_bar
                    save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name T '_Scatter_Time_Max_Change'],output_folder,savePlot)
            
                end
            end
        end
        
        %% time of max - for all responses
        for groupNum = 1:length(cfg.spk.compare_conditions)
            cond1_num = cfg.spk.compare_conditions{groupNum}(1);
            cond2_num = cfg.spk.compare_conditions{groupNum}(2);
            
            figure,
            set(gcf, 'Position', [723 152 799 832])
            
            colors = [1 0.53 0;
                0.05 0.65 0.7;
                1 0 0.6];
            
            groups                             = [];
            all_valid_cond1_times              = [];
            all_valid_cond2_times              = [];
            all_unit_counts                    = [];
            dt.(cfg.condition(cond1_num).name) = struct('cond2_cond1_bins', [], 'cond1_bins', [], 'cond2_bins', [], 'nonsig_bins', []);
            dt.(cfg.condition(cond2_num).name) = struct('cond2_cond1_bins', [], 'cond1_bins', [], 'cond2_bins', [], 'nonsig_bins', []);
            
            for a = 1: N_Areas
                T = Ana_TargetBrainArea{a};
                
                sign_consistent = ...
                    eq(Out.(T).(cfg.condition(cond1_num).name).sig_sign, Out.(T).(cfg.condition(cond2_num).name).sig_sign);
                
                unit_count = ~isnan([Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_SubtrSDP_signed, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_SubtrSDP_signed, ...
                    Out.(T).(cfg.condition(cond1_num).name).sig_time, Out.(T).(cfg.condition(cond2_num).name).sig_time]');
                valid_unit_ids = all(unit_count,1);
                unit_count = sum(valid_unit_ids);
                
                % choose only valid units with consistent sign of response
                valid_cond1_times      = Out.(T).(cfg.condition(cond1_num).name).sig_time(valid_unit_ids);
                valid_cond2_times      = Out.(T).(cfg.condition(cond2_num).name).sig_time(valid_unit_ids);
                
                valid_cond1_bins       = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins(valid_unit_ids);
                valid_cond2_bins       = Out.(T).(cfg.condition(cond2_num).name).sig_n_bins(valid_unit_ids);
                
                % count units
                non_sig       = valid_cond1_bins <= n_sig_bins & valid_cond2_bins <= n_sig_bins & ...
                    ~isnan(valid_cond1_times) & ~isnan(valid_cond2_times); % I have to search for non-significant explicitly as there are some nan units
                sig_cond2_cond1 = valid_cond1_bins > n_sig_bins & valid_cond2_bins > n_sig_bins;
                sig_cond2      = valid_cond1_bins <= n_sig_bins & valid_cond2_bins > n_sig_bins;
                sig_cond1      = valid_cond1_bins > n_sig_bins & valid_cond2_bins <= n_sig_bins;
                
                % split data into groups by area and response (temporary)
                groups_tmp                  = zeros(unit_count,1); % 0 - nonsig
                groups_tmp(sig_cond1)       = 1;                   % 1 - sig cond1
                groups_tmp(sig_cond2)       = 2;                   % 2 - sig cond2
                groups_tmp(sig_cond2_cond1) = 3;                   % 3 - sig cond2 cond1
                
                % units counts
                non_sig_count         = sum(non_sig);
                sig_cond2_cond1_count = sum(sig_cond2_cond1);
                sig_cond2_count       = sum(sig_cond2);
                sig_cond1_count       = sum(sig_cond1);
                
                % compute histograms - rest
                out.(T).bin_counts_cond1 = hist(valid_cond1_times, -250:50:250);
                out.(T).bin_counts_cond2 = hist(valid_cond2_times, -250:50:250);
                
                % gather the necessary data together
                groups               = [groups; groups_tmp+4*(a-1)];
                all_valid_cond1_times = [all_valid_cond1_times; valid_cond1_times];
                all_valid_cond2_times = [all_valid_cond2_times; valid_cond2_times];
                all_unit_counts      = [all_unit_counts; unit_count];
                
                if sum(sig_cond2_cond1) < 2
                    cc(a) = NaN;
                    pp(a) = NaN;
                    continue
                end
                
                [cc_tmp, pp_tmp] = corrcoef(valid_cond1_times(sig_cond2_cond1), valid_cond2_times(sig_cond2_cond1));
                cc(a) = cc_tmp(2,1);
                pp(a) = pp_tmp(2,1);
                clear cc_tmp pp_tmp
                
            end
            
            s = scatterhist(all_valid_cond1_times, all_valid_cond2_times, ...
                'Group', groups, ...
                ...%'Style', 'bar', ...
                'Marker', '.v^s', ...
                'MarkerSize', 6, ...
                'Location', 'SouthWest', ...
                'Direction', 'out', ...
                'Color', [repmat(colors(1,:), 4, 1); repmat(colors(2,:), 4, 1); repmat(colors(3,:), 4, 1)]);
            hold on
            
            legend('Non Significant', [cfg.condition(cond1_num).name ' Sig.'], [cfg.condition(cond2_num).name ' Sig'], [cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'])
%             legend('VPL: Non Significant', ['VPL: ' cfg.condition(cond1_num).name ' Sig.'], ['VPL: ' cfg.condition(cond2_num).name ' Sig'], ['VPL: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'], ...
%                 'dPul: Non Significant', ['dPul: ' cfg.condition(cond1_num).name ' Sig.'], ['dPul: ' cfg.condition(cond2_num).name ' Sig'], ['dPul: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'], ...
%                 'MD: Non Significant', ['MD: ' cfg.condition(cond1_num).name ' Sig.'], ['MD: ' cfg.condition(cond2_num).name ' Sig'], ['MD: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'])
            xlabel([cfg.condition(cond1_num).name ': Time of Max \DeltaFR, ms from R-peak'])
            ylabel([cfg.condition(cond2_num).name ': Time of Max \DeltaFR, ms from R-peak'])
            title({'Only Units with Consistent Modulation Sign', ...
                ...
                ['\color[rgb]{1.0000 0.5300 0}VPL: N = ' num2str(all_unit_counts(1)) ...
                '; \color[rgb]{0.0500 0.6500 0.7000}dPul: N = ' num2str(all_unit_counts(2)) ...
                '; \color[rgb]{1.0000 0 0.6000}MD: N = ' num2str(all_unit_counts(3))], ...
                ...
                ['\color[rgb]{0 0 0}Only with Significant Modulations in both Rest and Task:'], ...
                ...
                ['\color[rgb]{1.0000 0.5300 0}cc = ' num2str(cc(1)) ' ; p = ' num2str(pp(1)) ...
                '; \color[rgb]{0.0500 0.6500 0.7000}cc = ' num2str(cc(2)) ' ; p = ' num2str(pp(2)) ...
                '; \color[rgb]{1.0000 0 0.6000}cc = ' num2str(cc(3)) ' ; p = ' num2str(pp(3))]})
            
            xlim([-250 250])
            ylim([-250 250])
            axis square
            box on
            
            % plot the bottom marginal histogram
            bar_data = ...
                [out.(unqTargets{1}).bin_counts_cond1 / sum(out.(unqTargets{1}).bin_counts_cond1); ...
                out.(unqTargets{2}).bin_counts_cond1 / sum(out.(unqTargets{2}).bin_counts_cond1); ...
                out.(unqTargets{3}).bin_counts_cond1 / sum(out.(unqTargets{3}).bin_counts_cond1)]';
            st_bar1 = bar(s(2), centers, bar_data);
            ylabel(s(2), 'Unit Fraction')
            ylim(s(2), [0 0.6])
            % xlim(s(2), [-250 250])
            for ii = 1:3
                st_bar1(ii).FaceColor = colors(ii,:);
            end
            
            % plot the side marginal histogram
            bar_data = [out.(unqTargets{1}).bin_counts_cond2 / sum(out.(unqTargets{1}).bin_counts_cond2); ...
                out.(unqTargets{2}).bin_counts_cond2 / sum(out.(unqTargets{2}).bin_counts_cond2); ...
                out.(unqTargets{3}).bin_counts_cond2 / sum(out.(unqTargets{3}).bin_counts_cond2)]';
            st_bar2 = barh(s(3), centers, bar_data);
            xlabel(s(3), 'Unit Fraction')
            xlim(s(3), [0 0.6])
            % ylim(s(3), [-250 250])
            for ii = 1:3
                st_bar2(ii).FaceColor = colors(ii,:);
            end
            save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name '_All_Areas_Scatter_Time_Max_Change'],output_folder,savePlot)
        end
        
        %% time of max - only for consistent sign of responses
        for groupNum = 1:length(cfg.spk.compare_conditions)
            cond1_num = cfg.spk.compare_conditions{groupNum}(1);
            cond2_num = cfg.spk.compare_conditions{groupNum}(2);
            
            figure,
            set(gcf, 'Position', [723 152 799 832])
            
            colors = [1 0.53 0;
                0.05 0.65 0.7;
                1 0 0.6];
            
            groups                             = [];
            all_valid_cond1_times              = [];
            all_valid_cond2_times              = [];
            all_unit_counts                    = [];
            dt.(cfg.condition(cond1_num).name) = struct('cond2_cond1_bins', [], 'cond1_bins', [], 'cond2_bins', [], 'nonsig_bins', []);
            dt.(cfg.condition(cond2_num).name) = struct('cond2_cond1_bins', [], 'cond1_bins', [], 'cond2_bins', [], 'nonsig_bins', []);
            
            for a = 1: N_Areas
                T = Ana_TargetBrainArea{a};
                
                sign_consistent = ...
                    eq(Out.(T).(cfg.condition(cond1_num).name).sig_sign, Out.(T).(cfg.condition(cond2_num).name).sig_sign);
                
                unit_count = ~isnan([Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_SubtrSDP_signed, Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_SubtrSDP_signed, ...
                    Out.(T).(cfg.condition(cond1_num).name).sig_time, Out.(T).(cfg.condition(cond2_num).name).sig_time]') & sign_consistent';
                valid_unit_ids = all(unit_count,1);
                unit_count = sum(valid_unit_ids);
                
                % choose only valid units with consistent sign of response
                valid_cond1_times      = Out.(T).(cfg.condition(cond1_num).name).sig_time(valid_unit_ids & sign_consistent');
                valid_cond2_times      = Out.(T).(cfg.condition(cond2_num).name).sig_time(valid_unit_ids & sign_consistent');
                
                valid_cond1_bins       = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins(valid_unit_ids & sign_consistent');
                valid_cond2_bins       = Out.(T).(cfg.condition(cond2_num).name).sig_n_bins(valid_unit_ids & sign_consistent');
                
                % count units
                non_sig       = valid_cond1_bins <= n_sig_bins & valid_cond2_bins <= n_sig_bins & ...
                    ~isnan(valid_cond1_times) & ~isnan(valid_cond2_times); % I have to search for non-significant explicitly as there are some nan units
                sig_cond2_cond1 = valid_cond1_bins > n_sig_bins & valid_cond2_bins > n_sig_bins;
                sig_cond2      = valid_cond1_bins <= n_sig_bins & valid_cond2_bins > n_sig_bins;
                sig_cond1      = valid_cond1_bins > n_sig_bins & valid_cond2_bins <= n_sig_bins;
                
                % split data into groups by area and response (temporary)
                groups_tmp                  = zeros(unit_count,1); % 0 - nonsig
                groups_tmp(sig_cond1)       = 1;                   % 1 - sig cond1
                groups_tmp(sig_cond2)       = 2;                   % 2 - sig cond2
                groups_tmp(sig_cond2_cond1) = 3;                   % 3 - sig cond2 cond1
                
                % units counts
                non_sig_count         = sum(non_sig);
                sig_cond2_cond1_count = sum(sig_cond2_cond1);
                sig_cond2_count       = sum(sig_cond2);
                sig_cond1_count       = sum(sig_cond1);
                
                % compute histograms - rest
                out.(T).bin_counts_cond1 = hist(valid_cond1_times, -250:50:250);
                out.(T).bin_counts_cond2 = hist(valid_cond2_times, -250:50:250);
                
                % gather the necessary data together
                groups               = [groups; groups_tmp+4*(a-1)];
                all_valid_cond1_times = [all_valid_cond1_times; valid_cond1_times];
                all_valid_cond2_times = [all_valid_cond2_times; valid_cond2_times];
                all_unit_counts      = [all_unit_counts; unit_count];
                
                if sum(sig_cond2_cond1) < 2
                    cc(a) = NaN;
                    pp(a) = NaN;
                    continue
                end
                
                [cc_tmp, pp_tmp] = corrcoef(valid_cond1_times(sig_cond2_cond1), valid_cond2_times(sig_cond2_cond1));
                cc(a) = cc_tmp(2,1);
                pp(a) = pp_tmp(2,1);
                clear cc_tmp pp_tmp
                
            end
            
            s = scatterhist(all_valid_cond1_times, all_valid_cond2_times, ...
                'Group', groups, ...
                ...%'Style', 'bar', ...
                'Marker', '.v^s', ...
                'MarkerSize', 6, ...
                'Location', 'SouthWest', ...
                'Direction', 'out', ...
                'Color', [repmat(colors(1,:), 4, 1); repmat(colors(2,:), 4, 1); repmat(colors(3,:), 4, 1)]);
            hold on
            
            legend('Non Significant', [cfg.condition(cond1_num).name ' Sig.'], [cfg.condition(cond2_num).name ' Sig'], [cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'])
%             legend('VPL: Non Significant', ['VPL: ' cfg.condition(cond1_num).name ' Sig.'], ['VPL: ' cfg.condition(cond2_num).name ' Sig'], ['VPL: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'], ...
%                 'dPul: Non Significant', ['dPul: ' cfg.condition(cond1_num).name ' Sig.'], ['dPul: ' cfg.condition(cond2_num).name ' Sig'], ['dPul: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'], ...
%                 'MD: Non Significant', ['MD: ' cfg.condition(cond1_num).name ' Sig.'], ['MD: ' cfg.condition(cond2_num).name ' Sig'], ['MD: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'])
            xlabel([cfg.condition(cond1_num).name ': Time of Max \DeltaFR, ms from R-peak'])
            ylabel([cfg.condition(cond2_num).name ': Time of Max \DeltaFR, ms from R-peak'])
            title({'Only Units with Consistent Modulation Sign', ...
                ...
                ['\color[rgb]{1.0000 0.5300 0}VPL: N = ' num2str(all_unit_counts(1)) ...
                '; \color[rgb]{0.0500 0.6500 0.7000}dPul: N = ' num2str(all_unit_counts(2)) ...
                '; \color[rgb]{1.0000 0 0.6000}MD: N = ' num2str(all_unit_counts(3))], ...
                ...
                ['\color[rgb]{0 0 0}Only with Significant Modulations in both Rest and Task:'], ...
                ...
                ['\color[rgb]{1.0000 0.5300 0}cc = ' num2str(cc(1)) ' ; p = ' num2str(pp(1)) ...
                '; \color[rgb]{0.0500 0.6500 0.7000}cc = ' num2str(cc(2)) ' ; p = ' num2str(pp(2)) ...
                '; \color[rgb]{1.0000 0 0.6000}cc = ' num2str(cc(3)) ' ; p = ' num2str(pp(3))]})
            
            xlim([-250 250])
            ylim([-250 250])
            axis square
            box on
            
            % plot the bottom marginal histogram
            bar_data = ...
                [out.(unqTargets{1}).bin_counts_cond1 / sum(out.(unqTargets{1}).bin_counts_cond1); ...
                out.(unqTargets{2}).bin_counts_cond1 / sum(out.(unqTargets{2}).bin_counts_cond1); ...
                out.(unqTargets{3}).bin_counts_cond1 / sum(out.(unqTargets{3}).bin_counts_cond1)]';
            st_bar1 = bar(s(2), centers, bar_data);
            ylabel(s(2), 'Unit Fraction')
            ylim(s(2), [0 0.6])
            % xlim(s(2), [-250 250])
            for ii = 1:3
                st_bar1(ii).FaceColor = colors(ii,:);
            end
            
            % plot the side marginal histogram
            bar_data = [out.(unqTargets{1}).bin_counts_cond2 / sum(out.(unqTargets{1}).bin_counts_cond2); ...
                out.(unqTargets{2}).bin_counts_cond2 / sum(out.(unqTargets{2}).bin_counts_cond2); ...
                out.(unqTargets{3}).bin_counts_cond2 / sum(out.(unqTargets{3}).bin_counts_cond2)]';
            st_bar2 = barh(s(3), centers, bar_data);
            xlabel(s(3), 'Unit Fraction')
            xlim(s(3), [0 0.6])
            % ylim(s(3), [-250 250])
            for ii = 1:3
                st_bar2(ii).FaceColor = colors(ii,:);
            end
            save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name '_All_Areas_Scatter_Time_Max_Change_SignConsistent'],output_folder,savePlot)
        end
        
        %% histogram of modulation indices - % signal change
        
        figure
        set(gcf, 'Position', [2 381 1914 553])
        
        for a = 1: N_Areas
            T = Ana_TargetBrainArea{a};
            
            for c=1:N_conditions
                L=cfg.condition(c).name;
                
                subplot(N_conditions,N_Areas,a + N_Areas*(c-1))
                colororder([0.8500 0.3250 0.0980; 0 0.4470 0.7410; 1 1 1])
                MIs_Pc_all  = Out.(T).(L).FR_ModIndex_PcS_signed;
                pos_idx     = Out.(T).(L).sig_sign > 0;
                neg_idx     = Out.(T).(L).sig_sign < 0;
                sign_all    = Out.(T).(L).sig_n_bins > n_sig_bins;
                
                non_sig = histcounts(MIs_Pc_all(~sign_all), -60:10:60);
                sig_pos = histcounts(MIs_Pc_all(pos_idx & sign_all), -60:10:60);
                sig_neg = histcounts(MIs_Pc_all(neg_idx & sign_all), -60:10:60);
                
                bar(-55:10:60, [sig_pos; sig_neg; non_sig], 'stacked')
                
                if c == 1
                    title(T)
                end
                
                title([L ': ' T])
                
                if a == 1 && c == 1
                    legend({'Sig. Pos.', 'Sig. Neg.', 'Non-Significant'}, 'Location', 'Best')
                    xlabel('% signal change')
                    ylabel('Unit Counts')
                end
                
            end
            
        end
        
        sgtitle('Modulation Magnitude in Rest and Task')
        save_figure_as('Histograms_Magnitude_Signal_Change_Significant',output_folder,savePlot)
        
        %% plot correlation coefficients between rest and task
        hist_bins        = -1:0.1:1;
        hist_bin_centers = hist_bins + 0.05;
        b_colors         = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4645 0.5705 0.4330; 1 1 1]; % blue, white, yellow, green
        
        figure
        set(gcf, 'Position', [2 556 1914 378])
        
        for a = 1: N_Areas
            T = Ana_TargetBrainArea{a};
            
            % organize data
            cc         = Out.(T).cc_rest_vs_task;
            ids_cc_sig = Out.(T).pp_rest_vs_task < 0.05;
            sig        = (~isnan(Out.(T).Rest.sig_FR_diff) & ~isnan(Out.(T).Task.sig_FR_diff) ) & ...
                ( (Out.(T).Rest.sig_n_bins > n_sig_bins) | (Out.(T).Task.sig_n_bins > n_sig_bins) );
            
            % ids for unit groups
            noCorr_nonResp  = ~ids_cc_sig(:) & ~sig(:);
            noCorr_sigResp  = ~ids_cc_sig(:) & sig(:);
            sigCorr_sigResp = ids_cc_sig(:) & sig(:);
            sigCorr_noResp  = ids_cc_sig(:) & ~sig(:);
            
            % unit counts
            counts_noCorr_sigResp  = histc(cc(noCorr_sigResp),hist_bins);  % blue
            counts_sigCorr_sigResp = histc(cc(sigCorr_sigResp),hist_bins); % green
            counts_sigCorr_noResp  = histc(cc(sigCorr_noResp),hist_bins);  % yellow
            counts_noCorr_nonResp  = histc(cc(noCorr_nonResp),hist_bins);  % white
            
            counts_noCorr_sigResp = counts_noCorr_sigResp(:);
            counts_sigCorr_sigResp = counts_sigCorr_sigResp(:);
            counts_sigCorr_noResp = counts_sigCorr_noResp(:);
            counts_noCorr_nonResp = counts_noCorr_nonResp(:);
            
            subplot(1,3,a)
            b = bar(hist_bin_centers,[counts_noCorr_sigResp counts_sigCorr_sigResp counts_sigCorr_noResp counts_noCorr_nonResp],'stacked'); %
            
            for ii = 1:4
                b(ii).FaceColor = b_colors(ii,:);
            end
            
            xlim([-1 1])
            
            if a == 1
                xlabel('Corr.coef. Rest vs. Task')
                ylabel('Unit Counts')
            end
            legend('no corr + sig resp', 'sig corr + sig resp', 'sig corr + no resp', 'no corr + no resp')
            title(T)
            
        end
        
        save_figure_as('Histograms_Corr_Coef_Rest_vs_Task',output_folder,savePlot)
        
    end
    clear Out
end
end

function rounded=roundto(unrounded,n)
factor= 10^n;
rounded=round(unrounded*factor)/factor;
end

function make_correlation_plot(X,Y,sig,onlysignificant,ccol,ttl)
x=X(:);y=Y(:);
axis square; box on;hold on;
% all not significant & NaN units
if ~isempty(x) && ~isempty(y)
    if ~ onlysignificant
        scatter(x(~sig) , y(~sig), 'filled', 'MarkerFaceColor',ccol/2) % make nonsignificant circles empty maybe
    end
    scatter(x(sig) , y(sig), 'filled', 'MarkerFaceColor',ccol)
    
    if onlysignificant && sum(sig)>2
        [coef, pval] = corr(x(sig),y(sig), 'rows','complete') ;
        [p,S] = polyfit(x(~isnan(y) & sig),y(~isnan(y)  & sig),1); % fit all?
        [y_fit,delta] = polyval(p,x(~isnan(y) & sig),S);
        plot(x(~isnan(y)& sig), y_fit,'LineWidth', 2, 'Color', ccol);
    elseif ~onlysignificant && numel(x)>2
        [coef, pval] = corr(x,y, 'rows','complete') ;
        [p,S] = polyfit(x(~isnan(y)),y(~isnan(y)),1); % fit all?
        [y_fit,delta] = polyval(p,x(~isnan(y)),S);
        plot(x(~isnan(y)), y_fit,'LineWidth', 2, 'Color', ccol);
    else
        coef=NaN;pval=NaN;
    end
    ttl_set = get(gca, 'Title');
    if isempty(ttl_set.String)
        title({ttl})
        ttl_set = get(gca, 'Title');
    end
    ttl_str = [ttl_set.String; {['\color[rgb]{' num2str(ccol) '}cc = ' num2str(roundto(coef,2)) '; p = ' num2str(roundto(pval,3))]}];
    title(ttl_str)
end
end

function save_figure_as(filename,basepath_to_save,savePlot)
if savePlot;
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(gcf);
end
end
