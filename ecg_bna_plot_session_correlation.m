function ecg_bna_plot_session_correlation(session_info,cfg)

matlab_year = version('-release');
matlab_year = str2double(matlab_year(1:end-1));

dataFolder = [cfg.SPK_root_results_fldr filesep 'correlation_analysis' filesep];

condition_colors={cfg.condition.color};

fileList = dir([dataFolder session_info.Monkey(1:3) '_' session_info.Date '*correlation.mat']);

for untNum = 1:length(fileList)
    
    load([dataFolder filesep fileList(untNum).name], 'data')

    %% title with all the spike analysis parameters
    sgtitleText = {[data.unitId '_' data.target '; ch ' num2str(data.channel) ';'  ' unit ' data.unit], ... %  
        ['FR, Hz: ' num2str(data.FR) '; SNR: ' num2str(data.quantSNR) '; Stability ratio: ' num2str(data.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * data.Single_rating) '%'], ...
        ['Rest: Pearson r (lag 0) = ' num2str(data.Rest.pearson_r(4)) ' p = ' num2str(data.Rest.pearson_p(4)) ...
            '; Task: Pearson r (lag 0) = ' num2str(data.Task.pearson_r(4)) '; p = ' num2str(data.Task.pearson_p(4))]};
    
    %% rasters with sorting by RR duration
    f1 = figure;
    set(f1, 'Position', [376 64 1140 928])
    
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        col=cfg.condition(c).color;
        
        % compute raster
        AT        = data.(L).AT_one_stream;
        RR_starts = data.(L).valid_RRinterval_starts;
        RR_ends   = data.(L).valid_RRinterval_ends;
        
        % drop RRs before the 1st spike or after the last one
        drop_RR_ids = RR_starts < AT(1) | RR_ends > AT(end);
        RR_starts(drop_RR_ids) = [];
        RR_ends(drop_RR_ids)   = [];
        
        RR_starts              = RR_starts(~isnan(RR_starts));
        RR_ends                = RR_ends(~isnan(RR_ends));
        
        % select spikes per cardiac cycle
        raster = zeros(length(RR_starts),801);
        for cycleNum = 1:length(RR_starts)
            
            currCycleSpikeIds  = AT > RR_starts(cycleNum) & AT < RR_ends(cycleNum);
            raster(cycleNum,:) = histc(AT(currCycleSpikeIds)-RR_starts(cycleNum), 0:0.001:0.8);
            
        end
        
        % compute scatter input
        subplot(2,numel(cfg.condition),2*(c-1)+1)
        [raster_row, raster_col] = find(raster);
        scatter(raster_col,raster_row,2,col);
        hold on
        scatter(1000*(RR_ends-RR_starts),1:length(RR_starts), 'Marker','^')
        box on
        xlim([0 550])
        if ~isempty(raster)
            ylim([0 size(raster,1)])
        end
        xlabel('Time after R-peak, ms')
        ylabel('R-peak #')
        title([L ': Original R-peak Order'])
        
        subplot(2,numel(cfg.condition),2*c)
        [RRdur_sorted,ids] = sort(1000*(RR_ends-RR_starts));
        raster_sorted = raster(ids,:);
        [raster_row, raster_col] = find(raster_sorted);
        %scatter(raster_col,raster_row,2,'MarkerFaceColor',col,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',col,'MarkerEdgeAlpha',0.2);
        scatter(raster_col,raster_row,2,col);
        hold on
        scatter(RRdur_sorted,1:length(RR_starts), 'Marker','^')
        box on
        xlim([0 550])
        if ~isempty(raster)
            ylim([0 size(raster,1)])
        end
        title([L ': Sorted by RR duration from Low of High'])
        
    end
    
    filename= ['RastersSorted__' data.unitId, '_' data.target];
    saveas(f1,[dataFolder, filesep, filename], 'fig')
    saveas(f1,[dataFolder, filesep, filename], 'tiff')
    close(f1);
    
    
    %% results of correlation analysis
    f5 = figure;
    set(f5, 'Position', [381 424 1450 398])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        subplot(1,numel(cfg.condition)+2,c)
        scatter(data.(L).FRbyRR_Hz, data.(L).cycleDurations_s, [], condition_colors{c}(1:3), 'Marker', '.')
        hold on
        linTrend = lsline(gca);
        linTrend.Color = condition_colors{c}(1:3);
        xlabel('Firing Rate per Heart Cycle, Hz')
        ylabel('Heart Cycle Duration, s')
        box on
        title([L ': cc = ' num2str(data.(L).pearson_r(4)) '; p = ' num2str(data.(L).permuted_p(4))])
        legend({'Real Data', 'Least-Square Fit'}, 'Location', 'Best')
    end
    
    subplot(1,numel(cfg.condition)+2,numel(cfg.condition)+2)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        plot(data.cc_lag_list, data.(L).pearson_r, '-o', 'Color', condition_colors{c}(1:3))
    end
    ylim([-0.4 0.4])
    xlabel('Lag: Number of Heart Cycles Shifted')
    ylabel('Correlation Coefficient')
    legend({cfg.condition.name})
    
    filename= ['FR_RR_Correlations__' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% unit FR per heart cycle as a function of time
    if 1
        figure;
        set(gcf, 'Position', [1 41 1920 482])
        
        for c=1:numel(cfg.condition)
            L=cfg.condition(c).name;
%             yyaxis left
%             set(gca,'YColor',cfg.condition(c).color);
            subplot(2,1,1)
            hold on
            plot(data.(L).timeRRstart, data.(L).FRbyRR_Hz, 'o-', 'Color', cfg.condition(c).color)
            ylabel('FR per Heart Cycle, Hz')
%             yyaxis right
%             set(gca,'YColor',[0 0 0]);
%             hold on
            subplot(2,1,2)
            plot(data.(L).timeRRstart, data.(L).cycleDurations_s, 'kx--')
            ylabel('Heart-Cycle Duration, s')
            box on
            
            xlim_lows(c)  = data.(L).timeRRstart(1);
            xlim_highs(c) = data.(L).timeRRstart(end);
            
            %legend({'FR rest', 'cycle duration rest','FR task', 'cycle duration task'})
        end
        xlabel('Time from Session Start, s')
        if ~isnan(min(xlim_lows)) && ~isnan(max(xlim_highs))
            xlim([min(xlim_lows) max(xlim_highs)])
        end
        if matlab_year >= 2020
            sgtitle(sgtitleText, 'interpreter', 'none')
        end
        
        filename= ['Time_vs_FRperCycle_' data.unitId, '_' data.target];
        export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
        close(gcf);
        
    end
    
    %% unit FR by RR-cycle quintiles
    figure,
    set(gcf, 'Position', [381 201 1450 795])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        
        RRdur_quintiles = prctile(data.(L).cycleDurations_s, 0:20:100);
        
        avgFR_per_quintile = nan(length(RRdur_quintiles)-1,1);
        SEM_per_quintile   = nan(length(RRdur_quintiles)-1,1);
        for quintNum = 1:length(RRdur_quintiles)-1
            curr_RR_ids     = data.(L).cycleDurations_s > RRdur_quintiles(quintNum) & data.(L).cycleDurations_s < RRdur_quintiles(quintNum+1);
            avgFR_per_quintile(quintNum) = nanmean(data.(L).FRbyRR_Hz(curr_RR_ids));
            SEM_per_quintile(quintNum) = sterr(data.(L).FRbyRR_Hz(curr_RR_ids));
        end
        
        subplot(2,numel(cfg.condition), c)
        errorbar(1:5, avgFR_per_quintile, SEM_per_quintile, 'Color', cfg.condition(c).color)
        xlim([0 6])
        box on
        title([cfg.condition(c).name ': Pearson r (lag 0) = ' num2str(data.(L).pearson_r(4)) ' p = ' num2str(data.(L).pearson_p(4))])
        xlabel('RR duration Quintile Number')
        ylabel('Firing Rate, Hz')
        
        % RR durations sorted in ascending order vs. unit FR
        subplot(2,numel(cfg.condition), c + 2)
        [currRRdur, sort_ids] = sort(data.(L).cycleDurations_s);
        currFR                = data.(L).FRbyRR_Hz(sort_ids);
        RRbins                = min(data.(L).cycleDurations_s):0.01:max(data.(L).cycleDurations_s);
        avgFR_per_bin = NaN;
        for binNum = 1:length(RRbins)-1
            curr_RR_ids = currRRdur > RRbins(binNum) & currRRdur < RRbins(binNum+1);
            avgFR_per_bin(binNum) = mean(currFR(curr_RR_ids));
            semFR_per_bin(binNum) = sterr(currFR(curr_RR_ids));
        end
        
        plot(currRRdur, currFR, 'Color', cfg.condition(c).color/2)
        hold on
        if ~isnan(avgFR_per_bin)
            errorbar(RRbins(1:end-1)+0.005, avgFR_per_bin, semFR_per_bin, 'Color', cfg.condition(c).color, 'Marker', 'o', 'LineWidth', 2)
        end
        xlabel('RR durations (s) in ascending order')
        ylabel('FR in single cardiac cycles, Hz')
        legend('For Individual Cardiac Cycles', 'Averaged by Bin')
        clear curr_RR_ids avgFR_per_bin semFR_per_bin
    end
    
    filename= ['RRquintile_vs_FRperCycle_' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
end
