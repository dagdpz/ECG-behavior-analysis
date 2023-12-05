function ecg_bna_plot_session_ECG_related_spikePhase(session_info,cfg)

dataFolder = [cfg.SPK_root_results_fldr filesep 'cardioballistic' filesep];

Fs = 2.441406250000000e+04; % sampling frequency of BB signal, Hz
wf_times_ms = 1000 * (1/Fs:1/Fs:32/Fs); % in ms
wf_times_interp_ms = 1000 * (1/4/Fs:1/4/Fs:32/Fs); % in ms
peak_id = 10; % sample number of the trough in the spike waveform
bins = linspace(0, 2*pi, cfg.spk.N_phase_bins);
bin_centers   = 2*pi/cfg.spk.N_phase_bins:2*pi/cfg.spk.N_phase_bins:2*pi;

condition_colors={cfg.condition.color};

fileList = dir([dataFolder session_info.Monkey(1:3) '_' session_info.Date '*spikes_ECGphase.mat']);

for untNum = 1:length(fileList)
    
    load([dataFolder filesep fileList(untNum).name], 'data')
    
    % with the current Fs and number of samples per spike I have 1.4 ms per
    % spike
    
    %% AMP - absolute magnitude at the trough
    %     %% HW - in ms corresponding to the total time that the EAP trough is below half AMP
    %     HAMP = AMP_abs/2; % half of max amplitude
    %     beyond_halfWidth_idx = bsxfun(@(x,y) sign_trough*x>y, Y(15:40,:), HAMP); % multiply by signum: if spike is negative it will flip, if positiv - not
    % %     figure, imagesc(beyond_halfWidth_idx')
    %     HW = sum(beyond_halfWidth_idx,1)*(1/3/Fs); % in seconds
    %
    %     HWbyPhase      = arrayfun(@(x) mean(HW(bin == x)), unique(bin));
    %     HWbyPhase_std  = arrayfun(@(x) std(HW(bin == x)), unique(bin));
    
    %%
    sgtitleText = {[data.unitId '_' data.target '; ch ' num2str(data.channel) ';'  ' unit ' data.unit], ... %  
        ['FR, Hz: ' num2str(data.FR) '; SNR: ' num2str(data.quantSNR) '; Fano Factor: ' num2str(data.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * data.Single_rating) '%'], ...
        ['Task AMP_MI = ' num2str(data.Task.AMP_MI(1)) '; Task p = ' num2str(data.Task.AMP_MI(2))]};
    
    %% Overall picture of waveforms and PSTHs
    f1 = figure;
    set(f1, 'Position', [274 148 1452 797])
    %sgtitle(sgtitleText, 'interpreter', 'none')
    
    ax_sp = subplot(2,3,1);
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        line(wf_times_interp_ms, data.(L).waveforms_upsampled_microvolts(1:300:end, :), 'Color', [condition_colors{c} 0.1])
        line(repmat([wf_times_interp_ms(1) wf_times_interp_ms(end)],4,1)', repmat(data.thresholds_microV,1,2)', 'Color', 'r')
    end
    xlim([wf_times_interp_ms(1) wf_times_interp_ms(end)])
    xlabel('Time, ms')
    ylabel('Voltage, μV')
    title({'Example Waveforms: ', ...
        ['\color{blue}Rest: ' num2str(size(data.Rest.waveforms_upsampled_microvolts(1:300:end,:),1)) ' out of ' num2str(size(data.Rest.waveforms_upsampled_microvolts,1))], ...
        ['\color{red}Task: ' num2str(size(data.Task.waveforms_upsampled_microvolts(1:300:end,:),1)) ' out of ' num2str(size(data.Task.waveforms_upsampled_microvolts,1))]})
    
    subplot(2,3,2)
    box on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        h = hist(data.(L).spike_phases_radians, bin_centers);
        h = h / mean(h);
        %polarplot([phase_bin_centers phase_bin_centers(1)], [h h(1)], 'Color', condition_colors{c}(1:3))
        p=polar([bin_centers bin_centers(1)],[h h(1)]);%, 'Color', condition_colors{c}(1:3))
        set(p, 'Color', condition_colors{c}(1:3));
        hold on
    end
    hold off
    title({'ECG-triggered polar PSTH (spike counts / mean spike counts per bin)', ...
        '\color{blue}Rest', '\color{red}Task'})
    axis tight
    
    subplot(2,3,4)
    title('Average WFs by bin: \color{blue}Rest and \color{red}Task')
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        if sum(isnan(data.(L).waveforms_byBin_microvolts))==0
            line(wf_times_interp_ms, data.(L).waveforms_byBin_microvolts, 'Color', [condition_colors{c} 0.1])
        end
    end
    hold off
    box on
    xlabel('Spike Time, ms')
    ylabel('Voltage, μV')
    set(gca, 'XLim', get(ax_sp,'XLim'), 'YLim', get(ax_sp,'YLim'))
    
    subplot(2,3,5)
    imagesc(bsxfun(@minus,data.Rest.waveforms_byBin_microvolts,mean(data.Rest.waveforms_byBin_microvolts,1)))
    cbar = colorbar;
    cbar.Title.String = '\Delta Spike Amplitude, μV';
    caxis([-3 3])
    title({'{\color{blue}[Rest]}', 'WF Voltage Change over Heart Cycle'})
    xlabel('Spike Time')
    ylabel('Heart-Cycle Phase')
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
    
    subplot(2,3,6)
    imagesc(bsxfun(@minus,data.Task.waveforms_byBin_microvolts,mean(data.Task.waveforms_byBin_microvolts,1)))
    cbar = colorbar;
    cbar.Title.String = '\Delta Spike Amplitude, μV';
    caxis([-3 3])
    title({'{\color{red}[Task]}', 'WF Voltage Change over Heart Cycle'})
    xlabel('Spike Time')
    ylabel('Heart-Cycle Phase')
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
    
    filename= ['Spikes_PolarPSTH__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Distributions of features
    f2 = figure;
    set(f2, 'Position', [784   148   942   797])
    %sgtitle(sgtitleText, 'interpreter', 'none')
    
    subplot(2,2,1)
    box on
%     b = [histcounts(abs(data.Rest.AMP_microV), [0:15:350]); ...
%         histcounts(abs(data.Task.AMP_microV), [0:15:350])];
%     b = ( b ./ max(b,[],2) )';
%     barplot = bar([7.5:15:350], b);
%     for c = 1:2
%         barplot(c).FaceColor = condition_colors{c}(1:3);
%     end
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        h=hist(abs(data.(L).AMP_microV), [0:5:350]);
        plot(0:5:350,h,'Color', condition_colors{c})
    end
    hold off
    line([data.thresholds_microV(1) data.thresholds_microV(1)],ylim, 'Color', 'k')
    line([data.thresholds_microV(2) data.thresholds_microV(2)],ylim, 'Color', 'k')
    xlabel('AMP, μV');
    ylabel('Spike Counts')
    xlim([0 350])
    hold off
    title('AMP: Same X-axis for All Units')
    legend({cfg.condition.name}, 'Location', 'Best')
    
    subplot(2,2,2)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        h=hist(abs(data.(L).HW_ms), [0:0.005:0.5]);
        plot([0:0.005:0.5],h,'Color', condition_colors{c})
    end
    hold off
    xlim([0 0.5])
    xlabel('HW, ms');
    ylabel('Spike Counts')
    title('HW: Same X-axis for All Units')
    
    subplot(2,2,3)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        h=hist(abs(data.(L).TPW_ms), [0.2:0.005:1]);
        plot([0.2:0.005:1],h,'Color', condition_colors{c})
    end
    hold off
    xlabel('TPW, ms');
    ylabel('Spike Counts')
    title('TPW: X-axis Adjusts for Each Unit')
    
    subplot(2,2,4)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        [h,bins]=hist(abs(data.(L).REP_ms),[0:0.005:0.5]);
        plot(bins,h,'Color', condition_colors{c})
    end
    hold off
    xlabel('REP, ms');
    ylabel('Spike Counts')
    title('REP: X-axis Adjusts for Each Unit')
    
    filename= ['Distributions_AMP_HW_TPW_REP__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Real and reshuffled data
    f3 = figure;
    set(f3, 'Position', [1 41 1920 963])
    %sgtitle(sgtitleText, 'interpreter', 'none')
    feature_list = {'AMP', 'HW', 'TPW', 'REP'};
    feature_bin_list = {'AMP_microV_byBin', 'HW_ms_byBin', 'TPW_ms_byBin', 'REP_ms_byBin'};
    smoothed_feature_list = {'AMP_microV_byBin_smoothed', 'HW_ms_byBin_smoothed', 'TPW_ms_byBin_smoothed', 'REP_ms_byBin_smoothed'};
    subplot_numbers = {[1 2], [4 5], [11 12], [14 15]};
    
    for featureNum = 1:4
        for c=1:numel(cfg.condition)
            L=cfg.condition(c).name;
            subplot(3,5,subplot_numbers{featureNum}(c))
            %yyaxis left
            histogram(abs(data.(L).spike_phases_radians), bins, 'FaceColor', condition_colors{c}(1:3))
            ylabel('Spike Counts')
            currYLims = get(gca, 'YLim');
            currYLims(2) = 2 * currYLims(2);
            set(gca, 'YLim', currYLims)
            %yyaxis right
            filledArea = fill([bin_centers fliplr(bin_centers) bin_centers(1)], ...
                ...
                [data.(L).([feature_list{featureNum} '_upperPrctile_97_5']) ...
                fliplr(data.(L).([feature_list{featureNum} '_lowerPrctile_2_5'])) ...
                data.(L).([feature_list{featureNum} '_upperPrctile_97_5'])(1)], ...
                [0 0 0], 'FaceALpha', 0.15, 'EdgeColor', 'none');
            hold on;
            p1 = plot(bin_centers, data.(L).(feature_bin_list{featureNum}),'-k','LineWidth',2);
            p2 = plot(bin_centers, data.(L).(smoothed_feature_list{featureNum}), '-', 'Color', condition_colors{c}(1:3), 'LineWidth',2);
            currYLims = get(gca, 'YLim');
            currYLims(1) = currYLims(1) - diff(currYLims);
            set(gca, 'YLim', currYLims)
            title([feature_list{featureNum} ': MI = ' num2str(data.(L).([feature_list{featureNum} '_MI'])(1)) '; p = ' num2str(data.(L).([feature_list{featureNum} '_MI'])(2))])
            xlim([0 2*pi])
            xlabel('Heart-cycle Phase (0-2pi)')
            ylabel([feature_list{featureNum} ', microvolts'])
            if c == 1
                %             legend([filledArea p1 p2], {'95% Confidence Interval', 'Average by Bin', 'Rlowess-Smoothed'}, 'Location', 'southoutside')
            end
        end
    end
    
    filename= ['PSTH_overlaid_Feature_Dynamics__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% cosine fitted plots
    f4 = figure;
    set(f4, 'Position', [1 41 1920 963])
    %sgtitle(sgtitleText, 'interpreter', 'none')
    
    for c=1:numel(cfg.condition)
        col=condition_colors{c};
        L=cfg.condition(c).name;
        subplot(3,5,c)
        plot(bin_centers, data.(L).AMP_microV_byBin'/nanmean(data.(L).AMP_microV_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(bin_centers, ...
            data.(L).AMP_MI(5)+data.(L).AMP_MI(1)*cos(bin_centers-data.(L).AMP_MI(3)), 'Color', col,'LineWidth',2);
        plot(bin_centers, ...
            data.(L).AMP_microV_byBin_smoothed / nanmean(data.(L).AMP_microV_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).AMP_MI(1)) ...
            ', p = ' num2str(data.(L).AMP_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel('AMP Signal Change, signal / average(smoothed signal))')
        
        subplot(3,5,3+c)
        plot(bin_centers, data.(L).HW_ms_byBin'/nanmean(data.(L).HW_ms_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(bin_centers, ...
            data.(L).HW_MI(5)+data.(L).HW_MI(1)*cos(bin_centers-data.(L).HW_MI(3)), ...
            'Color', col,'LineWidth',2);
        plot(bin_centers, ...
            data.(L).HW_ms_byBin_smoothed / nanmean(data.(L).HW_ms_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).HW_MI(1)) ...
            ', p = ' num2str(data.(L).HW_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel('HW Signal Change, signal / average(smoothed signal))')
        
        subplot(3,5,10+c)
        plot(bin_centers, data.(L).TPW_ms_byBin'/nanmean(data.(L).TPW_ms_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(bin_centers, ...
            data.(L).TPW_MI(5)+data.(L).TPW_MI(1)*cos(bin_centers-data.(L).TPW_MI(3)), ...
            'Color', col,'LineWidth',2);
        plot(bin_centers, ...
            data.(L).TPW_ms_byBin_smoothed / nanmean(data.(L).TPW_ms_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).TPW_MI(1)) ...
            ', p = ' num2str(data.(L).TPW_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel('TPW Signal Change, signal / average(smoothed signal))')
        
        subplot(3,5,13+c)
        plot(bin_centers, data.(L).REP_ms_byBin'/nanmean(data.(L).REP_ms_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(bin_centers, ...
            data.(L).REP_MI(5)+data.(L).REP_MI(1)*cos(bin_centers-data.(L).REP_MI(3)), ...
            'Color', col,'LineWidth',2);
        plot(bin_centers, ...
            data.(L).REP_ms_byBin_smoothed / nanmean(data.(L).REP_ms_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).REP_MI(1)) ...
            ', p = ' num2str(data.(L).REP_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel('REP Signal Change, signal / average(smoothed signal))')
    end
    
    filename= ['Cosine_Fitted__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% results of correlation analysis
    f5 = figure;
    set(f5, 'Position', [381 424 1450 398])
    %sgtitle(sgtitleText, 'interpreter', 'none')
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        subplot(1,4,c)
        scatter(data.(L).FRbyRR_Hz, data.(L).cycleDurations_s, [], condition_colors{c}(1:3), 'Marker', '.')
        hold on
        linTrend = lsline(gca);
        linTrend.Color = condition_colors{c}(1:3);
        xlabel('Firing Rate per Heart Cycle, Hz')
        ylabel('Heart Cycle Duration, s')
        box on
        title([L ': cc = ' num2str(data.(L).pearson_r(4)) '; p = ' num2str(data.(L).pearson_p(4))])
        legend({'Real Data', 'Least-Square Fit'}, 'Location', 'Best')
    end
    
    subplot(1,4,4)
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
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    close all
    
end
