clear all, close all

% This script creates the figure on correlation analysis summary for the
% manuscript Vasileva, Kaduk et al., 2024

unqTargets    = {'VPL', 'dPul', 'MD'};
unqConditions = {'Rest', 'Task'};

bin_resolution = 0.05;
bar_colors = [244 149 173; ... % pink  - pos
            126 221 95; ...    % green - neg
            127 127 127]/255;  % grey  - non-significant

%% load data
figure,
set(gcf, 'Position', [80 45 1683 949])

%% sequences of FR and RRs in the rest
load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Bacchus_TaskRest\correlation_analysis\Bac_20220106_02_VPL_R_correlation.mat')

x_lim = [4878 4950];

% RR durations in the example data
subplot(7,3,1:2)
plot(data.Rest.timeRRstart, data.Rest.FRbyRR_Hz, 'b.-')
ylim([0 40])
xlim(x_lim)
set(gca,'XTickLabel',[])
ylabel('Firing Rate, Hz')

% FR per cardiac cycle in the example data
subplot(7,3,4:5)
plot(data.Rest.timeRRstart, 1000*data.Rest.cycleDurations_s, 'ko-')
ylim(1000*[0.35 0.58])
xlim(x_lim)
xlabel('Time from the recording onset, s')
ylabel('RR Duration, ms')

%% linear regression between FR and RR
fr_bins = 0:5:100;

% FR vs. RR + linear regression
subplot(3,3,3)
hold on
box on
scatter(data.Rest.FRbyRR_Hz, 1000*data.Rest.cycleDurations_s, 10, 'b', 'filled', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2)
scatter(data.Task.FRbyRR_Hz, 1000*data.Task.cycleDurations_s, 10, 'r', 'filled', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2)
rmnan_rest = isnan(data.Rest.FRbyRR_Hz) | isnan(data.Rest.cycleDurations_s);
rmnan_task = isnan(data.Task.FRbyRR_Hz) | isnan(data.Task.cycleDurations_s);

p_rest = polyfit(data.Rest.FRbyRR_Hz(~rmnan_rest),1000*data.Rest.cycleDurations_s(~rmnan_rest),1);
p_task = polyfit(data.Task.FRbyRR_Hz(~rmnan_task),1000*data.Task.cycleDurations_s(~rmnan_task),1);
line(fr_bins,p_rest(1)*fr_bins + p_rest(2),'Color','blue')
line(fr_bins,p_task(1)*fr_bins + p_task(2),'Color','red')

xlim([0 40])
xlabel('FR per Cardiac Cycle, Hz')
ylabel('RR Durations, ms')
ylim([300 650])
set(gca,'YTick',300:50:650)
legend({['Rest: Pearson''s r = ' num2str(round(data.Rest.pearson_r(ceil(end/2)),3)) '; p < 10^-^4'], ...
    ['Task: Pearson''s r = ' num2str(round(data.Task.pearson_r(ceil(end/2)),3)) '; p = ' num2str(round(data.Task.permuted_p(ceil(end/2)),3))]},'Location','best')
% legend({['Rest: r = ' num2str(round(data.Rest.pearson_r(ceil(end/2)),3)) '; p = ' num2str(round(data.Rest.permuted_p(ceil(end/2)),3))], ...
%     ['Task: r = ' num2str(round(data.Task.pearson_r(ceil(end/2)),3)) '; p = ' num2str(round(data.Task.permuted_p(ceil(end/2)),3))]},'Location','best')

%% histograms of unit FR-RR correlation signs
% loop through monkeys
unqMonkey = {'Bacchus', 'Magnus'};
monkeyInitials = {'B', 'M'};
for curr_monk = 1:2
    load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' unqMonkey{curr_monk} '_TaskRest\Population_correlation_analysis_after_SNR_exclusion_stable_noLow_amplitude_ccs_any\data.mat'])
    for a = 1:3
        
        T = unqTargets{a};
        
        for c=1:2
            L = unqConditions{c};
            
            curr_cc = dt.(T).(L).pearson_r(ceil(end/2), :);
            curr_pp = dt.(T).(L).pearson_p(ceil(end/2), :);
            [~, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
            [sig_idx,~] = fdr_bky(curr_pp, 0.05);
            pos_idx = curr_cc > 0;
            neg_idx = curr_cc < 0;
            nan_idx = isnan(curr_cc);
            
            counts_pos    = histc(curr_cc(sig_idx & pos_idx), -1:bin_resolution:1); % significant and positive
            counts_neg    = histc(curr_cc(sig_idx & neg_idx), -1:bin_resolution:1); % significant and negative
            counts_nonsig = histc(curr_cc(~sig_idx & ~nan_idx), -1:bin_resolution:1);
            
            counts_pos    = counts_pos(:);
            counts_neg    = counts_neg(:);
            counts_nonsig = counts_nonsig(:);
            
            unit_counts_by_area(c,a,:) = permute([sum(counts_pos) sum(counts_neg) sum(counts_nonsig)], [1 3 2]);
            
        end
    end
    
    for c = 1:2
        
        L = unqConditions{c};
        
        % proportion histograms
        subplot(3,4,4+c+(curr_monk-1)*2)
        tmp = permute(unit_counts_by_area(c,:,:), [2 3 1]);
        unit_percentages = 100 * tmp ./ sum(tmp,2);
        b = bar(unit_percentages,'stacked');
        for ii = 1:size(bar_colors,1)
            b(ii).FaceColor = bar_colors(ii,:);
        end
        % add percentages to bars
        xbarCnt = vertcat(b.XEndPoints);
        ybarTop = vertcat(b.YEndPoints);
        ybarCnt = ybarTop - unit_percentages'/2;
        
        % Create text strings
        txt = compose('%.1f%%',unit_percentages');
        th = text(xbarCnt(:), ybarCnt(:), txt(:), ...
            'HorizontalAlignment', 'center', ....
            'VerticalAlignment', 'middle', ...
            'Color', 'k',....
            'FontSize', 8);
        
        title(['Monkey ' monkeyInitials{curr_monk} ': ' L])
        % create x tick labels
        for ii = 1:3
            xtick_labels{ii} = [unqTargets{ii} ' (' num2str(length(dt.(unqTargets{ii}).unitId)) ')'];
        end
        set(gca,'XTickLabel',xtick_labels)
        if curr_monk == 1 & c == 1
            ylabel('Unit fractions, %')
        end
        
    end
    
end

%% histogram of ccs in example case at lag 0
unqConditions = {'Rest','Task'};
unqMonkey = {'Bacchus', 'Magnus'};
monkeyInitials = {'B', 'M'};

for monkNum = 1:2

    load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' unqMonkey{monkNum} '_TaskRest\Population_correlation_analysis_after_SNR_exclusion_stable_noLow_amplitude_ccs_any\data.mat'])
    
    for c = 1:2
%         subplot(3,4,9)
        s(c,monkNum) = subplot(3,4,8+c+(monkNum-1)*2);
%         s(c,monkNum) = subplot(4,4,12+c+(monkNum-1)*2);
        curr_cc = dt.VPL.(unqConditions{c}).pearson_r(ceil(end/2), :);
        curr_pp = dt.VPL.(unqConditions{c}).pearson_p(ceil(end/2), :);
        [sig_idx,~] = fdr_bky(curr_pp, 0.05);
        pos_idx = curr_cc > 0;
        neg_idx = curr_cc < 0;
        nan_idx = isnan(curr_cc);
        
        mean_cc = mean(curr_cc);
        
        % prepare data for plotting
        counts_pos    = histc(curr_cc(sig_idx & pos_idx),   -1+bin_resolution/2:bin_resolution:1-bin_resolution/2) / length(curr_cc); % significant and positive
        counts_neg    = histc(curr_cc(sig_idx & neg_idx),   -1+bin_resolution/2:bin_resolution:1-bin_resolution/2) / length(curr_cc); % significant and negative
        counts_nonsig = histc(curr_cc(~sig_idx & ~nan_idx), -1+bin_resolution/2:bin_resolution:1-bin_resolution/2) / length(curr_cc);
        
        counts_pos    = counts_pos(:);
        counts_neg    = counts_neg(:);
        counts_nonsig = counts_nonsig(:);
        
        plot_data = [counts_pos counts_neg counts_nonsig];
        
        hold on
        b = bar([-1+bin_resolution/2:bin_resolution:1-bin_resolution/2]+bin_resolution/2, plot_data, 'stacked');
        plot(mean_cc, 0, '^', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
        set(b, 'FaceColor', 'Flat')
        for ii = 1:size(bar_colors,1)
            b(ii).FaceColor = bar_colors(ii,:);
        end
        xlim([-0.4 0.4])
        set(gca,'XTick',-0.4:0.2:0.4,'XTickLabel',-0.4:0.2:0.4, 'TickDir', 'both')
        box on
        title(['Monkey ' monkeyInitials{monkNum} ' - VPL - ' unqConditions{c}])
        if monkNum == 1 && c == 1
            legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.', 'Mean CC'})
            xlabel('CC between FR and RR duration')
            ylabel('Proportion of Units')
        end
    end
    linkaxes(s)
end
