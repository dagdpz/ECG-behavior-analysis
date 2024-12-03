clear all, close all

bar_colors = [244 149 173; ... % pink  - pos
    126 221 95; ...    % green - neg
    127 127 127]/255;  % grey  - non-significant
deep_green = [6,64,43]/255;

pos_neg_lag_colors = [170 60 0;         % pos
    255 153 85;         % neg
    255 255 255]/255;   % zero

unqTargets    = {'VPL', 'dPul', 'MD'};
unqConditions = {'Rest', 'Task'};
unqMonkey     = {'Bacchus', 'Magnus'};
monkeyInitials = {'B', 'M'};

load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Bacchus_TaskRest\correlation_analysis\Bac_20210716_11_VPL_R_correlation.mat')

%% histograms with pos and neg lags

figure;
set(gcf,'Position',[141 480 1212 483])

for monkNum = 1:2
    
    load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' unqMonkey{monkNum} '_TaskRest\Population_correlation_analysis_after_SNR_exclusion_stable_noLow_amplitude_ccs_any\data.mat'])
    
    laglist = data.cc_lag_list;
    
    for a = 1:3
        T=unqTargets{a};
        
        for c=1:2
            
            L=unqConditions{c};
            
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
            
            % compute unit counts and fractions with lag <> 0
            lag_zero_counts(a,c) = sum(laglist(ids) == 0 | ~sig_idx);
            pos_lag_counts(a,c)  = sum(laglist(ids) > 0 & sig_idx);
            neg_lag_counts(a,c)  = sum(laglist(ids) < 0 & sig_idx);
            
            lag_zero_prc(a,c) = lag_zero_counts(a,c) / (lag_zero_counts(a,c) + pos_lag_counts(a,c) + neg_lag_counts(a,c));
            pos_lag_prc(a,c)  = pos_lag_counts(a,c) / (lag_zero_counts(a,c) + pos_lag_counts(a,c) + neg_lag_counts(a,c));
            neg_lag_prc(a,c)  = neg_lag_counts(a,c) / (lag_zero_counts(a,c) + pos_lag_counts(a,c) + neg_lag_counts(a,c));
            
        end
        
    end
    
    for c = 1:2
        
        
        subplot(3, 4, c + (monkNum - 1)*2)
        hold on
        box on
        b = bar(100 * [pos_lag_prc(:,c) neg_lag_prc(:,c) lag_zero_prc(:,c)], 'stacked');
        for ii = 1:size(bar_colors,1)
            b(ii).FaceColor = pos_neg_lag_colors(ii,:);
        end
        clear b
        set(gca, 'XTick', [1 2 3], 'XTickLabel', unqTargets)
%         legend({'Lag > 0', 'Lag < 0', 'Lag = 0 or non-sig.'},'Location','southoutside')
        title(['Monkey ' monkeyInitials{monkNum} ': ' unqConditions{c}])
        ylim([0 100])
        
        if c == 1 & monkNum == 1
            ylabel('Proportion of Units, %')
        end
        
    end
    
end

%% plot for lag of most effective correlation

for monkNum = 1:2

load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' unqMonkey{monkNum} '_TaskRest\Population_correlation_analysis_after_SNR_exclusion_stable_noLow_amplitude_ccs_any\data.mat'])

for c=1:2
            
    L=unqConditions{c};
    
    laglist=dt.dPul.cc_lag_list(1,:);
    R=dt.VPL.(L).pearson_r;
    P=dt.VPL.(L).pearson_p;
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
    
    % comptute histogram counts
    h0=hist(laglist(ids(sig==0)),laglist);
    h1=hist(laglist(ids(sig==-1)),laglist);
    h2=hist(laglist(ids(sig==1)),laglist);
    
    % compute medians by subgroup
    m0=median(laglist(ids(sig==0)));
    m1=median(laglist(ids(sig==-1)));
    m2=median(laglist(ids(sig==1)));
    % compute lag medians only sig.
    M = median(laglist(ids(sig_idx)));
    
    s(c) = subplot(2,4,4+c+(monkNum-1)*2);
    % x         = -20:0.125:20;
    % h1_h2_smo = smooth(interp1(laglist,h1+h2,x,'pchip'),50);
    % h1_smo    = smooth(interp1(laglist,h1,x,'pchip'),50);
    % h2_smo    = smooth(interp1(laglist,h2,x,'pchip'),50);
    x         = -20:20;
    h1_h2_smo = h1+h2;
    h1_smo    = h1;
    h2_smo    = h2;
    hold on
    
    plot(x,h1_h2_smo,'Color',[0.5 0.5 0.5],'LineWidth',0.25)
    % fill([x fliplr(x) x(1)], [h1_h2_smo zeros(size(h1_h2_smo)) h1_h2_smo(1)], deep_green, 'EdgeColor','none')
    fill([x fliplr(x) x(1)], [h1_smo zeros(size(h1_smo)) h1_smo(1)], bar_colors(2,:), 'EdgeColor','none','FaceAlpha',0.7)
    fill([x fliplr(x) x(1)], [h2_smo zeros(size(h2_smo)) h2_smo(1)], bar_colors(1,:), 'EdgeColor','none','FaceAlpha',0.7)
    % plot(M(a,c),  0, '^', 'MarkerFaceColor', deep_green, 'MarkerEdgeColor', deep_green)         % median overall
    plot(m1, 0, '^', 'MarkerFaceColor', bar_colors(2,:), 'MarkerEdgeColor', bar_colors(2,:)/2) % median sig pos
    plot(m2, 0, '^', 'MarkerFaceColor', bar_colors(1,:), 'MarkerEdgeColor', bar_colors(1,:)/2) % median sig neg
    box on
    
    if monkNum == 1 & c == 1
        xlabel('Lag for Most Effective Correlation, # cardiac cycles')
        ylabel('N units')
    end
    xlim([-20 20])
    % ylim([-0.5 0.5])
    box on
    title(['Monkey ' monkeyInitials{monkNum} ' - ' L ' - VPL'])
%     legend({'Sig. Neg. + Sig. Pos.', 'Sig. Neg.', 'Sig. Pos.', 'Median Neg.', 'Median Pos.'}, 'Location', 'southoutside')

end

linkaxes(s)

end

%% example unit with ccs for all the lags
% [~,idx_max_rest] = max(abs(data.Rest.pearson_r));
% [~,idx_max_task] = max(abs(data.Task.pearson_r));
% 
% subplot(3,4,5)
% hold on
% plot(data.cc_lag_list, data.Rest.pearson_r, 'bo-')
% plot(data.cc_lag_list, data.Task.pearson_r, 'rs-')
% plot(data.cc_lag_list(idx_max_rest), data.Rest.pearson_r(idx_max_rest), 'Marker', 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
% plot(data.cc_lag_list(idx_max_task), data.Task.pearson_r(idx_max_task), 'Marker', 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
% 
% box on
% xlabel('Lag, # Cardiac Cycles')
% ylabel('Correlation Coefficient')
% legend({'Rest', 'Task'}, 'Location', 'best')
