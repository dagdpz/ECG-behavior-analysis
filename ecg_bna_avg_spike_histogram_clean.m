function ecg_bna_avg_spike_histogram_clean(SPK_PSTH, output_folder, cfg)
% Here comes some sort of across population plot i assume?

n_sig_bins = 4;

bar_colors = [0.8500 0.3250 0.0980; 0 0.4470 0.7410; 1 1 1];
area_colors = [];

savePlot = 1;
saveTable= 0;
OnlyUnits_withRestANDTask = 1;
Graph_SelectionCriterion = 0;
colors = distinguishable_colors(25);
PSTH_bins = (cfg.analyse_states{1,3}:cfg.spk.PSTH_binwidth:cfg.analyse_states{1,4})*1000 ;

if ~exist(output_folder,'dir')
    mkdir(output_folder);
end

ConNames = {cfg.condition(:).name};

lineProps = cell(length(cfg.condition), 1);
for conNum = 1:length(cfg.condition)
    lineProps{conNum} = {'color',cfg.condition(conNum).color,'linewidth',4};
end

TargetBrainArea = {SPK_PSTH.target};
Ana_TargetBrainArea = cfg.targets;
if cfg.combine_hemispheres
    TargetBrainArea=cellfun(@(x) x(1:end-2),TargetBrainArea,'UniformOutput',false);
    Ana_TargetBrainArea=cellfun(@(x) x(1:end-2),Ana_TargetBrainArea,'UniformOutput',false);
end
Ana_TargetBrainArea=unique(Ana_TargetBrainArea);
Ana_TargetBrainArea=Ana_TargetBrainArea(ismember(Ana_TargetBrainArea,TargetBrainArea));
Ana_TargetBrainArea = {'VPL', 'dPul', 'MD'};

N_Areas=numel(Ana_TargetBrainArea);
N_conditions=numel(cfg.condition);


% Color_BrainArea = [[0 0 0];  colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]
Color_BrainArea = [colors(7,:);     colors(13,:); colors(21,:)  ];
Color_BrainArea = distinguishable_colors(6);
Color_BrainArea = [1.0000         0         0;
    1.0000    0.5300         0;
    0.5500    0.9000    0.1000;
    0.0500    0.6500    0.3000;
    0.0500    0.6500    0.7000;
    0.5500    0.2000    0.7500];
Color_BrainArea = [1 0.53 0;
                0.05 0.65 0.7;
                1 0 0.6];

%%  Calculations
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for c=1:N_conditions
        L=cfg.condition(c).name;
        dat=[SPK_PSTH(ismember(TargetBrainArea,T)).(L)];
        dat_fieldnames=fieldnames(rmfield(dat,{'raster','Rds','Rds_perm'}));% ,'Rts'
        for fn=1:numel(dat_fieldnames)
            N=dat_fieldnames{fn};
            Out.(T).(L).(N)=vertcat(dat.(N));
        end
        dat_fieldnames={'unit_ID','quantSNR','Single_rating','stability_rating'};
        for fn=1:numel(dat_fieldnames)
            N=dat_fieldnames{fn};
            if strcmp(N, 'unit_ID')
                Out.(T).(L).(N)=vertcat({SPK_PSTH(ismember(TargetBrainArea,T)).(N)}');
            else
                Out.(T).(L).(N)=vertcat(SPK_PSTH(ismember(TargetBrainArea,T)).(N));
            end
        end
        
        Out.(T).(L).FR_perECGTriggeredAverage = nanmean(vertcat(dat.SD),2);
        
        %% compute histograms - unit fractions and counts
        out = [Out.(T).(L)];
        
        Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
        Idx_Units_NaN =  sum(~Idx_Units_NonNaN);
        sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
        
        % increase, decrease, non-sign
        Out.(T).(L).Pc_SignFR = ([sum(out.sig_sign(sig) == 1) ,sum(out.sig_sign(sig) == -1),(sum(~sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
        Out.(T).(L).Nb_SignFR = ([sum(out.sig_sign(sig) == 1) ,sum(out.sig_sign(sig) == -1),(sum(~sig) -Idx_Units_NaN),] ) ;
        
        %% compute modulation indices with sign
        curr_sign = Out.(T).(L).sig_sign;
        curr_sign(curr_sign == 0) = 1; % to define real sign but not to lose non-significant ones
        Out.(T).(L).FR_ModIndex_SubtrSDP_signed = Out.(T).(L).FR_ModIndex_SubtrSDP .* curr_sign;
        Out.(T).(L).FR_ModIndex_PcS_signed      = Out.(T).(L).FR_ModIndex_PcS .* curr_sign;
        
    end
end

%% percentages of units with different responsivenes - pos. / neg. / non-sig.
figure('Name',sprintf('BarPlot_Pc'),'Position',[846   552   600   239],'PaperPositionMode', 'auto');
for c=1:N_conditions
    L=cfg.condition(c).name;
    % prepare data matrix for plotting
    data_mat = [];
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        data_mat = vertcat(data_mat, Out.(T).(L).Pc_SignFR);
    end
    subplot(1,N_conditions,c);
    b = bar(data_mat,'stacked');
    for ii = 1:3
        b(ii).FaceColor = bar_colors(ii,:);
    end
    title(L,'interpreter','none');
    set(gca,'XTickLabel',Ana_TargetBrainArea,'fontsize',10);
    ylim([0 100])
    ylabel('Percentage of Units, %')
    % legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')
    ax = gca;
    ax.FontSize = 12;
end
save_figure_as('Pc_CardiacRelatedUnits',output_folder,savePlot)

%% Fisher's exact test to compare unit fractions between task and rest
for groupNum = 1:length(cfg.spk.compare_conditions)
    
    cond1_num = cfg.spk.compare_conditions{groupNum}(1);
    cond2_num = cfg.spk.compare_conditions{groupNum}(2);
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        
        disp(T)
        
        p(a,:) = ecg_bna_fisher_test(Out.(T).(cfg.condition(cond2_num).name).Nb_SignFR, Out.(T).(cfg.condition(cond1_num).name).Nb_SignFR);
        
    end
    
    % find corrected p-values
    [p_corr, h] = bonf_holm(p, 0.05);
    Tab = table(Ana_TargetBrainArea', p_corr(:,1), p_corr(:,2), h(:,1), h(:,2), 'VariableNames', {'Brain Area', 'p_corr inc', 'p_corr dec', 'h inc', 'h dec'});
    savename = [output_folder filesep (cfg.condition(cond2_num).name) '_vs_' (cfg.condition(cond1_num).name) 'Table_Prevalences_pvalues_corrected.xlsx'];
    writetable(Tab, savename)
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
            scatter(out.FR, nanmean(out.SDP,2), 'filled', 'MarkerFaceColor',ccol)  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
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
            scatter(nanmean(out.SD,2) , nanmean(out.SDP,2), 'filled', 'MarkerFaceColor',ccol)  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
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
        if ~isempty(out.quantSNR)
            scatter(out.quantSNR(sig)  , out.FR_ModIndex_PcS(sig),  'filled', 'MarkerFaceColor', ccol);
            scatter(out.quantSNR(~sig) , out.FR_ModIndex_PcS(~sig), 'filled', 'MarkerFaceColor', ccol/2);
            axis square;
            
            %xf = [min(out.quantSNR), max(out.quantSNR)];
            [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_PcS)),out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_PcS)),1); %
            [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_PcS)),S);
            [coef, pval]  = corr(out.quantSNR,out.FR_ModIndex_PcS, 'rows','complete') ;
            plot(out.quantSNR(~isnan(out.FR_ModIndex_PcS)), y_fit,'LineWidth', 2, 'Color', ccol);
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
        text(-400,-15, out.unit_ID(UnitSig_Rest),'Color', ccol);
        if sum(UnitSig_Rest)
            line(PSTH_bins, out.SDsubstractedSDP_normalized(UnitSig_Rest,:), 'color', ccol, 'LineWidth', 1);
        end
        text(300,20, out.unit_ID(UnitNotSign_Rest),'Color','k');
        xlabel('Time from R-peak, ms')
        ylabel('% signal change')
        title('Units with Modulation strength > 30%');
        
        subplot(2,4,7:8);
        box on; hold on;
        text(300,20, out.unit_ID(UnitNotSign_Rest),'Color','k');
        if sum(UnitNotSign_Rest)
            line(PSTH_bins, out.SDsubstractedSDP_normalized(UnitNotSign_Rest,:), 'color', ccol,'LineWidth', 4);
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
            SDmean_SEM = nanstd(out.SDsubstractedSDP(sig, :))/ sqrt(length(nanmean(out.SDsubstractedSDP(sig, :)))) ;
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP(sig, :)) ,SDmean_SEM, lineProps{c}, 1);
            ylabel('Signal Change, spikes/s')
            title(['Pop:  (all significant) Cal Per Unit:SD-SDP ' T ' units'],'interpreter','none');
            legend(ConNames)
            
            subplot(2,2,3:4);
            hold on
            box on
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(sig, :),[], 1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(sig, :), 1))) ;
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(sig, :), 1) ,SDmean_SEM ,lineProps{c},1);
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
                Y1 = out.SD( idx1 & idx2,:)  ;
                Y2 = out.SDP( idx1 & idx2,:)  ;
                
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
            scatter(out.sig_n_bins , out.quantSNR, 'filled', 'MarkerFaceColor', ccol); hold on;
            scatter(out.sig_n_bins(idx_ex) , out.quantSNR(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
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
    for c=1:N_conditions
        L=cfg.condition(c).name;
        ccol=cfg.condition(c).color;
        
        out = [Out.(T).(L)];
        % from all not NAN units - how many were significant?
        Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
        Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
        sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
        dec = (out.sig_sign == -1);
        inc = (out.sig_sign == 1);
        
        % bar plot how many are significant & positive and negative FR?
        subplot(2,4,1:2);
        hold on;
        if sum(sig)
            line(PSTH_bins, out.SDsubstractedSDP_normalized(sig,:), 'color', ccol, 'LineWidth', 1);
        end
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(sig, :), [], 1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(sig, :), 1))) ;
        shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(sig, :), 1) ,SDmean_SEM,lineProps{c},1);
        title({[L ': units = ' ,num2str(sum(sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN)) ], ['Population:  (all significant)' (T) ' units']},'interpreter','none');
        ylabel('normalized Firing rate (%)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );
        
        %Bar graph how many
        subplot(2,4,3);
        toplot=arrayfun(@(x) Out.(T).(x.name).Pc_SignFR, cfg.condition, 'uniformoutput', false);
        barpairs =  vertcat(toplot{:});
        b = bar(barpairs,'stacked', 'Facecolor','flat' );
        b(3).FaceColor = [1 1 1];
        legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')
        set(gca,'XTickLabel',ConNames,'fontsize',10);
        
        % display only significant units showing a increase in FR
        subplot(2,4,5); %
        hold on;
        text(PSTH_bins(2),c*2, [L ': units = ' ,num2str(sum(inc & sig)) ],'Color',ccol)
        
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(inc & sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(inc & sig,:), 1))) ;
        shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(inc & sig,:),1), SDmean_SEM ,lineProps{c},1);
        set(gca,'ylim',[-10, 10]);
        title('units showing a sig. INCREASE in FR','interpreter','none');
        ylabel('normalized Firing rate (%)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );
        
        % display only significant units showing a decrease in FR
        subplot(2,4,6);
        text(PSTH_bins(2),c*2, [L ': units = ' ,num2str(sum(dec & sig)) ],'Color',ccol)
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(dec & sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(dec & sig,:), 1))) ;
        shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(dec & sig,:),1), SDmean_SEM ,lineProps{c},1);
        set(gca,'ylim',[-10, 10]);
        title('units showing a sig. DECREASE in FR','interpreter','none');
        ylabel('normalized Firing rate (spike/s)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );
        
        
        %% MODULATION INDEX
        subplot(2,4,4);
        hold on;
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );
        title('raw');
        axis square;
        make_correlation_plot(out,'quantSNR','FR_ModIndex_PcS',sig,0,ccol,30+5*c) %% or FR_ModIndex_SubtrSDP ??
        
        subplot(2,4,7);
        hold on;
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );
        title('-Surrogate');
        axis square;
        make_correlation_plot(out,'quantSNR','FR_ModIndex_SubtrSDP',sig,0,ccol,30+5*c)
        
        subplot(2,4,8);
        hold on
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('mean firing rate','fontsize',14);
        title('-Surrogate, sig correlated');
        axis square;
        make_correlation_plot(out,'FR','FR_ModIndex_SubtrSDP',sig,1,ccol,30+5*c)
        
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
        xlabel('Signal-to-Noise','fontsize',14 );
        title([T 'all units'])
        make_correlation_plot(out,'quantSNR','FR_ModIndex_PcS',sig,0,ccol,30+5*c)
        
        subplot(2,N_Areas,a +N_Areas);
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );
        title([T ' signficant only'])
        make_correlation_plot(out,'quantSNR','FR_ModIndex_PcS',sig,1,ccol,30+5*c)
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
        title([T 'all units'])
        make_correlation_plot(out,'quantSNR','FR_ModIndex_SubtrSDP',sig,0,ccol,10+2*c) %% was FR_ModIndex_PcS...
        
        subplot(2,N_Areas,a +N_Areas);
        ylabel('Modulation index (spike/s)','fontsize',14 );
        xlabel('Signal-to-Noise ratio (mV)','fontsize',14 );
        title([T ' signficicant only'])
        make_correlation_plot(out,'quantSNR','FR_ModIndex_SubtrSDP',sig,1,ccol,10+2*c)
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
        title([T 'all units'])
        make_correlation_plot(out,'FR_ModIndex_PcS','FR_ModIndex_SubtrSDP',sig,0,ccol,10+2*c) %% was FR_ModIndex_PcS...
        
        subplot(2,N_Areas,a +N_Areas);
        ylabel('Modulation index (spike/s)','fontsize',14 );
        xlabel('Modulation index(%pSc)','fontsize',14 );
        title([T ' signficicant only'])
        make_correlation_plot(out,'FR_ModIndex_PcS','FR_ModIndex_SubtrSDP',sig,1,ccol,10+2*c)
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
            scatter(out.quantSNR(sig),out.FR_ModIndex_PcS(sig) , 30, out.FR_perECGTriggeredAverage(sig)/max(out.FR_perECGTriggeredAverage(sig)) , 'filled')
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
        title([T ' all units'])
        make_correlation_plot(out,'FR','FR_ModIndex_PcS',sig,0,ccol,30+5*c)
        
        subplot(2,N_Areas,a +N_Areas);
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );
        title([T ' only significant'])
        make_correlation_plot(out,'FR','FR_ModIndex_PcS',sig,1,ccol,30+5*c)
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
        title([T ' all units'])
        make_correlation_plot(out,'quantSNR','FR_perECGTriggeredAverage',sig,0,ccol,15+2*c)
        
        subplot(2,N_Areas,a +N_Areas);
        xlabel('Signal-to-Noise ratio','fontsize',14 );
        ylabel('average Firing rate','fontsize',14 );
        title([T ' only significant'])
        make_correlation_plot(out,'quantSNR','FR_perECGTriggeredAverage',sig,1,ccol,15+2*c)
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
        title([T ' all units'])
        make_correlation_plot(out,'FR','FR_ModIndex_SubtrSDP',sig,0,ccol,10+2*c)
        
        subplot(2,N_Areas,a +N_Areas);
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );
        title([T ' only significant'])
        make_correlation_plot(out,'FR','FR_ModIndex_SubtrSDP',sig,1,ccol,10+2*c)
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
        title([T ' all units'])
        make_correlation_plot(out,'sig_n_bins','FR_ModIndex_PcS',sig,0,ccol,35+5*c)
        
        subplot(2,N_Areas,a +N_Areas);
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Bin size of sig. Interval','fontsize',14 );
        title([T ' only significant'])
        make_correlation_plot(out,'sig_n_bins','FR_ModIndex_PcS',sig,1,ccol,35+5*c)
    end
end
save_figure_as('ModulationIndex_Nspc_Binsize',output_folder,savePlot)

%% Decrease and Increase grouped for brain region
figure('Name',sprintf('CardiacRelated_Change_FR'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
for c=1:N_conditions
    L=cfg.condition(c).name;
    for a = 1: N_Areas
        T=Ana_TargetBrainArea{a};
        acol=Color_BrainArea(a,:);
        out = [Out.(T).(L)];
        
        sig      =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > n_sig_bins) ;
        dec      = (out.sig_sign == -1) & sig;
        inc      = (out.sig_sign == 1) & sig;
        %         idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
        %         idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
        %         idx_SigTime_After50       = (out.sig_time > 50);
        
        subplot(N_conditions,4,1+(c-1)*4);
        hold on; axis square; box on;
        lineProps={'color',acol,'linewidth',3};
        text(-400,-1* a, [(T), ' ' L ': units = ' ,num2str(sum(inc)) ],'Color',acol)
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(inc,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(inc,:),1))) ;
        shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(inc,:),1), SDmean_SEM ,lineProps,1);
        %set(gca,'ylim',[-10, 10]);
        title([L ' Increase' ],'interpreter','none')
        ylabel('normalized Firing rate (%)','fontsize',12 );
        xlabel('Time relative to R-peak (ms)','fontsize',12 );
        % set(gca, 'XTick', PSTH_bins)
        
        % display only significant units showing a decrease in FR
        subplot(N_conditions,4,3+(c-1)*4);
        hold on; axis square; box on;
        lineProps={'color',acol,'linewidth',3};
        text(-400,1* a, [T ' ' L ': units = ' ,num2str(sum(dec)) ],'Color',acol)
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(dec,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(dec,:),1))) ;
        shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(dec,:),1), SDmean_SEM ,lineProps,1);
        %set(gca,'ylim',[-10, 10]);
        title([L ' Decrease' ],'interpreter','none')
        ylabel('normalized Firing rate (spike/s)','fontsize',12 );
        xlabel('Time relative to R-peak (ms)','fontsize',12 );
        % set(gca, 'XTick', PSTH_bins)
        
        subplot(N_conditions,4,2+(c-1)*4);
        hold on; axis square; box on;
        scatter(out.sig_time(inc) , out.FR_ModIndex_SubtrSDP(inc), 'filled', 'MarkerFaceColor',acol);
        ylabel('Modulation strength (% pSc)','fontsize',12);
        xlabel('TinePoint of sig. highest diff in FR','fontsize',12 );
        title([L ' Increase' ],'interpreter','none')
        
        subplot(N_conditions,4,4+(c-1)*4);
        hold on; axis square; box on;
        scatter(out.sig_time(dec) , out.FR_ModIndex_SubtrSDP(dec), 'filled', 'MarkerFaceColor',acol);
        ylabel('Modulation strength (% pSc)','fontsize',12);
        xlabel('TinePoint of sig. highest diff in FR','fontsize',12 );
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
                acol=Color_BrainArea(a,:);
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
                title([L ' INCREASE ' ThreeTiming{i_Time}],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',11 );
                xlabel('Time relative to R-peak (ms)','fontsize',11 );
                %   text(-400,-4* a, [(T), ' R: n = ' ,num2str(sum(inc & idx_sig & idx_Time)) ],'Color',acol)
                
                %             % not sure about this weird distinction here
                %             if size(out.SDsubstractedSDP_normalized(inc & tim,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized(inc & tim,:),1) == 0
                %                 SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                %                 shadedErrorBar(PSTH_bins,out.SDsubstractedSDP_normalized(inc & tim,:), SDmean_SEM ,lineProps,1);
                %                 MI_groups =  max(out.SDsubstractedSDP_normalized(inc & tim,WindowIdx)) - min(out.SDsubstractedSDP_normalized(inc & tim,:)) ;
                %
                %             else
                SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(inc & tim,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(inc & tim,:)))) ;
                shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(inc & tim,:),1), SDmean_SEM ,lineProps,1);
                MI_groups =  max(nanmean(out.SDsubstractedSDP_normalized(inc & tim,WindowIdx),1)) - min(nanmean(out.SDsubstractedSDP_normalized(inc & tim,:),1)) ;
                %             end
                % set(gca,'ylim',[-15, 15]);
                % vline(50); vline(-50); vline(250); vline(-250);
                % set(gca, 'XTick', PSTH_bins)
                
                
                %% Table stuff
                %             Table_Units = table(T,L,{'sigINCREASE'}, ThreeTiming(i_Time),  sum(inc & sig & tim), roundto(MI_groups,2) );
                %             Table_NrUnits = [Table_NrUnits; Table_Units];
                
                %% DECREASE
                subplot(length(cfg.spk.compare_conditions{groupNum}),6,2*(i_Time-1)+2 +(c-1)*6);
                hold on; axis square; box on;
                title([L ' DECREASE ', ThreeTiming{i_Time}],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',11 );
                xlabel('Time relative to R-peak (ms)','fontsize',11 );
                %   text(-400,4* a, [(T), 'R: n = ' ,num2str(sum(idx_SigDec & idx_sig & idx_Time)) ],'Color',acol)
                %             if size(out.SDsubstractedSDP_normalized(dec & tim,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized(dec & tim,:),1) == 0
                %                 SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                %                 shadedErrorBar(PSTH_bins,out.SDsubstractedSDP_normalized(dec & tim,:), SDmean_SEM ,lineProps,1);
                %                 MI_groups =   max(out.SDsubstractedSDP_normalized(dec & tim,:)) - min(out.SDsubstractedSDP_normalized(dec & tim,WindowIdx));
                %            else
                SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(dec & tim,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(dec & tim,:),1))) ;
                shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(dec & tim,:),1), SDmean_SEM ,lineProps,1);
                MI_groups =   max(nanmean(out.SDsubstractedSDP_normalized(dec & tim,:),1)) - min(nanmean(out.SDsubstractedSDP_normalized(dec & tim,WindowIdx),1));
                %end
                %set(gca,'ylim',[-15, 15]); %please do not do this
                %vline(50); vline(-50); vline(250); vline(-250);
                %set(gca, 'XTick', PSTH_bins)
                
                
                %% Table stuff
                %             Table_Units = table(T,L,{'sigDECREASE'}, ThreeTiming(i_Time),  sum(idx_SigDec & idx_sig & idx_Time), roundto(MI_groups,2)  );
                %             Table_NrUnits = [Table_NrUnits; Table_Units];
                
                subplot(length(cfg.spk.compare_conditions{groupNum}),6,5 +(c-1)*6);
                hold on; axis square; box on;
                scatter(out.sig_time(dec) , out.FR_ModIndex_SubtrSDP(dec), 'filled', 'MarkerFaceColor',acol);
                ylabel('Modulation strength (% pSc)','fontsize',11);
                xlabel('TinePoint of sig. highest diff in FR','fontsize',11 );
                title([L ' Decrease' ])
                
                subplot(length(cfg.spk.compare_conditions{groupNum}),6,6 +(c-1)*6);
                hold on; axis square; box on;
                scatter(out.sig_time(inc) , out.FR_ModIndex_SubtrSDP(inc), 'filled', 'MarkerFaceColor',acol);
                ylabel('Modulation strength (% pSc)','fontsize',11);
                xlabel('TinePoint of sig. highest diff in FR','fontsize',11 );
                title([L ' Increase' ])
                
                %% guessing this can be deleted now?
                %             lineProps={'color',acol,'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(inc & idx_SigTime_BeforeMinus50 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(inc & idx_SigTime_BeforeMinus50 & idx_sig,:)))) ;
                %             shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(inc & idx_SigTime_BeforeMinus50 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                %
                %             lineProps={'color',acol,'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(inc & idx_SigTime_Around0 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(inc & idx_SigTime_Around0 & idx_sig,:)))) ;
                %             shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(inc & idx_SigTime_Around0 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                %
                %             lineProps={'color',acol,'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(inc & idx_SigTime_After50 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(inc & idx_SigTime_After50 & idx_sig,:)))) ;
                %             shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(inc & idx_SigTime_After50 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                
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
for c=1:N_conditions
    L=cfg.condition(c).name;
    
    sig_vpl  = ~isnan(Out.VPL.(L).sig_FR_diff) & (Out.VPL.(L).sig_n_bins > n_sig_bins);
    sig_dpul = ~isnan(Out.dPul.(L).sig_FR_diff) & (Out.dPul.(L).sig_n_bins > n_sig_bins);
    sig_md   = ~isnan(Out.MD.(L).sig_FR_diff) & (Out.MD.(L).sig_n_bins > n_sig_bins);
    
    inc_vpl  = Out.VPL.(L).sig_time((Out.VPL.(L).sig_sign == 1) & sig_vpl);
    mm(1+6*(c-1)) = median(inc_vpl);
    ci(1+6*(c-1)) = iqr(inc_vpl);
    inc_dpul = Out.dPul.(L).sig_time((Out.dPul.(L).sig_sign == 1) & sig_dpul);
    mm(2+6*(c-1)) = median(inc_dpul);
    ci(2+6*(c-1)) = iqr(inc_dpul);
    inc_md   = Out.MD.(L).sig_time((Out.MD.(L).sig_sign == 1) & sig_md);
    mm(3+6*(c-1)) = median(sig_md);
    ci(3+6*(c-1)) = iqr(inc_md);
    
    dec_vpl  = Out.VPL.(L).sig_time((Out.VPL.(L).sig_sign == -1) & sig_vpl);
    mm(4+6*(c-1)) = median(dec_vpl);
    ci(4+6*(c-1)) = iqr(dec_vpl);
    dec_dpul = Out.dPul.(L).sig_time((Out.dPul.(L).sig_sign == 1) & sig_dpul);
    mm(5+6*(c-1)) = median(dec_dpul);
    ci(5+6*(c-1)) = iqr(dec_dpul);
    dec_md   = Out.MD.(L).sig_time((Out.MD.(L).sig_sign == 1) & sig_md);
    mm(6+6*(c-1)) = median(dec_md);
    ci(6+6*(c-1)) = iqr(dec_md);
    
    disp('Increase: VPL vs. dPul')
    [pp(1+6*(c-1)), hh(1+6*(c-1))] = ranksum(inc_vpl, inc_dpul)
    disp('Increase: VPL vs. MD')
    [pp(2+6*(c-1)), hh(2+6*(c-1))] = ranksum(inc_vpl, inc_md)
    disp('Increase: dPul vs. MD')
    [pp(3+6*(c-1)), hh(3+6*(c-1))] = ranksum(inc_dpul, inc_md)
    
    disp('Decrease: VPL vs. dPul')
    [pp(4+6*(c-1)), hh(4+6*(c-1))] = ranksum(dec_vpl, dec_dpul)
    disp('Decrease: VPL vs. MD')
    [pp(5+6*(c-1)), hh(5+6*(c-1))] = ranksum(dec_vpl, dec_md)
    disp('Decrease: dPul vs. MD')
    [pp(6+6*(c-1)), hh(6+6*(c-1))] = ranksum(dec_dpul, dec_md)
end

%% GENERATE A TABLE FOR ALL THE NUMBER OF UNITS PER CATEGORY
figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
ThreeTiming = {'inc T<-50 or dec T>50', '-50>T<50', 'inc T>50 or dec T<-50'} ;
for c=1:N_conditions
    L=cfg.condition(c).name;
    for i_Time = 1: length(ThreeTiming)
        for a = 1: N_Areas
            T=Ana_TargetBrainArea{a};
            acol=Color_BrainArea(a,:);
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
            text(-400,-4*a, [T, ' R: n = ' ,num2str(sum(tim)) ],'Color',acol)
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(tim,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(tim,:),1))) ;
            shadedErrorBar(PSTH_bins,nanmean(out.SDsubstractedSDP_normalized(tim,:),1), SDmean_SEM ,lineProps,1);
            
            %set(gca,'ylim',[-15, 15]);
            title([L ' ' ThreeTiming(i_Time)],'interpreter','none');
            ylabel('normalized Firing rate (% pSc)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            xlim(1000 * [cfg.analyse_states{3} cfg.analyse_states{4}])
            ylim([-20 20])
%             vline(50); vline(-50); vline(250); vline(-250);
            vline(0, 'k')
            % set(gca, 'XTick', PSTH_bins)
            ax = gca;
            ax.FontSize = 14;
        end
    end
end
save_figure_as('Suppression_Enhancement_SeparatedForTime_GroupedAccordingTo',output_folder,savePlot)

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
            sig_cond2      = Out.(T).(cfg.condition(cond2_num).name).sig_n_bins <= n_sig_bins & Out.(T).(cfg.condition(cond2_num).name).sig_n_bins > n_sig_bins;
            sig_cond1      = Out.(T).(cfg.condition(cond1_num).name).sig_n_bins > n_sig_bins & Out.(T).(cfg.condition(cond1_num).name).sig_n_bins <= n_sig_bins;
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
        set(gcf, 'Position', [2 381 1914 553])
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
            [cc, pp] = ...
                corrcoef(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                'Rows','complete');
            [p,~] = ... % p(1) is linear slope, p(2) is b-member
                polyfit(Out.(T).(cfg.condition(cond1_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range), ...
                Out.(T).(cfg.condition(cond2_num).name).FR_ModIndex_PcS((sig.([cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name]) | sig.(cfg.condition(cond1_num).name) | sig.(cfg.condition(cond2_num).name)) & within_range),1);
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
            title({[T ': N = ' num2str(unit_count)], ...
                ['cc = ' num2str(cc(cond2_num)) '; p = ' num2str(pp(cond2_num))], ...
                ['y = ' num2str(p(1)) '*x + ' num2str(p(2))]})
            if a == 1
                xlabel(['% signal change in ' cfg.condition(cond1_num).name])
                ylabel(['% signal change in ' cfg.condition(cond2_num).name])
                legend([ls_line s1 s2 s3 s4], {'linear fit', ['non-significant', 'sig. ' cfg.condition(cond2_num).name ' & sig. ' cfg.condition(cond1_num).name], ['sig. only ' cfg.condition(cond2_num).name], ['sig. only ' cfg.condition(cond1_num).name]})
            end
            ax = gca;
            ax.FontSize = 16;
        end
        save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name '_Scatter_Pc_signal_change'],output_folder,savePlot)
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
            
            if isequal(sig_cond2_cond1, logical(zeros(5,1))) | isequal(sig_cond2_cond1, logical(sig_cond2_cond1))
                continue
            end
            
            [cc(a), pp(a)] = corr(valid_cond1_times(sig_cond2_cond1), valid_cond2_times(sig_cond2_cond1));
            
        end
        
        s = scatterhist(all_valid_cond1_times, all_valid_cond2_times, ...
            'Group', groups, ...
            ...%'Style', 'bar', ...
            'Marker', '.o+x', ...
            'MarkerSize', 6, ...
            'Location', 'SouthWest', ...
            'Direction', 'out', ...
            'Color', [repmat(colors(1,:), 4, 1); repmat(colors(2,:), 4, 1); repmat(colors(3,:), 4, 1)]);
        hold on
        x = -250:250;
        for a = 1:3
            line(x, cc(a)*x, 'Color', colors(a,:))
        end
        
        legend('VPL: Non Significant', ['VPL: ' cfg.condition(cond1_num).name ' Sig.'], ['VPL: ' cfg.condition(cond2_num).name ' Sig'], ['VPL: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'], ...
            'dPul: Non Significant', ['dPul: ' cfg.condition(cond1_num).name ' Sig.'], ['dPul: ' cfg.condition(cond2_num).name ' Sig'], ['dPul: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'], ...
            'MD: Non Significant', ['MD: ' cfg.condition(cond1_num).name ' Sig.'], ['MD: ' cfg.condition(cond2_num).name ' Sig'], ['MD: ' cfg.condition(cond2_num).name ' & ' cfg.condition(cond1_num).name ' Sig.'])
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
        bar_data = [out.VPL.bin_counts_cond1; out.dPul.bin_counts_cond1; out.MD.bin_counts_cond1]';
        st_bar = bar(s(2), centers, bar_data);
        ylabel(s(2), 'Unit Count')
        % xlim(s(2), [-250 250])
        for ii = 1:3
            st_bar(ii).FaceColor = colors(ii,:);
        end
        clear st_bar
        
        % plot the side marginal histogram
        bar_data = [out.VPL.bin_counts_cond2; out.dPul.bin_counts_cond2; out.MD.bin_counts_cond2]';
        st_bar = barh(s(3), centers, bar_data);
        xlabel(s(3), 'Unit Count')
        % xlim(s(3), [0 max(bar_data, [], 'all')])
        % ylim(s(3), [-250 250])
        for ii = 1:3
            st_bar(ii).FaceColor = colors(ii,:);
        end
        clear st_bar
        save_figure_as([cfg.condition(cond1_num).name '_vs_' cfg.condition(cond2_num).name '_All_Areas_Scatter_Time_Max_Change'],output_folder,savePlot)
    end
    
    %% histogram of modulation indices - % signal change

    figure
    set(gcf, 'Position', [2 381 1914 553])
    
    for a = 1: N_Areas
        T = Ana_TargetBrainArea{a};
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            subplot(N_conditions,N_Areas,a + N_Areas*(c-1))
            colororder([1 1 1; 0.8500 0.3250 0.0980; 0 0.4470 0.7410])
            MIs_Pc_all  = Out.(T).(L).FR_ModIndex_PcS_signed;
            pos_idx     = Out.(T).(L).sig_sign > 0;
            neg_idx     = Out.(T).(L).sig_sign < 0;
            sign_all    = Out.(T).(L).sig_n_bins > n_sig_bins;
            
            non_sig = histcounts(MIs_Pc_all(~sign_all), -60:10:60);
            sig_pos = histcounts(MIs_Pc_all(pos_idx & sign_all), -60:10:60);
            sig_neg = histcounts(MIs_Pc_all(neg_idx & sign_all), -60:10:60);
            
            bar(-55:10:60, [non_sig; sig_pos; sig_neg], 'stacked')
            
            if c == 1
                title(T)
            end
            
            title([L ': ' T])
            
            if a == 1 && c == 1
                legend({'Non-Significant', 'Sig. Pos.', 'Sig. Neg.'}, 'Location', 'Best')
                xlabel('% signal change')
                ylabel('Unit Counts')
            end
            
        end
        
    end
    sgtitle('Modulation Magnitude in Rest and Task')
    save_figure_as('Histograms_Magnitude_Signal_Change_Significant',output_folder,savePlot)
    
end
end

function rounded=roundto(unrounded,n)
factor= 10^n;
rounded=round(unrounded*factor)/factor;
end

function make_correlation_plot(out,X,Y,sig,onlysignificant,ccol,texty)
x=out.(X);y=out.(Y);
axis square; box on;hold on;
% all not significant & NaN units
if ~isempty(x) && ~isempty(y)
    if ~ onlysignificant
        scatter(x(~sig) , y(~sig), 'filled', 'MarkerFaceColor',ccol/2) % make nonsignificant circles empty maybe
    end
    scatter(x(sig) , y(sig), 'filled', 'MarkerFaceColor',ccol)
    
    % do we really want hardcoded limits here?
    %             set(gca,'xlim',[0, 80]);
    %             set(gca,'ylim',[0, 15]);
    
    %xf = [min(out.quantSNR), max(out.quantSNR)];
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
    text(8,texty, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', ccol,'fontsize',14)
    % plot(SNR,y_fit+2*delta,'color',ccol,SNR,y_fit-2*delta,'color',ccol)
    
    %legend({'Non-Significant Units', 'Significant Units', 'Linear Fit (Overall)'}, 'Location', 'Best')
end
end

function save_figure_as(filename,basepath_to_save,savePlot)
if savePlot;
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(gcf);
end
end
