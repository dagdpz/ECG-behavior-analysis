function ecg_bna_plot_session_spike_histogram(session_info, cfg, cfg_spec)

% check the Matlab version
v = version('-release');
v = str2double(v(1:end-1));

N_conditions=numel(cfg.condition);
condition_colors={cfg.condition.color};

for numTiming = 1:length(cfg_spec.analyse_states)
    
    curr_analyse_states = cfg_spec.analyse_states{numTiming};
    
    basepath_to_save=[cfg.SPK_root_results_fldr filesep 'per_unit_' num2str(curr_analyse_states{3}) '-' num2str(curr_analyse_states{4}) 's'];
    if ~exist(basepath_to_save,'dir')
        mkdir(basepath_to_save);
    end
    
    fileList = dir([basepath_to_save, filesep, session_info.Monkey(1:3) '_' session_info.Date '*.mat']);
    
    for flNum = 1:length(fileList)
        load([basepath_to_save filesep fileList(flNum).name], 'Output')
        
        % figure out the present units and check that they match in Task and Rest
        if ~isequal(Output.target, Output.target)
            error('Numbers and names of units don''t match for task and rest')
        end
        
        BINS=(curr_analyse_states{1,3}:cfg_spec.PSTH_binwidth:curr_analyse_states{1,4})*1000;
        BINS_1ms = (curr_analyse_states{1,3}:0.001:curr_analyse_states{1,4})*1000;
        
        unit_ID = Output.unit_ID;
        target = Output.target;
        
        sgtitleText = {[Output.unit_ID '_' Output.target], ... %
            ['SNR: ' num2str(Output.quantSNR) '; Fano Factor: ' num2str(Output.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * Output.Single_rating) '%']};
        
        %% R-peak triggered spike response
        figure; % raster
        set(gcf, 'Position', [50    40   450*N_conditions   956])
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            col=condition_colors{c};
            
            subplot(3,N_conditions,c)
            [raster_row, raster_col] = find(Output.(L).raster);
            scatter(raster_col+curr_analyse_states{3}*1000,raster_row,2,'MarkerFaceColor',col,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',col,'MarkerEdgeAlpha',0.2);
            hold on
            line([0 0],ylim,'color','k');
            set(gca, 'XTick', cfg_spec.x_ticks, 'XTickLabel', cfg_spec.x_tick_labels)
            if ~isempty(Output.(L).raster)
                set(gca, 'YTick', [0 size(Output.(L).raster,1)], 'YTickLabel', [0 size(Output.(L).raster,1)])
                xlim([curr_analyse_states{3} curr_analyse_states{4}]*1000)
                ylim([0 size(Output.(L).raster,1)])
            end
            xlabel('Time from R-peak, ms')
            ylabel('Number of R-peaks')
            box on
            if ~isnan(Output.(L).raster)
%                 % figure out how many R-peaks we have and if < 100, plot
%                 % everything; if > 100, choose only 100
%                 if size(Output.(L).raster,1) <= 100
%                     stepSize = 1; % just step through all the Rpeaks one by one
%                 else
%                     % choose min interval of linearly spaced intervals
%                     % between 1 and the number of Rpeaks
%                     stepSize = min(diff(floor(linspace(1, size(Output.(L).raster,1), 100))));
%                 end
%                 % plot spikes only for 100 Rpeaks or everything if we have
%                 % fewer
%                 a = 1; % introduce row counter
%                 for RpeakNum = 1:stepSize:size(Output.(L).raster,1)
%                     x = Output.(L).raster(RpeakNum,:);
%                     line([BINS_1ms; BINS_1ms], [x; 2*x]+a-1, 'Color', col, 'LineWidth', 1)
%                     a = a + 1;
%                 end
                
%                 title({['Raster for ' L ': 1-ms bins, '], ['100 R-peak intervals out of ' num2str(size(Output.(L).raster,1))], '(linearly spaced)'})
            end
            
        end
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            subplot(3,N_conditions,c+N_conditions)
            
            plot_PRTHs(Output.(L), curr_analyse_states, BINS, c, cfg)
            
            set(gca, 'XTick', cfg_spec.x_ticks, 'XTickLabel', cfg_spec.x_tick_labels)
            ylabel('Firing Rate, Hz');
            xlabel('Time from R-peak, ms');
            
        end
        
        
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            col=condition_colors{c};
            subplot(3,N_conditions,c+2*N_conditions);
            hold on
            box on
            % yyaxis left % dont even know what this does - but it doesnt work
            % in matlab 2015
            hstR=Output.(L).Rds;
            hstP=Output.(L).Rds_perm;
            hstR=hstR/sum(hstR);
            hstP=hstP/sum(hstP);
            if ~isnan(hstR)
                stairs(cfg_spec.histbins, hstR, 'Color', col, 'LineWidth', 2)
            end
            if ~isnan(hstP)
                stairs(cfg_spec.histbins, hstP, 'Color', [0.5,0.5,0.5], 'LineWidth', 2)
            end
            ylabel('fraction of intervals');
            xlabel('RR Duration, s');
            title('RR durations');
            legend({L, 'jittered'})
        end
        filename= ['Raster_PSTH_Rintervals_' unit_ID, '__' target];
        if v < 2016
            sgtitleText = [Output.unit_ID '_' Output.target,'SNR: ' num2str(Output.quantSNR) '; Fano Factor: ' num2str(Output.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * Output.Single_rating) '%'];
            mtit(sgtitleText,'xoff', 0, 'yoff', 0.2,'interpreter','none'); %% feel free to make this work with a cell, but suptitle is not a thing in matlab 2015
        else
            sgtitle(sgtitleText,'interpreter','none')
        end
        
        saveas(gcf,[basepath_to_save, filesep, filename], 'fig')
        saveas(gcf,[basepath_to_save, filesep, filename], 'tiff')
%         export_fig(gcf,[basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
        close(gcf);
        
        %% lowIBI and highIBI: R-peak responses
        figure; % raster
        set(gcf, 'Position', [50    40   450*N_conditions   600])
        
        for c=1:N_conditions
            L=cfg.condition(c).name;
            
            subplot(2,N_conditions,c)
            plot_PRTHs(Output.(L).lowIBI, curr_analyse_states, BINS, c, cfg)
            title([L ': lowIBI'])
            set(gca, 'XTick', cfg_spec.x_ticks, 'XTickLabel', cfg_spec.x_tick_labels)
            ylabel('Firing Rate, Hz');
            xlabel(['Time from ' cfg_spec.analyse_states{numTiming}{2} ', ms']);
            
            subplot(2,N_conditions,2+c)
            plot_PRTHs(Output.(L).highIBI, curr_analyse_states, BINS, c, cfg)
            title([L ': highIBI'])
            set(gca, 'XTick', cfg_spec.x_ticks, 'XTickLabel', cfg_spec.x_tick_labels)
            ylabel('Firing Rate, Hz');
            xlabel(['Time from ' cfg_spec.analyse_states{numTiming}{2} ', ms']);
            
        end
        
        if v < 2016
            sgtitleText = [Output.unit_ID '_' Output.target,'SNR: ' num2str(Output.quantSNR) '; Fano Factor: ' num2str(Output.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * Output.Single_rating) '%'];
            mtit(sgtitleText,'xoff', 0, 'yoff', 0.2,'interpreter','none'); %% feel free to make this work with a cell, but suptitle is not a thing in matlab 2015
        else
            sgtitle(sgtitleText,'interpreter','none')
        end
        
        filename = ['lowIBI_highIBI_Raster_PSTH_Rintervals_' Output.unit_ID '_' Output.target];
        export_fig(gcf,[basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
        close(gcf);
        
%         %% cosine and von Mises fits
%         figure;
%         set(gcf, 'Position', [2 38 500*numel(cfg.condition) 958])
%         if v >= 2020
%             sgtitle(sgtitleText, 'interpreter', 'none')
%         end
%         for c=1:numel(cfg.condition)
%             L=cfg.condition(c).name;
%             
%             distList = {'linear', 'cosine', 'vonMisesPos', 'vonMisesNeg'};
%             for distNum = 1:4
%                 currFit = distList{distNum};
%                 subplot(4,length(cfg.condition),(distNum-1)*length(cfg.condition)+c)
%                 yyaxis left
%                 line(Output.phase_bin_centers, Output.(L).SD(ismember(BINS, Output.bin_centers)), 'Color', cfg.condition(c).color, 'LineWidth', 2)
%                 yyaxis right
% %                 line(Output.phase_bin_centers, Output.(L).cosine.average, 'Color', cfg.condition(c).color, 'LineWidth', 2) % plot average data
% %                 hold on
%                 plot(Output.phase_bin_centers, Output.(L).(currFit).yfit, 'k--')
%                 if strcmp(currFit, 'linear')
%                     title({[currFit ': slope = ' num2str(Output.(L).linear.coefs(2)) '; p = ' num2str(Output.(L).linear.pvalue(2))], ...
%                         ['R^2 = ' num2str(Output.Task.linear.rsquared)]})
%                 elseif strcmp(currFit, 'cosine')
%                     if ~isnan(Output.(L).(currFit).coefs(2))
%                         xline(Output.(L).(currFit).coefs(2), ':', 'Color', cfg.condition(c).color)
%                     end
%                     title({[currFit 'Mean Phase = ' num2str(Output.(L).(currFit).coefs(2)) '; R^2 = ' num2str(Output.(L).(currFit).rsquared) '; p = ' num2str(Output.(L).(currFit).pvalue)]})
%                 elseif strcmp(currFit, 'vonMisesPos') || strcmp(currFit, 'vonMisesNeg')
%                     if ~isnan(Output.(L).(currFit).coefs(4))
%                         xline(Output.(L).(currFit).coefs(4), ':', 'Color', cfg.condition(c).color)
%                     end
%                     title({['Mean Phase = ' num2str(Output.(L).(currFit).coefs(4)) '; R^2 = ' num2str(Output.(L).(currFit).rsquared) '; p = ' num2str(Output.(L).(currFit).pvalue)], ...
%                         ['a1 = ' num2str(Output.(L).(currFit).coefs(1)) '; d1 = ' num2str(Output.(L).(currFit).coefs(2)) '; \kappa = ' num2str(Output.(L).(currFit).coefs(3))]})
%                 end
%                 
%                 xlim([min(Output.phase_bin_centers) max(Output.phase_bin_centers)])
%                 ylabel('Spikes per bin across heart-cycles')
%                 legend({'Average Smoothed PSTH', 'Circular Fit'}, 'Location', 'best')
%                 box on
%             end
%             
%         end
%         
%         filename = ['FittedPSTH__' Output.unit_ID, '_' Output.target];
%         export_fig(gcf, [basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
%         close(gcf);
        
%         %% plot fitted data by low/high IBIs
%         figure;
%         set(gcf, 'Position', [1 38 1919 958])
%         if v >= 2020
%             sgtitle(sgtitleText, 'interpreter', 'none')
%         end
%         
%         distList = {'linear', 'cosine', 'vonMisesPos', 'vonMisesNeg'};
%         intList  = {'lowIBI_', 'highIBI_'};
%         
%         for distNum = 1:4
%             currFit = distList{distNum};
%             
%             for IBIgroupNum = 1:length(intList)
%                 
%                 for c=1:numel(cfg.condition)
%                     L=cfg.condition(c).name;
%                     
%                     subplot(4, length(cfg.condition)*length(intList), (distNum-1)*(length(cfg.condition)+ length(intList)) + (c-1)*2 + IBIgroupNum)
%                     yyaxis left
%                     line(Output.phase_bin_centers, Output.(L).([intList{IBIgroupNum} 'SD'])(ismember(BINS, Output.bin_centers)), 'Color', cfg.condition(c).color, 'LineWidth', 2)
%                     yyaxis right
% %                     line(Output.phase_bin_centers, Output.(L).([intList{IBIgroupNum} 'cosine']).average, 'Color', cfg.condition(c).color, 'LineWidth', 2) % plot average data
% %                     hold on
%                     plot(Output.phase_bin_centers, Output.(L).([intList{IBIgroupNum} currFit]).yfit, 'k--')
%                     if strcmp(currFit, 'linear')
%                         title({[intList{IBIgroupNum}(1:end-1) ' ' currFit ': slope = ' num2str(Output.(L).([intList{IBIgroupNum} 'linear']).coefs(2)) '; p = ' num2str(Output.(L).([intList{IBIgroupNum} 'linear']).pvalue(2))], ...
%                             ['R^2 = ' num2str(Output.(L).([intList{IBIgroupNum} 'linear']).rsquared)]})
%                     elseif strcmp(currFit, 'cosine')
%                         if ~isnan(Output.(L).([intList{IBIgroupNum} currFit]).coefs(2))
%                             xline(Output.(L).([intList{IBIgroupNum} currFit]).coefs(2), ':', 'Color', cfg.condition(c).color)
%                         end
%                         title({[intList{IBIgroupNum}(1:end-1) ' ' currFit 'Mean Phase = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).coefs(2)) '; R^2 = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).rsquared) '; p = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).pvalue)]})
%                     elseif strcmp(currFit, 'vonMisesPos') || strcmp(currFit, 'vonMisesNeg')
%                         if ~isnan(Output.(L).([intList{IBIgroupNum} currFit]).coefs(4))
%                             xline(Output.(L).([intList{IBIgroupNum} currFit]).coefs(4), ':', 'Color', cfg.condition(c).color)
%                         end
%                         title({[intList{IBIgroupNum}(1:end-1) ' Mean Phase = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).coefs(4)) '; R^2 = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).rsquared) '; p = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).pvalue)], ...
%                             ['a1 = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).coefs(1)) '; d1 = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).coefs(2)) '; \kappa = ' num2str(Output.(L).([intList{IBIgroupNum} currFit]).coefs(3))]})
%                     end
%                     
%                     xlim([min(Output.phase_bin_centers) max(Output.phase_bin_centers)])
%                     ylabel('Spikes per bin across heart-cycles')
%                     legend({'Average PSTH', 'Circular Fit'}, 'Location', 'best')
%                     box on
%                 end
%                 
%             end
%         end
%         
%         filename = ['LowHighIBI_FittedPSTH__' Output.unit_ID, '_' Output.target];
%         export_fig(gcf, [basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
%         close(gcf);
        
%         %% AIC and BIC
%         figure,
%         set(gcf, 'Position', [381 607 821 389])
%         if v >= 2020
%             sgtitle(sgtitleText, 'interpreter', 'none')
%         end
%         
%         for c=1:numel(cfg.condition)
%             L=cfg.condition(c).name;
%             
%             [maxAIC, maxAICid] = max([Output.(L).linear.aic Output.(L).cosine.aic Output.(L).vonMisesPos.aic Output.(L).vonMisesNeg.aic]);
%             [maxBIC, maxBICid] = max([Output.(L).linear.bic Output.(L).cosine.bic Output.(L).vonMisesPos.bic Output.(L).vonMisesNeg.bic]);
%             
%             subplot(1,2,c)
%             bar([Output.(L).linear.aic Output.(L).cosine.aic Output.(L).vonMisesPos.aic Output.(L).vonMisesNeg.aic; ...
%                 Output.(L).linear.bic Output.(L).cosine.bic Output.(L).vonMisesPos.bic Output.(L).vonMisesNeg.bic]')
%             hold on
%             if ~isnan(maxAIC)
%                 plot(maxAICid, 0, 'or')
%             end
%             if ~isnan(maxBIC)
%                 plot(maxBICid, 0, 'xb')
%             end
%             set(gca, 'XTickLabel', {'linear', 'cosine', 'pos. von Mises', 'neg. von Mises'})
%             xlabel('Type of Fit')
%             ylabel('AIC or BIC value')
%         end
%         
%         filename= ['AIC_BIC_' Output.unit_ID, '_' Output.target];
%         export_fig(gcf,[basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
%         close(gcf);
    end
end
end

function plot_PRTHs(data, curr_analyse_states, BINS, c, cfg)
hold on
box on
xlim([curr_analyse_states{3} curr_analyse_states{4}]*1000)
set(gca, 'XTick', BINS(1:5:end), 'XTickLabel', BINS(1:5:end))
col=cfg.condition(c).color;
lineProps={'color',col,'linewidth',1};
shadedErrorBar(BINS, data.SD, data.SD_SEM, lineProps,1);
lineProps={'color', col, 'linewidth', 1, 'linestyle', ':'};
shadedErrorBar(BINS, data.SDP, [data.SDPCu;data.SDPCL], lineProps,1);
line([0 0],ylim,'color','k');
ypos=NaN;
if data.sig_sign==-1;
    ypos=min(data.SD)*-1;
elseif data.sig_sign==1;
    ypos=max(data.SD);
end
to_plot=data.sig;
to_plot(to_plot==0)=NaN;
plot(BINS,to_plot*ypos,'color',col,'linewidth',5);

y_lims=get(gca,'ylim');
if isfield(data,'NrTrials')
    text(BINS(10),y_lims(2)-diff(y_lims)*c*1/20, [cfg.condition(c).name ': trials = ' ,num2str(data.NrTrials), ' events = ' ,num2str(data.NrEvents) ],'Color',col, 'Interpreter', 'none');
end
end
