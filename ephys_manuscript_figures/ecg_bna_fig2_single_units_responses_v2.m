clear all, close all

% This script plots figure 2 for the Vasileva, Kaduk et al., 2024
% manuscript

%% folder
dir2save = 'Y:\Manuscripts\2024_Thalamus_ephys_heart_brain\';

%% set up parameters
N_monkeys        = 2;
N_conditions     = 2;
N_areas          = 3;
monkey_names     = {'Bacchus','Magnus'};
condition_names  = {'Rest', 'Task'};
condition_colors = {[0 0 1], [1 0 0]};
area_list        = {'VPL', 'dPul', 'MD'};
area_colors      = {[1 0.53 0]; % VPL
                [0.05 0.65 0.7]; % dPul
                [1 0 0.6]}; % MD
figure_font      = 'Arial';
font_size        = 10;
BINS             = -250:5:250;

nRpeaks          = 400; % number of R-peaks for rasters
max_Rpeaks       = 8942;

ECG_fs = 2034.5; % Hz

%% load files
dataFolder = ...
    {'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Bacchus_TaskRest_state4\per_unit_-0.25-0.25s\', ...
    'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest_state4\per_unit_-0.25-0.25s\'};
dataFiles = ...
    {{'Bac_20210720_16_VPL_R.mat', 'Bac_20220322_06_dPul_L.mat', 'Bac_20220322_23_MD_R.mat'}, ...
    {'Mag_20230511_18_VPL_L.mat', 'Mag_20230524_29_dPul_L.mat', 'Mag_20230621_38_MD_L.mat'}};

% dataFiles = ...
%     {{'Bac_20210720_16_VPL_R.mat', '', ''}, ...
%     {'', '', ''}};

load('Y:\Data\BodySignals\ECG_CAP\Bacchus\20210720_ecg_cap.mat','out')
x_times   = 1000*[-0.5:1/ECG_fs:0.5];
ECG_1     = mean(out(1).ECG_Rpeaks_valid,1);
ECG_1_SD  = std(out(1).ECG_Rpeaks_valid,[],1);
clear out

load('Y:\Data\BodySignals\ECG_CAP\Magnus\20230511_ecg_cap.mat','out')
ECG_2    = mean(out(1).ECG_Rpeaks_valid,1);
ECG_2_SD = std(out(1).ECG_Rpeaks_valid,[],1);
clear out

f1 = figure;
set(gcf, 'Position',[541    42   776   954], ...
    'defaultUicontrolFontName',figure_font, ...
    'defaultUitableFontName',figure_font, ...
    'defaultAxesFontName',figure_font, ...
    'defaultTextFontName',figure_font, ...
    'defaultUipanelFontName',figure_font)
tiledlayout(6, 3, 'TileSpacing', 'none', 'Padding', 'none');

for m = 1:2 % loop through monkeys
    
    for a = 1:N_areas
        T = area_list{a};
        
        % load neuronal data
        if exist([dataFolder{m} dataFiles{m}{a}],'file') == 2
            data = load([dataFolder{m} dataFiles{m}{a}]);
        end
        
        nexttile
        for c=1:N_conditions
            L=condition_names{c};
            col=condition_colors{c};
            %% plot R-peak triggered PETHs
            %             subplot(6,N_areas,9*(m-1)+6+a)
            hold on
            box on
            xlim([-0.25 0.25]*1000)
            lineProps={'color',col,'linewidth',1};
            s1 = shadedErrorBar(BINS, data.Output.(L).SD, data.Output.(L).SD_SEM,lineProps,1);
            lineProps={'color',col,'linewidth',1,'linestyle',':'};
            s2 = shadedErrorBar(BINS, data.Output.(L).SDP, [data.Output.(L).SDPCu; data.Output.(L).SDPCL], lineProps,1);
            ypos = max(abs(data.Output.(L).SDsubstractedSDP));
            to_plot=abs(data.Output.(L).sig);
            to_plot(to_plot==0)=NaN;
            plot(BINS,to_plot*(mean(data.Output.(L).SD)+data.Output.(L).sig_sign*ypos),'color',col,'linewidth',5);
            set(gca, 'XTick', [-200 -100 0 100 200], 'XTickLabel', [-200 -100 0 100 200], 'FontSize', font_size)
            
            if c == 2
                
                title(['\color[rgb]{' num2str(area_colors{a}) '}' T],'FontSize',14)
                
                y_lims=get(gca,'ylim');
                
                % extend ylims by 1 on both sides
                y_lims(1) = y_lims(1)-1;
                y_lims(2) = y_lims(2)+1;
                set(gca,'YLim',y_lims)
                
                line([0 0],ylim,'color','k');
                
            end
            
            if m == 1 && a == 1
                
                xlabel('Time from R-peak, ms')
                ylabel('Firing Rate, Hz')
                
            elseif m == 2 && a == 3
                
                legend([s1.mainLine s2.mainLine s2.patch], ...
                    {'Average PETH','Average Surrogate Data','95% Confidence Interval'}, ...
                    'Location','northwest')
                
            end
            
            if m == 1 && a == 1 && c == 2
                
                yyaxis right
                hold on
                fill([x_times fliplr(x_times) x_times(1)],[ECG_1-ECG_1_SD fliplr(ECG_1+ECG_1_SD) ECG_1_SD(1)],[0 0 0],'FaceAlpha',0.3,'EdgeColor','none')
                plot(x_times,ECG_1,'k-')
                set(gca,'YTick',[],'YTickLabel',[],'YColor','k')
                
            end
            
            if m == 2 && a == 1 && c == 2
                
                yyaxis right
                hold on
                fill([x_times fliplr(x_times) x_times(1)],[ECG_2-ECG_2_SD fliplr(ECG_2+ECG_2_SD) ECG_2_SD(1)],[0 0 0],'FaceAlpha',0.3,'EdgeColor','none')
                plot(x_times,ECG_2,'k-')
                set(gca,'YTick',[],'YTickLabel',[],'YColor','k')
                
            end
        end
    end
    
    for c=1:N_conditions
        L=condition_names{c};
        col=condition_colors{c};
        
        for a = 1:N_areas
            T = area_list{a};
            
            % load neuronal data
            if exist([dataFolder{m} dataFiles{m}{a}],'file') == 2
                data = load([dataFolder{m} dataFiles{m}{a}]);
            end
            

            %% plot raster with 1 ms resolution
            figure(f1);
            nexttile
            % take 400 equally spaced R-peaks
            Rpeak_ids = round(1 : data.Output.(L).NrEvents/nRpeaks : data.Output.(L).NrEvents);
            [raster_row, raster_col] = find(data.Output.(L).raster(Rpeak_ids,:)); % row - R-peak number, col - bin number
            scatter(raster_col-250,raster_row,3,'filled','MarkerFaceColor',col,'MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
            hold on
            line([0 0],ylim,'color','k');
            set(gca, 'XTick', [-200 -100 0 100 200], 'XTickLabel', [-200 -100 0 100 200], 'FontSize', font_size)
            set(gca, 'YTick', [], 'YTickLabel', [])
            xlim([-250 250])
            
            box on
            if c == 1
                
                text(-235,330,['Neuron ' data.Output.unit_ID(5:end)],'Interpreter','none','BackgroundColor','w','EdgeColor','k','FontSize',font_size)
                
            end
            
            if m == 1 && a == 1
                set(gca, 'YTick', [1 400], 'YTickLabel', [1 400])
                ylabel('Number of R-peaks')
            end
            
        end
        
    end
    
end
save_figure_as('Fig2_',dir2save,1)

%% create a legend to paste it separately
figure,
set(gcf, 'Position', [1258 606 198 101], ...
    'defaultUicontrolFontName',figure_font, ...
    'defaultUitableFontName',figure_font, ...
    'defaultAxesFontName',figure_font, ...
    'defaultTextFontName',figure_font, ...
    'defaultUipanelFontName',figure_font)
p = plot(NaN, NaN,'b', NaN, NaN, 'r', 'LineWidth', 1);
axis off
legend({'Rest','Task'},'FontSize',font_size,'Orientation','horizontal')

save_figure_as('Fig2_legend',dir2save,1)

%% plot rasters for all R-peaks
for m = 1:2 % loop through monkeys

    for a = 1:N_areas
        T = area_list{a};
    
        %% plot rasters with all collected R-peaks
        figure;
        set(gcf, ...
            'defaultUicontrolFontName',figure_font, ...
            'defaultUitableFontName',figure_font, ...
            'defaultAxesFontName',figure_font, ...
            'defaultTextFontName',figure_font, ...
            'defaultUipanelFontName',figure_font)
        
        % load neuronal data
        if exist([dataFolder{m} dataFiles{m}{a}],'file') == 2
            data = load([dataFolder{m} dataFiles{m}{a}]);
        end
        
        currMaxRpeakNum = max([data.Output.Rest.NrEvents data.Output.Task.NrEvents]);
        set(gcf,'Position',[683    41   560   currMaxRpeakNum*955/max_Rpeaks]) % 955 - screen height
        for c=1:N_conditions
            L=condition_names{c};
            col=condition_colors{c};
            
            subplot(1,2,c)
            [raster_row, raster_col] = find(data.Output.(L).raster); % row - R-peak number, col - bin number
            scatter(raster_col-250,raster_row,3,'filled','MarkerFaceColor',col,'MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
            hold on
            line([0 0],ylim,'color','k');
            set(gca, 'XTick', [-200 -100 0 100 200], 'XTickLabel', [-200 -100 0 100 200])
            set(gca, 'YTick', [0 data.Output.(L).NrEvents])
            xlim([-250 250])
            ylim([0 data.Output.(L).NrEvents])
            title(condition_names{c},'Color',condition_colors{c})
            box on
            
            if c == 1
                
                xlabel('Time from R-peak, ms')
                ylabel('Number of Detected R-peaks')
                
            end
            
        end
        saveas(gcf,[dir2save filesep 'FigS2_' monkey_names{m} '_' T '.tif'])
        save_figure_as(['FigS2_' monkey_names{m} '_' T],dir2save,1)
        
    end
end

function save_figure_as(filename,basepath_to_save,savePlot)
if savePlot
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(gcf);
end
end
