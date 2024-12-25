clear all, close all

%% set up parameters
N_monkeys        = 2;
N_conditions     = 2;
N_areas          = 3;
subset_names     = {'stable', 'selected'};
monkey_names     = {'Bacchus','Magnus'};
condition_names  = {'Rest', 'Task'};
condition_colors = {[0 0 1], [1 0 0]};
nuclei_list      = {'VPL', 'dPul', 'MD'};
nuclei_nums      = [1 2 3];
figure_font      = 'Arial';
font_size        = 10;
output_folder    = 'Y:\Manuscripts\2024_Thalamus_ephys_heart_brain\';

for s = 1:2
    
    figure('Name',sprintf(['BarPlot_Pc_' subset_names{s}]),'Position',[728 347 600 466],'PaperPositionMode', 'auto', ...
        'defaultUicontrolFontName',figure_font, ...
        'defaultUitableFontName',figure_font, ...
        'defaultAxesFontName',figure_font, ...
        'defaultTextFontName',figure_font, ...
        'defaultUipanelFontName',figure_font);
    
    for m = 1:2 % loop through monkeys
        
        % load data
        load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' monkey_names{m} ...
            '_TaskRest_state4\Population_time_domain_per_unit_-0.25-0.25s__after_SNR_exclusion_' subset_names{s} '_noLow_amplitude_ccs_any_VPL_dPul_MD\Output.mat'])
        
        %     load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' monkey_names{m} ...
        %         '_TaskRest\Population_time_domain_after_SNR_exclusion_selected_noLow_amplitude_ccs_any_VPL_dPul_MD\Output.mat'])
        nb_grand_mat = [];
        for c=1:N_conditions
            L=condition_names{c};
            
            % prepare color matrix
            bar_colors_merged = [condition_colors{c}; 0.7 0.7 0.7];
            
            % prepare data matrix for plotting
            pc_mat = [];
            nb_mat = [];
%             nb_grand_mat = [];
            for a = 1: N_areas
                T=nuclei_list{a};
                pc_mat = vertcat(pc_mat, [Out.(T).(L).Pc_SignFR(1)+Out.(T).(L).Pc_SignFR(2) Out.(T).(L).Pc_SignFR(3)]);
                nb_mat = vertcat(nb_mat, [Out.(T).(L).Nb_SignFR(1)+Out.(T).(L).Nb_SignFR(2) Out.(T).(L).Nb_SignFR(3)]);
            end
            
            nb_grand_mat = cat(3,nb_grand_mat,nb_mat);
            
            % I. compute chi-square test and save statistics
            if c == 2
                [chi2_stat, p_val, stats] = DAG_chi_square_test(nb_grand_mat);
                
                filename = [monkey_names{m} '_Grand_chi-square_between_areas_' subset_names{s}];
                TBL = table({'Rest'; 'Task'}, chi2_stat', p_val', ...
                    'VariableNames', {'Condition', 'Chi^2 Stat.', 'P-values'});
                writetable(TBL, [output_folder filesep filename '.xls'])
                clear TBL
            end
            
            % II. Compute exact Fisher's test between nuclei
            nuclei_combinations = {};
            p_values            = [];
            odds_ratios         = [];
            for i = 1:length(nuclei_list)
                for j = i+1:length(nuclei_list)
                    
                    % Extract data for the two nuclei
                    data = [nb_mat(i,:); nb_mat(j,:)];
                    
                    [~, p, stats] = fishertest(data);
                    
                    % store data
                    nuclei_combinations = [nuclei_combinations; [nuclei_list{i} ' vs. ' nuclei_list{j}]];
                    p_values            = [p_values; p];
                    odds_ratios         = [odds_ratios; stats.OddsRatio];
                    
                end
            end
            
            % correct for multiple comparisons
            h = fdr_bky(round(p_values,10), 0.05);
            
            % put stuff into the table
            filename = [monkey_names{m} '_Fisher_test_between_areas_' subset_names{s} '_' condition_names{c}];
            TBL = table(nuclei_combinations, p_values, h, odds_ratios, ...
                'VariableNames', {'Nuclei Conbinations', 'Fisher''s p', 'Fisher''s h (BKY-corrected)', 'Odds Ratios'});
            writetable(TBL, [output_folder filesep filename '.xls'])
            clear TBL
            
            % plot bar plots
            subplot(2,N_conditions, 2*(m-1) + c);
            b = bar(pc_mat,'stacked');
            for ii = 1:length(b)
                b(ii).FaceColor = bar_colors_merged(ii,:);
            end
            
            % add percentages to bars
            xbarCnt = vertcat(b.XEndPoints);
            ybarTop = vertcat(b.YEndPoints);
            ybarCnt = ybarTop - pc_mat'/2;
            
            % create text strings
            data_mat(1:2:5,:) = pc_mat;
            data_mat(2:2:6,:) = nb_mat;
            txt = compose('%.1f%%\n(%d)',data_mat');
            
            th = text(xbarCnt(:), ybarCnt(:), txt(:), ...
                'HorizontalAlignment', 'center', ....
                'VerticalAlignment', 'middle', ...
                'Color', 'w',...
                'FontName', figure_font, ...
                'FontSize', 8, ...
                'FontWeight', 'normal');
            
            title(['Monkey ' monkey_names{m}(1) ': ' L])
            set(gca,'XTickLabel',nuclei_list,'fontsize',font_size);
            ylim([0 100])
            if c == 1
                ylabel('Percentage of Units', 'FontSize', font_size)
            end
            ax = gca;
            ax.FontSize = font_size;
            
        end
    end
    save_figure_as(['Fig3_' subset_names{s}],output_folder,1)
end

%% plot legend
legend_colors = [condition_colors{1}; condition_colors{2}; 0.7 0.7 0.7];
figure,
set(gcf, 'Position', [1258 606 198 101], ...
    'defaultUicontrolFontName',figure_font, ...
    'defaultUitableFontName',figure_font, ...
    'defaultAxesFontName',figure_font, ...
    'defaultTextFontName',figure_font, ...
    'defaultUipanelFontName',figure_font)
dummy_data = nan(3);
b = bar(dummy_data,'stacked');
for ii = 1:3
    b(ii).FaceColor = legend_colors(ii,:);
end
axis off
leg = legend({'Rest: Responsive', 'Task: Responsive', 'Non-Responsive'}, 'FontSize', font_size);
save_figure_as('Fig3_legend',output_folder,1)


function save_figure_as(filename,basepath_to_save,savePlot)
if savePlot
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(gcf);
end
end
