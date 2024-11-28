clear all, close all

%% set up parameters
N_conditions     = 2;
N_areas          = 3;
condition_names  = {'Rest', 'Task'};
condition_colors = {'b', 'r'};
area_list        = {'VPL', 'dPul', 'MD'};

results_folder = ...
    {'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Bacchus_TaskRest\Population_time_domain_after_SNR_exclusion_stable_noLow_amplitude_ccs_any_VPL_dPul_MD\', ...
    'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest\Population_time_domain_after_SNR_exclusion_stable_noLow_amplitude_ccs_any_VPL_dPul_MD\'};

filenames = ...
    {'R-peakResponsiveness_Condition_Contingency_VPL.xls', 'R-peakResponsiveness_Condition_Contingency_dPul.xls', 'R-peakResponsiveness_Condition_Contingency_MD.xls'};

for m = 1:2 % loop through monkeys

    for a = 1:N_areas
        T = area_list{a};
        
        tbl = readmatrix([results_folder{m} filenames{a}],'Range','B2:C3');
        
        spl_num = 3*(m - 1) + a;
        subplot(2,3,spl_num)
        
        % set up correct labels for x-axis
        hm = heatmap(tbl);
        axis square
        hm.XData        = {'Rest', 'Task'};
        hm.XDisplayData = {'Rest', 'Task'}';
        hm.XLimits      = {'Rest', 'Task'};
        hm.XLabel       = {'Condition'};
        
        % set up correct labels for y-axis
        hm.YData        = {'Responsive', 'Non-Responsive'};
        hm.YDisplayData = {'Responsive', 'Non-Responsive'};
        hm.YLimits      = {'Responsive', 'Non-Responsive'};
        hm.YLabel       = {'Responsiveness'};
        
    end
    
end

