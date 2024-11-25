function keys=ecg_bna_get_unit_list(cfg,compute_unit_subsets)
% apply exclusion criteria, save lists of included and excluded units
% requires the processed part of cardioballistic analysis

keys=struct;
keys=ph_general_settings(cfg.project,keys);
project_specific_settings=[keys.db_folder 'ph_project_settings.m'];
run(project_specific_settings)
version_specific_settings=[keys.db_folder cfg.spikes_version filesep 'ph_project_version_settings.m'];
run(version_specific_settings)
keys.anova_table_file=[keys.basepath_to_save cfg.spikes_version filesep 'tuning_table_combined_CI.mat'];
keys.tuning_table = ph_load_tuning_table(keys); %% load tuning table

% find column indices
unit_column_index   = DAG_find_column_index(keys.tuning_table, 'unit_ID');
target_column_index = DAG_find_column_index(keys.tuning_table, 'target');
site_column_index   = DAG_find_column_index(keys.tuning_table, 'site_ID');
grid_x_column_index = DAG_find_column_index(keys.tuning_table, 'grid_x');
grid_y_column_index = DAG_find_column_index(keys.tuning_table, 'grid_y');
depths_column_index = DAG_find_column_index(keys.tuning_table, 'electrode_depth');

% create 'hemisphere' column in the tuning table
keys.tuning_table(2:end, end+1) = cellfun(@(x) x(end), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);
keys.tuning_table(1,end)        = {'hemisphere'};
hemisphere_column_index = DAG_find_column_index(keys.tuning_table, 'hemisphere');

% drop hemisphere indices in the 'target' column
keys.tuning_table(2:end, target_column_index) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);

% I. Unit list before exclusion - this one doesn't take number of
% R-peaks into account
% create unit list before exclusion criteria and save it
unit_ids    = keys.tuning_table(2:end, unit_column_index);
targets     = keys.tuning_table(2:end, target_column_index);
sites       = keys.tuning_table(2:end, site_column_index);
grid_x      = keys.tuning_table(2:end, grid_x_column_index);
grid_y      = keys.tuning_table(2:end, grid_y_column_index);
depths      = keys.tuning_table(2:end, depths_column_index);
hemispheres = keys.tuning_table(2:end, hemisphere_column_index);

filename = [cfg.unit_lists filesep 'unitInfo_before_exclusion'];
save(filename, 'unit_ids', 'targets', 'sites', 'depths', 'grid_x', 'grid_y', 'hemispheres')

% keep unit list before exclusion for later
unit_ids_before_exclusion    = unit_ids;
targets_before_exclusion     = targets;
sites_before_exclusion       = sites;
grid_x_before_exclusion      = grid_x;
grid_y_before_exclusion      = grid_y;
depths_before_exclusion      = depths;
hemispheres_before_exclusion = hemispheres;

% create table of unit numbers by area
write_unit_table(targets, filename);
clear unit_ids targets sites grid_x grid_y depths hemispheres

if compute_unit_subsets
    % II. Unit list for BOTH task and rest
    % 0. Set up exclusion crtieria
    SNR       = keys.tt.avg_SNR;
    % return to default keys
    keys.tt.avg_stability = [-Inf Inf]; % drop stability thresholds, we'll implement it from 'criteria' in the neuronal data
    keys.tt.avg_SNR       = [-Inf Inf]; % drop SNR thresholds, we'll implement it from 'criteria' in the neuronal data
    keys.tt.tasktypes     = {''};       % drop conditions here, we'll implement exclusion by those in the next step
    keys.tt.FR            = [0 Inf];    % drop FR, we'll use it on the next step
    keys.tt.n_spikes      = keys.cal.n_spikes;
    keys.monkey           = '';         % empty because we ignore which monkey it is basically
    % 1. Start with units with at least 1 spike (non-empty clusters)
    keys.tt.tasktypes     = {'Fsac_opt','Vsac_opt'}; % WE NEED BOTH CONDITIONS HERE
    keys                  = ph_tuning_table_correction(keys);
    % 1.1. create unit list
    unit_ids    = keys.tuning_table(2:end, unit_column_index);
    targets     = keys.tuning_table(2:end, target_column_index);
    sites       = keys.tuning_table(2:end, site_column_index);
    grid_x      = keys.tuning_table(2:end, grid_x_column_index);
    grid_y      = keys.tuning_table(2:end, grid_y_column_index);
    depths      = keys.tuning_table(2:end, depths_column_index);
    hemispheres = keys.tuning_table(2:end, hemisphere_column_index);
    filename = [cfg.unit_lists filesep 'unitInfo_after_condition_exclusion_stable'];
    save(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres')
    % 1.2. [passed] create and save unit count table
    write_unit_table(targets, filename);
    % 1.3. figure out units excluded by number of spikes and save those
    unit_ids_after_condition_exclusion_stable    = unit_ids;
    targets_after_condition_exclusion_stable     = targets;
    sites_after_condition_exclusion_stable       = sites;
    grid_x_after_condition_exclusion_stable      = grid_x;
    grid_y_after_condition_exclusion_stable      = grid_y;
    depths_after_condition_exclusion_stable      = depths;
    hemispheres_after_condition_exclusion_stable = hemispheres;
    % find excluded units
    [unit_ids, ia] = setdiff(unit_ids_before_exclusion, unit_ids_after_condition_exclusion_stable); % 1 - before current exclusion step, 2 - after
    targets     = targets_before_exclusion(ia);
    sites       = sites_before_exclusion(ia);
    grid_x      = grid_x_before_exclusion(ia);
    grid_y      = grid_y_before_exclusion(ia);
    depths      = depths_before_exclusion(ia);
    hemispheres = hemispheres_before_exclusion(ia);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_condition_stable'];
    save(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres')
    % 1.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites grid_x grid_y depths hemispheres ia
    
    % 2. Exclude by SNR
    % 2.0. load 'criteria' structure for units that passed FR cirterion
    curr_unit_ids = unit_ids_after_condition_exclusion_stable;
    dt = ecg_bna_load_variables(cfg, curr_unit_ids, 'per_unit_-0.25-0.25s', 'Output', {}, 0);
    % 2.1. [passed] exclude by SNR
    task_SNR_ids  = dt.criteria.SNR_F >= SNR(1) & dt.criteria.SNR_F <= SNR(2);
    rest_SNR_ids  = dt.criteria.SNR_V >= SNR(1) & dt.criteria.SNR_V <= SNR(2);
    SNR_ids       = task_SNR_ids & rest_SNR_ids;
    unit_ids      = unit_ids_after_condition_exclusion_stable(SNR_ids);
    targets       = targets_after_condition_exclusion_stable(SNR_ids);
    sites         = sites_after_condition_exclusion_stable(SNR_ids);
    grid_x        = grid_x_after_condition_exclusion_stable(SNR_ids);
    grid_y        = grid_y_after_condition_exclusion_stable(SNR_ids);
    depths        = depths_after_condition_exclusion_stable(SNR_ids);
    hemispheres   = hemispheres_after_condition_exclusion_stable(SNR_ids);
    filename = [cfg.unit_lists filesep 'unitInfo_after_SNR_exclusion_stable'];
    save(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres')
    % 2.2. [passed] save table
    write_unit_table(targets, filename);
    % 2.3. figure out units excluded by SNR and save those
    unit_ids_after_SNR_exclusion_stable    = unit_ids;
    targets_after_SNR_exclusion_stable     = targets;
    sites_after_SNR_exclusion_stable       = sites;
    grid_x_after_SNR_exclusion_stable      = grid_x;
    grid_y_after_SNR_exclusion_stable      = grid_y;
    depths_after_SNR_exclusion_stable      = depths;
    hemispheres_after_SNR_exclusion_stable = hemispheres;
    
    % find excluded units
    [unit_ids, ia] = setdiff(unit_ids_after_condition_exclusion_stable, unit_ids_after_SNR_exclusion_stable);
    targets     = targets_after_condition_exclusion_stable(ia);
    sites       = sites_after_condition_exclusion_stable(ia);
    grid_x      = grid_x_after_condition_exclusion_stable(ia);
    grid_y      = grid_y_after_condition_exclusion_stable(ia);
    depths      = depths_after_condition_exclusion_stable(ia);
    hemispheres = hemispheres_after_condition_exclusion_stable(ia);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_SNR_stable'];
    save(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres')
    % 2.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites grid_x grid_y depths hemispheres
    
    % III. Unit list for either task or rest
    % apply exclusion criteria - choose cells that had either a task or a rest
    % block
    % 1.0. choose units that had a rest block
    keys.tt.tasktypes = {'Fsac_opt'}; % fixation saccade - for rest
    keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
    
    % create 'hemisphere' column in the tuning table
    keys.tuning_table(2:end, end+1) = cellfun(@(x) x(end), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);
    keys.tuning_table(1,end)        = {'hemisphere'};
    hemisphere_column_index         = DAG_find_column_index(keys.tuning_table, 'hemisphere');
    
    % drop hemisphere indices in the 'target' column
    keys.tuning_table(2:end, target_column_index) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);
    keys = ph_tuning_table_correction(keys);
    unit_ids_after_exclusion_rest    = keys.tuning_table(2:end, unit_column_index);
    targets_after_exclusion_rest     = keys.tuning_table(2:end, target_column_index);
    sites_after_exclusion_rest       = keys.tuning_table(2:end, site_column_index);
    grid_x_after_exclusion_rest      = keys.tuning_table(2:end, grid_x_column_index);
    grid_y_after_exclusion_rest      = keys.tuning_table(2:end, grid_y_column_index);
    depths_after_exclusion_rest      = keys.tuning_table(2:end, depths_column_index);
    hemispheres_after_exclusion_rest = keys.tuning_table(2:end, hemisphere_column_index);
    
    % 1.1. exclude by SNR threshold
    % 1.1.0. load criteria from the data
    dt = ecg_bna_load_variables(cfg, unit_ids_after_exclusion_rest, 'per_unit_-0.25-0.25s', 'Output', {}, 0);
    % 1.1.1. [passed] exclude by SNR
    task_SNR_ids  = dt.criteria.SNR_F >= SNR(1) & dt.criteria.SNR_F <= SNR(2);
%     rest_SNR_ids  = dt.criteria.SNR_V >= SNR(1) & dt.criteria.SNR_V <= SNR(2);
%     SNR_ids       = task_SNR_ids & rest_SNR_ids;
    unit_ids_task    = unit_ids_after_exclusion_rest(task_SNR_ids);
    targets_task     = targets_after_exclusion_rest(task_SNR_ids);
    sites_task       = sites_after_exclusion_rest(task_SNR_ids);
    grid_x_task      = grid_x_after_exclusion_rest(task_SNR_ids);
    grid_y_task      = grid_y_after_exclusion_rest(task_SNR_ids);
    depths_task      = depths_after_exclusion_rest(task_SNR_ids);
    hemispheres_task = hemispheres_after_exclusion_rest(task_SNR_ids);
%     filename = [cfg.unit_lists filesep 'unitInfo_after_SNR_exclusion_stable'];
%     save(filename, 'unit_ids', 'targets', 'sites')
    % 2.2. [passed] save table
%     write_unit_table(targets, filename);
    % 2.3. figure out units excluded by SNR and save those
%     unit_ids_after_SNR_exclusion_stable = unit_ids;
%     targets_after_SNR_exclusion_stable  = targets;
%     sites_after_SNR_exclusion_stable    = sites;
    
    % 2.0. choose units that had a rest block
    keys.tt.tasktypes= {'Vsac_opt'}; % fixation saccade - for rest
    keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
    % create 'hemisphere' column in the tuning table
    keys.tuning_table(2:end, end+1) = cellfun(@(x) x(end), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);
    keys.tuning_table(1,end)        = {'hemisphere'};
    hemisphere_column_index = DAG_find_column_index(keys.tuning_table, 'hemisphere');
    % drop hemisphere indices in the 'target' column
    keys.tuning_table(2:end, target_column_index) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);
    keys = ph_tuning_table_correction(keys);
    unit_ids_after_exclusion_task    = keys.tuning_table(2:end, unit_column_index);
    targets_after_exclusion_task     = keys.tuning_table(2:end, target_column_index);
    sites_after_exclusion_task       = keys.tuning_table(2:end, site_column_index);
    grid_x_after_exclusion_task      = keys.tuning_table(2:end, grid_x_column_index);
    grid_y_after_exclusion_task      = keys.tuning_table(2:end, grid_y_column_index);
    depths_after_exclusion_task      = keys.tuning_table(2:end, depths_column_index);
    hemispheres_after_exclusion_task = keys.tuning_table(2:end, hemisphere_column_index);
    % 2.1. exclude by SNR threshold
    % 2.1.0. load criteria from the data
    dt = ecg_bna_load_variables(cfg, unit_ids_after_exclusion_task, 'per_unit_-0.25-0.25s', 'Output', {}, 0);
    % 2.1.1. [passed] exclude by SNR
    rest_SNR_ids  = dt.criteria.SNR_V >= SNR(1) & dt.criteria.SNR_V <= SNR(2);
    unit_ids_rest    = unit_ids_after_exclusion_task(rest_SNR_ids);
    targets_rest     = targets_after_exclusion_task(rest_SNR_ids);
    sites_rest       = sites_after_exclusion_task(rest_SNR_ids);
    grid_x_rest      = grid_x_after_exclusion_task(rest_SNR_ids);
    grid_y_rest      = grid_y_after_exclusion_task(rest_SNR_ids);
    depths_rest      = depths_after_exclusion_task(rest_SNR_ids);
    hemispheres_rest = hemispheres_after_exclusion_task(rest_SNR_ids);
    
    % 3. find units that had at least 1 rest block or 1 task block, figure out
    % targets for those
    unit_ids_after_exclusion = union(unit_ids_rest, unit_ids_task);
    ids_rest = ismember(unit_ids_after_exclusion, unit_ids_rest);
    ids_task = ismember(unit_ids_after_exclusion, unit_ids_task);
    ids_both = ids_rest + 2*ids_task; % 1 - rest, 2 - task, 3 - both
    targets_after_exclusion           = cell(length(unit_ids_after_exclusion), 1);
    targets_after_exclusion(ids_rest) = targets_rest;
    targets_after_exclusion(ids_task) = targets_task;
    sites_after_exclusion             = cell(length(unit_ids_after_exclusion), 1);
    sites_after_exclusion(ids_rest)   = sites_rest;
    sites_after_exclusion(ids_task)   = sites_task;
    grid_x_after_exclusion            = cell(length(unit_ids_after_exclusion), 1);
    grid_x_after_exclusion(ids_rest)  = grid_x_rest;
    grid_x_after_exclusion(ids_task)  = grid_x_task;
    grid_y_after_exclusion            = cell(length(unit_ids_after_exclusion), 1);
    grid_y_after_exclusion(ids_rest)  = grid_y_rest;
    grid_y_after_exclusion(ids_task)  = grid_y_task;
    depths_after_exclusion            = cell(length(unit_ids_after_exclusion), 1);
    depths_after_exclusion(ids_rest)  = depths_rest;
    depths_after_exclusion(ids_task)  = depths_task;
    hemispheres_after_exclusion       = cell(length(unit_ids_after_exclusion), 1);
    hemispheres_after_exclusion(ids_rest) = hemispheres_rest;
    hemispheres_after_exclusion(ids_task) = hemispheres_task;
    
    filename = [cfg.unit_lists filesep 'unitInfo_after_SNR_exclusion_selected'];
    unit_ids    = unit_ids_after_exclusion;
    targets     = targets_after_exclusion;
    sites       = sites_after_exclusion;
    grid_x      = grid_x_after_exclusion;
    grid_y      = grid_y_after_exclusion;
    depths      = depths_after_exclusion;
    hemispheres = hemispheres_after_exclusion;
    save(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    write_unit_table(targets, filename);
    
    % IV. List of units existing in at least in one condition
    unit_ids_after_condition_exclusion_selected = union(unit_ids_after_exclusion_rest, unit_ids_after_exclusion_task);
    ids_either_rest = ismember(unit_ids_after_condition_exclusion_selected, unit_ids_after_exclusion_rest);
    ids_either_task = ismember(unit_ids_after_condition_exclusion_selected, unit_ids_after_exclusion_task);
    ids_both = ids_either_rest + 2*ids_either_task; % 1 - rest, 2 - task, 3 - both
    targets_after_condition_exclusion_selected           = cell(length(unit_ids_after_condition_exclusion_selected), 1);
    targets_after_condition_exclusion_selected(ids_either_rest) = targets_after_exclusion_rest;
    targets_after_condition_exclusion_selected(ids_either_task) = targets_after_exclusion_task;
    sites_after_condition_exclusion_selected             = cell(length(unit_ids_after_condition_exclusion_selected), 1);
    sites_after_condition_exclusion_selected(ids_either_rest)   = sites_after_exclusion_rest;
    sites_after_condition_exclusion_selected(ids_either_task)   = sites_after_exclusion_task;
    grid_x_after_condition_exclusion_selected            = cell(length(unit_ids_after_condition_exclusion_selected), 1);
    grid_x_after_condition_exclusion_selected(ids_either_rest)  = grid_x_after_exclusion_rest;
    grid_x_after_condition_exclusion_selected(ids_either_task)  = grid_x_after_exclusion_task;
    grid_y_after_condition_exclusion_selected            = cell(length(unit_ids_after_condition_exclusion_selected), 1);
    grid_y_after_condition_exclusion_selected(ids_either_rest)  = grid_y_after_exclusion_rest;
    grid_y_after_condition_exclusion_selected(ids_either_task)  = grid_y_after_exclusion_task;
    depths_after_condition_exclusion_selected            = cell(length(unit_ids_after_condition_exclusion_selected), 1);
    depths_after_condition_exclusion_selected(ids_either_rest)  = depths_after_exclusion_rest;
    depths_after_condition_exclusion_selected(ids_either_task)  = depths_after_exclusion_task;
    hemispheres_after_condition_exclusion_selected       = cell(length(unit_ids_after_condition_exclusion_selected), 1);
    hemispheres_after_condition_exclusion_selected(ids_either_rest) = hemispheres_after_exclusion_rest;
    hemispheres_after_condition_exclusion_selected(ids_either_task) = hemispheres_after_exclusion_task;
    
    unit_ids    = unit_ids_after_condition_exclusion_selected;
    targets     = targets_after_condition_exclusion_selected;
    sites       = sites_after_condition_exclusion_selected;
    grid_x      = grid_x_after_condition_exclusion_selected;
    grid_y      = grid_y_after_condition_exclusion_selected;
    depths      = depths_after_condition_exclusion_selected;
    hemispheres = hemispheres_after_condition_exclusion_selected;
    filename = [cfg.unit_lists filesep 'unitInfo_after_condition_exclusion_selected'];
    save(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres', 'ids_both')
    % 1.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites grid_x grid_y depths hemispheres
    
    % V. Unit list of excluded units (that aren't units existing either task or rest)
    ids_excluded = ~ismember(unit_ids_before_exclusion, unit_ids);
    unit_ids_excluded = unit_ids_before_exclusion(ids_excluded);
    targets_excluded  = targets_before_exclusion(ids_excluded);
    sites_excluded    = sites_before_exclusion(ids_excluded);
    grid_x_excluded   = grid_x_before_exclusion(ids_excluded);
    grid_y_excluded   = grid_y_before_exclusion(ids_excluded);
    depths_excluded   = depths_before_exclusion(ids_excluded);
    hemispheres_excluded = hemispheres_before_exclusion(ids_excluded);
    unit_ids    = unit_ids_excluded;
    targets     = targets_excluded;
    sites       = sites_excluded;
    grid_x      = grid_x_excluded;
    grid_y      = grid_y_excluded;
    depths      = depths_excluded;
    hemispheres = hemispheres_excluded;
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_SNR_selected'];
    save(filename, 'unit_ids', 'targets', 'sites', 'grid_x', 'grid_y', 'depths', 'hemispheres')
    write_unit_table(targets, filename);
    %     clear T unique_areas ic unit_counts filename unit_ids_excluded targets_excluded unit_ids_excluded targets_excluded unit_ids_before_exclusion targets_before_exclusion
end

end

function write_unit_table(targets, filename)
[unique_areas, ~, ic] = unique(targets);
unit_counts = accumarray(ic,1);
T = table(unique_areas, unit_counts);
writetable(T, [filename '.xls'])
end
