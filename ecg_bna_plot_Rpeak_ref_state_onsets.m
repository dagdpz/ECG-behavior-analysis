function ecg_bna_plot_Rpeak_ref_state_onsets( Rpeak_evoked_states, ecg_bna_cfg, plottitle, results_file, varargin )
%ecg_bna_plot_evoked_lfp  - Plots the Rpeak onset probability for given
%time windows and different conditions to be compared
%
% USAGE:
%   ecg_bna_plot_Rpeak_ref_state_onsets( Rpeak_evoked_states, ecg_bna_cfg,
%   plottitle)
%
% INPUTS:
%       evoked_lfp       - average LFP power spectrum for different
%       hand-space conditions to be compared
%		ecg_bna_cfg      - struct containing the required settings
%           Required Fields: see ecg_bna_settings
%               compare.reach_hands     - reach hands
%               compare.reach_spaces    - reach spaces               
%       plottitle        - title for the plot
%       results_file     - path to filename to store the resulting image
%
% See also ecg_bna_compute_session_Rpeak_evoked_state_onsets,
% ecg_bna_avg_sessions_Rpeak_evoked_state_onsets

    h = figure;
    set(h, 'position', [100, 100,900, 675]);
    
    %set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    % number of subplots required
    nstates = size(ecg_bna_cfg.analyse_Rpeak_states, 1);
    nhandlabels = length(ecg_bna_cfg.compare.reach_hands);
    nspacelabels = length(ecg_bna_cfg.compare.reach_spaces);
    nhandspacelabels = nhandlabels * nspacelabels;
    colors = ['b'; 'r'; 'g'; 'y'; 'm'; 'c'; 'k'];
    
    % loop through handspace
    for hs = 1:size(Rpeak_evoked_states, 2)
        if ~isempty([Rpeak_evoked_states(:,hs).abs_timeprob]) && ...
                ~isempty([Rpeak_evoked_states(:,hs).rel_timeprob])
            for st = 1:size(Rpeak_evoked_states, 1)
                if isempty(Rpeak_evoked_states(st,hs).abs_timeprob)
                    continue;
                end
                if isempty(Rpeak_evoked_states(st,hs).rel_timeprob)
                    continue;
                end
                for cn = 1:length(Rpeak_evoked_states(st, hs).rel_timeprob)
                    if isempty(Rpeak_evoked_states(st, hs).rel_timeprob(cn).prob)
                        continue;
                    end
    %               % state onset probability vs. Rpeak phase
                    timebinedges = Rpeak_evoked_states(st, hs).rel_timeprob(cn).timebins * 2*pi;
                    timebinmiddle = (timebinedges(1:end-1) + timebinedges(2:end))/2;
                    mean_prob = nanmean(Rpeak_evoked_states(st, hs).rel_timeprob(cn).prob, 1);
                    null_hyp_prob = ones(size(mean_prob)) * (1/length(mean_prob));
                    std_prob = nanstd(Rpeak_evoked_states(st, hs).rel_timeprob(cn).prob, 1);
                    
                    if isfield(Rpeak_evoked_states(st, hs), 'colors') && ...
                            ~isempty(Rpeak_evoked_states(st, hs).colors)
                        plot_color = Rpeak_evoked_states(st, hs).colors(cn, :);
                    else
                        plot_color = colors(cn, :);
                    end

                    % state onset probability vs relative time from Rpeak
                    subplot(nhandspacelabels*2, nstates, (hs-1)*nstates*2 + st);
                    % join the mean
                    prob_plot = polar([timebinmiddle timebinmiddle(1)], ...
                        [mean_prob mean_prob(1)], '-');
                    prob_plot.LineWidth = 2;
                    prob_plot.Color = plot_color;
                    % SN(14.10.2019) - not showing polar plot when hold on is written before
                    % first polar plot
                    hold on;
                    polar([timebinmiddle timebinmiddle(1)], ...
                        [null_hyp_prob null_hyp_prob(1)], 'k-');
                    view([90 -90]);
                    %phase_plot.LineSize = 2;
                    % standard deviation
                    if size(Rpeak_evoked_states(st, hs).rel_timeprob(cn).prob, 1) > 1
                        meanplusstd = polar([timebinmiddle timebinmiddle(1)], ...
                            [mean_prob mean_prob(1)] + ...
                            [std_prob std_prob(1)], '--');
                        meanplusstd.Color = plot_color;
                        meanminusstd = polar([timebinmiddle timebinmiddle(1)], ...
                            [mean_prob mean_prob(1)] - ...
                            [std_prob std_prob(1)], '--');
                        meanminusstd.Color = plot_color;
                    end
                    %set(gca, 'XTick', timebinedges);
                    subplottitle = Rpeak_evoked_states(st, hs).state_name;
                    if isfield(Rpeak_evoked_states(st, hs), 'ntrials') && ...
                            ~isempty(Rpeak_evoked_states(st, hs).ntrials)
                        subplottitle = [subplottitle, sprintf(' (%g)', ...
                            Rpeak_evoked_states(st, hs).ntrials)];
                    end
                    if isfield(Rpeak_evoked_states, 'nsessions') && ...
                            ~isempty(Rpeak_evoked_states(st, hs).nsessions)
                        subplottitle = [subplottitle, ' (nsessions = ' ...
                            num2str(Rpeak_evoked_states(st, hs).nsessions) ')'];
                    end
                    title(subplottitle); %box on;
                    %xlabel('Rel. Time from Rpeak')
                    ylabel('Rpeak phase[�]');
                end
                
                for cn = 1:length(Rpeak_evoked_states(st, hs).abs_timeprob)    
                    if isempty(Rpeak_evoked_states(st, hs).abs_timeprob(cn).prob)
                        continue;
                    end
                    timebinedges = Rpeak_evoked_states(st, hs).abs_timeprob(cn).timebins;
                    timebinmiddle = (timebinedges(1:end-1) + timebinedges(2:end))/2;
    %                 if size(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1) > 1
    %                     plot(timebinmiddle, Rpeak_evoked_states(st, hs).abs_timeprob.prob, ...
    %                         'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.5);
    %                 end
    
%                     if isfield(Rpeak_evoked_states(st, hs), 'colors') && ...
%                             ~isempty(Rpeak_evoked_states(st, hs).colors)
%                         plot_color = Rpeak_evoked_states(st, hs).colors(cn, :);
%                     else
%                         plot_color = colors(cn, :);
%                     end
                    
                    % state onset probability vs absolute time around Rpeak
                    subplot(nhandspacelabels*2, nstates, (hs-1)*nstates*2 + nstates + st);
                    if isfield(Rpeak_evoked_states(st, hs), 'colors') && ...
                            ~isempty(Rpeak_evoked_states(st, hs).colors)
                        colors = Rpeak_evoked_states(st, hs).colors;
                    end
                    hold on;
                    % join the mean
                    plot(timebinmiddle, ...
                        nanmean(Rpeak_evoked_states(st, hs).abs_timeprob(cn).prob, 1), '-', ...
                        'Color', colors(cn, :), 'LineWidth', 2);
                end
                
                if isfield(Rpeak_evoked_states(st, hs), 'legend')
                    legend(Rpeak_evoked_states(st, hs).legend);
                end
                set(gca, 'XLim', [ecg_bna_cfg.analyse_Rpeak_states{st, 3:4}]);
                box on;
                % line at Rpeak onset
                line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
                %set(gca, 'XTick', timebinedges);
                %grid on;
                subplottitle = Rpeak_evoked_states(st, hs).state_name;
                if isfield(Rpeak_evoked_states(st, hs), 'ntrials') && ...
                        ~isempty(Rpeak_evoked_states(st, hs).ntrials)
                    subplottitle = [subplottitle, sprintf(' (%g)', ...
                        Rpeak_evoked_states(st, hs).ntrials)];
                elseif isfield(Rpeak_evoked_states, 'nsessions') && ...
                        ~isempty(Rpeak_evoked_states(st, hs).nsessions)
                    subplottitle = [subplottitle, ' (nsessions = ' ...
                        num2str(Rpeak_evoked_states(st, hs).nsessions) ')'];
                end
                title(subplottitle); box on;
                xlabel('Time from Rpeak (s)')
                ylabel('P(onset)');
                
                % standard deviation
                for cn = 1:length(Rpeak_evoked_states(st, hs).abs_timeprob)
                    if size(Rpeak_evoked_states(st, hs).abs_timeprob(cn).prob, 1) > 1
                        plot(timebinmiddle, ...
                            nanmean(Rpeak_evoked_states(st, hs).abs_timeprob(cn).prob, 1) + ...
                            nanstd(Rpeak_evoked_states(st, hs).abs_timeprob(cn).prob, 1), '--', ...
                            'Color', plot_color);
                        plot(timebinmiddle, ...
                            nanmean(Rpeak_evoked_states(st, hs).abs_timeprob(cn).prob, 1) - ...
                            nanstd(Rpeak_evoked_states(st, hs).abs_timeprob(cn).prob, 1), '--', ...
                            'Color', plot_color);
                    end
                end
            end
        end
    end
    
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h, results_file, '-pdf');
    %print(h, '-depsc', [results_file '.ai']);

end

