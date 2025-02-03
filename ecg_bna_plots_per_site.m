function ecg_bna_plots_per_site( data, cfg, PlotMethod,varargin )


%lfp_tfa_plot_hs_tuned_tfr_multiple_img  - Plots the LFP time frequency spectrogram
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, cfg, plottitle, results_file )
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, cfg, plottitle, results_file, cm )
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, cfg, plottitle, results_file, cm, plot_significant )
%
%
% INPUTS:
%       avg_tfr         - average LFP time frequency response for different
%       hand-space conditions to be compared
%		cfg     - struct containing the required settings
%           Required Fields: see settings/lfp_tfa_settings_example
%               1. baseline_method             - method used for baseline
%               normalization
%               2. compare.reach_hands          - hand labels to compare
%               3. compare.reach_spaces         - space labels to compare
%       plottitle       - title for the plot
%       results_file    - path to filename to store the resulting image
%       varargin        - colormap to be used (default = 'jet', can be any
%                       standard colormap additionally supported is 'bluewhitered')
%                       - flag to indicate if only significant difference
%                       bins should be plotted
%
% REQUIRES:	bluewhitered, export_fig
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_plot_site_average_tfr,
% lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_tfr_across_sites,
% bluewhitered, colormap, lfp_tfa_compute_difference_condition_tfr
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%


plot_significant = 0;
if nargin > 3
    plot_significant = varargin{2};
end

% colorbar title
if strcmp(cfg.lfp.normalization, 'zscore') % can change to baseline_method
    % 'zscore' - P_norm(t,f) = ( mean(real) - mean(shuffled) ) / std(shuffled)
    %     cbtitle = 'Z-score'; % for baseline normalization
    cbtitle = '(P - \mu) / std'; % for shuffle predictor normalization
    %     imscale = [-1, 1];
elseif strcmp(cfg.lfp.normalization, 'division')
    cbtitle = 'P / \mu';
elseif strcmp(cfg.lfp.normalization, 'subtraction')
    cbtitle = 'P - \mu';
elseif strcmp(cfg.lfp.normalization, 'relchange')
    cbtitle = '(P - \mu) / \mu';
elseif strcmp(cfg.lfp.normalization, 'none')
    cbtitle = 'not-Normalized';
end

% number of subplots required
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked','average phase'};
nsubplots=numel(plot_names);
ncolumns=2;
nrows=ceil(nsubplots/ncolumns);
results_folder=[cfg.sites_lfp_fldr filesep];

%% Smoothing Kernel here:
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

%% loop through conditions
for cn= 1:numel(data.condition)
    if isempty(fieldnames(data.condition(cn))) || isempty(data.condition(cn).event) % &&~isempty([sites_data(i).condition(cn).event.lfp])
        continue
    end    
    % this here is eventually indicating number of triggers for each alignment
    con_data=data.condition(cn).event;
    
    % create figure
    h = figure('units','normalized','position',[0 0 1 1]);
    
%    states_valid=[];
    collim{1}=[];
    collim{2}=[];
    collim{6}=[];
    
    %for hs = 1:size(con_data, 2)
%     if isempty(cat(3, con_data(:, hs).(PlotMethod).pow_mean)) % this is a strange break condition to be honest
%         continue;
%     end
    % concatenate tfs for different state windows for plotting
    concat.pow = [];
    concat.pha = [];
    concat.itpc = [];
    concat.itpcbp = [];
    concat.powbp = [];
    concat.lfp = [];
    concat.pow_sgnf = [];
    concat.itpc_sgnf = [];
    concat.itpcbp_sgnf = [];
    concat.powbp_sgnf = [];
    concat.lfp_sgnf = [];
    concat.tfr_time = [];
    concat.lfp_time = [];    
    concat.lfp_shufmean = [];
    concat.lfp_shufstd = [];
    concat.freq  = cfg.lfp.foi; %% this is actually in the settings...
    
    %event_info = struct();
    ticksamples_tfr=[];
    ticksamples_lfp=[];
    for e = 1:size(con_data, 2)
        shufmean=con_data(e).shuffled.lfp.mean;
        shufstd=con_data(e).shuffled.lfp.std;
        
        pow_sgnf        =con_data(e).significance.pow;
        itpc_sgnf       =con_data(e).significance.itpc;
        itpcbp_sgnf     =con_data(e).significance.itpcbp;
        powbp_sgnf      =con_data(e).significance.powbp;
        lfp_sgnf        =con_data(e).significance.lfp;
        
        pow_mean        =con_data(e).(PlotMethod).pow.mean;
        pha_mean        =con_data(e).(PlotMethod).pha.mean;
        itpc_mean       =con_data(e).(PlotMethod).itpc.mean;
        itpcbp_mean     =con_data(e).(PlotMethod).itpcbp.mean;
        powbp_mean      =con_data(e).(PlotMethod).powbp.mean;
        lfp_mean        =con_data(e).(PlotMethod).lfp.mean;
        
        tfr_time=con_data(e).tfr_time;
        lfp_time=con_data(e).time;
        
        onset_s = find(tfr_time <= 0, 1, 'last'); % state onset time
        start_s = 1; % state start time
        finish_s = length(tfr_time); % start end sample        
        ticksamples_tfr=[ticksamples_tfr length(concat.tfr_time)+[start_s onset_s finish_s]];
        
        onset_s = find(lfp_time <= 0, 1, 'last'); % state onset time
        start_s = 1; % state start time
        finish_s = length(lfp_time); % start end sample        
        ticksamples_lfp=[ticksamples_lfp length(concat.lfp_time)+[start_s onset_s finish_s]];
        
        %% smooth here !
        itpcbp_mean =smooth_bp(itpcbp_mean, gaussian_kernel,half_win);
        powbp_mean  =smooth_bp(powbp_mean,  gaussian_kernel,half_win);
        lfp_mean    =smooth_bp(lfp_mean,    gaussian_kernel,half_win);
        
        % concatenate across events with a NaN separation in between
        NaNseparator=100/25;
        concat.pow          = cat(3, concat.pow,        pow_mean,   nan(size(pow_mean, 1),      size(pow_mean, 2),   NaNseparator));
        concat.pha          = cat(3, concat.pha,        pha_mean,   nan(size(pha_mean, 1),      size(pha_mean, 2),   NaNseparator));
        concat.itpc         = cat(3, concat.itpc,       itpc_mean,  nan(size(itpc_mean, 1),     size(itpc_mean, 2),  NaNseparator));
        concat.itpcbp       = cat(3, concat.itpcbp,     itpcbp_mean,nan(size(itpcbp_mean, 1),   size(itpcbp_mean, 2),NaNseparator));
        concat.powbp        = cat(3, concat.powbp,      powbp_mean, nan(size(powbp_mean, 1),    size(powbp_mean, 2), NaNseparator));
        concat.lfp          = cat(3, concat.lfp,        lfp_mean,   nan(size(lfp_mean, 1),      size(lfp_mean, 2),   NaNseparator));
        concat.pow_sgnf     = cat(3, concat.pow_sgnf,   pow_sgnf,   nan(size(pow_sgnf, 1),      size(pow_sgnf, 2),   NaNseparator));
        concat.itpc_sgnf    = cat(3, concat.itpc_sgnf,  itpc_sgnf,  nan(size(itpc_sgnf, 1),   	size(itpc_sgnf, 2),  NaNseparator));
        concat.itpcbp_sgnf  = cat(3, concat.itpcbp_sgnf,itpcbp_sgnf,nan(size(itpcbp_sgnf, 1),   size(itpcbp_sgnf, 2),NaNseparator));
        concat.powbp_sgnf   = cat(3, concat.powbp_sgnf, powbp_sgnf, nan(size(powbp_sgnf, 1),    size(powbp_sgnf, 2), NaNseparator));
        concat.lfp_sgnf     = cat(3, concat.lfp_sgnf,   lfp_sgnf,   nan(size(lfp_sgnf, 1),      size(lfp_sgnf, 2),   NaNseparator));
        
        
        concat.lfp_shufmean = cat(3, concat.lfp_shufmean,  shufmean,  nan(size(shufmean, 1),     size(shufmean, 2),  NaNseparator));
        concat.lfp_shufstd  = cat(3, concat.lfp_shufstd,   shufstd,   nan(size(shufstd, 1),      size(shufstd, 2),   NaNseparator));
        
        concat.tfr_time     = [concat.tfr_time, tfr_time, nan(1, NaNseparator)];
        concat.lfp_time     = [concat.lfp_time, lfp_time,     nan(1, NaNseparator)];
        
%         % somehow needed for (not) labelling not existing alignments
%         if ~all(isnan(tfr_time))
%             states_valid=[states_valid e];
%         end
    end
    
    %% plot
    tfr_events.onset        = find(concat.tfr_time == 0);
    %tfr_events.name         ={con_data(states_valid).event_name};    
    tfr_events.name         ={con_data.event_name};    
    tfr_events.ticksamples  = sort(ticksamples_tfr);
    tfr_events.startsamples = ticksamples_tfr(1:3:end);
    tfr_events.endsamples   = ticksamples_tfr(3:3:end);
    tfr_events.ticks        =round(concat.tfr_time(tfr_events.ticksamples)*10)/10;
    
    lfp_events.onset        = find(concat.lfp_time == 0);
    %lfp_events.name         ={con_data(states_valid).event_name};    
    lfp_events.name         ={con_data.event_name};   
    lfp_events.ticksamples  = sort(ticksamples_lfp);
    lfp_events.ticks        =round(concat.lfp_time(lfp_events.ticksamples)*10)/10;
    
%     
%     if isfield(con_data, 'nsessions')
%         plottitle = [plottitle ' (nsessions = ' num2str(con_data(1, hs).nsessions) ')'];
%     elseif isfield(con_data, 'nsites')
%         plottitle = [plottitle ' (nsites = ' num2str(con_data(1, hs).nsites) ')'];
%     elseif isfield(con_data, 'ntrials') && ~isempty(con_data(1, hs).ntrials)
%         plottitle = [plottitle ' (ntrials = ' num2str(con_data(1, hs).ntrials) ')'];
%     end
%     
    
    %% POW and ITPC
    toplot={concat.pow,concat.itpc,concat.powbp,concat.itpcbp};
    sigplot = {concat.pow_sgnf, concat.itpc_sgnf,concat.powbp_sgnf,concat.itpcbp_sgnf};
    for sp=1:2 % frequency spectra
        sph(sp)=subplot(nrows, ncolumns, sp);
        image(1:size(toplot{sp},3), 1:numel(concat.freq), squeeze(toplot{sp}),'CDataMapping','scaled');
        set(gca,'YDir','normal');
        hold on;
        % calculate the significance here:
        significance = double(squeeze(sigplot{sp})); %+repmat(randi([0 1], 60,19),[1,11]); checking the plot results
        %significance(significance==0)=NaN;
        % plotting a contour around the significant parts:
        contour(1:size(toplot{sp},3),1:numel(concat.freq),significance,1,'linecolor','k')
        nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
        collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
        
        % horizontal lines to separate frequency bands
        fbandstart = unique(cfg.lfp.frequency_bands(:))';
        fbandstart_idx = zeros(size(fbandstart));
        for f = fbandstart
            f_idx = find(abs(concat.freq - f) == min(abs(concat.freq - f)), 1, 'first');
            line([tfr_events.startsamples' tfr_events.endsamples']-1/2, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
            fbandstart_idx(fbandstart == f) = f_idx;
        end
        
        set(gca,'TickDir','out')
        set(gca, 'ytick', fbandstart_idx);
        set(gca, 'yticklabel', fbandstart);
        % add 0.5 at end since the time value is the center of the bin
        % add 0 at beginning to make x-axis visible
        set(gca, 'ylim', [0.5,numel(concat.freq) + 0.5]);
        add_ticks_and_labels(tfr_events,[0.5,numel(concat.freq) + 0.5],8)
        
        set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5);
        ylabel('Frequency (Hz)');
        title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);            
    end
    
    %% Bandpassed POW and ITPC
    for sp=3:4
    sph(sp)=subplot(nrows, ncolumns, sp);             %% change color order to something nicer
    hold on;
    set(gca,'ColorOrder',jet(size(toplot{sp},2)));
    %plot(repmat(concat.lfp_time,size(concat.powbp,2),1)', squeeze(concat.powbp)')
    plot(squeeze(toplot{sp})')
    xlabel('Time(s)'); ylabel('Power (W)');
    
    % adding the signifiance horizontal lines:
    ylm = get(gca,'Ylim');
    stp = (ylm(2)-ylm(1))/20;
    ylm(1)=ylm(1)-size(toplot{sp},2)*stp;
    set(gca,'Ylim',ylm);
    
    significance = double(squeeze(sigplot{sp}));         % i needed to create concat.itpcbp_sgnf, it basically appends Nans for a (potential) separator with a second alignment
    significance(significance==0)=NaN;                          % replacing zeros with Nans means once we plot, lines will be discontinoous there
    multiplicator= ylm(1)+(1:size(significance,1))*stp;              % multiplicator basically defines position of significance line
    significance=significance.*repmat(multiplicator',1,size(significance,2));
    X=~all(isnan(diff(significance,1,2)),2);
    plot(significance(X,:)','linewidth',3);
    
    add_ticks_and_labels(lfp_events,ylm,stp)
    
    legend({strcat(num2str(round(cfg.lfp.frequency_bands(:,1))), '-',num2str(round(cfg.lfp.frequency_bands(:,2))), ' Hz')},'fontsize',3);
    title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);
    set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5); %%should be from lfp_events   
        
    end
    
   %% POW and ITPC
   toplot={concat.pow,concat.itpc,concat.powbp,concat.itpcbp,0,concat.pha};
   sp=6; % frequency spectra
   sph(sp)=subplot(nrows, ncolumns, sp);
   image(1:size(toplot{sp},3), 1:numel(concat.freq), squeeze(toplot{sp}),'CDataMapping','scaled');
   set(gca,'YDir','normal');
   hold on;
   
   nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
   collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
   
   % horizontal lines to separate frequency bands
   fbandstart = unique(cfg.lfp.frequency_bands(:))';
   fbandstart_idx = zeros(size(fbandstart));
   for f = fbandstart
       f_idx = find(abs(concat.freq - f) == min(abs(concat.freq - f)), 1, 'first');
       line([tfr_events.startsamples' tfr_events.endsamples']-1/2, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
       fbandstart_idx(fbandstart == f) = f_idx;
   end
   
   set(gca,'TickDir','out')
   set(gca, 'ytick', fbandstart_idx);
   set(gca, 'yticklabel', fbandstart);
   % add 0.5 at end since the time value is the center of the bin
   % add 0 at beginning to make x-axis visible
   set(gca, 'ylim', [0.5,numel(concat.freq) + 0.5]);
   add_ticks_and_labels(tfr_events,[0.5,numel(concat.freq) + 0.5],8)
   
   set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5);
   ylabel('Frequency (Hz)');
   title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);
   
    
    
    %% Evoked LFP
    sp=5;
    sph(sp)=subplot(nrows, ncolumns, sp);
    hold on;
    plot(squeeze(concat.lfp)','linewidth',1.5)
    line([0 0], ylim, 'color', 'k');
    xlabel('Time(s)'); ylabel('Voltage (V)');
    title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);
    
    if strcmp(PlotMethod,'real')
        lineprops={};
        shadedErrorBar(1:size(concat.lfp_shufmean,3), squeeze(concat.lfp_shufmean),squeeze(concat.lfp_shufstd),lineprops,1);
    end
    ylm = get(gca,'Ylim');
    significance = double(squeeze(concat.lfp_sgnf));
    significance(significance==0)=NaN;
    significance=significance.*ylm(1);
    % adding the signifiance horizontal lines:
    %if ~all(isnan(diff(significance,1,1)));
        plot(significance','linewidth',3);
    %end
    add_ticks_and_labels(lfp_events,ylm,diff(ylim)/10)
    set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5); %%should be from lfp_events 
    % end
    
    %% format spectra colors
    collim{3}=collim{1};
    collim{4}=collim{1};
    collim{5}=collim{1};
    for sp=1:6
        subplot(sph(sp));
        set(gca,'CLim',collim{sp})
        cm = colormap('jet');
        if nargin > 4
            cm = colormap(varargin{1});
        end
        cb = colorbar('EastOutside');%('North');
        if strcmp(PlotMethod,'normalized')
            set(get(cb,'title'),'string', cbtitle, 'fontsize',8);
        else
            if sp==1
                set(get(cb,'title'),'string', 'power(W)', 'fontsize',8);
            elseif sp==2
                set(get(cb,'title'),'string', 'ITPC', 'fontsize',8);
            end
        end
        switch sp
            case {1,2,6}
                colormap(cm);
            case {3,4,5}
                set(cb,'Visible','off');
        end
        
    end
    
    
    %% plot title...
    R=[con_data(:).real];R=[R(:).ntriggers];
    S=[con_data(:).shuffled];S=[S(:).ntriggers];
    
    realtriggers=[num2str(R') repmat('/',size(R'))]';realtriggers=realtriggers(:)';realtriggers=strrep(realtriggers,' ','');
    shuffledtriggers=[num2str(S') repmat('/',size(S'))]';shuffledtriggers=shuffledtriggers(:)';shuffledtriggers=strrep(shuffledtriggers,' ','');
    plottitle = [data.site_ID ' - ' data.target ', ' data.condition(cn).label ', ' num2str(cfg.lfp.n_permutations) ' shuffles, ntriggers:' realtriggers ' real, ' shuffledtriggers ' shuffled'];
    ntriggers = [ realtriggers ' real, ' shuffledtriggers ' shuffled'];
    
    results_file = fullfile(results_folder, [data.site_ID '_' data.condition(cn).label ' ' PlotMethod]);
    if strcmp(PlotMethod,'normalized')
        mtit([plottitle ' ' PlotMethod ' (' cfg.lfp.normalization ')'],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    else
        mtit([plottitle ' ' PlotMethod],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    end
    
    
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
    
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end
end


function add_ticks_and_labels(events,ylm,stp)


%state_ticks=round(concat.tfr_time(event_samples)*10)/10;

event_onsets    =events.onset;
event_names     =events.name;
event_samples   =events.ticksamples;
state_ticks     =events.ticks;

for so = event_onsets
    line([so so], ylm, 'color', 'k');
    %if isfield(con_data(event_onsets == so), 'event_name') && ~isempty(event_names(event_onsets == so))
    event_name = event_names{event_onsets == so};
    event_name=strrep(event_name,'_',' ');
    text(so+1, ylm(1)+stp, event_name, 'fontsize', 8, 'fontweight', 'bold');
    %end
end

% mark state onsets
set(gca,'xtick',event_samples(~isnan(state_ticks)))
set(gca,'xticklabels', state_ticks(~isnan(state_ticks)), 'fontsize', 6)
% add 0.5 since the time value is the center of the bin
% add 0 at the beginning to make the y-axis visib
end

function smoothed=smooth_bp(input,gaussian_kernel,half_win)

jnk = [];
concat_input = cat(3,(input(:,:,half_win:-1:1)),(input(:,:,:)));
concat_input = cat(3,concat_input, (input(:,:,end:-1:end-half_win+1)));
for m = 1: size(input,1)
    for k=1:size(input,2)
        jnk(m,k,:)= conv(squeeze(concat_input(m,k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
    end
end
clear concat_powbp
smoothed = jnk(:,:,half_win+1:end-half_win);
end


