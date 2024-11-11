function ecg_bna_compute_IBI_raw_diff(cfg,withunits,type)


if cfg.combine_hemispheres
    cfg.targets=unique(cellfun(@(x) x(1:strfind(x,'_')-1),cfg.targets,'uniformoutput',false));
end
targets = cfg.targets;%unique({out.target}); %%cfg.monkey

if strcmp(type , 'normalized')
    
    IBIhigh_filepath = fullfile([cfg.LFP_root_results_fldr filesep 'grand_average_IBIhigh']);
    IBIlow_filepath = fullfile([cfg.LFP_root_results_fldr filesep 'grand_average_IBIlow' ]);
    filename = [cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Condition_diff_grand_avg_sessions_sites',withunits,'.mat'];
    
    cd(IBIhigh_filepath)
    % read the files
    all_IBIhigh_data = dir(['*Grand_grand_avg_sessions_sites',withunits,'.mat']);
    load(all_IBIhigh_data.name);
    IBIhigh = grand_avg;
    
    cd(IBIlow_filepath)
    % read the files
    all_IBIlow_data = dir(['*Grand_grand_avg_sessions_sites',withunits,'.mat']);
    load(all_IBIlow_data.name);
    IBIlow = grand_avg;
    clear grand_avg
    %%
    % storing all the Triggered parameters of the sites in the results folder
    condName_idx = find(~[IBIhigh.nSites]==0);
    condname = {IBIhigh(condName_idx(1)).avg.cond_name};
    for tr = 1:length(IBIhigh)
        for cn = 1:2
            IBIdiff(tr).target = IBIhigh(tr).target;
            IBIdiff(tr).nSites = IBIhigh(tr).nSites;
            IBIdiff(tr).avg(cn).condition = condname{cn};
            if (IBIdiff(tr).nSites>0)
                IBIdiff(tr).avg(cn).lfp = IBIhigh(tr).avg(cn).lfp - IBIlow(tr).avg(cn).lfp;
                IBIdiff(tr).avg(cn).itpc = IBIhigh(tr).avg(cn).itpc - IBIlow(tr).avg(cn).itpc;
                IBIdiff(tr).avg(cn).pow = IBIhigh(tr).avg(cn).pow - IBIlow(tr).avg(cn).pow;
                IBIdiff(tr).avg(cn).itpcbp = IBIhigh(tr).avg(cn).itpcbp - IBIlow(tr).avg(cn).itpcbp;
                IBIdiff(tr).avg(cn).powbp = IBIhigh(tr).avg(cn).powbp - IBIlow(tr).avg(cn).powbp;
                
                IBIdiff(tr).avg(cn).lfp_avg = IBIhigh(tr).avg(cn).lfp_avg - IBIlow(tr).avg(cn).lfp_avg;
                IBIdiff(tr).avg(cn).itpc_avg = IBIhigh(tr).avg(cn).itpc_avg - IBIlow(tr).avg(cn).itpc_avg;
                IBIdiff(tr).avg(cn).pow_avg = IBIhigh(tr).avg(cn).pow_avg - IBIlow(tr).avg(cn).pow_avg;
                IBIdiff(tr).avg(cn).itpcbp_avg = IBIhigh(tr).avg(cn).itpcbp_avg - IBIlow(tr).avg(cn).itpcbp_avg;
                IBIdiff(tr).avg(cn).powbp_avg = IBIhigh(tr).avg(cn).powbp_avg - IBIlow(tr).avg(cn).powbp_avg;
            else
                IBIdiff(tr).avg(cn).lfp = [];
                IBIdiff(tr).avg(cn).itpc = [];
                IBIdiff(tr).avg(cn).pow = [];
                IBIdiff(tr).avg(cn).itpcbp = [];
                IBIdiff(tr).avg(cn).powbp = [];
                
                IBIdiff(tr).avg(cn).lfp_avg = [];
                IBIdiff(tr).avg(cn).itpc_avg = [];
                IBIdiff(tr).avg(cn).pow_avg = [];
                IBIdiff(tr).avg(cn).itpcbp_avg = [];
                IBIdiff(tr).avg(cn).powbp_avg = [];
            end
        end
    end
    fileName2save = fullfile([cfg.LFP_root_results_fldr filesep, cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise-','IBIdifference_grand_avg_sites',withunits,'_',cfg.lfp.IBIdiff_type]);
    save(fileName2save,'IBIdiff', 'tfr_time', 'time');
    
elseif strcmp(type, 'raw')
    
    IBIhigh_filepath = fullfile([cfg.LFP_root_results_fldr filesep 'Per_Site_IBIhigh']);
    IBIlow_filepath = fullfile([cfg.LFP_root_results_fldr filesep 'Per_Site_IBIlow' ]);
    
    cd(IBIhigh_filepath)
    % read the files
    all_IBIhigh_data = dir('*.mat');
    
    cd(IBIlow_filepath)
    % read the files
    all_IBIlow_data = dir('*.mat');
    
    %%
    % storing all the Triggered parameters of the sites in the results folder:
    out = struct;
    i=0;
    for f = 1:length(all_IBIhigh_data)
        cd(IBIhigh_filepath)
        IBIhigh_data = load(all_IBIhigh_data(f).name);
        cd(IBIlow_filepath)
        IBIlow_data = load(all_IBIlow_data(f).name);
        condname = {IBIlow_data.triggered_site_data.condition.label};
        
        session = IBIhigh_data.triggered_site_data.session;
        site_ID = IBIhigh_data.triggered_site_data.site_ID;
        switch withunits
            case 'w_units'
                if ~ismember(site_ID,cfg.site_IDS)
                    continue
                end
            case 'wo_units'
                if ismember(site_ID,cfg.site_IDS)
                    continue
                end
            case 'all_sites'
        end
        i=i+1;
        target = IBIhigh_data.triggered_site_data.target;
        out(i).all_conditions_present = 1;
        
        for cn = 1:length(IBIhigh_data.triggered_site_data.condition)
            con_high=IBIhigh_data.triggered_site_data.condition(cn);
            con_low=IBIlow_data.triggered_site_data.condition(cn);
            if isempty(con_high.event)|| isempty(con_low.event)
                out(i).all_conditions_present = 0;
                continue;
            end
            condition = con_high.label;
            nTriggers = con_high.event.real.ntriggers;
            time      = con_high.event.time;
            tfr_time  = con_high.event.tfr_time;
            lfp       = squeeze(con_high.event.real.lfp.mean)'  - squeeze(con_low.event.real.lfp.mean)';
            itpc      = squeeze(con_high.event.real.itpc.mean)  - squeeze(con_low.event.real.itpc.mean);
            power     = squeeze(con_high.event.real.pow.mean)   - squeeze(con_low.event.real.pow.mean);
            itpcbp    = squeeze(con_high.event.real.itpcbp.mean)- squeeze(con_low.event.real.itpcbp.mean);
            powerbp   = squeeze(con_high.event.real.powbp.mean) - squeeze(con_low.event.real.powbp.mean);
            
            out(i).site_ID = site_ID;
            out(i).target = target;
            out(i).condition(cn).condition_name = condition;
            out(i).condition(cn).nTriggers = nTriggers;
            out(i).condition(cn).lfp = lfp;
            out(i).condition(cn).itpc = itpc;
            out(i).condition(cn).pow = power;
            out(i).condition(cn).itpcbp = itpcbp;
            out(i).condition(cn).powbp = powerbp;
            
            out(i).condition(cn).time = time;
            out(i).condition(cn).tfr_time = tfr_time;
        end
        
    end
    
    %% this part to remove sites with not all conditions?
    % out_mask = cellfun(@isempty, {out.site_ID});
    % out = out(~out_mask);
    out = out([out.all_conditions_present]==1);
    
    %%
    % Computing the Target-wise averaging of total available sites
    IBI_raw_diff = struct;
    for tr = 1: length(targets)
        sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{tr})),out);
        target_sites = out(sites_for_this_target);
        avg = struct('lfp',[],'itpc',[],'pow',[],'itpcbp',[],'powbp',[],...
            'lfp_avg',[],'itpc_avg',[],'pow_avg',[],'itpcbp_avg',[],'powbp_avg',[]);
        if ~isempty(target_sites)
            for cn = 1:2
                lfp_concat = [];
                itpc_concat = [];
                pow_concat = [];
                itpcbp_concat = [];
                powbp_concat = [];
                for trc = 1:length(target_sites)
                    lfp_concat = cat(3,lfp_concat,target_sites(trc).condition(cn).lfp);
                    itpc_concat = cat(3,itpc_concat,target_sites(trc).condition(cn).itpc);
                    pow_concat = cat(3,pow_concat,target_sites(trc).condition(cn).pow);
                    itpcbp_concat = cat(3,itpcbp_concat,target_sites(trc).condition(cn).itpcbp);
                    powbp_concat = cat(3,powbp_concat,target_sites(trc).condition(cn).powbp);
                end
                avg(cn).condition = target_sites(1).condition(cn).condition_name;
                avg(cn).lfp = lfp_concat;
                avg(cn).itpc = itpc_concat;
                avg(cn).pow = pow_concat;
                avg(cn).itpcbp = itpcbp_concat;
                avg(cn).powbp = powbp_concat;
                avg(cn).itpc_avg = mean(avg(cn).itpc,3);
                avg(cn).pow_avg = mean(avg(cn).pow,3);
                avg(cn).lfp_avg = mean(avg(cn).lfp,3);
                avg(cn).itpcbp_avg = mean(avg(cn).itpcbp,3);
                avg(cn).powbp_avg = mean(avg(cn).powbp,3);
            end
            
            IBI_raw_diff(tr).target = targets{tr};
            IBI_raw_diff(tr).nSites = length(target_sites);
            IBI_raw_diff(tr).avg = avg;
            
        else
            IBI_raw_diff(tr).target = targets{tr};
            IBI_raw_diff(tr).nSites = 0;
            IBI_raw_diff(tr).avg = [];
            continue;
        end
    end
    
    fileName2save = fullfile([cfg.LFP_root_results_fldr filesep, cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise-','IBIdifference_grand_avg_sites',withunits,'_',cfg.lfp.IBIdiff_type]);
    save(fileName2save,'IBI_raw_diff', 'tfr_time', 'time');
    IBIdiff = IBI_raw_diff;
    clear avg
end

%%
% plotting the results:
cond = condname;
freq = cfg.lfp.foi;
frequency_bands=cfg.lfp.frequency_bands;
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked'};

% Smoothing Kernel here:
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

% Task_Rest_raw_diff(tr).avg = IBIdiff(tr)
% if strcmp(type , 'normalized')


%% tfr_time not defined well!
for tr = 1: length(targets)
    %
    if (IBIdiff(tr).nSites >0)
        for cn = 1:2
            % create figure
            h(cn) = figure('units','normalized','position',[0 0 1 1]);
            toplot={IBIdiff(tr).avg(cn).pow_avg,IBIdiff(tr).avg(cn).itpc_avg};
            % =========================== Power ============================= %
            sp=1;
            sph(cn,sp)=subplot(3,2,sp);
            image(tfr_time, 1:numel(freq), squeeze(IBIdiff(tr).avg(cn).pow_avg),'CDataMapping','scaled');
            set(gca,'YDir','normal');
            line([0 0], ylim, 'color', 'k');
            
            % horizontal lines to separate frequency bands
            fbandstart = unique(cfg.lfp.frequency_bands(:))';
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                line([tfr_time(1) tfr_time(end)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            set(gca,'TickDir','out')
            set(gca, 'ytick', fbandstart_idx);
            set(gca, 'yticklabel', fbandstart);
            set(gca, 'ylim', [0.5,numel(freq) + 0.5]);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
            collim{cn,sp}=[min(nonnan(:)) max(nonnan(:))];
            %         collim{cn,sp}=[-3 3];% from Bacchus
            set(gca,'Xlim',[-.25 .25]);
            
            % =========================== ITPC ============================= %
            sp=2;
            sph(cn,sp)=subplot(3,2,sp);
            image(tfr_time, 1:numel(freq), squeeze(IBIdiff(tr).avg(cn).itpc_avg),'CDataMapping','scaled');
            set(gca,'YDir','normal');
            line([0 0], ylim, 'color', 'k');
            
            % horizontal lines to separate frequency bands
            fbandstart = unique(cfg.lfp.frequency_bands(:))';
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                line([tfr_time(1) tfr_time(end)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            set(gca,'TickDir','out')
            set(gca, 'ytick', fbandstart_idx);
            set(gca, 'yticklabel', fbandstart);
            set(gca, 'ylim', [0.5,numel(freq) + 0.5]);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
            collim{cn,sp}=[min(nonnan(:)) max(nonnan(:))];
            %         collim{cn,sp}=[0 .3];% from Bacchus
            set(gca,'Xlim',[-.25 .25]);
            
            %========================== Bandpassed POWER ==================== %
            % Smoothing of the POWbp here:
            jnk = [];
            concat_input = cat(2,(IBIdiff(tr).avg(cn).powbp_avg(:,half_win:-1:1)),(IBIdiff(tr).avg(cn).powbp_avg(:,:)));
            concat_input = cat(2,concat_input, (IBIdiff(tr).avg(cn).powbp_avg(:,end:-1:end-half_win+1)));
            for k=1:size(IBIdiff(tr).avg(cn).powbp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=3;
            sph(cn,sp)=subplot(3,2,sp);
            hold on;
            set(gca,'ColorOrder',jet(size(IBIdiff(tr).avg(cn).powbp_avg,1)),'xlim',[time(1) time(end)]);
            plot(repmat(time,size(IBIdiff(tr).avg(cn).powbp_avg,1),1)', squeeze(smoothed)')
            line([0 0], ylim, 'color', 'k');
            xlabel('Time(s)'); ylabel('Power (\mu V^2)');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            set(gca,'Xlim',[-.25 .25]);
            %         set(gca,'Ylim',[-3 3]);
            
            %========================== Bandpassed ITPC ==================== %
            % Smoothing of the ITPCbp here:
            jnk = [];
            concat_input = cat(2,(IBIdiff(tr).avg(cn).itpcbp_avg(:,half_win:-1:1)),(IBIdiff(tr).avg(cn).itpcbp_avg(:,:)));
            concat_input = cat(2,concat_input, (IBIdiff(tr).avg(cn).itpcbp_avg(:,end:-1:end-half_win+1)));
            for k=1:size(IBIdiff(tr).avg(cn).itpcbp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=4;
            sph(cn,sp)=subplot(3,2,sp);
            hold on;
            set(gca,'ColorOrder',jet(size(IBIdiff(tr).avg(cn).itpcbp_avg,1)),'xlim',[time(1) time(end)]);
            plot(repmat(time,size(IBIdiff(tr).avg(cn).itpcbp_avg,1),1)', squeeze(smoothed)')
            line([0 0], ylim, 'color', 'k');
            xlabel('Time(s)'); ylabel('ITPC');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            set(gca,'Xlim',[-.25 .25]);
            %         set(gca,'Ylim',[0 0.3]);
            
            %========================== LFP evoked Potential ==================== %
            % Smoothing of the  LFP evoked Potential here:
            jnk = [];
            lfp_se=sterr(IBIdiff(tr).avg(cn).lfp,3,0);
            lfp_std=std(IBIdiff(tr).avg(cn).lfp,0,3);
            percentile25=prctile(IBIdiff(tr).avg(cn).lfp,25,3);
            percentile75=prctile(IBIdiff(tr).avg(cn).lfp,75,3);
            concat_input = cat(2,(IBIdiff(tr).avg(cn).lfp_avg(:,half_win:-1:1)),(IBIdiff(tr).avg(cn).lfp_avg(:,:)));
            concat_input = cat(2,concat_input, (IBIdiff(tr).avg(cn).lfp_avg(:,end:-1:end-half_win+1)));
            for k=1:size(IBIdiff(tr).avg(cn).lfp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_mean = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(lfp_std(:,half_win:-1:1)),(lfp_std(:,:)));
            concat_input = cat(2,concat_input, (lfp_std(:,end:-1:end-half_win+1)));
            for k=1:size(IBIdiff(tr).avg(cn).lfp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_std = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
            concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
            for k=1:size(IBIdiff(tr).avg(cn).lfp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_25 = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
            concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
            for k=1:size(IBIdiff(tr).avg(cn).lfp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_75 = jnk(:,half_win+1:end-half_win);
            
            sp=5;
            sph(cn,sp)=subplot(3,2,sp);
            hold on;
            set(gca,'ColorOrder',jet(size(IBIdiff(tr).avg(cn).lfp_avg,1)),'xlim',[time(1) time(end)]);
            lineProps={'color',[0 0 1]};
            %shadedErrorBar(time,smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);
            shadedErrorBar(time,smoothed_mean,smoothed_std,lineProps,1);
            %plot(repmat(time,size(grand_avg(tr).avg.lfp_avg,1),1)', squeeze(smoothed_mean)')
            line([0 0], ylim, 'color', 'k');
            line([time(1) time(end)], [0 0], 'color', 'k');
            xlabel('Time(s)'); ylabel(' LFP evoked Potential');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            set(gca,'Xlim',[-.25 .25]);
            %         set(gca,'Ylim',[-60 40]);% for Bacchus
            %         set(gca,'Ylim',[-6 6]);% for Magnus
            
            %=================================================================%
            %% format spectra colors
            
            if strcmp(type , 'normalized')
                cbtitle = {'(P - \mu) / std','P - \mu'};
            elseif strcmp(type , 'raw')
                cbtitle = {'raw','raw'};
            end
            for sp=1:2
                subplot(3,2,sp);
                cm = colormap('jet');
                cb = colorbar;
                set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                colormap(cm);
            end
            
            results_file{cn} = fullfile([cfg.analyse_lfp_folder filesep, cfg.monkey,' - ',targets{tr},'-','IBIdifference_grand_avg_sites',withunits, '-',cond{cn},'-',cfg.lfp.IBIdiff_type]);
            mtit([ cfg.monkey,'-',targets{tr},'-avg of ',num2str(IBIdiff(tr).nSites),' sites ' ,withunits, ' - ',cond{cn},'-',cfg.lfp.IBIdiff_type],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
            
        end
        
        
        for cn = 1:length(cond)
            figure(h(cn));
            for sp=1:2
                subplot(sph(cn,sp));
                set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
            end
            export_fig(h(cn),[results_file{cn},'.pdf']);
            %         saveas(h,results_file{cn});
        end
        
    end
end
close all,
clc


end