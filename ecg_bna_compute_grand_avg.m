function grand_avg = ecg_bna_compute_grand_avg(cfg,withunits)

reprocess=1;
fileName = fullfile([cfg.analyse_lfp_folder filesep cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Grand_grand_avg_sessions_sites',withunits,'.mat']);

if cfg.combine_hemispheres
    cfg.targets=unique(cellfun(@(x) x(1:strfind(x,'_')-1),cfg.targets,'uniformoutput',false));
end
targets = cfg.targets;%unique({sites.target}); %%cfg.monkey

if reprocess
    data_path = cfg.sites_lfp_fldr;
    cd(data_path)
    % read the files
    all_lfp_data = dir('*.mat');
    % storing all the Triggered parameters of the sites in the results folder:
    sites = struct;
    s=0;
    for f = 1:length(all_lfp_data)
        load(all_lfp_data(f).name)
        session = triggered_site_data.session;
        site_ID = triggered_site_data.site_ID;
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
        s=s+1;
        sites(s).all_conditions_present = 1;        
        sites(s).site_ID = site_ID;
        sites(s).target = triggered_site_data.target;
        
        for c = 1:length(triggered_site_data.condition)
            con=triggered_site_data.condition(c);
            if isempty(con.event)
                sites(s).all_conditions_present = 0;
                continue;
            end
            sites(s).condition(c).condition_name = con.label;
            for e=1:size(cfg.analyse_states,1)
                event=con.event(e);
                %% think about event present?
                
                nTriggers = event.real.ntriggers;
                time      = event.time;
                tfr_time  = event.tfr_time;
                lfp       = squeeze(event.normalized.lfp.mean)';
                itpc      = squeeze(event.real.itpc.mean)-squeeze(event.shuffled.itpc.mean);
                power     = squeeze(event.normalized.pow.mean);
                itpcbp    = squeeze(event.real.itpcbp.mean)-squeeze(event.shuffled.itpcbp.mean);
                powerbp   = squeeze(event.normalized.powbp.mean);
                
                add_index=[0:size(itpcbp,1)-1]'*size(itpcbp,2);
                [~,max_itpcbp_idx] = max(abs(itpcbp),[],2);
                max_itpcbp_time=time(max_itpcbp_idx)';
                tmp=itpcbp';
                max_itpcbp=tmp(max_itpcbp_idx+add_index);
                [~,max_powbp_idx] = max(abs(powerbp),[],2);
                tmp=powerbp';
                max_powbp=tmp(max_powbp_idx+add_index);
                max_powbp_time=time(max_powbp_idx)';
                
                lfp_sig=squeeze(event.significance.lfp)'; %% weird dimension thing
                itpc_sig=squeeze(event.significance.itpc);
                pow_sig=squeeze(event.significance.pow);
                itpcbp_sig=squeeze(event.significance.itpcbp);
                powbp_sig=squeeze(event.significance.powbp);
                                
                sites(s).condition(c).event(e).nTriggers = nTriggers;
                sites(s).condition(c).event(e).lfp = lfp;
                sites(s).condition(c).event(e).itpc = itpc;
                sites(s).condition(c).event(e).pow = power;
                sites(s).condition(c).event(e).itpcbp = itpcbp;
                sites(s).condition(c).event(e).powbp = powerbp;
                sites(s).condition(c).event(e).max_itpcbp_time = max_itpcbp_time;
                sites(s).condition(c).event(e).max_powbp_time = max_powbp_time;
                sites(s).condition(c).event(e).max_itpcbp = max_itpcbp;
                sites(s).condition(c).event(e).max_powbp = max_powbp;
                sites(s).condition(c).event(e).lfp_sig = lfp_sig;
                sites(s).condition(c).event(e).itpc_sig = itpc_sig;
                sites(s).condition(c).event(e).pow_sig = pow_sig;                
                
                sites(s).condition(c).event(e).itpcbp_sig = itpcbp_sig;
                sites(s).condition(c).event(e).powbp_sig = powbp_sig;
                sites(s).condition(c).event(e).time = time;
                sites(s).condition(c).event(e).tfr_time = tfr_time;
            end
        end
    end
    
    %% this part to remove sites with not all conditions?
    % out_mask = cellfun(@isempty, {sites.site_ID});
    % sites = sites(~out_mask);
    sites = sites([sites.all_conditions_present]==1);
    save(fileName,'sites')
    
else
    load(fileName);
end


%% Computing the Target-wise averaging of total available sites
for tr = 1: length(targets)
    sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{tr})),sites);
    target_sites = sites(sites_for_this_target);
    condition=vertcat(target_sites.condition);
    %avg=struct();
    for c = 1:size(condition,2)
        events=cat(1,condition(:,c).event);
        
        concat.cond_name = cfg.condition(c).name;
        concat.nTriggers=round(mean(reshape([events.nTriggers],size(events)),1));
        concat.max_itpcbp=reshape([events.max_itpcbp],[size(cfg.lfp.frequency_bands,1),size(events)]);
        concat.max_powbp=reshape([events.max_powbp],[size(cfg.lfp.frequency_bands,1),size(events)]);
        concat.max_itpcbp_time=reshape([events.max_itpcbp_time],[size(cfg.lfp.frequency_bands,1),size(events)]);
        concat.max_powbp_time=reshape([events.max_powbp_time],[size(cfg.lfp.frequency_bands,1),size(events)]);
        
        concat.pow          = [];
        concat.itpc         = [];
        concat.lfp          = [];
        concat.itpcbp       = [];
        concat.powbp        = [];
        concat.powsig       = [];
        concat.itpcsig      = [];
        concat.powbpsig     = [];
        concat.itpcbpsig    = [];
        concat.lfpsig       = [];
        concat.tfr_time     = [];
        concat.lfp_time     = [];
        concat.lfp_std      = [];
        
        for e=1:size(cfg.analyse_states,1)
            pow_mean    = mean(cat(3,events(:,e).pow),3);
            itpc_mean   = mean(cat(3,events(:,e).itpc),3);
            lfp_mean    = mean(cat(3,events(:,e).lfp),3);
            itpcbp_mean = mean(cat(3,events(:,e).itpcbp),3);
            powbp_mean  = mean(cat(3,events(:,e).powbp),3);
            pow_sig     = mean(cat(3,events(:,e).pow_sig),3);
            itpc_sig    = mean(cat(3,events(:,e).itpc_sig),3);
            powbp_sig   = mean(cat(3,events(:,e).powbp_sig),3);
            itpcbp_sig  = mean(cat(3,events(:,e).itpcbp_sig),3);
            lfp_sig     = mean(cat(3,events(:,e).lfp_sig),3);
            tfr_time    = events(1,e).tfr_time;
            lfp_time    = events(1,e).time;
            lfp_std     = std(cat(3,events(:,e).lfp),0,3);
            
            NaNseparator=100/25;
            concat.pow          = cat(2, concat.pow,        pow_mean,   nan(size(pow_mean, 1),   NaNseparator));
            concat.itpc         = cat(2, concat.itpc,       itpc_mean,  nan(size(itpc_mean,1),   NaNseparator));
            concat.lfp          = cat(2, concat.lfp,        lfp_mean,	nan(size(lfp_mean, 1),   NaNseparator));
            concat.itpcbp       = cat(2, concat.itpcbp,     itpcbp_mean,nan(size(itpcbp_mean, 1),NaNseparator));
            concat.powbp        = cat(2, concat.powbp,      powbp_mean, nan(size(powbp_mean, 1), NaNseparator));
            concat.powsig       = cat(2, concat.powsig,     pow_sig,    nan(size(pow_sig, 1),    NaNseparator));
            concat.itpcsig      = cat(2, concat.itpcsig,    itpc_sig,   nan(size(itpc_sig, 1),   NaNseparator));
            concat.powbpsig     = cat(2, concat.powbpsig,   powbp_sig,  nan(size(powbp_sig, 1),  NaNseparator));
            concat.itpcbpsig    = cat(2, concat.itpcbpsig,  itpcbp_sig, nan(size(itpcbp_sig, 1), NaNseparator));
            concat.lfpsig       = cat(2, concat.lfpsig,     lfp_sig,    nan(size(lfp_sig, 1),    NaNseparator));
            concat.lfp_std      = cat(2, concat.lfp_std,    lfp_std,	nan(size(lfp_std, 1),   NaNseparator));
            concat.tfr_time     = [concat.tfr_time, tfr_time, nan(1, NaNseparator)];
            concat.lfp_time     = [concat.lfp_time, lfp_time, nan(1, NaNseparator)];
        end
        grand_avg(tr).avg(c)=concat;
    end
    grand_avg(tr).target = targets{tr};
    grand_avg(tr).nSites = length(target_sites);
end


%% plotting the results:
cond = {cfg.condition.name};
freq = cfg.lfp.foi;
frequency_bands=cfg.lfp.frequency_bands;
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked'};

% Smoothing Kernel here:
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

% figure 1 - overview
for tr = 1: length(targets)
    if grand_avg(tr).nSites == 0
        continue;
    end
    for c = 1:length(cond)
        % tfr_time not defined well!
        tfr_time=grand_avg(tr).avg(c).tfr_time;
        tfr_events.onset        = find(tfr_time == 0);
        tfr_events.name         = cfg.analyse_states(:,1);
        tfr_events.ticksamples  = sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
        tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
        tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
        tfr_events.ticks        = round(tfr_time(tfr_events.ticksamples)*100)/100;
        
        % create figure
        h(c) = figure('units','normalized','position',[0 0 1 1]);
        toplot={squeeze(grand_avg(tr).avg(c).pow),squeeze(grand_avg(tr).avg(c).itpc)};
        for sp=1:2
            % =========================== Power & ITPC ============================= %
            sph(c,sp)=subplot(3,2,sp);
            image(toplot{sp},'CDataMapping','scaled');
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
            collim{c,sp}=[min(nonnan(:)) max(nonnan(:))];
            set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5);
            add_ticks_and_labels(tfr_events,[0.5,numel(freq) + 0.5],8)
        end
                
        %========================== Bandpassed POWER ==================== %
        % Smoothing of the POWbp here:
        jnk = [];
        concat_input = cat(2,(grand_avg(tr).avg(c).powbp(:,half_win:-1:1)),(grand_avg(tr).avg(c).powbp(:,:)));
        concat_input = cat(2,concat_input, (grand_avg(tr).avg(c).powbp(:,end:-1:end-half_win+1)));
        for k=1:size(grand_avg(tr).avg(c).powbp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=3;
        sph(c,sp)=subplot(3,2,sp);
        hold on;
        set(gca,'ColorOrder',jet(size(grand_avg(tr).avg(c).powbp,1)));
        plot(squeeze(smoothed)')
        xlabel('Time(s)'); ylabel('Power (W)');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        add_ticks_and_labels(tfr_events,get(gca,'ylim'),8)
        
        %========================== Bandpassed ITPC ==================== %
        % Smoothing of the ITPCbp here:
        jnk = [];
        concat_input = cat(2,(grand_avg(tr).avg(c).itpcbp(:,half_win:-1:1)),(grand_avg(tr).avg(c).itpcbp(:,:)));
        concat_input = cat(2,concat_input, (grand_avg(tr).avg(c).itpcbp(:,end:-1:end-half_win+1)));
        for k=1:size(grand_avg(tr).avg(c).itpcbp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=4;
        sph(c,sp)=subplot(3,2,sp);
        hold on;
        set(gca,'ColorOrder',jet(size(grand_avg(tr).avg(c).itpcbp,1)));
        plot(squeeze(smoothed)')
        xlabel('Time(s)'); ylabel('ITPC');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        add_ticks_and_labels(tfr_events,get(gca,'ylim'),8)
        
        %========================== LFP evoked Potential ==================== %
        % Smoothing of the  LFP evoked Potential here:
        jnk = [];
        %lfp_se=sterr(grand_avg(tr).avg(c).lfp,3,0);
        lfp_std=grand_avg(tr).avg(c).lfp_std;
        %percentile25=prctile(grand_avg(tr).avg(c).lfp,25,3);
        %percentile75=prctile(grand_avg(tr).avg(c).lfp,75,3);
        concat_input = cat(2,(grand_avg(tr).avg(c).lfp(:,half_win:-1:1)),(grand_avg(tr).avg(c).lfp(:,:)));
        concat_input = cat(2,concat_input, (grand_avg(tr).avg(c).lfp(:,end:-1:end-half_win+1)));
        for k=1:size(grand_avg(tr).avg(c).lfp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_mean = jnk(:,half_win+1:end-half_win);
        concat_input = cat(2,(lfp_std(:,half_win:-1:1)),(lfp_std(:,:)));
        concat_input = cat(2,concat_input, (lfp_std(:,end:-1:end-half_win+1)));
        for k=1:size(grand_avg(tr).avg(c).lfp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_std = jnk(:,half_win+1:end-half_win);
        %         concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
        %         concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
        %         for k=1:size(grand_avg(tr).avg(c).lfp,1)
        %             jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        %         end
        %         clear concat_input
        %         smoothed_25 = jnk(:,half_win+1:end-half_win);
        %         concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
        %         concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
        %         for k=1:size(grand_avg(tr).avg(c).lfp,1)
        %             jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        %         end
        %         clear concat_input
        %smoothed_75 = jnk(:,half_win+1:end-half_win);
        
        sp=5;
        sph(c,sp)=subplot(3,2,sp);
        hold on;
        set(gca,'ColorOrder',jet(size(grand_avg(tr).avg(c).lfp,1)));
        lineProps={'color',[0 0 1]};
        %shadedErrorBar(time,smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);
        shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
        %plot(repmat(time,size(grand_avg(tr).avg(c).lfp_avg,1),1)', squeeze(smoothed_mean)')
        xlabel('Time(s)'); ylabel(' LFP evoked Potential');
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        add_ticks_and_labels(tfr_events,get(gca,'ylim'),8)
        
        
        %         set(gca,'Ylim',[-60 40]);% for Bacchus
        %         set(gca,'Ylim',[-6 6]);% for Magnus
        
        %=================================================================%
        %% format spectra colors
        cbtitle = {'(P - \mu) / std','P - \mu'};
        for sp=1:2
            subplot(3,2,sp);
            cm = colormap('jet');
            cb = colorbar;
            set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
            colormap(cm);
        end
        
        results_file{c} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','avg of ',num2str(grand_avg(tr).nSites),' sites ',withunits, '- ',grand_avg(tr).avg(c).cond_name]);
        mtit([ cfg.monkey,'-',targets{tr},'-avg of ',num2str(grand_avg(tr).nSites),' sites ' ,withunits, ' - ',grand_avg(tr).avg(c).cond_name],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
        
    end
    
    for c = 1:length(cond)
        figure(h(c));
        for sp=1:2
            subplot(sph(c,sp));
            set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
        end
        export_fig(h(c),[results_file{c},'.pdf']);
    end
end
close all,
clc

%%
freqb = {'delta 2-4 Hz','theta 4-8 Hz','alpha 8-14 Hz','beta 14-30 Hz','lowgamma 30-50 Hz','highgamma 70-150 Hz'};
freqName = {'delta','theta','alpha','beta','lowGamma','highGamma'};
color = jet(length(frequency_bands));

for tr = 1:length(targets)
    if grand_avg(tr).nSites == 0
        continue;
    end
    h = figure('Name',['Max ITPCbp/POWbp of Target=',grand_avg(tr).target],'NumberTitle','off');
    for c = 1:length(cond)
        for fb = 1: length(freqb)
            % itpcbp
            subplot(2,2,2*c-1)
            scatter(grand_avg(tr).avg(c).max_itpcbp_time(fb,:),grand_avg(tr).avg(c).max_itpcbp(fb,:),15,color(fb,:))
            hold on
            title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(grand_avg(tr).nSites),...
                ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('time (s)','Interpreter','latex'), ylabel('Max ITPC value','Interpreter','latex');
%             set(gca,'xlim',[-0.5,0.5],'ylim',[0,1])
            set(gca,'xlim',[-0.25,0.25],'ylim',[0,1])
            
            % powbp
            subplot(2,2,2*c)
            scatter(grand_avg(tr).avg(c).max_powbp_time(fb,:),grand_avg(tr).avg(c).max_powbp(fb,:),15,color(fb,:))
            hold on
            title([' max power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(grand_avg(tr).nSites),...
                ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('time (s)','Interpreter','latex'), ylabel('Max POWER value','Interpreter','latex');
%             set(gca,'xlim',[-0.5,0.5])
            set(gca,'xlim',[-0.25,0.25])
            
        end
        
        subplot(2,2,2*c-1)
        plot([0,0],get(gca,'ylim'),'k--')
        %
        subplot(2,2,2*c)
        plot([0,0],get(gca,'ylim'),'k--')
        legend(freqb,'FontSize',5)
        legend('boxoff')
        
    end
    mtit(['Target=',strrep(grand_avg(tr).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
    results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','Max_ITPCbp_POWbp_of ',num2str(grand_avg(tr).nSites),' sites ', withunits]);
    %     mtit([ ecg_bna_cfg.monkey,'-',targets{tr},'-avg of ',num2str(grand_avg(tr).nSites),' sites in all sessions - ',grand_avg(tr).avg(c).cond_name],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    export_fig(h,[results_file,'.pdf']);
%     saveas(h,[results_file]);
end
close all,
clc

%%
% ploting the avgogram of the timing of the Max ITPC/POW:
bins=cfg.analyse_states{1,3}:cfg.lfp.timestep*10:cfg.analyse_states{1,4};
for tr = 1: length(targets)
    if grand_avg(tr).nSites == 0
        continue;
    end
    h = figure;
    for c = 1:length(cond)
        for fb = 1: length(freqb)
            % itpcbp
            sp1=subplot(2,2,2*c-1);
            [nelements,centers]= hist(squeeze(grand_avg(tr).avg(c).max_itpcbp_time(fb,:)),bins);
            plot(centers,nelements,'-','Color',color(fb,:));
            hold on
            title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).max_itpcbp,2)),...
                ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Max ITPC times','Interpreter','latex');
            set(gca,'xlim',[-0.25,0.25])
            
            % powbp
            sp2=subplot(2,2,2*c);
            [nelements,centers]= hist(squeeze(grand_avg(tr).avg(c).max_powbp_time(fb,:)),bins);
            plot(centers,nelements,'-','Color',color(fb,:));
            hold on
            title([' max power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).max_powbp,2)),...
                ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Max POWER times','Interpreter','latex');
            set(gca,'xlim',[-0.25,0.25])
            
        end
        legend(freqb,'FontSize',5)
        legend('boxoff')
        
        plot(sp1,[0,0],get(sp1,'ylim'),'k--')
        plot(sp2,[0,0],get(sp2,'ylim'),'k--')
    end
    mtit(['Target=',strrep(grand_avg(tr).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
    results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','Time_of_Max_ITPCbp_POWbp_for ',num2str(grand_avg(tr).nSites),' sites ', withunits]);
    export_fig(h,[results_file,'.pdf']);
%     saveas(h,[results_file]);
end

close all,
clc


% ploting number of significant sites in each bin for ITPC/POW/evoked lfp:
bins=tfr_time;
for tr = 1: length(targets)
    if grand_avg(tr).nSites == 0
        continue;
    end
    h = figure;
    for c = 1:length(cond)
        for fb = 1: length(freqb)
            % itpcbp
            sp1=subplot(3,2,c);
            plot(bins,squeeze(grand_avg(tr).avg(c).itpcbp_sig_avg(fb,:)),'-','Color',color(fb,:));
            hold on
            title([' fraction significant itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).max_itpcbp,2)),...
                ' sites,',num2str(grand_avg(tr).avg(c).nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Significant ITPC times','Interpreter','latex');
            set(gca,'xlim',[-0.25,0.25])
%             xlim([bins(1) bins(end)]);
            
            % powbp
            sp2=subplot(3,2,c+2);
            plot(bins,squeeze(grand_avg(tr).avg(c).powbp_sig_avg(fb,:)),'-','Color',color(fb,:));
            hold on
            title([' fraction significant power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).max_powbp,2)),...
                ' sites,',num2str(grand_avg(tr).avg(c).nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Significant POWER times','Interpreter','latex');
            set(gca,'xlim',[-0.25,0.25])
%             xlim([bins(1) bins(end)]);            
        end
        
        legend(freqb,'FontSize',5)
        legend('boxoff')
        
        % evoked lfp
        sp3=subplot(3,2,c+4);
        plot(bins,grand_avg(tr).avg(c).lfp_sig_avg,'-','Color',[0 0 1]);
        hold on
        title([' fraction significant evoked lfp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).lfp_sig_avg,2)),...
            ' sites,',num2str(grand_avg(tr).avg(c).nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
        xlabel('Significant evoked LFP times','Interpreter','latex');
        set(gca,'xlim',[-0.25,0.25])
%         xlim([bins(1) bins(end)]);
        
        
        plot(sp1,[0,0],get(sp1,'ylim'),'k--')
        plot(sp2,[0,0],get(sp2,'ylim'),'k--')
    end
    mtit(['Target=',strrep(grand_avg(tr).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
    results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','Significant_bins_for ',num2str(grand_avg(tr).nSites),' sites ', withunits]);
    export_fig(h,[results_file,'.pdf']);
%     saveas(h,[results_file]);
end

close all,
clc


end


function add_ticks_and_labels(events,ylm,stp)
event_onsets    =events.onset;
event_names     =events.name;
event_samples   =events.ticksamples;
state_ticks     =events.ticks;

for so = event_onsets
    line([so so], ylm, 'color', 'k');
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
