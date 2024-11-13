function Task_Rest_raw_diff = ecg_bna_compute_TaskRest_raw_diff(cfg,withunits)

reprocess=1;
fileName = fullfile([cfg.grand_avg_fldr filesep cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Condition_diff_grand_avg_sessions_sites',withunits,'.mat']);

if cfg.combine_hemispheres
    cfg.targets=unique(cellfun(@(x) x(1:strfind(x,'_')-1),cfg.targets,'uniformoutput',false));
end
targets = cfg.targets;%unique({out.target}); %%cfg.monkey

if reprocess
    data_path = cfg.sites_lfp_fldr;
    cd(data_path)
    
    
    % read the files
    all_lfp_data = dir('*.mat');
    %%
    % storing all the Triggered parameters of the sites in the results folder:
    out = struct;
    i=0;
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
        i=i+1;
        target = triggered_site_data.target;
        out(i).all_conditions_present = 1;
        
        for cn = 1:length(triggered_site_data.condition)
            con=triggered_site_data.condition(cn);
            if isempty(con.event)
                out(i).all_conditions_present = 0;
                continue;
            end
            condition = con.label;
            nTriggers = con.event.real.ntriggers;
            time      = con.event.time;
            tfr_time  = con.event.tfr_time;
            lfp       = squeeze(con.event.real.lfp.mean)';
            itpc      = squeeze(con.event.real.itpc.mean);
            power     = squeeze(con.event.real.pow.mean);
            itpcbp    = squeeze(con.event.real.itpcbp.mean);
            powerbp   = squeeze(con.event.real.powbp.mean);
            
%             add_index=[0:size(itpcbp,1)-1]'*size(itpcbp,2);
%             [~,max_itpcbp_idx] = max(abs(itpcbp),[],2);
%             max_itpcbp_time=time(max_itpcbp_idx)';
%             tmp=itpcbp';
%             max_itpcbp=tmp(max_itpcbp_idx+add_index);
%             [~,max_powbp_idx] = max(abs(powerbp),[],2);
%             tmp=powerbp';
%             max_powbp=tmp(max_powbp_idx+add_index);
%             max_powbp_time=time(max_powbp_idx)';
%             
%             lfp_sig=squeeze(con.event.significance.lfp);
%             itpcbp_sig=squeeze(con.event.significance.itpcbp);
%             powerbp_sig=squeeze(con.event.significance.powbp);
%             
            out(i).site_ID = site_ID;
            out(i).target = target;
            out(i).condition(cn).condition_name = condition;
            out(i).condition(cn).nTriggers = nTriggers;
            out(i).condition(cn).lfp = lfp;
            out(i).condition(cn).itpc = itpc;
            out(i).condition(cn).pow = power;
            out(i).condition(cn).itpcbp = itpcbp;
            out(i).condition(cn).powbp = powerbp;
%             out(i).condition(cn).max_itpcbp_time = max_itpcbp_time;
%             out(i).condition(cn).max_powbp_time = max_powbp_time;
%             out(i).condition(cn).max_itpcbp = max_itpcbp;
%             out(i).condition(cn).max_powbp = max_powbp;
%             out(i).condition(cn).lfp_sig = lfp_sig;
%             out(i).condition(cn).itpcbp_sig = itpcbp_sig;
%             out(i).condition(cn).powerbp_sig = powerbp_sig;
            out(i).condition(cn).time = time;
            out(i).condition(cn).tfr_time = tfr_time;
        end
        
        if (out(i).all_conditions_present == 1)
            fields = fieldnames(out(i).condition(1));
            fields = fields(3:end-2);
            % Loop over each field name and apply minus function
            for fil = 1:length(fields)
                field = fields{fil};  % Get current field name
                out(i).condition_diff.(field) = minus(out(i).condition(2).(field), out(i).condition(1).(field));  % Apply minus
            end
        else
            out(i).condition_diff = [];
        end

    end
    
    %% this part to remove sites with not all conditions?
    % out_mask = cellfun(@isempty, {out.site_ID});
    % out = out(~out_mask);
    listOfExcluded_sites = {out(find([out.all_conditions_present]==0)).site_ID};
    save(fullfile([cfg.results_folder,filesep, 'listOfExcluded_sites_',cfg.monkey ]),"listOfExcluded_sites",'-v7.3');
    out = out([out.all_conditions_present]==1);
    
    %%
    % Computing the Target-wise averaging of total available sites
    Task_Rest_raw_diff = struct;
    for tr = 1: length(targets)
        sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{tr})),out);
        target_sites = out(sites_for_this_target);
        avg = struct('lfp',[],'itpc',[],'pow',[],'itpcbp',[],'powbp',[]);
        
        for trc = 1:length(target_sites)
            avg.lfp = cat(3,avg.lfp,target_sites(trc).condition_diff.lfp);
            avg.itpc = cat(3,avg.itpc,target_sites(trc).condition_diff.itpc);
            avg.pow = cat(3,avg.pow,target_sites(trc).condition_diff.pow);
            avg.itpcbp = cat(3,avg.itpcbp,target_sites(trc).condition_diff.itpcbp);
            avg.powbp = cat(3,avg.powbp,target_sites(trc).condition_diff.powbp);
        end
            
        avg.itpc_avg = mean(avg.itpc,3);
        avg.pow_avg = mean(avg.pow,3);
        avg.lfp_avg = mean(avg.lfp,3);
        avg.itpcbp_avg = mean(avg.itpcbp,3);
        avg.powbp_avg = mean(avg.powbp,3);

        Task_Rest_raw_diff(tr).target = targets{tr};
        Task_Rest_raw_diff(tr).nSites = length(target_sites);
        Task_Rest_raw_diff(tr).avg = avg;
        clear avg
        
    end
    save(fileName,'Task_Rest_raw_diff','tfr_time','time')
    
else
    load(fileName);
end


%%
% plotting the results:
cond = {[cfg.condition(2).name,'-',cfg.condition(1).name]};
freq = cfg.lfp.foi;
frequency_bands=cfg.lfp.frequency_bands;
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked'};

% Smoothing Kernel here:
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

%% tfr_time not defined well!
for tr = 1: length(targets)
    %
    if Task_Rest_raw_diff(tr).nSites == 0
        continue;
    end
    for cn = 1
        % create figure
        h = figure('units','normalized','position',[0 0 1 1]);
        toplot={Task_Rest_raw_diff(tr).avg.pow_avg,Task_Rest_raw_diff(tr).avg.itpc_avg};
        % =========================== Power ============================= %
        sp=1;
        sph(cn,sp)=subplot(3,2,sp);
        image(tfr_time, 1:numel(freq), squeeze(Task_Rest_raw_diff(tr).avg.pow_avg),'CDataMapping','scaled');
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
        image(tfr_time, 1:numel(freq), squeeze(Task_Rest_raw_diff(tr).avg.itpc_avg),'CDataMapping','scaled');
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
        concat_input = cat(2,(Task_Rest_raw_diff(tr).avg.powbp_avg(:,half_win:-1:1)),(Task_Rest_raw_diff(tr).avg.powbp_avg(:,:)));
        concat_input = cat(2,concat_input, (Task_Rest_raw_diff(tr).avg.powbp_avg(:,end:-1:end-half_win+1)));
        for k=1:size(Task_Rest_raw_diff(tr).avg.powbp_avg,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=3;
        sph(cn,sp)=subplot(3,2,sp);
        hold on;
        set(gca,'ColorOrder',jet(size(Task_Rest_raw_diff(tr).avg.powbp_avg,1)),'xlim',[time(1) time(end)]);
        plot(repmat(time,size(Task_Rest_raw_diff(tr).avg.powbp_avg,1),1)', squeeze(smoothed)')
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel('Power (\mu V^2)');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        set(gca,'Xlim',[-.25 .25]);
%         set(gca,'Ylim',[-3 3]);
        
        %========================== Bandpassed ITPC ==================== %
        % Smoothing of the ITPCbp here:
        jnk = [];
        concat_input = cat(2,(Task_Rest_raw_diff(tr).avg.itpcbp_avg(:,half_win:-1:1)),(Task_Rest_raw_diff(tr).avg.itpcbp_avg(:,:)));
        concat_input = cat(2,concat_input, (Task_Rest_raw_diff(tr).avg.itpcbp_avg(:,end:-1:end-half_win+1)));
        for k=1:size(Task_Rest_raw_diff(tr).avg.itpcbp_avg,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=4;
        sph(cn,sp)=subplot(3,2,sp);
        hold on;
        set(gca,'ColorOrder',jet(size(Task_Rest_raw_diff(tr).avg.itpcbp_avg,1)),'xlim',[time(1) time(end)]);
        plot(repmat(time,size(Task_Rest_raw_diff(tr).avg.itpcbp_avg,1),1)', squeeze(smoothed)')
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel('ITPC');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        set(gca,'Xlim',[-.25 .25]);
%         set(gca,'Ylim',[0 0.3]);
        
        %========================== LFP evoked Potential ==================== %
        % Smoothing of the  LFP evoked Potential here:
        jnk = [];
        lfp_se=sterr(Task_Rest_raw_diff(tr).avg.lfp,3,0);
        lfp_std=std(Task_Rest_raw_diff(tr).avg.lfp,0,3);
        percentile25=prctile(Task_Rest_raw_diff(tr).avg.lfp,25,3);
        percentile75=prctile(Task_Rest_raw_diff(tr).avg.lfp,75,3);
        concat_input = cat(2,(Task_Rest_raw_diff(tr).avg.lfp_avg(:,half_win:-1:1)),(Task_Rest_raw_diff(tr).avg.lfp_avg(:,:)));
        concat_input = cat(2,concat_input, (Task_Rest_raw_diff(tr).avg.lfp_avg(:,end:-1:end-half_win+1)));
        for k=1:size(Task_Rest_raw_diff(tr).avg.lfp_avg,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_mean = jnk(:,half_win+1:end-half_win);
        concat_input = cat(2,(lfp_std(:,half_win:-1:1)),(lfp_std(:,:)));
        concat_input = cat(2,concat_input, (lfp_std(:,end:-1:end-half_win+1)));
        for k=1:size(Task_Rest_raw_diff(tr).avg.lfp_avg,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_std = jnk(:,half_win+1:end-half_win);
        concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
        concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
        for k=1:size(Task_Rest_raw_diff(tr).avg.lfp_avg,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_25 = jnk(:,half_win+1:end-half_win);
        concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
        concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
        for k=1:size(Task_Rest_raw_diff(tr).avg.lfp_avg,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_75 = jnk(:,half_win+1:end-half_win);
        
        sp=5;
        sph(cn,sp)=subplot(3,2,sp);
        hold on;
        set(gca,'ColorOrder',jet(size(Task_Rest_raw_diff(tr).avg.lfp_avg,1)),'xlim',[time(1) time(end)]);
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
%         cbtitle = {'(P - \mu) / std','P - \mu'};
        cbtitle = {'raw','raw'};
        for sp=1:2
            subplot(3,2,sp);
            cm = colormap('jet');
            cb = colorbar;
            set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
            colormap(cm);
        end
        
        results_file{cn} = fullfile(cfg.grand_avg_fldr, [cfg.monkey,'-',targets{tr},' - ','avg of ',num2str(Task_Rest_raw_diff(tr).nSites),' sites ',withunits, '- ',cond{cn}]);
        mtit([ cfg.monkey,'-',targets{tr},'-avg of ',num2str(Task_Rest_raw_diff(tr).nSites),' sites ' ,withunits, ' - ',cond{cn}],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
        
    end
    
    
    for cn = 1:length(cond)
        figure(h);
        for sp=1:2
            subplot(sph(cn,sp));
            set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
        end
        export_fig(h,[results_file{cn},'_raw_diff','.pdf']);
%         saveas(h,results_file{cn});
    end
    
end
close all,
clc

% %%
% freqb = {'delta 2-4 Hz','theta 4-8 Hz','alpha 8-14 Hz','beta 14-30 Hz','lowgamma 30-50 Hz','highgamma 70-150 Hz'};
% freqName = {'delta','theta','alpha','beta','lowGamma','highGamma'};
% color = jet(length(frequency_bands));
% 
% for tr = 1:length(targets)
%     if Task_Rest_raw_diff(tr).nSites == 0
%         continue;
%     end
%     h = figure('Name',['Max ITPCbp/POWbp of Target=',Task_Rest_raw_diff(tr).target],'NumberTitle','off');
%     for cn = 1:length(cond)
%         for fb = 1: length(freqb)
%             % itpcbp
%             subplot(2,2,2*cn-1)
%             scatter(Task_Rest_raw_diff(tr).avg.max_itpcbp_time(fb,:),Task_Rest_raw_diff(tr).avg.max_itpcbp(fb,:),15,color(fb,:))
%             hold on
%             title([' max itpc-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(Task_Rest_raw_diff(tr).nSites),...
%                 ' sites,',num2str(Task_Rest_raw_diff(tr).avg.nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('time (s)','Interpreter','latex'), ylabel('Max ITPC value','Interpreter','latex');
% %             set(gca,'xlim',[-0.5,0.5],'ylim',[0,1])
%             set(gca,'xlim',[-0.25,0.25],'ylim',[0,1])
%             
%             % powbp
%             subplot(2,2,2*cn)
%             scatter(Task_Rest_raw_diff(tr).avg.max_powbp_time(fb,:),Task_Rest_raw_diff(tr).avg.max_powbp(fb,:),15,color(fb,:))
%             hold on
%             title([' max power-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(Task_Rest_raw_diff(tr).nSites),...
%                 ' sites,',num2str(Task_Rest_raw_diff(tr).avg.nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('time (s)','Interpreter','latex'), ylabel('Max POWER value','Interpreter','latex');
% %             set(gca,'xlim',[-0.5,0.5])
%             set(gca,'xlim',[-0.25,0.25])
%             
%         end
%         
%         subplot(2,2,2*cn-1)
%         plot([0,0],get(gca,'ylim'),'k--')
%         %
%         subplot(2,2,2*cn)
%         plot([0,0],get(gca,'ylim'),'k--')
%         legend(freqb,'FontSize',5)
%         legend('boxoff')
%         
%     end
%     mtit(['Target=',strrep(Task_Rest_raw_diff(tr).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
%     results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','Max_ITPCbp_POWbp_of ',num2str(Task_Rest_raw_diff(tr).nSites),' sites ', withunits]);
%     %     mtit([ ecg_bna_cfg.monkey,'-',targets{tr},'-avg of ',num2str(grand_avg(tr).nSites),' sites in all sessions - ',grand_avg(tr).avg.cond_name],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
%     export_fig(h,[results_file,'.pdf']);
% %     saveas(h,[results_file]);
% end
% close all,
% clc
% 
% %%
% % ploting the avgogram of the timing of the Max ITPC/POW:
% bins=cfg.analyse_states{1,3}:cfg.lfp.timestep*10:cfg.analyse_states{1,4};
% for tr = 1: length(targets)
%     if Task_Rest_raw_diff(tr).nSites == 0
%         continue;
%     end
%     h = figure;
%     for cn = 1:length(cond)
%         for fb = 1: length(freqb)
%             % itpcbp
%             sp1=subplot(2,2,2*cn-1);
%             [nelements,centers]= hist(squeeze(Task_Rest_raw_diff(tr).avg.max_itpcbp_time(fb,:)),bins);
%             plot(centers,nelements,'-','Color',color(fb,:));
%             hold on
%             title([' max itpc-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(size(Task_Rest_raw_diff(tr).avg.max_itpcbp,2)),...
%                 ' sites,',num2str(Task_Rest_raw_diff(tr).avg.nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Max ITPC times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
%             
%             %             % powbp
%             sp2=subplot(2,2,2*cn);
%             [nelements,centers]= hist(squeeze(Task_Rest_raw_diff(tr).avg.max_powbp_time(fb,:)),bins);
%             plot(centers,nelements,'-','Color',color(fb,:));
%             hold on
%             title([' max power-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(size(Task_Rest_raw_diff(tr).avg.max_powbp,2)),...
%                 ' sites,',num2str(Task_Rest_raw_diff(tr).avg.nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Max POWER times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
%             
%         end
%         legend(freqb,'FontSize',5)
%         legend('boxoff')
%         
%         plot(sp1,[0,0],get(sp1,'ylim'),'k--')
%         plot(sp2,[0,0],get(sp2,'ylim'),'k--')
%     end
%     mtit(['Target=',strrep(Task_Rest_raw_diff(tr).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
%     results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','Time_of_Max_ITPCbp_POWbp_for ',num2str(Task_Rest_raw_diff(tr).nSites),' sites ', withunits]);
%     export_fig(h,[results_file,'.pdf']);
% %     saveas(h,[results_file]);
% end
% 
% close all,
% clc
% 
% 
% % ploting number of significant sites in each bin for ITPC/POW/evoked lfp:
% bins=tfr_time;
% for tr = 1: length(targets)
%     if Task_Rest_raw_diff(tr).nSites == 0
%         continue;
%     end
%     h = figure;
%     for cn = 1:length(cond)
%         for fb = 1: length(freqb)
%             % itpcbp
%             sp1=subplot(3,2,cn);
%             plot(bins,squeeze(Task_Rest_raw_diff(tr).avg.itpcbp_sig_avg(fb,:)),'-','Color',color(fb,:));
%             hold on
%             title([' fraction significant itpc-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(size(Task_Rest_raw_diff(tr).avg.max_itpcbp,2)),...
%                 ' sites,',num2str(Task_Rest_raw_diff(tr).avg.nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Significant ITPC times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
% %             xlim([bins(1) bins(end)]);
%             
%             % powbp
%             sp2=subplot(3,2,cn+2);
%             plot(bins,squeeze(Task_Rest_raw_diff(tr).avg.powbp_sig_avg(fb,:)),'-','Color',color(fb,:));
%             hold on
%             title([' fraction significant power-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(size(Task_Rest_raw_diff(tr).avg.max_powbp,2)),...
%                 ' sites,',num2str(Task_Rest_raw_diff(tr).avg.nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Significant POWER times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
% %             xlim([bins(1) bins(end)]);
%             
%             
%             
%             
%         end
%         
%         legend(freqb,'FontSize',5)
%         legend('boxoff')
%         
%         % evoked lfp
%         sp3=subplot(3,2,cn+4);
%         plot(bins,Task_Rest_raw_diff(tr).avg.lfp_sig_avg,'-','Color',[0 0 1]);
%         hold on
%         title([' fraction significant evoked lfp in ', strrep(cond{cn},'_',' '),' for ',num2str(size(Task_Rest_raw_diff(tr).avg.lfp_sig_avg,2)),...
%             ' sites,',num2str(Task_Rest_raw_diff(tr).avg.nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%         xlabel('Significant evoked LFP times','Interpreter','latex');
%         set(gca,'xlim',[-0.25,0.25])
% %         xlim([bins(1) bins(end)]);
%         
%         
%         plot(sp1,[0,0],get(sp1,'ylim'),'k--')
%         plot(sp2,[0,0],get(sp2,'ylim'),'k--')
%     end
%     mtit(['Target=',strrep(Task_Rest_raw_diff(tr).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
%     results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','Significant_bins_for ',num2str(Task_Rest_raw_diff(tr).nSites),' sites ', withunits]);
%     export_fig(h,[results_file,'.pdf']);
% %     saveas(h,[results_file]);
% end
% 
% close all,
% clc


end