function ecg_bna_rawLFP_butterfly_plots(cfg,allSitesData,sr)

session = cfg.session_info.Date;
path2save = cfg.results_folder;
time = 1/sr:1/sr:5;
nsamples = 5*sr;

clear legend_list
h = figure; 
hold on

for site = 1: length(allSitesData)
    sites = allSitesData(site).site;
    legend_list(site) = {['Channel-', num2str(sites.channel),'-',sites.target, '- depth =',num2str(sites.electrode_depth)]};
    plot(time,sites.LFP(1:nsamples))
end
hold off

legend(legend_list, 'Location','eastoutside','Interpreter','none','FontSize',6);
legend_list = legend_list';

set(gcf, 'WindowState', 'maximized');
figName = fullfile([path2save,filesep,'rawLFP channel Check',filesep,'session_',session,'_RawLFPsample']);
saveas(h, figName,'fig');
saveas(h, figName,'png');
%%
close all

end