function ecg_bna_plot_venns(cfg)

ccs = load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG\unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_any.mat']);
CBE = load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG\unitInfo_after_SNR_exclusion_stable_noCB_corr.mat']);

both_sets = intersect(CBE.unit_ids, ccs.unit_ids);
set_CBE = ~ismember(CBE.unit_ids,ccs.unit_ids);
set_ccs = ~ismember(ccs.unit_ids,CBE.unit_ids);

figure,
venn([length(CBE.unit_ids) length(ccs.unit_ids)],[length(both_sets), sum(set_CBE), sum(set_ccs)])
legend('cosine fit based','corr. coef. based')
axis off

unqAreas = {'VPL', 'dPul', 'MD'};

figure,

for a = 1:length(unqAreas)
    
    currAreaIds_CBE = cellfun(@(x) strcmp(x,unqAreas{a}), CBE.targets, 'UniformOutput',true);
    currAreaIds_ccs = cellfun(@(x) strcmp(x,unqAreas{a}), ccs.targets, 'UniformOutput',true);
    currArea_both = intersect(CBE.unit_ids(currAreaIds_CBE), ccs.unit_ids(currAreaIds_ccs));
    currArea_CBE  = ~ismember(CBE.unit_ids(currAreaIds_CBE),ccs.unit_ids(currAreaIds_ccs));
    currArea_ccs  = ~ismember(ccs.unit_ids(currAreaIds_ccs),CBE.unit_ids(currAreaIds_CBE));
    
    subplot(1,3,a)
    venn([length(CBE.unit_ids(currAreaIds_CBE)) length(ccs.unit_ids(currAreaIds_ccs))],[length(currArea_both), sum(currAreaIds_CBE), sum(currAreaIds_ccs)])
    axis off
    title(unqAreas{a})
    if a == 1
        legend('cosine fit based','corr. coef. based')
    end
end
end

