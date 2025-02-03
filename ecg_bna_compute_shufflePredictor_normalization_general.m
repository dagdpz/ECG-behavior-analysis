function [normalized] = ecg_bna_compute_shufflePredictor_normalization_general(real,shuffled,ecg_bna_cfg)
% ecg_bna_compute_shufflePredictor_normalization_general - normalizing the real tfs
% and evoked data based on the shuffle predictor results
%
% USAGE:
%	[normalized] = ecg_bna_compute_shufflePredictor_normalization(variable,ecg_bna_cfg)
%
% INPUTS:
%		real            - 1xN struct containing the real data
%		shuffled  	    - 1xN struct containing the shuffled data
%       ecg_bna_cfg     - struct containing configuration settings
%
% OUTPUTS: 
%       normalized        - 1xN struct containing the normalized results of
%       pow, itpc, lfp, or itpcbp data
% ======================================================================= %


method = ecg_bna_cfg.lfp.normalization;
parameters={'pow','itpc','lfp','itpcbp','powbp','pha'};
for p=1:numel(parameters)
    parameter=parameters{p};
    realmean=real.(parameter).mean;
    realstd=real.(parameter).std;
    shuffledmean=shuffled.(parameter).mean;
    shuffledstd=shuffled.(parameter).std;
    if strcmp(method , 'subtraction')
        normalized.(parameter).mean    = realmean-shuffledmean;
        normalized.(parameter).std    = shuffledstd;  %%??
    elseif strcmp(method , 'division')
        normalized.(parameter).mean    = realmean./shuffledmean;
        normalized.(parameter).std    = realstd;%./shuffledmean;  %%??
    elseif strcmp(method , 'zscore')
        normalized.(parameter).mean    = (realmean-shuffledmean)./shuffledstd;
        normalized.(parameter).std    = realstd;%./shuffledstd; %% ??
    elseif strcmp(method , 'not normalized')
        normalized.(parameter).mean    = realmean;
        normalized.(parameter).std    = realstd;
    end
end

