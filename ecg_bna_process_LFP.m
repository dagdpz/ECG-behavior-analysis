function site_lfp = ecg_bna_process_LFP(sites,cfg,ts_original)

% lfp_tfa_process_LFP - function to read in the trial-wise LFP data for all
% sites recorded in a session, compute the LFP time frequency spectrogram,
% detect the noisy trials, and compute site-wise baseline pow
%
% USAGE:
%	session_info = ecg_bna_process_LFP(sites,cfg,ts_original)
%
% INPUTS:
%       session_info        - structure containing information about the
%       session to be processed, see settings/lfp_tfa_settings_example
%       Required fields:
%           Input               - path to the file which contains the
%                               trial-wise LFP data for all sites recorded
%                               in a session
%           proc_results_folder - folder where the results has to be
%           stored, see lfp_tfa_define_settings
%
% OUTPUTS:
%		session_info            - same as input structure session_info
%
% REQUIRES: lfp_tfa_compute_site_tfr, lfp_tfa_reject_noisy_lfp_trials,
% lfp_tfa_compute_site_baseline, lfp_tfa_compute_sitepair_csd
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings,
% lfp_tfa_global_states, lfp_tfa_reject_noisy_lfp_trials,
% lfp_tfa_compute_site_baseline

% struct to save data for a site
site_lfp = rmfield(sites,'LFP');
fprintf('=============================================================\n');
fprintf('Processing site, %s\n', sites.site_ID);
site_lfp.session = sites.site_ID(1:12);
site_lfp.recorded_hemisphere = upper(sites.target(end));


N_cycles=cfg.lfp.n_cycles;
frequencies = cfg.lfp.foi;
frequency_bands=cfg.lfp.frequency_bands;
morlet_borders=1/min(frequencies)*N_cycles/2;
s = N_cycles./(2*pi*frequencies);

ts=round(cfg.lfp.timestep/ts_original);

% fT parameters (use next-pow-of-2)
time = -morlet_borders:1*ts_original:morlet_borders;
n_wavelet     = length(time);


site_lfp.tfs.resampling_factor=1/ts;
site_lfp.tfs.sr=1/ts_original/ts;
sizepreallocator=[size(frequency_bands,1),floor(numel(sites.LFP)/ts)];
site_lfp.tfs.phabp=NaN(sizepreallocator);
site_lfp.tfs.powbp=NaN(sizepreallocator);
sizepreallocator=[numel(frequencies),floor(numel(sites.LFP)/ts)];
site_lfp.tfs.pha=NaN(sizepreallocator);
site_lfp.tfs.pow=NaN(sizepreallocator);
site_lfp.tfs.lfp=NaN([1,floor(numel(sites.LFP)/ts)]);

trials_block=[sites.block];
blocks_with_LFP=unique(trials_block);
samples_past=0;
samples_past_resampled=0;


for b=1:numel(blocks_with_LFP) 
    B=blocks_with_LFP(b);
    b_samples=trials_block==B;
    bs=samples_past_resampled+1;
    be=floor(sum(sites.LFP_samples(b_samples))/ts)+samples_past_resampled;  
    bs_original=samples_past+1;
    be_original=sum(sites.LFP_samples(b_samples))+samples_past;  
    samples_past=be_original;
    samples_past_resampled=be;    
    
    concat_raw = double(sites.LFP(bs_original:be_original));
    
    [concat_raw, noisy_trials_lfp_mean , noisy_trials_lfp_zscore] = ecg_bna_noisy_LFP_detection(concat_raw);
    
    
    n_data        = size(concat_raw,2);
    n_convolution = n_wavelet+n_data;
    n_conv_pow2   = pow2(nextpow2(n_convolution));
    half_wavelet_len = ceil(length(time)/2);
    
    % get fT of data
    dataft = fft(concat_raw,n_conv_pow2);
    site_lfp.tfs.freq             = frequencies;
    
    site_lfp.tfs.n_samples_per_block(:,b)=[B,be-bs+1];
    dat = nanmean(reshape(concat_raw(1:end-mod(size(concat_raw,2), ts)),ts,[]),1);
    site_lfp.tfs.lfp(1,bs:be)= dat;
    
    
    for f=1:length(frequencies)
        % create wavelet        
        wavelet = (ts_original/(s(f)*sqrt(2*pi))) * exp(2*1i*pi*frequencies(f).*time) .* exp(-time.^2./(2*(s(f)^2)));   
        % convolution
        datconv = ifft(fft(wavelet,n_conv_pow2).*dataft);
        
        dat = datconv(1:n_convolution);
        dat = dat(half_wavelet_len+1:end-half_wavelet_len); % +1 is for one overlap point at the start
        % resample here:
        dat = nanmean(reshape(dat(1:end-mod(size(dat,2), ts)),ts,[]),1);
        
        % extract pha values of reshaped data:
        site_lfp.tfs.pha(f,bs:be)= dat./abs(dat);
        % extracted Power of each trial
        site_lfp.tfs.pow(f,bs:be) = abs(dat).^2;        
    end
    
    for f=1:size(frequency_bands,1)
        %% think about how to use good filters without causing errors for short periods
        %  fltered_data = eegfilt(concat_LFP, round(1/ts),frequency_bands(f,1), []);
        %  fltered_data = eegfilt(fltered_data, round(1/ts), [], frequency_bands(f,2));
        [u, v]=butter(3, 2*frequency_bands(f,:)*ts_original); % band-pass filter
        dat = filtfilt(u,v,concat_raw);
        H=hilbert(dat);
        H = mean(reshape(H(1:end-mod(size(H,2), ts)),ts,[]),1);
        absH=abs(H);
        site_lfp.tfs.phabp(f,bs:be) = H./absH;
        site_lfp.tfs.powbp(f,bs:be) = absH.^2;               
    end
end
site_lfp.tfs.pha=site_lfp.tfs.pha(:,1:be);
site_lfp.tfs.pow=site_lfp.tfs.pow(:,1:be);
site_lfp.tfs.phabp=site_lfp.tfs.phabp(:,1:be);
site_lfp.tfs.powbp=site_lfp.tfs.powbp(:,1:be);
site_lfp.tfs.lfp=site_lfp.tfs.lfp(:,1:be);

% Noise rejection - is this even still feasable?
% tic
% site_lfp = ecg_bna_reject_noisy_lfp_trials( site_lfp, lfp_tfa_cfg.noise );
% toc
end

