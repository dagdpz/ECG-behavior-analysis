function [ site_lfp ] = ecg_bna_compute_site_lfp_tfr( site_lfp, cfg )
% lfp_tfa_compute_site_tfr - Computes the LFP time frequency spectrogram
% for all trials of a site. This function calls the ft_freqanalysis routine
% of feldtrip toolbox for the calculation of LFP power spectrogram
%
% USAGE:
%	site_lfp = lfp_tfa_compute_site_tfr( site_lfp, cfg )
%
% INPUTS:
%		site_lfp        - struct containing LFP signal for each trial for a
%		single site
%       cfg      - settings for TFR computation, see
%       settings/lfp_tfa_settings_example
%           Required felds:
%               1. tfr.method   - method to be used for calculating the LFP
%               power spectra. Can be 'mtmconvol' or 'wavelet'
%               2. tfr.width    - width of the wavelets in number of
%               cycles (For method = 'wavelet', Ignored for 'mtmconvol')
%               3. tfr.taper    - taper (single or multiple) to be used.
%               Can be 'dpss', 'hanning' or many others (Used when
%               tfr.method = 'mtmconvol', ignored for 'wavelet'
%               4. tfr.foi      - frequencies of interest (in Hz), a vector
%               of freq values or a 1x2 vector with start and end freq
%               5. tfr.t_ftimwin - length of the sliding time-window in
%               seconds, should be vector of length 1 x numfoi (Only used
%               when tfr.method = 'mtmconvol')
%               6. tfr.timestep - number of lfp samples to step for the
%               sliding time window
% OUTPUTS:
%		site_lfp         - structure containing trial-wise LFP power
%		spectrograms
%
% REQUIRES: ft_freqanalysis
%
% See also lfp_tfa_process_lfp, settings/lfp_tfa_settings_example,
% ft_freqanalysis
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	frst Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfle header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

%loop through each run (actually good choice not to loop through runs)
for r = (unique([site_lfp.trials.run]))
    fprintf('Computing TFR for run %g\n-----------------------\n', r);
    % concatenate all trials for this run
    trials_idx = find([site_lfp.trials.run] == r);
    
    concat_raw = double([site_lfp.trials(trials_idx).lfp_data]);
    concat_TFR = single(resample(concat_raw,1,cfg.tfr.timestep));   % this one is resampled! --> think of a better way
    
    ts_original = site_lfp.trials(trials_idx(1)).tsample;           % original time step
    ts = ts_original*cfg.tfr.timestep;                              % resampled time step

        
    N_cycles=cfg.tfr.n_cycles;
    frequencies = cfg.tfr.foi; 
    morlet_borders=1/min(frequencies)*N_cycles/2;
    time = -morlet_borders:1*ts:morlet_borders;
    s = N_cycles./(2*pi*frequencies);
    
    % fT parameters (use next-power-of-2)
    n_wavelet     = length(time); 
    n_data        = size(concat_TFR,2);
    n_convolution = n_wavelet+n_data;
    n_conv_pow2   = pow2(nextpow2(n_convolution));
    half_wavelet_len = ceil(length(time)/2);
        
    % get ft of data
    dataft = fft(concat_TFR,n_conv_pow2);
    
    for f=1:length(frequencies)
        % create wavelet
        wavelet = (1/(s(f)*sqrt(2*pi))) * exp(2*1i*pi*frequencies(f).*time) .* exp(-time.^2./(2*(s(f)^2)));
        % convolution
        datconv = ifft(fft(wavelet,n_conv_pow2).*dataft);
        datconv = datconv(1:n_convolution);
        datconv = datconv(half_wavelet_len:end-half_wavelet_len);
        
        % now reshape to the original number of trials:
        n_sample_start=1;
        remainder=0;
        for t = trials_idx
            %% for some reason, last trial is 2 samples short (thats 50 ms!)
            n_smpls=numel(site_lfp.trials(t).lfp_data)/cfg.tfr.timestep+remainder;
            n_samples=round(n_smpls);
            remainder=n_smpls-n_samples;
            t_start_sample=ceil(cfg.tfr.timestep*(1/2-remainder));
            n_sample_end =n_sample_start + n_samples-1;
            % extract phase values of reshaped data:
            site_lfp.trials(t).tfs.phase(1,f,:)= angle(datconv(n_sample_start:n_sample_end));
            % extracted Power of each trial
            site_lfp.trials(t).tfs.powspctrm(1,f,:)= abs(datconv(n_sample_start:n_sample_end));%.^2;
            % time ( :(, only cause needed in this format in other bits)
            site_lfp.trials(t).tfs.time= site_lfp.trials(t).time(t_start_sample:cfg.tfr.timestep:end);
            site_lfp.trials(t).tfs.freq=frequencies;
            
            n_sample_start = n_sample_end +1;
        end
    end
    
    
    frequency_bands=cfg.tfr.frequency_bands;
    for f=1:size(frequency_bands,1)

%% think about how to use good filters without causing errors for short periods
%         fltered_data = eegfilt(concat_LFP, round(1/ts),frequency_bands(f,1), []);
%         fltered_data = eegfilt(fltered_data, round(1/ts), [], frequency_bands(f,2));
%         
        [b, a]=butter(3, frequency_bands(f,1)*ts_original); % low-pass filter
        fltered_data = filtfilt(b,a,concat_raw);
        [b, a]=butter(3, frequency_bands(f,2)*ts_original,'high'); % hgh-pass filter
        fltered_data = filtfilt(b,a,fltered_data);
        
        phase_data = angle(hilbert(fltered_data));
        n_sample_start=1;
        
        for t = trials_idx
            n_sample_end = n_sample_start + numel(site_lfp.trials(t).lfp_data)-1;
            % extract phase values of reshaped data:
            site_lfp.trials(t).phase_bandpassed(1,f,:) = phase_data(n_sample_start:n_sample_end);
            n_sample_start = n_sample_end +1;
        end
    end
end

end
