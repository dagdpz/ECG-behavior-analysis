function Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(session_info,cfg)
%% load seed/reset RGN seed
load(session_info.Input_ECG);

% get the number of permutations from settings
global N
N=cfg.n_permutations;

% % condition for using CAP data instead of ECG 
% if (ecg_bna_cfg.outNameCap)  
%     oldnames = fieldnames(out_cap);
%     newnames = strrep(oldnames,'B2B','R2R');
%     idx = find(~strcmpi(newnames,oldnames));
%     for n = 1: numel(newnames)
%         [out_cap.(newnames{n})] = deal(out_cap.(oldnames{n}));
%     end
%     for n = 1:numel(idx)
%         out_cap = rmfield(out_cap,oldnames{idx(n)});
%     end
%     out = out_cap;
%     nrblockscell=num2cell(ses.nrblock_combinedFiles);
%     [out.nrblock_combinedFiles] = deal(nrblockscell{:});
% end

offset_blocks_Rpeak=0;
for b=1:numel(out)
    Rpeaks(b).block=out(b).nrblock_combinedFiles;
    Rpeaks(b).offset=offset_blocks_Rpeak(b);
    if isempty(out(b).nrblock_combinedFiles) || isempty(out(b).Rpeak_t) || isempty(out(b).R2R_t)
        Rpeaks(b).block=NaN;
        Rpeaks(b).RPEAK_ts=[];
        Rpeaks(b).RR_durations = [];
        Rpeaks(b).shuffled_ts=[];
        offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b);
        continue
    end
    
    [RPEAK_ts, RPEAK_dur, RPEAK_ts_p, RPEAK_dur_p, allowed_jitter_range] = jittering_core_function(out(b).Rpeak_t, out(b).R2R_valid, out(b).R2R_t, out(b).idx_valid_R2R_consec);
    
    %% put the data together
    Rpeaks(b).RPEAK_ts=RPEAK_ts+offset_blocks_Rpeak(b);         % this offset is just a trick to be able to append Rpeaks across blocks easily
    Rpeaks(b).RPEAK_dur=RPEAK_dur; % durations of RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).RPEAK_ts)
    Rpeaks(b).shuffled_ts=RPEAK_ts_p+offset_blocks_Rpeak(b);
    Rpeaks(b).shuffled_ts(Rpeaks(b).shuffled_ts==offset_blocks_Rpeak(b))=0;
    
    Rpeaks(b).shuffled_dur = RPEAK_dur_p; % durations of reshuffled RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).shuffled_ts)
    offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b)+max(RPEAK_ts)+allowed_jitter_range*2;
    
    clear RPEAK_ts RPEAK_dur RPEAK_ts_p RPEAK_dur_p allowed_jitter_range
    
    %% jitter data for exhalation and inhalation
    if isfield(out, 'is_RR_insp') && sum(out(b).is_RR_insp)>1 && sum(out(b).is_RR_exp)>1
        % inhalation first
        RPEAK_ts             = out(b).Rpeak_t(out(b).is_R_peak_insp);
        R2R_valid            = out(b).R2R_valid(out(b).is_RR_insp);
        R2R_t                = out(b).R2R_t(out(b).is_RR_insp);
        
        % prepare idx_valid_R2R_consec
        R2R_consec     = out(b).R2R_t(out(b).idx_valid_R2R_consec);
        R2R_intersect  = intersect(R2R_consec,R2R_t);
        idx_valid_R2R_consec = find(ismember(R2R_intersect, R2R_t));
        if idx_valid_R2R_consec(1) == 1
            % drop the first Rpeak if it's consecutive as the jittering 
            % function complains about such a situation
            idx_valid_R2R_consec(1) = [];
        end
        
        [RPEAK_ts, RPEAK_dur, RPEAK_ts_p, RPEAK_dur_p] = ...
            jittering_core_function(RPEAK_ts, R2R_valid, R2R_t, idx_valid_R2R_consec);
        
        Rpeaks(b).RPEAK_ts_insp=RPEAK_ts+offset_blocks_Rpeak(b);         % this offset is just a trick to be able to append Rpeaks across blocks easily
        Rpeaks(b).RPEAK_dur_insp=RPEAK_dur; % durations of RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).RPEAK_ts)
        Rpeaks(b).shuffled_ts_insp=RPEAK_ts_p+offset_blocks_Rpeak(b);
        Rpeaks(b).shuffled_ts_insp(Rpeaks(b).shuffled_ts_insp==offset_blocks_Rpeak(b))=0;
    
        Rpeaks(b).shuffled_dur_insp = RPEAK_dur_p; % durations of reshuffled RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).shuffled_ts)
        
        clear RPEAK_ts RPEAK_dur RPEAK_ts_p RPEAK_dur_p allowed_jitter_range
        
        % exhalation second
        RPEAK_ts             = out(b).Rpeak_t(out(b).is_R_peak_exp);
        R2R_valid            = out(b).R2R_valid(out(b).is_RR_exp);
        R2R_t                = out(b).R2R_t(out(b).is_RR_exp);
        
        % prepare idx_valid_R2R_consec
        R2R_consec     = out(b).R2R_t(out(b).idx_valid_R2R_consec);
        R2R_intersect  = intersect(R2R_consec,R2R_t);
        idx_valid_R2R_consec = find(ismember(R2R_intersect, R2R_t));
        if idx_valid_R2R_consec(1) == 1
            % drop the first Rpeak if it's consecutive as the jittering 
            % function complains about such a situation
            idx_valid_R2R_consec(1) = [];
        end
        
        [RPEAK_ts, RPEAK_dur, RPEAK_ts_p, RPEAK_dur_p] = ...
            jittering_core_function(RPEAK_ts, R2R_valid, R2R_t, idx_valid_R2R_consec);
        
        Rpeaks(b).RPEAK_ts_exp=RPEAK_ts+offset_blocks_Rpeak(b);         % this offset is just a trick to be able to append Rpeaks across blocks easily
        Rpeaks(b).RPEAK_dur_exp=RPEAK_dur; % durations of RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).RPEAK_ts)
        Rpeaks(b).shuffled_ts_exp=RPEAK_ts_p+offset_blocks_Rpeak(b);
        Rpeaks(b).shuffled_ts_exp(Rpeaks(b).shuffled_ts_exp==offset_blocks_Rpeak(b))=0;
    
        Rpeaks(b).shuffled_dur_exp = RPEAK_dur_p; % durations of reshuffled RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).shuffled_ts)
        
        clear RPEAK_ts RPEAK_dur RPEAK_ts_p RPEAK_dur_p allowed_jitter_range
        
    end
    
    %% This figure plots all R-peaks and consecutive R-peaks to make sure our selection procedure does the right thing
%     figure
%     hold on
%     for this_index_wont_be_used = 1:length(iv_starts)
%         f= fill([iv_starts(this_index_wont_be_used) iv_starts(this_index_wont_be_used) ...
%             iv_ends(this_index_wont_be_used) iv_ends(this_index_wont_be_used) ...
%             iv_starts(this_index_wont_be_used)], ...
%             [0 1 1 0 0], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
%     end
%     st1 = stem(out(b).Rpeak_t, ones(length(out(b).Rpeak_t),1)); % all the R-peaks from out.Rpeak_t
%     st2 = stem(out(b).R2R_t, 0.75*ones(length(out(b).R2R_t),1));
%     st3 = stem(RPEAK_ts, 0.5*ones(length(RPEAK_ts), 1),'g'); % plot consecutive R-peaks
%     legend([f st1 st2 st3], {'Invalid Time Window', 'All R-peaks', 'Ends of valid RRs', 'Consecutive R-peaks: 1 valid RR before & 1 valid RR after'})
%     xlabel('Time, s')
end

end

function [RPEAK_ts, RPEAK_dur, RPEAK_ts_p, RPEAK_dur_p, allowed_jitter_range] = jittering_core_function(RPEAK_ts, R2R_valid, R2R_t, idx_valid_R2R_consec)
global N
% arguments: Rpeak_t, R2R_valid, R2R_t, idx_valid_R2R_consec
RPEAKS_intervals=diff(RPEAK_ts);                               % These are the intervals of our valid Rpeaks, the ones we are going to jitter

%% shuffling the VALID R2R intervals to append at the end of our jittered intervals,
%  just in case we randomly end up not covering the entire time window
%  Don't worry, most (typically all) of it is going to be removed
[~,ix]=sort(rand(N,length(R2R_valid)),2);
perm=R2R_valid(ix);

%% here we jitter all Rpeak intervals
allowed_jitter_range = std(R2R_valid);                    % one-sided std of normally distributed jitter values
RPEAKS_intervals_p   = [repmat(RPEAK_ts(1),N,1) repmat(RPEAKS_intervals,N,1)+randn(N,length(RPEAKS_intervals))*allowed_jitter_range perm]; % jittering every interval!
RPEAK_ts_p           = cumsum(RPEAKS_intervals_p,2);
RPEAK_dur_p          = [zeros(N,1) diff(RPEAK_ts_p,1,2)];

%% figure out consecutive RR-intervals - we do it only after jittering because we need to know durations of jittered intervals corresponding to the consecutive ones
[~,consecutive_idx]=ismember(R2R_t(idx_valid_R2R_consec),RPEAK_ts); % indexing of consecutive_idx corrsponds to out(b).R2R_t
valid_idx=consecutive_idx-1;
next_invalid=diff(valid_idx)~=1;                               % is the next Rpeak invalid (i.e. followed by invalid R2R interval)
iv_starts  =[0  RPEAK_ts(valid_idx([next_invalid true]))];     % start of invalid intervals: Timestamps of valid Rpeaks followed by invalid ones
% First Segment (for 0 to first valid Rpeak) and last segment
% (everything after last valid Rpeak) are always invalid
iv_ends    =[RPEAK_ts(valid_idx([true next_invalid]))   inf];  % end of invalid intervals: Timestamps of valid Rpeaks PRECEDED by invalid ones
% we can get Rpeak-ts preceded by invalid R2R by shifting next_invalid
grace_window=mean(R2R_valid)/2;                         % +/- Range for shuffled Rpeaks to be allowed inside invalid segments

%% take data corresponding to consecutive R-peaks
RPEAK_ts     = RPEAK_ts(valid_idx);                                  % take only Rpeaks surrounded by valid R2R
RPEAK_dur    = R2R_valid(idx_valid_R2R_consec-1);      %

%% remove jittered Rpeaks and corresponding durations that fell into invalid segments
for iv=1:numel(iv_starts)
    idx2exclude_1 = RPEAK_ts_p>iv_starts(iv)+grace_window & RPEAK_ts_p<iv_ends(iv)-grace_window;
    RPEAK_ts_p   (idx2exclude_1) = 0;
    RPEAK_dur_p (idx2exclude_1) = NaN;
end
idx2exclude_2 = RPEAK_ts_p>max(RPEAK_ts)+allowed_jitter_range;
RPEAK_ts_p    (idx2exclude_2) = 0;
RPEAK_dur_p  (idx2exclude_2) = NaN;

RPEAK_ts_p   (:,all(RPEAK_ts_p==0,1)) = [];
RPEAK_dur_p (:,all(isnan(RPEAK_dur_p),1)) = [];

end