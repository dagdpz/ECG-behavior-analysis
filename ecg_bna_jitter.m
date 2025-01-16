function Triggers=ecg_bna_jitter(IN,cfg)

N=cfg.n_permutations;

ts=IN.ts;
intervals=IN.intervals;
valid_idx=IN.valid_ts;
blocks=IN.blocks;
td=diff(ts); % These are the ACTUAL distances of our Rpeaks

%% shuffling the VALID R2R intervals to append at the end of our jittered intervals,
switch cfg.jitter_method
    case 'interval_jitter'
        %  just in case we randomly end up not covering the entire time window
        %  Don't worry, most (typically all) of it is going to be removed
        [~,ix]=sort(rand(N,length(intervals)),2);
        perm=intervals(ix);
        %% here we jitter all Rpeak intervals
        jitter_range = 2*std(intervals);                    % two-sided std of normally distributed jitter values
        td_jit   = [repmat(ts(1),N,1) repmat(td,N,1)+(randn(N,length(td))-0.5)*jitter_range perm]; % jittering every interval!
        ts_jit   = cumsum(td_jit,2);
        %td_jit(1,:)=[];
        
        % now we define
        next_invalid=diff(valid_idx)~=1;                               % is the next Rpeak invalid (i.e. followed by invalid R2R interval)
        iv_starts  =[0  ts(valid_idx([next_invalid true]))];     % start of invalid intervals: Timestamps of valid Rpeaks followed by invalid ones
        % First Segment (for 0 to first valid Rpeak) and last segment
        % (everything after last valid Rpeak) are always invalid
        iv_ends    =[ts(valid_idx([true next_invalid]))   inf];  % end of invalid intervals: Timestamps of valid Rpeaks PRECEDED by invalid ones
        % we can get Rpeak-ts preceded by invalid R2R by shifting next_invalid
        grace_window=mean(intervals)/2;                         % +/- Range for shuffled Rpeaks to be allowed inside invalid segments
        
        
        %% remove jittered Rpeaks and corresponding durations that fell into invalid segments
        for iv=1:numel(iv_starts)
            idx2exclude_1 = ts_jit>iv_starts(iv)+grace_window & ts_jit<iv_ends(iv)-grace_window;
            ts_jit   (idx2exclude_1) = 0;
            td_jit (idx2exclude_1) = NaN;
        end
        idx2exclude_2 = ts_jit>max(ts)+jitter_range;
        ts_jit (idx2exclude_2) = 0;
        td_jit (idx2exclude_2) = NaN;
        %
        %         ts_jit (:,all(ts_jit==0,1)) = [];
        %         td_jit (:,all(isnan(td_jit),1)) = [];
        
    case 'jisi dithering'
        %% here we shuffle (valid!) interval PROPORTIONS (before/after), keeping interval sum equal
        I1=[ts(1) td];
        I2=[td td(end)+td(1)];
        interval_sums=(I2+I1);
        proportions=I1./(I2+I1);
        valid_proportions=proportions(valid_idx);
        valid_sums=interval_sums(valid_idx);
        real_distances=valid_proportions.*valid_sums;
        
        [~,IX]=sort(rand(N,numel(valid_idx)),2);
        displacement=valid_proportions(IX).*repmat(valid_sums,N,1)-repmat(real_distances,N,1);
        valid_mat=bsxfun(@plus,(valid_idx-1)*N,(1:N)');
        
        ts_jit=repmat(ts,N,1);
        ts_jit(valid_mat)=ts_jit(valid_mat)+displacement;
        td_jit          = [zeros(N,1) diff(ts_jit,1,2)];
        
    case 'pseudo jisi dithering'
        %% here we shuffle (valid!) interval PROPORTIONS (before/after), keeping interval sum equal
        I1=[ts(1) td];
        I2=[td td(end)+td(1)];
        interval_sums=(I2+I1);
        proportions=I1./(I2+I1);
        valid_proportions=proportions(valid_idx);
        valid_sums=interval_sums(valid_idx);
        real_distances=valid_proportions.*valid_sums;
        
        
        
        valid_mat=bsxfun(@plus,(valid_idx-1)*N,(1:N)');
        %shuffled_I1=repmat(I1,N,1);
        valid_I1   =I1(valid_idx);
        [~,IX]=sort(rand(N,numel(valid_idx)),2);
        shuffled_valid_I1=valid_I1(IX);
        %shuffled_I1(valid_mat)=shuffled_valid_I1;
        
        %shuffled_I2=repmat(I2,N,1);
        valid_I2   =I2(valid_idx);
        [~,IX]=sort(rand(N,numel(valid_idx)),2);
        shuffled_valid_I2=valid_I2(IX);
        %shuffled_I2(valid_mat)=shuffled_valid_I2;
        
        shuffled_proportions=shuffled_valid_I1./(shuffled_valid_I1+shuffled_valid_I2);
        
        
        displacement=shuffled_proportions.*repmat(valid_sums,N,1)-repmat(real_distances,N,1);
        
        %displacement=valid_proportions(IX).*repmat(valid_sums,N,1)-repmat(real_distances,N,1);
        
        ts_jit=repmat(ts,N,1);
        ts_jit(valid_mat)=ts_jit(valid_mat)+displacement;
        td_jit          = [zeros(N,1) diff(ts_jit,1,2)];
        
        
    case 'isi dithering'
        %  just in case we randomly end up not covering the entire time window
        %  Don't worry, most (typically all) of it is going to be removed
        [~,ix]=sort(rand(length(td),1)); % index of a randomly chosen event that we keep and jitter intervals around
        complete_sets=floor(N/numel(ix));
        rest=mod(N,numel(ix));
        ix=[repmat(ix,complete_sets); ix(1:rest)];
        tdf=fliplr(td);
        
        %% here we jitter all Rpeak intervals
        jitter_range = 2*std(intervals);                    % two-sided std of normally distributed jitter values
        td_jit   = repmat([tdf td tdf],N,1)+(randn(N,3*length(td))-0.5)*jitter_range; % jittering every interval!
        ts_jit   = cumsum(td_jit,2);
        adjustments=ts(ix)-ts_jit(sub2ind(size(ts_jit),1:N,ix'+numel(td)));  % could add jitter to this one too?
        ts_jit=ts_jit+repmat(adjustments',1,size(ts_jit,2));
        
        % now remove excess (timesteps
        grace_window=mean(intervals)/2;
        ts_jit(ts_jit<0)=0;
        ts_jit(ts_jit>ts(end)+grace_window)=0;
        to_keep=any(ts_jit,1);
        ts_jit=ts_jit(:,to_keep);
        td_jit=td_jit(:,to_keep);
        
        
        %now we define invalid intervals
        next_invalid=diff(valid_idx)~=1;                         % is the next Rpeak invalid (i.e. followed by invalid R2R interval)
        iv_starts  =[0  ts(valid_idx([next_invalid true]))];     % start of invalid intervals: Timestamps of valid Rpeaks followed by invalid ones
        % First Segment (for 0 to first valid Rpeak) and last segment(everything after last valid Rpeak) are always invalid
        iv_ends    =[ts(valid_idx([true next_invalid]))   inf];  % end of invalid intervals: Timestamps of valid Rpeaks PRECEDED by invalid ones
        % we can get Rpeak-ts preceded by invalid R2R by shifting next_invalid
        grace_window=mean(intervals)/2;                          % +/- Range for shuffled Rpeaks to be allowed inside invalid segments
        
        
        %% remove jittered Rpeaks and corresponding durations that fell into invalid segments
        
        idx2exclude_1          = ts_jit==0;
        td_jit (idx2exclude_1) = NaN;
        for iv=1:numel(iv_starts)
            idx2exclude_2 = ts_jit>iv_starts(iv)+grace_window & ts_jit<iv_ends(iv)-grace_window;
            ts_jit (idx2exclude_2) = 0;
            td_jit (idx2exclude_2) = NaN;
        end
        idx2exclude_3 = ts_jit>max(ts)+jitter_range;
        ts_jit (idx2exclude_3) = 0;
        td_jit (idx2exclude_3) = NaN;
        ts_jit (:,all(ts_jit==0,1)) = [];
        td_jit (:,all(isnan(td_jit),1)) = [];
        
    case 'train_jitter'
        %% here we jitter the entire train of all Rpeaks instead
        jitter_range    = mean(intervals)*2;                    %
        jitterthing     = repmat(rand(N,1)-0.5,1,numel(ts))*jitter_range;
        ts_jit          = repmat(ts,N,1)+jitterthing;    % jittered timestamps
        td_jit          = [zeros(N,1) diff(ts_jit,1,2)]; % durations of intervals (before corresponding trigger n)
    case 'uniform dithering'
        %% here we jitter all Rpeaks individually instead - within a box that spans half into both intervals
        % ignore those with long intervals around them - not valid because not consecutive
        interval_min    = IN.interval_starts; %[0 td/2];
        interval_ran    = IN.interval_ends-IN.interval_starts; %[0 td/2];
        %interval_min    = [0 td/2];
        %interval_ran    = interval_min+[td/2 0];
        
        %jitterthing     = rand([N,numel(ts)]).*repmat(interval_ran,N,1)-repmat(interval_min,N,1);
        jitterthing     = rand([N,numel(ts)]).*repmat(interval_ran,N,1)+repmat(interval_min,N,1);
        ts_jit          = repmat(ts,N,1)+jitterthing;
        td_jit          = [zeros(N,1) diff(ts_jit,1,2)];
        
end

%% take data corresponding to valid events
if ismember(cfg.jitter_method,{'uniform dithering','train_jitter','jisi dithering','pseudo jisi dithering'})
    ts  = ts(valid_idx);                                          % take only Rpeaks surrounded by valid R2R
    V   = repmat(ismember(1:size(ts_jit,2),valid_idx),N,1);       % logical index to reduce
    Triggers.shuffled_ts  =reshape(ts_jit(V),N,numel(ts));
    Triggers.shuffled_intervals = reshape(td_jit(V),N,numel(ts)); % durations of reshuffled RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).shuffled_ts)
else
    Triggers.shuffled_ts  =ts_jit;
    Triggers.shuffled_intervals = td_jit; % durations of reshuffled RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).shuffled_ts)
end

%% put the data together
Triggers.blocks       =blocks;
Triggers.ts           =ts;
Triggers.intervals    =intervals; %% all intervals deemed valid (see outside function)








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

