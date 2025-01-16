function Event=ecg_bna_get_events(session_info,trials,event,block_anchors,varargin)

%ecg_bna_get_triggers(session_info
eventtype=strrep(event{1,2}, ' ',''); %??

%% load respective file
switch eventtype
    case {'CAP'}
        load(session_info.Input_ECG);
        nrblockscell={out.nrblock_combinedFiles};
        oldnames = fieldnames(out_cap);
        newnames = strrep(oldnames,'B2B','R2R');
        idx = find(~strcmpi(newnames,oldnames));
        for n = 1: numel(newnames)
            [out_cap.(newnames{n})] = deal(out_cap.(oldnames{n}));
        end
        for n = 1:numel(idx)
            out_cap = rmfield(out_cap,oldnames{idx(n)});
        end
        out = out_cap;
        [out.nrblock_combinedFiles] = deal(nrblockscell{:});
        
    case {'Rpeak','Rpeak_insp','Rpeak_exp','Rpeak_lowIBI','Rpeak_highIBI'} %% why load it several times ... not that it matters right now
        load(session_info.Input_ECG);
end


Event.valid_ts     = [];
Event.ts           = [];
Event.intervals    = [];
Event.blocks       = [];
Event.interval_starts       = [];
Event.interval_ends       = [];
switch eventtype
    case {'Rpeak','Rpeak_insp','Rpeak_exp','CAP','Rpeak_lowIBI','Rpeak_highIBI'}
        %% treat block_merging as invalid?
        for b=1:numel(out)
                       
            %% adjustment for isnp/exp/ median split 
            th=median(out(b).R2R_valid); %% eventually I will replace this with per session threshold values
            IBIs=diff(out(b).Rpeak_t);
            valid_interval=out(b).idx_valid_R2R_consec;
            if strcmp(eventtype,'Rpeak_insp')
                valid_interval=intersect(find(out(b).is_RR_insp),valid_interval);
            elseif strcmp(eventtype,'Rpeak_exp')
                valid_interval=intersect(find(out(b).is_RR_exp),valid_interval);
            elseif strcmp(eventtype,'Rpeak_highIBI')
                valid_interval=intersect(find(IBIs(1:end-1)>=th & IBIs(2:end)>=th),valid_interval);
            elseif strcmp(eventtype,'Rpeak_lowIBI')
                valid_interval=intersect(find(IBIs(1:end-1)<=th & IBIs(2:end)<=th),valid_interval);
            end
            
            [~,consecutive_idx]=ismember(out(b).R2R_t(valid_interval),out(b).Rpeak_t); % indexing of consecutive_idx corrsponds to out(b).R2R_t
            valid_idx=consecutive_idx-1;
            intervals    = out(b).R2R_valid(valid_interval);
            
            if isempty(valid_idx) %isempty(block) || isempty(out(b).Rpeak_t) || isempty(out(b).R2R_t) || isempty(out(b).idx_valid_R2R_consec)
                continue
            end
            
            
            %% add to offset
            
            %% use trial information to get INI length
            %% for a given block, this is what it is (?) -> do we have a global offset (?)
            %trials
            block=out(b).nrblock_combinedFiles;
            %offset=offset+max(max(Event.ts))+max(intervals);%onset of recording relative to state 2
            offset=block_anchors(block); %%%%
            
            Event.valid_ts        = [Event.valid_ts valid_idx+numel(Event.ts) ]; % +numel(ts)
            Event.ts              = [Event.ts out(b).Rpeak_t+offset ]; %+offset
            Event.intervals       = [Event.intervals intervals];
            Event.blocks          = [Event.blocks repmat(block,size(intervals))];
            
            tmp_before=-diff([0 out(b).Rpeak_t]);
            tmp_after=diff([out(b).Rpeak_t inf]);
            Event.interval_starts  = [Event.interval_starts tmp_before/2];
            Event.interval_ends    = [Event.interval_ends tmp_after/2];
%         
%           % now we define
%         next_invalid=diff(valid_idx)~=1;                               % is the next Rpeak invalid (i.e. followed by invalid R2R interval)
%         iv_starts  =[0  ts(valid_idx([next_invalid true]))];     % start of invalid intervals: Timestamps of valid Rpeaks followed by invalid ones
%         % First Segment (for 0 to first valid Rpeak) and last segment
%         % (everything after last valid Rpeak) are always invalid
%         iv_ends    =[ts(valid_idx([true next_invalid]))   inf];  % end of invalid intervals: Timestamps of valid Rpeaks PRECEDED by invalid ones
      
        end
    case 'microstim'
    case 'state'
        state=event{3};
        %find state in each trial
        [O]=ecg_bna_synchronize_with_trigger(block_anchors,trials);
        Event.ts=O.state_onsets(O.states==state);
        Event.valid_ts     = 1:numel(Event.ts);
        Event.blocks       = O.state_blocks(O.states==state);
        Event.interval_starts  = O.state_trial_starts(O.states==state)-Event.ts;
        Event.interval_ends    = O.state_trial_ends(O.states==state)-Event.ts;
        Event.intervals    = Event.interval_ends-Event.interval_starts;
        
end