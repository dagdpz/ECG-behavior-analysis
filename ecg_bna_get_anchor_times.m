function blockstart=ecg_bna_get_anchor_times(monkey,session,blocks)

%u_block=unique([trials.block]);

for B=blocks
    BLOCK_PATH=['Y:\Data\TDTtanks\' monkey '_phys\' session '\Block-' num2str(B)];
    OPTIONS={'TYPE', {'EPOCS'}};
    data = TDTbin2mat_working(BLOCK_PATH, OPTIONS{:});
    blockstart(B)=seconds(datetime(data.info.utcStartTime)-datetime('00:00:00'));
end

% 
% runsstart=
% 
% trials
% 
% 
%             FIX_ACQ_start_time              =find(TDT_DATA.Trial(tr).purestates==2)./data.streams.stat.fs;
%             FIX_ACQ_start_time              =FIX_ACQ_start_time(1);
%         else
%             trial_time                  =[trialonsets(tr_block) trialonsets(tr_block+1)];
%             trial_states_indexes        =data.epocs.SVal.onset>=trial_time(1) & data.epocs.SVal.onset<=trial_time(2);
%             TDT_DATA.Trial(tr).states   =data.epocs.SVal.data(trial_states_indexes);
%             state_onsets_temp           =data.epocs.SVal.onset(trial_states_indexes);
%             FIX_ACQ_start_time          =data.epocs.SVal.onset(data.epocs.SVal.data==2)-data.epocs.Tnum.onset(1:end-1);

% OPTIONS={};
% data = TDTbin2mat_working(BLOCK_PATH, OPTIONS{:});
 blockstart= blockstart- blockstart(1);

end