function [One_stream]=ecg_bna_synchronize_with_trigger(Blockoffsets,trcell,popcell)
%ecg_bna_synchronize_with_trigger(trcell,popcell,Rpeaks)
%Bockoffsets=Rpeaks([Rpeaks.block] == z).offset
%Blockoffsets=Blockoffsets-Blockoffsets(1);
% 0. Prepare data variables
states_onset                = {trcell.states_onset};
states                      = {trcell.states};

%% does this NEED to be LFPx, not general t0 from TDT file start?
rec_start                   = {trcell.TDT_LFPx_t0_from_rec_start}; % this is the relative time (First_trial_INI is NOT included (?))
block_nums                  = {trcell.block};

state_onsets                = cellfun(@(x,y,z) x+y+Blockoffsets(z), states_onset,  rec_start, block_nums, 'Uniformoutput', false);
state_blocks                = cellfun(@(x,y) repmat(y,size(x)), states_onset,  block_nums, 'Uniformoutput', false);
state_trial_start           = cellfun(@(x,y) repmat(x(y==2),size(x)), state_onsets,  states, 'Uniformoutput', false);
state_trial_end             = cellfun(@(x,y) repmat(x(y==90),size(x)), state_onsets,  states, 'Uniformoutput', false);



One_stream.state_onsets     = [state_onsets{:}];
One_stream.states           = [states{:}];
One_stream.state_blocks     = [state_blocks{:}];
One_stream.trial_starts     = One_stream.state_onsets(One_stream.states==2);
One_stream.trial_ends       = One_stream.state_onsets(One_stream.states==90);
One_stream.state_trial_starts     = [state_trial_start{:}];
One_stream.state_trial_ends       = [state_trial_end{:}];

%% is there a way, starts and ends can be of different length? hope not
if nargin>2
    arrival_times              = {popcell.arrival_times};
    arrival_times              = cellfun(@(x,y,z) x+y+Blockoffsets(z), arrival_times, rec_start, block_nums, 'Uniformoutput', false);
    One_stream.AT              = vertcat(arrival_times{:});
    One_stream.AT(One_stream.AT>One_stream.trial_ends(end))=[];
end
end