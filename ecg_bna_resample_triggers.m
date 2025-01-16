
function                    triggers = ecg_bna_resample_triggers(Rpeaks,fn,t_offset_per_block,n_samples,ts)
b=0;
past_blocks=[];
for bl=1:numel(Rpeaks)
    B=Rpeaks(bl).block;
    
    BT=t_offset_per_block(1,:)==B;
    BS=n_samples(1,:)==B;   %this block's samples in the LFP
    
    if ~any(BS) || ~any(BT)
        continue
    end
    
    R=Rpeaks(bl).(fn)-Rpeaks(bl).offset;%+t_offset_per_block(2,BT);
    S=round(R/ts);
    S(S<=0)=0;
    S(S>n_samples(2,BS))=0;
    b=b+1;
    triggers(b).event_sample=S+sum(n_samples(2,ismember(n_samples(1,:),past_blocks)));
    triggers(b).block=B;
    
    %% suspect something's wrong here!!!!!!
    %past_blocks=[past_blocks triggers.block]; 
    past_blocks=[triggers.block]; 
end
end