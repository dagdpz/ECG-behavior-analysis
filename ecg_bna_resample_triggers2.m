function  Trigger_samples = ecg_bna_resample_triggers2(Triggers,offset_per_block,LFP_samples_per_block,sr)

LFP_blocks=LFP_samples_per_block(1,:);
LFP_samples=LFP_samples_per_block(2,:);
trigger_blocks=offset_per_block(1,:);
trigger_offset=offset_per_block(2,:);

events=fieldnames(Triggers)';

for f=events
    
    past_blocks=[];
    ts_real=Triggers.(f{:}).ts;
    ts_surrogate=Triggers.(f{:}).shuffled_ts;
    
    Real_ts=[];    
    Surrogate_ts=[];    
    for b=LFP_blocks
        %B=Rpeaks(bl).block;
        BS=LFP_blocks==b;   %this block's samples in the LFP
        BT=trigger_blocks==b;
        if ~any(BT)
            display(['Block' num2str(b) 'not present in Triggers, check for potential bugs!']);
        end
                
        idx=Triggers.(f{:}).blocks==b;
        
        R=ts_real(idx)-trigger_offset(BT);%+t_offset_per_block(2,BT);
        S=round(R*sr);
        S(S<=0)=0;
        S(S>LFP_samples(BS))=0;
        S=S+sum(LFP_samples(ismember(LFP_blocks,past_blocks)));
        S(S==LFP_samples(BS))=0;
        
        Real_ts=[Real_ts S];
        
        %% not shuffled, but surrogate!
        
        R=ts_surrogate(:,idx)-trigger_offset(BT);%+t_offset_per_block(2,BT);
        S=round(R*sr);
        S(S<=0)=0;
        S(S>LFP_samples(BS))=0;
        S=S+sum(LFP_samples(ismember(LFP_blocks,past_blocks)));
        S(S==LFP_samples(BS))=0;
        
        Surrogate_ts=[Surrogate_ts S];
%         
%         
%         Triggers.(f{:}).ts(idx)=Triggers.(f{:}).ts(idx)-trigger_offset(BT);
%         Triggers.(f{:}).shuffled_ts
                
        %% suspect something's wrong here!!!!!!
        past_blocks=[past_blocks b];
    end
    Trigger_samples.([f{:} '_real'])=Real_ts;
    Trigger_samples.([f{:} '_shuffled'])=Surrogate_ts;
end