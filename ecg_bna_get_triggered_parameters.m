function [triggered, significance] = ecg_bna_get_triggered_parameters(site,triggers, width_in_samples,realD,cfg)
significance=struct;
% ecg_bna_get_triggered_split_shuffled - gets the time-frequency
% pow,phs spectrogram and phsBP for a specified time window around all shuffled Rpeak onset for a single site for
% given trials (usually trials belonging to a condition) in a session
%
% USAGE:
%	triggered = ecg_bna_get_triggered_split_shuffled( site_lfp,
%	cond_ecg, state, ecg_bna_cfg )
%
% INPUTS:
%		site_lfp            - struct containing LFP data for all trials of
%		a session from a single site
%       cond_ecg         - ecg data of trials whose time freq spectrogram
%       are to be obtaied
%       state               - a cell array specifying time window during
%       which LFP tfs should be obtained
% OUTPUTS:
%		triggered   - struct containing LFP tfs for all given
%		trials
%
% REQUIRES:	lfp_tfa_baseline_normalization ???
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_Rpeak_evoked_LFP,
% ecg_bna_get_shuffled_Rpeak_evoked_LFP




triggered = struct;
w_samples=width_in_samples(1):width_in_samples(2);

%% get rid of some events to get a matrix (!!!!)
n_shuffles=size(triggers,1);
NNN=triggers==0;
AAA=all(NNN,1);
TT=triggers(:,~AAA);
NNN=TT==0;
nRpeaks=sum(~NNN,2);

%% reducing rpeaks - this can lead to 0 remaining shuffled Rpeaks quite easily...
nRpeaks=sum(~NNN,2);
minvalid=min(nRpeaks);
nrows=size(NNN,2);
todrawfrom=1:nrows;
for sh = 1: n_shuffles
    F=find(NNN(sh,:));
    todraw=todrawfrom;
    todraw(F)=[];
    randindex=rand(numel(todraw),1);
    [~, randindex]=sort(randindex);
    ndrawn=nrows-numel(F)-minvalid;
    toremove(sh,:)=   [F todraw(randindex(1:ndrawn))]+(sh-1)*nrows;
end

trigs=TT';
trigs(toremove)=[];
%trigs=reshape(trigs,[n_shuffles,minvalid]); %% check dimensions
trigs=reshape(trigs,[minvalid,n_shuffles]); %% check dimensions
trigs=trigs';

% %% alternative: remove indexes where more than half are invalid
% 
% trigs=zeros(n_shuffles,max(nRpeaks));
% 
% [~,sortix]=sort(-nRpeaks);
% for sh = sortix'
%     trigs(sh,1:nRpeaks(sh))=TT(sh,~TT(sh,:)==0);
% end
%     
% toremove=sum(trigs==0,1)>size(trigs,1)/2;
% TT(:,toremove)=[];
% for tr = 1:size(TT,2)
%     F=TT(:,tr)==0;
%     F
%     
%     todraw=todrawfrom;
%     todraw(F)=[];
%     randindex=rand(numel(todraw),1);
%     [~, randindex]=sort(randindex);
%     ndrawn=nrows-numel(F)-minvalid;
%     toremove(sh,:)=   [F todraw(randindex(1:ndrawn))]+(sh-1)*nrows;
% end
% 
% %%%

tfs=[site.tfs];
%% check pre-allocation time saving (?)
% 
    FN_in={'pow','powbp','lfp','pha','phabp'};
    FN={'pow','powbp','lfp','itpc','itpcbp'};

    for f=1:numel(FN_in)  
        fn=FN{f};      
        fni=FN_in{f};
        S2=size(tfs.(fni),1);
        S3=numel(w_samples);
        
        triggered.(fn).mean        = zeros(1,S2,S3);
        triggered.(fn).std         = zeros(1,S2,S3);
        triggered.(fn).conf95      = zeros(2,S2,S3);
    end

%% NOW: loop through window samples to compute shuffled parameters all at once!

for s = 1: numel(w_samples)
    w=w_samples(s);
    t=trigs+w;
    triggered.ntriggers=minvalid;
    
    for f=1:numel(FN_in)
        fn=FN{f};
        fni=FN_in{f};
        
        AA=      reshape(tfs.(fni)(:,t),size(tfs.(fni),1),size(t,1),size(t,2));
        if ismember(fni,{'pha','phabp'})
            BB=abs(mean(AA,3));
        else
            BB=mean(AA,3);
        end
        
        triggered.(fn).mean(1,:,s)        = mean(BB,2);
        triggered.(fn).std(1,:,s)         = std(BB,0,2);
        %% this one not sure what to put
        %triggered.(fn).std_mean(:,s)    = mean(abs(mean(AA,3)),2);
        %% check dimensions
        triggered.(fn).conf95(1,:,s)      = prctile(BB,97.5,2);
        triggered.(fn).conf95(2,:,s)      = prctile(BB,2.5,2);
        %triggered.(fn).complete(:,:,s)    = BB';
        
        complete.(fn)(:,:,s)    = BB';
    end
    
    
    %% differentiate between
    %%      a) mean of std and std of (shuffle) means
    %%      a) mean confidence interval and confidence interval of (shuffle) means
    
end

if nargin>3 %% surrogate data,add signficance test here !
    for f=1:numel(FN_in)
        fn=FN{f};
        dat=complete.(fn);
        s_ix=1:size(dat,1);
        
        
        tThreshold = abs(tinv(0.05, numel(s_ix)-1));
        
        % find (for each surrogate) the max cluster significance (somehow i
        % feel this should be done PER frequency, hmmmm....)
        for s=s_ix
            dat_s=squeeze(dat(s,:,:));
            est_mean=squeeze(mean(dat(s_ix~=s,:,:),1));
            est_var=squeeze(var(dat(s_ix~=s,:,:),1));
            T=abs((dat_s-est_mean)./sqrt(est_var));
            sig=T>tThreshold;
            CC = bwconncomp(sig);
            sumT = cellfun(@(x) sum(T(x)),CC.PixelIdxList);
            T_max(s) = max([sumT 0]);
        end
        
        %% now find significant clusters in the real data        
        %T_max_thr=prctile(T_max,95);          
        T_max_thr=prctile(T_max,80);    % sum of     
        est_mean=squeeze(mean(dat,1));
        est_var=squeeze(var(dat,1));
        T=abs((squeeze(realD.(fn).mean)-est_mean)./sqrt(est_var));
        sig=T>tThreshold;
        CC = bwconncomp(sig);
        sigC = cellfun(@(x) sum(T(x))>=T_max_thr,CC.PixelIdxList);
        sigpix=vertcat(CC.PixelIdxList{sigC});
        significance.(fn)=zeros(size(realD.(fn).mean));
        significance.(fn)(sigpix)=true;
    end
end


end

