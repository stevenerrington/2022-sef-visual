for session = 14:29
    session
    session_contacts_list = [];
    session_contacts_list = find(corticalLFPmap.session == session);
    n_contacts = length(session_contacts_list);
    
    csd_lfp_aligned = [];
    
    for ch_i=1:n_contacts
        
        filename_LFP = ['lfp_session' int2str(session) '_' ...
            corticalLFPmap.channelNames{session_contacts_list(ch_i)} '_target'];
        
        data = load(fullfile(dataDir,'lfp', filename_LFP)); 
        depth_i = corticalLFPmap.depth(session_contacts_list(ch_i));
        
        for trl_i = 1:size(data.filteredLFP,1)
            csd_lfp_aligned(depth_i,:,trl_i) = data.filteredLFP(trl_i,:);
        end
        
    end
    
    csd_analysis{session} = D_CSD_BASIC(csd_lfp_aligned, 'cndt', 0.0004, 'spc', 0.15);

end

test = nan(16,3000,length(14:29));

for session = 14:29
    inputCSD = []; inputCSD = csd_analysis{session};
    test(1:size(inputCSD,1),1:size(inputCSD,2),session-13) =...
        nanmean(inputCSD(:, :, executiveBeh.ttx.NC{session}),3);
end

f_h = figure; hold on;

ax1 = subplot(1, 1, 1);
P_CSD_BASIC(nanmean(test(2:end-1,:,:),3), [-999:2000], [-100 300], f_h, ax1)

    
    