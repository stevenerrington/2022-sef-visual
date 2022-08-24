function dirs = setup_dirs(user)

switch user
    case 'erringsp'
        if ispc()
            dirs.server = 'D:';
            dirs.root = fullfile(dirs.server,'projectCode','project_visualSEF\');
            dirs.data_proc = fullfile(dirs.root,'data');
            dirs.data_raw = fullfile(dirs.server,'data','2012_Cmand_EuX');
            dirs.data_spike = fullfile(dirs.data_raw,'spikingData');

        else
            dirs.server = '/Volumes/Alpha/';
            dirs.root = fullfile(dirs.server,'projectCode','project_visualSEF');
            dirs.data_proc = fullfile(dirs.root,'data');
            dirs.data_raw = fullfile(dirs.server,'data','2012_Cmand_EuX');
            dirs.data_spike = fullfile(dirs.data_raw,'spikingData');
        end
end

end