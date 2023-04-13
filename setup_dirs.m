function dirs = setup_dirs(user)

switch user
    case 'erringsp'
        if ispc()
            dirs.server = 'D:';
            dirs.root = fullfile(dirs.server,'projectCode','2022-sef-visual');
            dirs.data_proc = fullfile(dirs.root,'data');
            dirs.data_raw = fullfile(dirs.server,'data','2012_Cmand_EuX');
            dirs.data_spike = fullfile(dirs.data_raw,'spikingData');
            
        else
            dirs.server = '/Volumes/Alpha/';
            dirs.root = fullfile(dirs.server,'projectCode','2022-sef-visual');
            dirs.data_proc = fullfile(dirs.root,'data');
            dirs.data_raw = fullfile(dirs.server,'data','2012_Cmand_EuX');
            dirs.data_spike = fullfile(dirs.data_raw,'spikingData');
        end
        
    case 'erringsp_WH633'
        dirs.server = 'C:';
        dirs.root = fullfile(dirs.server,'Users','Steven','Desktop','Projects','2022-sef-visual');
        dirs.data_proc = fullfile(dirs.root,'data');
        dirs.data_raw = fullfile('T:/','Users','Steven','dataRepo','2012_Cmand_EuX');
        dirs.data_spike = fullfile(dirs.data_raw,'ePhys');
        
end

addpath(genpath(dirs.root));


end