function data = util_bandpassFilter( data_in, data_fs, t_lims, band, band_name, filt_order, varargin )

% c/o Jacob Westerberg


band_name_time = [band_name '_time'];
band_name_fs = [band_name '_fs'];

data.(band_name).hpc = band(1);
data.(band_name).lpc1 = band(2);
data.(band_name).lpc2 = band(1)/2;

data.(band_name).filt_order = filt_order;

do_power = true;
do_oscil = true;

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-p','power'}
            do_power = varargin{varStrInd(iv)+1};
        case {'-o','oscil'}
            do_oscil = varargin{varStrInd(iv)+1};
    end
end

if strcmp(band_name, 'mua')
    do_oscil = false;
end

hWn = data.(band_name).hpc / (data_fs/2);
[ bwb, bwa ] = butter( filt_order, hWn, 'high' );

hphga = filtfilt( bwb, bwa, data_in' );

lWn = data.(band_name).lpc1 / (data_fs/2);
[ bwb, bwa ] = butter( filt_order, lWn, 'low' );
hphga = filtfilt( bwb, bwa, hphga );

if do_oscil
    
    hphga_d = [];
    if data_fs > G_FS('slow')*1.5
        for i = 1 : size(hphga, 2)
            hphga_d = cat(2,hphga_d, decimate2( hphga(:,i), single(floor(data_fs / G_FS('slow'))) ));
        end
        new_fs = data_fs / (single(floor(data_fs / G_FS('slow'))));
    else
        hphga_d = hphga;
        new_fs = data_fs;
    end
    

    data.data = single(hphga_d');
    data.(band_name_fs) = new_fs;
    data.(band_name_time) = 0:(size(hphga_d, 1) - 1);
    data.(band_name_time) = data.(band_name_time) .* (1000 / data.(band_name_fs));
    data.(band_name_time) = round(data.(band_name_time) + t_lims(1));
end

if do_power
    
    if strcmp(band_name, 'mua')
        band_name_power_time = [band_name '_time'];
        band_name_power_fs = [band_name '_fs'];
        band_name_power = band_name;
    else
        band_name_power_time = [band_name '_pwr_time'];
        band_name_power_fs = [band_name '_pwr_fs'];
        band_name_power = [band_name '_pwr'];
    end
    
    hphga = abs( hphga );
    
    lWn = data.(band_name).lpc2 / (data_fs/2);
    [ bwb, bwa ] = butter( filt_order, lWn, 'low' );
    hphga = filtfilt( bwb, bwa, hphga );
    
    hphga_p = [];
    if data_fs > G_FS('slow')*1.5
        for i = 1 : size(hphga, 2)
            hphga_p = cat(2,hphga_p, decimate2( hphga(:,i), single(floor(data_fs / G_FS('slow'))) ));
        end
        new_fs = data_fs / (single(floor(data_fs / G_FS('slow'))));
    else
        hphga_p = hphga;
        new_fs = data_fs;
    end
    
    data.(band_name_power) = single(hphga_p');
    data.(band_name_power_fs) = new_fs;
    data.(band_name_power_time) = 0:(size(hphga_p, 1) - 1);
    data.(band_name_power_time) = data.(band_name_power_time) .* (1000 / data.(band_name_power_fs));
    data.(band_name_power_time) = round(data.(band_name_power_time) + t_lims(1));
end

end