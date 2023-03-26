%%this is the main code for processing the ESM data. 
% the code is a part of ESM_proc1m which contains even the functions that
% are being called. this code does not contain the funtions. dont medddle
% with ESM_proc1.m because that is the most important working version.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

tstart = tic;

direc = "../../../../electron_data/jet_data/20072022/"

imext = '.gel';
run_name = 'run1_';
run_name2 = '_sc1';

[scan1,scan_last,total_scan,save_name] = file_nam(direc,run_name,run_name2,imext);
folder =  direc;


%reader settings
%microns_per_pixel = 50; %resolution
microns_per_pixel = 25;
sensit = 4000; %sensitivity of the ip reader
latitude = 5; %latitude of the reader
dyn_range =  16; %16 bit dynamic range of the data

en_array = [35 50 75 100:50:700 800:100:1200 1350 1500];
calib_file = strcat(direc,'ph_jet_mag1.csv');

%name of the table headers
name_total_arr = ["Rand"];
ip_resp_file = "";

num_ESM = 1; 
is_ang_dist = 1;
%the space on the ip before the magnet begins (in cm)
%the calibration file has its zero at the start  of the magnet 
%the gap here is the gap between the aperture and the corner of the 
%magnet
gap = 1.5;


%here we have to enter the coordinates of the parts where the actual image exists within the scan
% 1 is the left top corner
% 2 is the bottom right corner
% there is a slight adjustment done later. i have taken notations from imageJ
%in imageJ x-axis is counted from left to right; y-axis is counted from top to bottom

topLeft = [186,774];
botright = [8220,1830];


%creating this list just to go as an input to the img_proc function
total_crop = [topLeft;botright];

%this is the width in pixels around the central point where summation is taken columnwise for the image 
%required for ESM
trend_width = 40;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp("variable loading done");
tim = toc(tstart);
tim = tic;
%reading the calibration files and doing performing the calibration
%this takes values from the comsol files
calib_cell = cell(2,num_ESM);
data_elec = readtable(calib_file);
data_en_xax = en_array;
[pos_val,en_val] = en_calib(data_elec,microns_per_pixel,data_en_xax);
calib_cell(1) = {en_val};
calib_cell(2) = {pos_val};

disp("Callibration files processed")
toc(tim);

scan_num = 1;
inv_arr = [0];
axis_cell = cell(length(total_scan),num_ESM);
trend_cell = cell(length(total_scan),num_ESM);
satidx_cell = cell(length(total_scan),num_ESM);
isSat_cell = cell(length(total_scan),num_ESM);


%% this is reading the last file and doing operations on it
for i = 1:length(total_scan)
    tim = tic;
    disp('reading file')
    if (strcmp(imext,".tif"))
        img_fil2 = Tiff(total_scan(i));
        %cl = class(img_fil)
        %getTag(img_fil2,'YResolution');
        img = read(img_fil2);
    end
    if(strcmp(imext,".gel"))
        total_scan(i)
        img = imread(total_scan(i));
    end
    disp("file read")
    toc(tim)
    [axis_arr,trend_arr,isSatarr,satidx_arr,inv_arr] = img_proc(scan_num,img,total_crop,microns_per_pixel,num_ESM,sensit,latitude,dyn_range,trend_width,gap,inv_arr); 
    axis_cell(i,:) = axis_arr(:);
    trend_cell(i,:) = trend_arr(:);
    isSat_cell(i,:) = isSatarr(:);
    satidx_cell(i,:) = satidx_arr(:);
    scan_num = scan_num + 1;
    disp("1 file has been processed")
    toc(tim)
end

tim = tic;
%removing saturation
for i = 1:num_ESM
    unsat_final_all = total_sat_rem(trend_cell(:,i),isSat_cell(:,i),satidx_cell(:,i));
    trend_cell(:,i) = unsat_final_all;
end
disp("removed all saturation");
toc(tim);

%gap = gap*(10^-4);
%gap = gap/(microns_per_pixel^(10^-6));

final_dat_cell = cell(3,num_ESM);
for i  = 1:num_ESM

    temp_trend_ar = cell2mat(trend_cell(1,i));
    temp_en_ar = cell2mat(calib_cell(1,i));
    temp_calib_pos_ar = cell2mat(calib_cell(2,i));
    calib_length = length(temp_en_ar);
    temp_trend_ar = temp_trend_ar(1:calib_length);
    %plot(temp_trend_ar)
    temp_trend_ar = jacobian(temp_en_ar,temp_calib_pos_ar,temp_trend_ar);
    %plot(temp_trend_ar)
    temp_en_ar = temp_en_ar(1:end-1);
    temp_trend_ar = eff_corr(temp_en_ar,temp_trend_ar);
    %plot(temp_en_ar,temp_trend_ar)
    temp_norm_ar = temp_trend_ar/(max(temp_trend_ar));

    
    plot(temp_en_ar,temp_trend_ar)
    final_dat_cell(2,i) = {temp_trend_ar};
    final_dat_cell(1,i) = {temp_en_ar};
    final_dat_cell(3,i) = {temp_norm_ar};
end

dat_len_max = max(max(cellfun('length',final_dat_cell)));
for i = 1:num_ESM
    for j = 1:3
        if length(final_dat_cell(j,i))<dat_len_max
            temp_ar = cell2mat(final_dat_cell(j,i));
            temp_ar = [temp_ar, zeros(1, dat_len_max - length(temp_ar))];
            final_dat_cell(j,i) = {temp_ar};
        end
    end
end

dat_len_max = length(cell2mat(final_dat_cell(1,1)));

final_dat_mat = zeros(dat_len_max,num_ESM*3);
csv_ind_names = string.empty;

%just naming the columns before printing
for i = 1:num_ESM
    for j = 1:3
        %3*(i-1)+j
        %plot(cell2mat(final_dat_cell(j,i)))
        final_dat_mat(:,3*(i-1)+j) = cell2mat(final_dat_cell(j,i)); 
    end
    en_nam = strcat(name_total_arr(i),"_en(keV)");
    csv_ind_names(end+1) = en_nam;
    t_nam = name_total_arr(i);
    csv_ind_names(end+1) = t_nam;
    t_norm_nam = strcat(t_nam,"_norm");
    csv_ind_names(end+1) = t_norm_nam;
end

savetable = array2table(final_dat_mat,'VariableNames',csv_ind_names);
writetable(savetable,save_name);