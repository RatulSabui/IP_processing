%%this code is for processing the ESM data
%taken from elec_dat_precode2.m which does only the angular distriburtion
%data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

tstart = tic;

direc = "../../particle_exp_round2/4Oct2021/";
%img = 'run1_50_4000_sc1-[Phosphor].tif';
run_name = 'run1_';

[scan1,scan_last,total_scan,save_name] = file_nam(direc,run_name);
%direc = fig_dir;
folder =  direc;
%direc = "elec_data_ip_esm_angdist"


%reader settings
microns_per_pixel = 50; %resolution
sensit = 4000; %sensitivity of the ip reader
latitude = 5; %latitude of the reader
dyn_range =  16; %16 bit dynamic range of the data


en_array1 = [50:50:700 800:100:1200 1350 1500];
calib_file1 = strcat(direc,'round2_mag1_final.csv');
calib_file2 = strcat(direc,'round2_mag2_final.csv');
%calib_file1 = strcat(direc,'height_mag1_1.csv');
%calib_file2 = strcat(direc,'height_mag1_1.csv');
coms_en_ax_total_arr = [en_array1;en_array1;en_array1];
coms_file_total = [calib_file1;calib_file2;calib_file2];
%front is 45 deg, side is 90 deg and back is 135 deg
%name_total_arr = ["back","side","front"];
name_total_arr = ["135","90","45"];
ip_resp_file = "";

num_ESM = 3;
is_ang_dist = 1;
%the space on the ip before the magnet begins (in cm)
%the calibration file has its zero at the start  of the magnet 
%the gap here is the gap between the aperture and the corner of the 
%magnet
%gap = [1.5,1,1];
gap = 1.5;


%here we have to enter the coordinates of the parts where the actual image exists within the scan
% 1 is the left top corner
% 2 is the bottom right corner
% there is a slight adjustment done later. i have taken notations from imageJ
%in imageJ x-axis is counted from left to right; y-axis is counted from top to bottom


back135_1 = [47,1779];
back135_2 = [3009,2132];

side90_1 = [44,2312];
side90_2 = [3068,2779];

front45_1 = [49,2911];
front45_2 = [2998,3360];
% 
% back135_1 = [49,2911];
% back135_2 = [2998,3360];
% 
% side90_1 = [49,2911];
% side90_2 = [2998,3360];

%creating this list just to go as an input to the img_proc function
total_crop = [back135_1;back135_2;side90_1;side90_2;front45_1;front45_2];

%this is the width in pixels around the central point where summation is taken columnwise for the image 
%required for ESM
trend_width = 40;

%this part is automated now
%make flip = 0 if the image need not be flipped
%make flip = 1 if the image needs to be flipped
%flip135 = 0;
%flip90 = 0;
%flip45 = 0;

%tot_flip = [flip135,flip90,flip45];

%len_135 = 100.0;
%len_90 = 100.0;
%len_45 = 100.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this is to find the circles and the correction angle
%img_fil = Tiff(scan1);
%cl = class(img_fil)
%getTag(img_fil,'YResolution');
%B = read(img_fil);
%class(B)
%[corn,length,breadth] = rect_find(B,num_ESM,is_ang_dist);



disp("variable loading done");
tim = toc(tstart);
tim = tic;
%reading the calibration files and doing performing the calibration
%this takes values from the comsol files
calib_cell = cell(2,num_ESM);
for i =1:num_ESM
    data_elec = readtable(coms_file_total(i,:));
    data_en_xax = coms_en_ax_total_arr(i,:);
    [pos_val,en_val] = en_calib(data_elec,microns_per_pixel,data_en_xax);
    %plot(pos_val,en_val)
    calib_cell(1,i) = {en_val};
    calib_cell(2,i) = {pos_val};
end

disp("Callibration files processed")
toc(tim);


%cd (folder);

scan_num = 1;
inv_arr = [0,0,0];
axis_cell = cell(length(total_scan),num_ESM);
trend_cell = cell(length(total_scan),num_ESM);
satidx_cell = cell(length(total_scan),num_ESM);
isSat_cell = cell(length(total_scan),num_ESM);


%% this is reading the last file and doing operations on it
for i = 1:length(total_scan)
    tim = tic;
    disp('reading file')
    img_fil2 = Tiff(total_scan(i));
    %cl = class(img_fil)
    %getTag(img_fil2,'YResolution');
    img = read(img_fil2);
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
    %final_dat_cell(1,i) = trend_cell(1,i);
    %final_dat_cell(2,i) = calib_cell(1,i);
    %final_dat_cell(3,i) = {cell2mat(final_dat_cell(1,i)) / max(cell2mat(final_dat_cell(1,i)))};
%     temp_trend_ar = cell2mat(final_dat_cell(1,i));
%     temp_trend_ar = temp_trend_ar(gap(i):end);
%     final_dat_cell(1,i) = {temp_trend_ar};
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
%csv_ind_names
%final_dat_mat = [csv_ind_names;final_dat_mat]
%savetable = array2table(final_dat_mat,'VariableNames',{'135en_keV','135_count','135_norm','90en_keV','90_count','90_norm','45en_keV','45_count','45_norm'});
savetable = array2table(final_dat_mat,'VariableNames',csv_ind_names);
writetable(savetable,save_name);
%csv_ind_names
%csv_ind_names = cellstr(csv_ind_names);
%savetable.Properties.VariableNames = {'en(keV)','135_count','135_norm','en(keV)','90_count','90_norm','en(keV)','45_count','45_norm'}

% 
% 
%     %C_crop = ang_dist_crop(C_rot,holes_rot,microns_per_pixel);
%     %figure('name','cropped image');
%     %imshow(C_crop);
%     Sat = isSat(C_crop);
%     %surf(C_crop);
%     %colorbar;
%     %C_psl = psl_calc(C_crop,microns_per_pixel,sensit,latitude,dyn_range);
%     %surf(C_psl);
%     %colorbar;
%     [x_ax,dist] = trendgen(C_crop);
%     %figure('name','trend');
%     %plot(x_ax,dist);
%     dist_norm = dist / max(dist);
%     ax_name = strcat('scan',int2str(i),'_ax');
%     dist_name = strcat('scan',int2str(i),'_dist');
%     dist_norm_name = strcat('scan',int2str(i),'_distnorm');
%     dat_table.(ax_name) = transpose(x_ax);
%     dat_table.(dist_name) = transpose(dist);
%     dat_table.(dist_norm_name) = transpose(dist_norm);
% 
% writetable(dat_table,save_name,'Delimiter',',');

function [corn,length,breadth] = rect_find(img,num_ESM,is_ang_dist)
    T = 0.00007;
    bin = imbinarize(img,T);
    %figure('name','binary');
    %imshow(bin)
    bin = medfilt2(bin, [6 6]);
    %imshow(bin);
    %se = strel('rectangle',[6 5]);
    %bin = imclose(bin,se);
    %imshow(bin);
    
    if(is_ang_dist==1)
        iter = num_ESM + 1;
        %class(bin)
        ang_rem = bwareafilt(bin,1);
        bin = bin - ang_rem;
        bin = logical(bin);
        %imshow(bin)
    else
        iter = num_ESM;
        
    end
    
    %class(bin)
    bin = bwareafilt(bin,3);
    %imshow(bin)
    s = regionprops(bin,'Area','BoundingBox');
    temp = struct2table(s)
    
%     [y,x] = find(bin==1);
%     idx = kmeans([x y],iter);
%     
%     for i = 1:iter
%         plot(x(idx==i),y(idx==i),'.')
%         % Get L and B of rectangles 
%         x0 = min(x(idx==i)) ; x1 = max(x(idx==i)) ;
%         y0 = min(y(idx==i)) ; y1 = max(y(idx==i)) ;
%         %
%         L = x1-x0 ;
%         B = y1-y0 ;
%         coor = [x0 y0 ; x0 y1 ; x1 y1 ; x1 y0] ;
%         patch(coor(:,1),coor(:,2),rand(1,3))
% 
%     end
%     
    corn=0;length=0;breadth = 0;
end

%% this function finds the holes on the angular distribution plate
function [cent,select_cent] = cent_find(img)

    %tagNames = img_fil.getTagNames();
    %imshow(B);
    %T = graythresh(B)
    T = 0.0005;
    bin = imbinarize(img,T);
    %figure('name','binary');
    imshow(bin)
    %bleh
    bin = bwareaopen(bin,1000);
    com = imcomplement(bin);
    rem = bwareaopen(com,1000);
    rem2 = rem - bwareaopen(rem,100000);
    %imshow(rem2);
    %rem3 = imrotate(rem2,-50);
    %imshow(rem3);
    [B,L] = bwboundaries(rem2,'noholes');
    %bwim = rem2;
    
    %%this part is just for plotting the circles and the boundaries
    imshow(label2rgb(L,@jet,[.5 .5 .5]))
    h = gca;
    h.Visible = 'On';
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
    end

    stats = regionprops(L,'Area','Centroid');
    
    struct2table(stats)
    cent = stats;
%    threshold = 0.94;

    
    area_array = [];
    % loop over the boundaries
    % selecting the 2 biggest regions
    for k = 1:length(B)

        % obtain (X,Y) boundary coordinates corresponding to label 'k'
        %boundary = B{k};

        % compute a simple estimate of the object's perimeter
        %delta_sq = diff(boundary).^2;
        %perimeter = sum(sqrt(sum(delta_sq,2)));

        % obtain the area calculation corresponding to label 'k'
        %area = stats(k).Area;

        % compute the roundness metric
        %metric = 4*pi*area/perimeter^2;

        % display the results
%         metric_string = sprintf('%2.2f',metric);

        % mark objects above the threshold with a black circle
%         if metric > threshold
            %centroid = stats(k).Centroid;
            %plot(centroid(1),centroid(2),'ko');
%         end
          area_temp = stats(k).Area;
          area_array(end+1) = area_temp;
%         text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y',...
%             'FontSize',14,'FontWeight','bold')

    end
    
    %the next part selects out the biggest 
    [max1, ind1] = max(area_array);
    area_array(ind1)      = -Inf;
    [max2, ind2] = max(area_array);
    
    max_area(1) = stats(ind1);
    max_area(2) = stats(ind2);
    %this var cocntains  the data about the two largest holes detected
    select_cent = max_area;
    struct2table(max_area)
    %[centers,radii] = imfindcircles(rgb,[20 25],'ObjectPolarity','dark')
    
    for k = 1:length(max_area)
        centroid = max_area(k).Centroid;
        plot(centroid(1),centroid(2),'ko');
    end
    
    
    %pixel_vals = impixel
    %imhist(B);
    %imshow(his);
    %bin = imbinarize(B);
    %imshow(bin);
    %bin = imquantize(B,50);
    %imshow(bin);
end

%% this function determines if the image needs to be rotated or not
% also it calculates the angle for which the image has to be corrected
function ang = angle_correction(cent)
    temp.Flip = 0;
    % rotation anticlockwise
    temp.Rotate = 0;
    if (cent(1).Area < cent(2).Area)
        if (cent(1).Centroid(1) > cent(2).Centroid(1))
            temp.Flip = 1;
        else

        end
    else 
        if (cent(1).Centroid(1) < cent(2).Centroid(1))
            temp.Flip = 1;
        end
    end
    x1 = cent(1).Centroid(1);
    x2 = cent(2).Centroid(1);
    y1 = cent(1).Centroid(2);
    y2 = cent(2).Centroid(2);

    slope = (y2 - y1)/(x2 - x1);
    % the rotate angle goes directly into the imrotate function thats why
    % the minus sign has been put
    temp.Rotate = rad2deg(atan(slope));
    ang = temp;
end

%% takes in the image and crops out the part with the angular distribution. 
% The img param is the image and the hole_dat is the struct with the data about the hole centers
function cropimg = ang_dist_crop(img,hole_dat,res)
    %these are the coordinates for cropping
    %input2crop = struct2table(hole_dat)
    y_ext_neg = ((240 * 2)/50)*res;
    y_ext_pos = ((240 * 2)/50)*res;
    x1 = 0;
    y1 = 0;
    x2 = 0;
    y2 = 0;
    %for imcrop the syntax is  imcrop(img,[xmin ymin width height])
    if hole_dat(1).Centroid(1) < hole_dat(2).Centroid(1)
        x1 = hole_dat(1).Centroid(1);
        x2 = hole_dat(2).Centroid(1);
        y1 = hole_dat(1).Centroid(2) - y_ext_neg;
    else
        x1 = hole_dat(2).Centroid(1);
        x2 = hole_dat(1).Centroid(1);
        y1 = hole_dat(2).Centroid(2) - y_ext_neg;
    end
    width = x2 - x1;
    height = y_ext_neg + y_ext_pos;
    cropimg = imcrop(img,[x1 y1 width height]);
end

%rec dat is the 2 corners of the rectangle to be cropped
%rec dat is given as [[x1,y1],[x2,y2]]
%res is the resolution of the image
function cropimg = ESM_crop(img,rec_dat,res)
    %these are the coordinates for cropping
    %input2crop = struct2table(hole_dat)
    %y_ext_neg = ((240 * 2)/50)*res;
    %y_ext_pos = ((240 * 2)/50)*res;
    x1 = rec_dat(1,1);
    y1 = rec_dat(1,2);
    x2 = rec_dat(2,1);
    y2 = rec_dat(2,2);
    %for imcrop the syntax is  imcrop(img,[xmin ymin width height])
%     if hole_dat(1).Centroid(1) < hole_dat(2).Centroid(1)
%         x1 = hole_dat(1).Centroid(1);
%         x2 = hole_dat(2).Centroid(1);
%         y1 = hole_dat(1).Centroid(2) - y_ext_neg;
%     else
%         x1 = hole_dat(2).Centroid(1);
%         x2 = hole_dat(1).Centroid(1);
%         y1 = hole_dat(2).Centroid(2) - y_ext_neg;
%     end
    width = x2 - x1;
    height = y2 - y1;
    cropimg = imcrop(img,[x1 y1 width height]);
end

%decides if image needs to to be inverted
%useful for ESM processing
function is_invert = Is_invert(img)
    colsum = sum(img,1);
    len = length(colsum);
    half_len = len/2;
    sum1 = sum(colsum(1:round(half_len)));
    sum2 = sum(colsum(round(half_len):end));
    if(sum1>sum2)
       is_invert = 0;
    else
        is_invert = 1;
    end
end

%returns the trend
%returns the axis in microns
%gap in cm gap is the length on the ip that lies before the magnet begins 
function [axis,trend] = trendgen_ESM(img,wid,invert,res,gap)
    %imshow(img,[])
    if invert==1
        img  = imcrop(img,180);
    end
    [rows columns] = size(img);
    row_sum = sum(img,2);
    [rowmax,rowmax_idx] = max(row_sum);
    y1 = rowmax_idx - round(wid/2);
    gap = round(gap* 0.01 / (res*(10^(-6))));
    cropimg = imcrop(img,[gap y1 columns wid]);
    %imshow(cropimg);
    trend = sum(cropimg,1);
    trend = trend/wid;
    %plot(trend);
    axis = linspace(1,columns-gap,columns-gap);
    axis = axis * (res*(10^(-6)));
    trend = trend(1:end-1);
    %length(trend)
    %length(axis)
    %plot(axis,trend)
end

%% determines the trend of the angular distribution
% takes in the cropped image nad finds the trend
function [axis,trend] = trendgen(img,res,wid)
    hole_rad1 = 150;
    hole_rad2 = 220;
    siz = size(img);
    length = siz(2);
    ang_res = 180.0 / length;
    ang_axis = linspace(0,180,length+1);
    ang_axis = ang_axis + (ang_res/2);
    ang_axis(end) = [];
    sum_trend = sum(img,1);
    size(sum_trend);
    %figure('name','trend');
    %plot(ang_axis,sum_trend);

    ang_axis = ang_axis(hole_rad1:end-hole_rad2);
    sum_trend = sum_trend(hole_rad1:end-hole_rad2);
    axis = ang_axis;
    trend = sum_trend;
end

%% Calculates the PSL value of each pixel in the image
function img_psl = psl_calc(img,res,sensit,lat,ran)
    %imshow(img,[])
    %imhist(img)
    %set(gca,'YScale','log')
    data = double(img);
    Res = res;
    L = lat;
    S = sensit;
    G = (2^ran) - 1;
    h = 3.5;
    row=size(data,1);
    col=size(data,2);
    %PSLtotal=0;
    %you guessed it right. I donno where i found this expression
    %psl = (h*(10^(L/2))*((data./G).^2)*(Res/100)^2);
    %here is teh correct expression
    psl = ((Res/100)^2) * (4000/S) * 10.^(L * ((data./G)-0.5));
    PSLnew = psl;
    %PSLnew = psl/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));
%     for i=1:row
%         for j =1:col
%             %signal(i,j)=data(i,j);
%             %signal(i,j)=data(i,j)-noise(i,j);
%             %gray value to PSL conversion FLA 7000
%             %PSL(i,j)=double(((Res/100)^2)*(10^(L*((signal(i,j)/G)-0.5))));
%             %gray value to PSL conversion for GE
%             PSL(i,j)=double(h*(10^(L/2))*((data(i,j)/G)^2)*(Res/100)^2);
%             PSLnew(i,j)= PSL(i,j)/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));
%             if data(i,j)>60000
%                PSL(i,j)
%                PSLnew(i,j) 
%                %fprintf("bleh")
%             end
%             %PSLtotal=PSLtotal+PSLnew(i,j);
%         end
%     end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    img_psl = PSLnew;
    %imshow(img_psl,[])
    %imhist(img_psl)
    %set(gca,'YScale','log')
    %fprintf("bleh")
end

function [first_scan,last_scan,file_list_run,sav_nam] = file_nam(directory,run_num)
    
    dir_list = dir(directory);
    file_list = string(zeros(1,length(dir_list)));
    %file_tif = [];

    for i = 1:length(dir_list)
        name = dir_list(i).name;
        file_list(i) = name;
    end
    % selecting the names with tiff
    TF = contains(file_list,'.tif');
    file_list = file_list(TF);
    %TF_sc1 = contains(file_list,'_sc1');
    %file_list_sc1 = file_list(TF_sc1);
    TF_run = contains(file_list,run_num);
    file_list_run = file_list(TF_run);
    for i = 1:length(file_list_run)
        file_list_run(i) = strcat(directory,file_list_run(i));
    end
    last_scan = file_list_run(end);
    first_scan = file_list_run(1);
    sav_nam = erase(first_scan,["_sc1",".tif"]);
    sav_nam = strcat(sav_nam,".csv");
end


% 
function satVal = isSat(img)
    [counts,binLocations] = imhist(img);
    %figure('name','img_hist');
    %plot(binLocations,counts);
    [count_max,ind_max] = max(counts);
    if ind_max == 256
        satVal = true;
    else
        satVal = false;
    end
 
end
% 

%energy callibration 
function [pos_val,en_val] = en_calib(cal_table,res_mic,en_ar)
    cal_table.Properties.VariableNames{2} = 'yPos';
    cal_table.Properties.VariableNames{1} = 'xPos';
    x = cal_table.xPos;
    y = cal_table.yPos;
    en_vals = length(en_ar);
    id = pdist(y);
    sqF = squareform(id);
    z = linkage(sqF);
    T = cluster(z,'maxclust',en_vals);
    Tu = unique(T);
    for i = 1:length(en_ar)
        xCurr  = x(T == Tu(i));
        yCurr =  y(T == Tu(i));
        %scatter(xCurr,yCurr);
        yCentre(i) = mean(yCurr);
        hold on
    end
    yCentre = sort(yCentre);
    %getting the calibration to pixels according to resolution
    %here we are using the convention y-axis because in the calibration
    %file the position is in the y_axis. the y axis of the calibration
    %file is the same as the x-axis of the imafge file
    %yCentre = yCentre/ (res_mic * (10^-6));
    y_axis_n_points = (max(yCentre)-min(yCentre))/(res_mic * (10^-6));
    pos_val = linspace(min(yCentre),max(yCentre),y_axis_n_points+1);
    %en_val = interp1(yCentre,en_ar,pos_val,'spline');
    p = polyfit(yCentre,en_ar,3);
    en_val = polyval(p,pos_val);
end

%decides whether trend is saturated 
function [isSat,idx] = isSat_ESM(trend,Res,sensit,lat,dyn_range)
    sat_val = (2^dyn_range)- 1000;
    %plot(trend);
    sat_val = psl_calc_point(sat_val,Res,sensit,lat,dyn_range);
    tf = islocalmax(trend,'FlatSelection','all','MinProminence',sat_val);
    idx = tf;
    %plot(trend(idx));
    num_trues = sum(tf,'all');
    if num_trues>10
        isSat = 1;
    else
        isSat = 0;
    end
end

%this is for removing the scanner lines
%please refer to the independent code sca_lin_trial_1.m
%the threshold value and the 
function filt_im = scan_line_rem(im)
    
    log_thresh = 16 ;
    %img = Tiff(file_name);
    %img = read(img);
    %imgdisp = imcrop(img,[2530 1044 159 135]);
    %subplot(2,3,1);
    %imshow(imgdisp);

    % Compute the 2D fft.
    frequencyImage = fftshift(fft2(im));
    % Take log magnitude so we can see it better in the display.
    amplitudeImage = log(abs(frequencyImage));
    minValue = min(min(amplitudeImage));
    maxValue = max(max(amplitudeImage));
    %subplot(2,3,2);
    %imshow(amplitudeImage, []);
    %title(caption, 'FontSize', fontSize);
    %axis on;

    amplitudeThreshold = log_thresh;
    brightSpikes = amplitudeImage > amplitudeThreshold; % Binary image.
    %subplot(2, 3, 3);
    %imshow(brightSpikes);
    %axis on;
    
    %please check the filtering values for every set of data
    %brightSpikes(1715:1785, :) = 0;
    %brightSpikes(1700:1800, :) = 0;
    brightSpikes(1650:1850, :) = 0;
    brightSpikes(:,1:1900) = 0;
    brightSpikes(:,2100:end) = 0;
    %subplot(2, 3, 4);
    %imshow(brightSpikes);
    %axis on;

    % Filter/mask the spectrum.
    frequencyImage(brightSpikes) = 0;
    %frequencyImage = frequencyImage - brightSpikes;
    % Take log magnitude so we can see it better in the display.
    %amplitudeImage2 = log(abs(frequencyImage));
    %subplot(2, 3, 5);
    %imshow(log(abs(frequencyImage)),[]);
    %imshow(amplitudeImage2);
    %bleh
    %axis on;

    filteredImage = ifft2(fftshift(frequencyImage));
    filteredImage = abs(filteredImage);
    %imgdisp2 = imcrop(filteredImage,[2530 1044 159 135]);
    %subplot(2,3,6);
    %imagesc(filteredImage);
    %imshow(filteredImage,[]);
    %bleh
    %imshow(imgdisp2,[])
    filt_im = filteredImage;
end

function dNdE = jacobian(en_vals,x_axis,dat)
    %figure;
    %plot(x_axis)
    %figure
    %plot(x_axis,dat);
    figure;
    %plot(x_axis,en_vals)
    dx = diff(x_axis);
    dE = diff(en_vals);
    dxdE = dx./dE;
    %figure;
    %plot(en_vals(1:end-1),dxdE);

    %dEdx = dE./dx;
    %plot(en_vals(1:end-1),dxdE);
    %size(dxdE)
    %size(dat)
    dNdE = dat(1:end-1).*dxdE;
    %plot(en_vals(1:end-1),dNdE)

end

%idx is the array of locations 
function recon_sat = sat_rem(satscan,unsatscan,idx)
    %plot(satscan - unsatscan);

    %this part selects the indices just after the region where the
    %saturation has been detected
    %plot(idx)
    num_trues = sum(idx,'all');
    num_selec = round(num_trues /2) ;
    %last 1 detected in idx
    last_idx = find(idx==1,1,'last');
    temp_idx = zeros(size(idx));
    %plot(temp_idx)
    temp_idx(last_idx+1:last_idx+num_selec+1) = 1;
    temp_idx = logical(temp_idx);
    %plot(temp_idx)
    %plot(unsatscan)
    sat_scan_temp = satscan(temp_idx);
    %plot(satscan)
    unsat_scan_temp = unsatscan(temp_idx);
    %plot(unsat_scan_temp)
    delta = sat_scan_temp - unsat_scan_temp;
    %plot(sat_scan_temp)
    length(unsat_scan_temp);
    %length(delta)
    %size(unsat_scan_temp)
    %size(delta)
    %size(delta)
    p = polyfit(unsat_scan_temp,delta,1);
    %p = fit(unsat_scan_temp,delta,'poly1');
    %plot(unsat_scan_temp,delta)
    m = p(1);
    c = p(2);
    %saturatedd part of the unsaturated scan
    %unsat_scan_sat = unsatscan(idx);
    %recon_sat = (unsatscan + c) / (1 + m);
    recon_sat = (unsatscan * (m+1)) + c;
    plot(recon_sat)
    hold on
    %recon = satscan;
    %recon_temp = recon;
    plot(satscan)
    hold off
    %recon(idx) = recon_sat;
    %fprintf("bleh");
end

function unsat_final = total_sat_rem(lin_arr,isSatarr,idx)
    lin_arr  = cell2mat(lin_arr);
    isSatarr = cell2mat(isSatarr);
    idx = cell2mat(idx);
    num_rows = size(lin_arr,1);
    num_cols = size(lin_arr,2);
    unsat = [];
    unsat(1,:) = lin_arr(end,:);
    for i = 1:num_rows-1
        sat_arr = lin_arr(num_rows-i,:);
        idx_temp = idx(num_rows-i,:);
%         size(sat_arr)
%         size(unsat)
        unsat_arr = unsat(1,:);
        if (isSatarr(num_rows-i)==1)
            res = sat_rem(sat_arr,unsat_arr,idx_temp);
        else
            res = sat_arr;
        end
        %size(res)
        %size(unsat)
        temp = [res;unsat];
        unsat = temp;
    end
    %unsat_size = size(unsat)
    unsat_final = mat2cell(unsat,ones(1,num_rows),[num_cols]);
end

%this takes one image and does all the proessing on that
%ESM_num is the number of ESM plates in a single image
%scan_num is the scan number that the image belongs to
function [axis_arr,trend_arr,isSatarr,satidx_arr,inv_arr] = img_proc(scan_num,img,crop_loc,res,ESM_num,sensit,lat,ran,wid,gap,inv_arr)
    corner_index = 1;
    fil_tim = tic;
    %disp("Starting processing for one file")
    fprintf('Starting processing for scan num %.15g',scan_num);
    %imshow(img,[]);
    %the removal of scan lines is causing issues with the saturation
    %so either we remove scan lines or the saturation coorection works
    %img = scan_line_rem(img);
    %imshow(img,[])
    toc(fil_tim);
    temp_fil_tic = tic;
    %disp("scan lines removed for one file");
    %inv_arr = zeros([1,ESM_num]);
    psl_im = psl_calc(img,res,sensit,lat,ran);
    %imshow(psl_im,[])
    %imshow(psl_im,[])
    toc(temp_fil_tic);
    disp("psl values calculated for one file");
    axis_arr = {};
    trend_arr = {};
    isSatarr = {};
    satidx_arr = {};
    for i = 1:ESM_num
        temp_fil_tic = tic;
        lef_corn = crop_loc(corner_index,:);
        corner_index = corner_index + 1;
        right_corn = crop_loc(corner_index,:);
        corner_index = corner_index + 1;
        rec_dat = [lef_corn;right_corn];
        cropim = ESM_crop(psl_im,rec_dat,res);
        imshow(cropim,[])
        toc(temp_fil_tic);
        disp("image cropped");
        temp_fil_tic = tic;
        if(scan_num == 1)
            inv_arr(i) = Is_invert(cropim);
        end
        toc(temp_fil_tic);
        disp("image inverted");
        temp_fil_tic = tic;
        [axis,trend] = trendgen_ESM(cropim,wid,inv_arr(i),res,gap);
        %plot(axis,trend)
        toc(temp_fil_tic);
        disp("trend_generated");
        [isSat,idx] = isSat_ESM(trend,res,sensit,lat,ran);
        axis_arr(i) = {axis};
        trend_arr(i) = {trend};
        isSatarr(i) = {isSat};
        satidx_arr(i) = {idx};
    end
end

%finds psl value for a single intensity
function point_psl = psl_calc_point(data,Res,sensit,lat,ran)
    L = lat;
    S = sensit;
    G = (2^ran) - 1;
    h = 3.5;
    %PSL=double(h*(10^(L/2))*((data/G)^2)*(Res/100)^2);
    %PSLnew= PSL/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));    
    PSLnew=double(((Res/100)^2) * (4000/S) * 10.^(L * ((data./G)-0.5)));
    point_psl = PSLnew;
end

%corrects the ESM trace for the ip quantum efficiency
%takes the energy axis and the trend
function trend_eff_corr = eff_corr(en_axis,trend)

    %now we neeed to fit the efficiency curves for our dataset
    %for all values above the max(at 273 keV) single exponential decay fitting
    %works
    %for vlaues below the max, we can do 1 - exponential decay
    %we have used origin to get the different values of the parameters and
    %fit the efficiency curve
    % we have used the function of the form (y = y0 + A*exp(-x/t))
    %this callibration works only till 4000keV
    %values are available for higher energies also but they have not been
    %included in the fitting

    %values for below 273keV
    y1 = 0.12486;
    A1 = -0.15699;
    t1 = 135.17822;

    %values for above 273keV
    y2 = 0.02631;
    A2 = 0.14433;
    t2 = 460.07596;

    eff_calc = zeros(1,length(en_axis));
    low_idx = find(en_axis<=273);
    high_idx = find(en_axis>273);
    eff_calc(low_idx) = y1 + A1*(exp(-en_axis(low_idx)/t1));
    eff_calc(high_idx) = y2 + A2*(exp(-en_axis(high_idx)/t2));

    trend_eff_corr = trend./eff_calc;
end





