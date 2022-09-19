%%this code is for processing the ang_dist data
%taken from elec_dat_precode2.m which does only the angular distriburtion
%data
%this was the last code created by Ratul on 4th Apr 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;


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

width = 2; %width of the cutout which is finally used for trend generation

trend_type = "slice";
%trend_type = "whole";


cd (folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this is to find the circles and the correction angle
img_fil = Tiff(scan1);
%cl = class(img_fil)
%getTag(img_fil,'YResolution');
B = read(img_fil);
%class(B)
%[corn,length,breadth] = rect_find(img);
%bleh
[holes,big_holes] = cent_find(B);
%hol_tab = struct2table(holes)
%bhol_tab = struct2table(big_holes)
if (length(holes) ~= 2)
    fprintf("more than 2 holes detected");
    %holes_table = struct2table(holes)
    %hole_mat = cell2mat(struct2cell(holes))
    %areas = holes.Area
    holes = big_holes;
end
%hol_tab2 = struct2table(holes)
an = angle_correction(holes);
if an.Flip == 1
    B_rot = imrotate(B,180 + an.Rotate);
else
    B_rot = imrotate(B,an.Rotate);
end
[holes_rot,big_holes_rot] = cent_find(B_rot);
%length(holes_rot)
%length(big_holes_rot)

if (length(holes_rot) ~= 2)
    holes_rot = big_holes_rot;
    %struct2table(holes_rot)
end

%bhol_tab3 = struct2table(holes_rot)

dat_table = table;
unsat_dat_table = table;
%% this is reading the last file and doing operations on it
for i = 1:length(total_scan)
    img_fil2 = Tiff(total_scan(i));
    %cl = class(img_fil)
    %getTag(img_fil2,'YResolution');
    C = read(img_fil2);
    if an.Flip == 1
        C_rot = imrotate(C,180 + an.Rotate);
    else
        C_rot = imrotate(C,an.Rotate);
    end

    %
    %C_rot = scan_line_rem(C_rot);
    C_crop = ang_dist_crop(C_rot,holes_rot,microns_per_pixel);
    figure('name','cropped image');
    imshow(C_crop);

    Sat = isSat(C_crop);
    %surf(C_crop);
    %colorbar;
    %C_psl = psl_calc(C_crop,microns_per_pixel,sensit,latitude,dyn_range);
    %surf(C_psl);
    %colorbar;
    
    %this tredgen sums the entire angular distribution plate columnwise and
    %then generted the trend
    if trend_type == "whole"
        [x_ax,dist] = trendgen_ang(C_crop);
    end

    if trend_type == "slice"    
        %this takes a thin slice from the ang dist cutout and then generates
        %teh trend
        [x_ax,dist] = trendgen_ang_slice(C_crop,width,microns_per_pixel);
        [isSat_val,idx] = isSat_ESM(dist,microns_per_pixel,sensit,latitude,dyn_range);

    end
    
    %figure('name','trend');
    %plot(x_ax,dist);
    dist_norm = dist / max(dist);
    ax_name = strcat('scan',int2str(i),'_ax');
    dist_name = strcat('scan',int2str(i),'_dist');
    dist_norm_name = strcat('scan',int2str(i),'_distnorm');
    dat_table.(ax_name) = transpose(x_ax);
    dat_table.(dist_name) = transpose(dist);
    dat_table.(dist_norm_name) = transpose(dist_norm);
end



writetable(dat_table,save_name,'Delimiter',',');



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

%% determines the trend of the angular distribution
% takes in the cropped image
function [axis,trend] = trendgen_angDist(img)
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

%theres a working version of this function later
% %% Calculates the PSL value of each pixel in the image
% function img_psl = psl_calc(img,res,sensit,lat,ran)
%     img_psl = zeros(size(img));
%     for i = 1:size(img,1)
%         for j = 1:size(img,2)
%             img_psl(i,j) = ((res/100)^2) * (4000/sensit) * (10 ^ (lat*((img(i,j)/((2^ran)-1))-0.5)));
%         end
%     end
% end

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
    psl = (h*(10^(L/2))*((data./G).^2)*(Res/100)^2);
    %PSLnew = psl;
    PSLnew = psl/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
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

%returns the trend
%returns the axis in degrees
%this takes a linecut of the angular distribution and then generates the trend 
function [axis,trend] = trendgen_ang_slice(img,wid,res)
    hole_rad1 = 150;
    hole_rad2 = 220;
    siz = size(img);
    length = siz(2);
    ang_res = 180.0 / length;
    ang_axis = linspace(0,180,length+1);
    ang_axis = ang_axis + (ang_res/2);
    ang_axis(end) = [];
    
    [rows,columns] = size(img);
    row_sum = sum(img,2);
    [rowmax,rowmax_idx] = max(row_sum);
    y1 = rowmax_idx - round(wid/2);
    %gap = round(gap* 0.01 / (res*(10^(-6))));

    cropimg = imcrop(img,[hole_rad1 y1 columns-hole_rad2 wid]);
    %imshow(cropimg);
    trend = sum(cropimg,1);
    trend = trend/wid;
    %plot(trend);
    ang_axis = ang_axis(hole_rad1:end-hole_rad2);
    %length(trend)
    %length(axis)
    %plot(axis,trend)
    axis = ang_axis;
    %trend = sum_trend;
end

%finds psl value for a single intensity value
function point_psl = psl_calc_point(data,Res,sensit,lat,ran)
    L = lat;
    S = sensit;
    G = (2^ran) - 1;
    h = 3.5;
    PSL=double(h*(10^(L/2))*((data/G)^2)*(Res/100)^2);
    PSLnew= PSL/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    point_psl = PSLnew;
end


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

function unsat_final = total_sat_rem_ang(lin_arr,isSatarr,idx)
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

