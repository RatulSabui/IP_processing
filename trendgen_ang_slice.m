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