%returns the trend
%returns the axis in microns
%gap in cm gap is the length on the ip that lies before the magnet begins 
function [axis,trend] = trendgen_ESM(img,wid,invert,res,gap)
    %imshow(img,[])
    if invert==1
        img  = imrotate(img,180);
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
