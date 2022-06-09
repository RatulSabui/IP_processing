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