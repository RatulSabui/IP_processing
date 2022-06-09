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