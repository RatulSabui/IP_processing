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