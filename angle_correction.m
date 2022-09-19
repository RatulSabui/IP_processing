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