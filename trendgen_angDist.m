%this sums up the entire angular distribution plate and give sthe trend
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
