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