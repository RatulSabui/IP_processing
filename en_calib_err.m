%energy callibration along with estimation of error
function [pos_val,en_val,en_err] = en_calib_err(cal_table,res_mic,en_ar)
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
        ystd(i) = std(yCurr);
        hold on
    end
    %yCentre = sort(yCentre);
    
    [~, sorted_indices] = sort(yCentre);
    yCentre = yCentre(sorted_indices);
    ystd = ystd(sorted_indices);

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
    %plot(pos_val,en_vals)

    p1 = polyfit(en_ar,ystd,3);
    er_val = polyval(p1,en_val);
    
    dx = diff(pos_val);
    dE = diff(en_val);
    dEdx = dE./dx;
    en_err = dEdx.*er_val;

end