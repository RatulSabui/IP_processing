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