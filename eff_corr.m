%corrects the ESM trace for the ip quantum efficiency
%takes the energy axis and the trend
function trend_eff_corr = eff_corr(en_axis,trend)

    %now we neeed to fit the efficiency curves for our dataset
    %for all values above the max(at 273 keV) single exponential decay fitting
    %works
    %for vlaues below the max, we can do 1 - exponential decay
    %we have used origin to get the different values of the parameters and
    %fit the efficiency curve
    % we have used the function of the form (y = y0 + A*exp(-x/t))
    %this callibration works only till 4000keV
    %values are available for higher energies also but they have not been
    %included in the fitting

    %values for below 273keV
    y1 = 0.12486;
    A1 = -0.15699;
    t1 = 135.17822;

    %values for above 273keV
    y2 = 0.02631;
    A2 = 0.14433;
    t2 = 460.07596;

    eff_calc = zeros(1,length(en_axis));
    low_idx = find(en_axis<=273);
    high_idx = find(en_axis>273);
    eff_calc(low_idx) = y1 + A1*(exp(-en_axis(low_idx)/t1));
    eff_calc(high_idx) = y2 + A2*(exp(-en_axis(high_idx)/t2));

    trend_eff_corr = trend./eff_calc;
end
