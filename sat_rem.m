%idx is the array of locations 
function recon_sat = sat_rem(satscan,unsatscan,idx)
    %plot(satscan - unsatscan);

    %this part selects the indices just after the region where the
    %saturation has been detected
    %plot(idx)
    num_trues = sum(idx,'all');
    num_selec = round(num_trues /2) ;
    %last 1 detected in idx
    last_idx = find(idx==1,1,'last');
    temp_idx = zeros(size(idx));
    %plot(temp_idx)
    temp_idx(last_idx+1:last_idx+num_selec+1) = 1;
    temp_idx = logical(temp_idx);
    %plot(temp_idx)
    %plot(unsatscan)
    sat_scan_temp = satscan(temp_idx);
    %plot(satscan)
    unsat_scan_temp = unsatscan(temp_idx);
    %plot(unsat_scan_temp)
    delta = sat_scan_temp - unsat_scan_temp;
    %plot(sat_scan_temp)
    length(unsat_scan_temp);
    %length(delta)
    %size(unsat_scan_temp)
    %size(delta)
    %size(delta)
    p = polyfit(unsat_scan_temp,delta,1);
    %p = fit(unsat_scan_temp,delta,'poly1');
    %plot(unsat_scan_temp,delta)
    m = p(1);
    c = p(2);
    %saturatedd part of the unsaturated scan
    %unsat_scan_sat = unsatscan(idx);
    %recon_sat = (unsatscan + c) / (1 + m);
    recon_sat = (unsatscan * (m+1)) + c;
    plot(recon_sat)
    hold on
    %recon = satscan;
    %recon_temp = recon;
    plot(satscan)
    hold off
    %recon(idx) = recon_sat;
    %fprintf("bleh");
end
