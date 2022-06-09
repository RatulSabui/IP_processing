function unsat_final = total_sat_rem(lin_arr,isSatarr,idx)
    lin_arr  = cell2mat(lin_arr);
    isSatarr = cell2mat(isSatarr);
    idx = cell2mat(idx);
    num_rows = size(lin_arr,1);
    num_cols = size(lin_arr,2);
    unsat = [];
    unsat(1,:) = lin_arr(end,:);
    for i = 1:num_rows-1
        sat_arr = lin_arr(num_rows-i,:);
        idx_temp = idx(num_rows-i,:);
%         size(sat_arr)
%         size(unsat)
        unsat_arr = unsat(1,:);
        if (isSatarr(num_rows-i)==1)
            res = sat_rem(sat_arr,unsat_arr,idx_temp);
        else
            res = sat_arr;
        end
        %size(res)
        %size(unsat)
        temp = [res;unsat];
        unsat = temp;
    end
    %unsat_size = size(unsat)
    unsat_final = mat2cell(unsat,ones(1,num_rows),[num_cols]);
end

