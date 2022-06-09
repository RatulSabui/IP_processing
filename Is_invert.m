%decides if image needs to to be inverted
%useful for ESM processing
function is_invert = Is_invert(img)
    colsum = sum(img,1);
    len = length(colsum);
    half_len = len/2;
    sum1 = sum(colsum(1:round(half_len)));
    sum2 = sum(colsum(round(half_len):end));
    if(sum1>sum2)
       is_invert = 0;
    else
        is_invert = 1;
    end
end
