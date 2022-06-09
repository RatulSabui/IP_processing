%finds psl value for a single intensity value
function point_psl = psl_calc_point(data,Res,sensit,lat,ran)
    L = lat;
    S = sensit;
    G = (2^ran) - 1;
    h = 3.5;
    PSL=double(h*(10^(L/2))*((data/G)^2)*(Res/100)^2);
    PSLnew= PSL/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    point_psl = PSLnew;
end


