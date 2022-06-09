function dNdE = jacobian(en_vals,x_axis,dat)
    %figure;
    %plot(x_axis)
    %figure
    %plot(x_axis,dat);
    figure;
    %plot(x_axis,en_vals)
    dx = diff(x_axis);
    dE = diff(en_vals);
    dxdE = dx./dE;
    %figure;
    %plot(en_vals(1:end-1),dxdE);

    %dEdx = dE./dx;
    %plot(en_vals(1:end-1),dxdE);
    %size(dxdE)
    %size(dat)
    dNdE = dat(1:end-1).*dxdE;
    %plot(en_vals(1:end-1),dNdE)

end