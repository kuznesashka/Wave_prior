function plot_topo(X, Y, data)
    % scatter3(meg_loc_3d(:,1), meg_loc_3d(:,2), meg_loc_3d(:,3));
    % data = ones(151,1);
    
    F = scatteredInterpolant(X,Y,data, 'natural');

    r = 1.01 * max(sqrt(X .^ 2 + Y .^ 2));

    [phi,RR] = meshgrid(0:0.001:2.01 * pi,0:0.001:r);
    ZZ = F(RR .* cos(phi),RR .* sin(phi));

    surf(RR .* cos(phi), RR .* sin(phi), zeros(size(ZZ)),  ZZ, 'EdgeColor', 'none');
    view(-90,90)
    hold on;
    scatter(X,Y,100,'k.');
    scatter(X(5), Y(5), 100, 'r')
    axis off
end


