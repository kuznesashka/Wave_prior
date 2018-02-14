function plttp(MRI, cortex, pos, G1, chan_loc)
    coord_mri = pos;
    coord_mri(3) = 256-coord_mri(3);
    coord_scs = cs_convert(MRI, 'voxel', 'scs', coord_mri)';
    
    cm_viridis = viridis(100);
    p4 = subplot(2,3,4);
    cla(p4)
    trisurf(cortex.Faces, cortex.Vertices(:,1), ...
        cortex.Vertices(:,2), cortex.Vertices(:,3), zeros(1,size(cortex.Vertices(:,1),1)),...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    light;
    lighting phong;
    axis equal
    grid off
    axis off
    view(360,360)
    hold on
    scatter3(coord_scs(1), coord_scs(2), coord_scs(3), 100, 'filled', 'r')
    colormap(p4, cm_viridis);

    dist = sum((bsxfun(@minus, cortex.Vertices, coord_scs').^2)');
    [val, ind] = min(dist);
    
    p5 = subplot(2,3,5);
    plot_topo(chan_loc, G1(:,ind))
    colormap(p5, cm_viridis);
end