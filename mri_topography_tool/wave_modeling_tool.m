function wave_modeling_tool(MRI, channel, cortex, G)
% script to draw the topographies of the source selected on MRI scan
% too specific for Jordy
% Launch brainstorm
% Inputs:
% 1. Import MRI a MRI from bst
% 2. Import Channel as channel from bst
% 3. Import cortex as cortex
% 4. Import G as G from bst

faces = cortex.Faces;
vertices = cortex.Vertices;

% mean locations of MEG coils
k = 1;
for i = 32:99
    chan_loc(k,:) = mean(channel.Channel(i).Loc, 2)';
    k = k+1;
end
[X,Y] = bst_project_2d(chan_loc(:,1), chan_loc(:,2), chan_loc(:,3), '2dlayout');


% forward model for restricted orientations
G1 = bst_gain_orient(G.Gain,G.GridOrient);
G1 = G1(32:99,:);

% graph
h = figure;
setappdata(h, 'pos', [128, 128, 128])
set(h, 'WindowButtonDownFcn', @(src,evnt)PlotMri(MRI, getappdata(h, 'pos'), h))
subplot(2,3,1)
imagesc(flipud(squeeze(MRI.Cube(:,128,:))'))
set(gca,'tag',num2str(1))
subplot(2,3,2)
imagesc(flipud(squeeze(MRI.Cube(128,:,:))'))
set(gca,'tag',num2str(2))
subplot(2,3,3)
imagesc(flipud(squeeze(MRI.Cube(:,:,128))'))
set(gca,'tag',num2str(3))
colormap('gray')
btn = uicontrol('Style', 'pushbutton', 'String', 'Pick points',...
        'Position', [150 100 150 50],...
        'Callback', @(src, evnt)btn(h, MRI));
btn2 = uicontrol('Style', 'pushbutton', 'String', 'Compute wave',...
        'Position', [300 100 150 50],...
        'Callback', @(src, evnt)btn2(h, MRI, cortex, G1, X, Y));
end

function PlotMri(MRI, pos, h)
    k = str2num(get(gca, 'tag'));
    click = get(gca,'CurrentPoint');
    click = round(click);

    if k == 1
        clx = round(click(1,1)); % sagittal
        cly = pos(2);        % coronal
        clz = round(click(1,2)); % axial

        subplot(2,3,1)
        imagesc(flipud(squeeze(MRI.Cube(:,pos(2),:))'))
        set(gca,'tag',num2str(1))
        hold on
        plot(clx,clz,'r+', 'MarkerSize', 15);

        subplot(2,3,2)
        imagesc(flipud(squeeze(MRI.Cube(256-clx,:,:))'))
        set(gca,'tag',num2str(2))
        hold on
        plot(pos(2),clz,'r+', 'MarkerSize', 15);

        subplot(2,3,3)
        imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
        set(gca,'tag',num2str(3))
        hold on
        plot(clx,256-pos(2),'r+', 'MarkerSize', 15);
        pos([1,3]) = [256-clx, 256-clz];
    else
        if k == 2
            clx = pos(1); % sagittal
            cly = round(click(1,1));        % coronal
            clz = round(click(1,2)); % axial

            subplot(2,3,2)
            imagesc(flipud(squeeze(MRI.Cube(256-pos(1),:,:))'))
            set(gca,'tag',num2str(2))
            hold on
            plot(cly,clz,'r+', 'MarkerSize', 15);

            subplot(2,3,1)
            imagesc(flipud(squeeze(MRI.Cube(:,cly,:))'))
            set(gca,'tag',num2str(1))
            hold on
            plot(pos(1),clz,'r+', 'MarkerSize', 15);

            subplot(2,3,3)
            imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
            set(gca,'tag',num2str(3))
            hold on
            plot(pos(1),256-cly,'r+', 'MarkerSize', 15);
            pos([2,3]) = [cly, clz];

        else
            clx = round(click(1,1)); % sagittal
            cly = round(click(1,2));        % coronal
            clz = pos(3); % axial

            subplot(2,3,3)
            imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
            set(gca,'tag',num2str(3))
            hold on
            plot(clx,cly,'r+', 'MarkerSize', 15);

            subplot(2,3,1)
            imagesc(flipud(squeeze(MRI.Cube(:,256-cly,:))'))
            set(gca,'tag',num2str(1))
            hold on
            plot(clx,clz,'r+', 'MarkerSize', 15);

            subplot(2,3,2)
            imagesc(flipud(squeeze(MRI.Cube(clx,:,:))'))
            set(gca,'tag',num2str(2))
            hold on
            plot(256-cly,clz,'r+', 'MarkerSize', 15);
            pos([1,2]) = [clx, 256-cly];
        end
    end

    setappdata(h, 'pos', pos)

end

function btn(h, MRI)
    set(h, 'WindowButtonDownFcn', @(src,evnt)pickpoints(h, MRI))
end

function pickpoints(h, MRI)
    k = str2num(get(gca, 'tag'));
    click = get(gca,'CurrentPoint');
    pos = getappdata(h, 'pos');
    if k == 1
        setappdata(h, 'coord', [getappdata(h, 'coord'), [click(1,1), pos(2), click(1,2)]]);
        coord = getappdata(h, 'coord');
        subplot(2,3,1)
        imagesc(flipud(squeeze(MRI.Cube(:,pos(2),:))'))
        set(gca,'tag',num2str(1))
        hold on
        plot(coord(1:3:length(coord)), coord(3:3:length(coord)), 'r*', 'MarkerSize', 2);
    elseif k == 2
        setappdata(h, 'coord', [getappdata(h, 'coord'), [256-pos(1), click(1,1), click(1,2)]]);
        coord = getappdata(h, 'coord');
        subplot(2,3,2)
        imagesc(flipud(squeeze(MRI.Cube(256-pos(1),:,:))'))
        set(gca,'tag',num2str(2))
        hold on
        plot(coord(2:3:length(coord)), coord(3:3:length(coord)), 'r*', 'MarkerSize', 2);
    else
        setappdata(h, 'coord', [getappdata(h, 'coord'), [click(1,1), click(1,2), 256-pos(3)]]);
        coord = getappdata(h, 'coord');
        subplot(2,3,3)
        imagesc(flipud(squeeze(MRI.Cube(:,:,256-pos(3)))'))
        set(gca,'tag',num2str(3))
        hold on
        plot(coord(1:3:length(coord)), coord(2:3:length(coord)), 'r*', 'MarkerSize', 2);
    end
    
end

function btn2(h, MRI, cortex, G1, X, Y)
    mri_coord = getappdata(h, 'coord');
    mri_coord = reshape(mri_coord, [3,length(mri_coord)/3])';
    mri_coord(:,3) = 257-mri_coord(:,3); 
    coord_scs = cs_convert(MRI, 'voxel', 'scs', mri_coord)';
    
    cm_viridis = viridis(100);
    p5 = subplot(2,3,5);
    cla(p5)
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
    scatter3(coord_scs(1,:), coord_scs(2,:), coord_scs(3,:), 100, 'filled', 'r')
    colormap(p5, cm_viridis);
    
    for i = 1:size(coord_scs,2)
        dist(:,i) = sum((bsxfun(@minus, cortex.Vertices, coord_scs(:,i)').^2)');
        [val(i), ind(i)] = min(dist(:,i));
    end
    
    npoints = size(coord_scs,2);
    n = 0:3;
    t = 0:(npoints-1);
    waves = zeros(npoints, npoints);
    for i = 1:npoints
       wave = (1 + cos(2*pi*(n-t(i))/npoints));
       if npoints >= (max(n)+1)
        waves(i, n(n<npoints)+1) = wave;
       else waves(i, n(n<npoints)+1) = wave(1:length(n)-(max(n)+1-npoints));
       end
       n = n+1;
    end
    
    FM = G1(:,ind);
    % waves on sensors  
    sensor_waves = FM*waves';
    setappdata(h, 'sensor_waves', sensor_waves)
    
    
    p6 = subplot(2,3,6)
    plot_topo(X, Y, sensor_waves(:,1))
    colormap(p6, cm_viridis);
    sld = uicontrol('Style', 'slider', 'Min',1,'Max',size(sensor_waves,2),'Value',1, ...
        'Position', [1400 20 250 20], 'Callback', ...
        @(src,evnt)changetopo(src, evnt, h, X, Y), 'SliderStep',[1/(size(sensor_waves,2)-1) 1]);
    
end

function changetopo(source, event, h, X, Y)
    cm_viridis = viridis(100);
    sensor_waves = getappdata(h, 'sensor_waves');
    p6 = subplot(2,3,6)
    plot_topo(X, Y, sensor_waves(:,source.Value))
    colormap(p6, cm_viridis);
end
