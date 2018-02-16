function mri_topography_tool(MRI, channel, cortex, G)
% script to draw the topographies of the source selected on MRI scan
% too specific for Jordy
% Launch brainstorm
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
btn = uicontrol('Style', 'pushbutton', 'String', 'Show topography',...
        'Position', [150 100 150 50],...
        'Callback', @(src, evnt)plttp(MRI, cortex, getappdata(h, 'pos'), G1, chan_loc));
end

function PlotMri(MRI, pos, h)
    k = str2num(get(gca, 'tag'));
    click = get(gca,'CurrentPoint');
        
    if k == 1
        clx = round(click(1,1)); % sagittal
        cly = pos(2);        % coronal
        clz = round(click(1,2)); % axial
        
        subplot(2,3,1)
        imagesc(flipud(squeeze(MRI.Cube(:,cly,:))'))
        set(gca,'tag',num2str(1))
        hold on
        plot(clx,clz,'r+', 'MarkerSize', 30);
        
        subplot(2,3,2)
        imagesc(flipud(squeeze(MRI.Cube(256-clx,:,:))'))
        set(gca,'tag',num2str(2))
        hold on
        plot(pos(2),clz,'r+', 'MarkerSize', 30);
        
        subplot(2,3,3)
        imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
        set(gca,'tag',num2str(3))
        hold on
        plot(clx,256-pos(2),'r+', 'MarkerSize', 30);
        pos([1,3]) = [clx, clz];
    else
        if k == 2
            clx = pos(1); % sagittal
            cly = round(click(1,1));        % coronal
            clz = round(click(1,2)); % axial
            
            subplot(2,3,2)
            imagesc(flipud(squeeze(MRI.Cube(256-pos(1),:,:))'))
            set(gca,'tag',num2str(2))
            hold on
            plot(cly,clz,'r+', 'MarkerSize', 30);
            
            subplot(2,3,1)
            imagesc(flipud(squeeze(MRI.Cube(:,cly,:))'))
            set(gca,'tag',num2str(1))
            hold on
            plot(256-pos(1),clz,'r+', 'MarkerSize', 30);
            
            subplot(2,3,3)
            imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
            set(gca,'tag',num2str(3))
            hold on
            plot(pos(1),256-cly,'r+', 'MarkerSize', 30);
            pos([2,3]) = [cly, clz];
            
        else
            clx = round(click(1,1)); % sagittal
            cly = round(click(1,2));        % coronal
            clz = pos(3); % axial
            
            subplot(2,3,3)
            imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
            set(gca,'tag',num2str(3))
            hold on
            plot(clx,cly,'r+', 'MarkerSize', 30);
            
            subplot(2,3,1)
            imagesc(flipud(squeeze(MRI.Cube(:,256-cly,:))'))
            set(gca,'tag',num2str(1))
            hold on
            plot(clx,clz,'r+', 'MarkerSize', 30);
            
            subplot(2,3,2)
            imagesc(flipud(squeeze(MRI.Cube(clx,:,:))'))
            set(gca,'tag',num2str(2))
            hold on
            plot(256-cly,clz,'r+', 'MarkerSize', 30);
            pos([1,2]) = [clx, 256-cly];
        end
    end

    setappdata(h, 'pos', pos)
 
end

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
    disp(ind)
    
    p5 = subplot(2,3,5);
    plot_topo(chan_loc, G1(:,ind))
    colormap(p5, cm_viridis);
end
