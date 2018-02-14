% script to draw the topographies of the source selected on MRI scan
% too specific for Jordy
% Launch brainstorm
% 1. Import MRI a MRI from bst
% 2. Import Channel as channel from bst
% 3  Import cortex as cortex

% import cortex from bst
% PARAMS = struct();
% PARAMS.protocol_name = 'Waves';
% PARAMS.subject_name = 'Jordy_3';
% PARAMS.study_name = 'jordy_spont03';
% PARAMS.forward.name = 'Overlapping spheres_hr';
% PARAMS.protocol_name = 'Spikes2';
% PARAMS.cortex.low.comment = 'cortex_247306V';
% PARAMS.visualize.flattening = false;
% cortex_comment = PARAMS.cortex.low.comment;
% cortex = from_bst_get_surface(cortex_comment, PARAMS);

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
    