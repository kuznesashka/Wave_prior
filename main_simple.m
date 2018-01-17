cd C:\Users\Administrator\Documents\brainstorm3
brainstorm
cd C:\Kuznetsova\WavePrior\git\waves

%% Initial parameters
PARAMS = struct();
PARAMS.protocol_name = 'Waves';
PARAMS.subject_name = 'Jordy_3';
PARAMS.study_name = 'jordy_spont03';
PARAMS.channels_idx = 32:99; 
PARAMS.sampling_rate = 250;
PARAMS.max_distance = 0.04; 
PARAMS.forward.name = 'Overlapping_spheres3';
PARAMS.protocol_name = 'Spikes2';

PARAMS.cortex.low.comment = 'cortex_247306V';
PARAMS.visualize.flattening = false;
cortex_comment = PARAMS.cortex.low.comment;
cortex = from_bst_get_surface(cortex_comment, PARAMS);
faces = cortex.Faces;
vertices = cortex.Vertices;
PARAMS.wave.half_width = 0.005; % (m) half width of one wave
PARAMS.wave.duration = 0.02; % (s) time of wave duration 
PARAMS.wave.speeds = [0.001, 0.0025, 0.005, 0.01, 0.03 , 0.05, 0.07, 0.1, 0.2, 0.3, 0.4]; % (mm/ms = m/s)1

spike_time = {[11.488, 11.584],[12.116, 12.212], [41.940, 42.036], [42.160, 42.256], [42.640, 42.736]};
[meas, measfull]  = from_bst_get_measurements('jordy_spont03 (#1) band(1-40Hz)', PARAMS);
Ff = meas{1}.F(:,:);

% Forward model matrices
[G G3] = from_bst_get_gain_matrix(PARAMS.forward.name, PARAMS); 
Ns = size(G3,2)/3; % number of sources

ind = find(meas{1}.Time > 0.456 & meas{1}.Time < 0.5520); % averaged spike

figure
plot(Ff(:,ind)')

%% First point
Data = Ff(:,ind(7:9));

[U,S,V] = svd(Data);
MaxC = zeros(1,Ns);
SC = G'*U(:,1:7); 
for i=1:Ns
    MaxC(i) = sum(SC(i,:).^2); 
end;
[ValMax, IndMax] = max(MaxC); % the best fitted source, starting point


colors = GetColors(3);
figure;
h_br = trisurf(cortex.Faces, ...
               cortex.Vertices(:,1), ...
               cortex.Vertices(:,2), ...
               cortex.Vertices(:,3), ...
               'FaceColor', colors(1).RGB, ...
               'EdgeColor','none', ...
               'FaceAlpha', 1);
light
lighting phong;
axis equal;
view(360,360);
axis off;
grid off;
%camlight('left');
hold on
sphere_marker(vertices(IndMax,1),vertices(IndMax,2), vertices(IndMax,3), 0.002, colors(7).RGB)


%% MSE
ind = find(meas{1}.Time > 0.468 & meas{1}.Time < 0.5240); % averaged spike
Data = Ff(:,ind);

figure
plot(Data')

%% Neighbourhood
DIST = zeros(1, size(vertices,1));
for d = 1:size(vertices,1)
    DIST(d) = norm(vertices(d,:)-vertices(idx,:));
end

% DIST = zeros(1, size(vertices,1));
% for d = 1:size(vertices,1)
%     DIST(d) = norm(vertices(d,:)-[-0.015,-0.06, 0.045]);
% end
% [val, idx] = min(DIST);

ind_n = find(DIST<0.015);

surf_hr = cortex;
surf_lr = cortex;
vert_inds_lr = [];
vertices = cortex.Vertices;
minmsecort = repmat(0,1,size(vertices,1));
minmsecort(ind_n) = 1;
data_lr = minmsecort';

h = plot_brain_cmap2(surf_hr, surf_lr, vert_inds_lr, data_lr)
hold on
sphere_marker(vertices(idx, 1),vertices(idx, 2),vertices(idx, 3), 0.002, [0,0,0])

for npoint = 1:length(ind_n)
    PARAMS.seed_vertices = ind_n(npoint);
    sensor_waves{npoint} = wave_on_sensors_simple(cortex, PARAMS);
    npoint
end

%% New part
load ind_n
load sensor_waves

Nsens = size(Ff,1);
wavelength = PARAMS.wave.duration/(1/PARAMS.sampling_rate) + 1;

clear bestcoef bestMSE bestshift
l = 1;
for s = 8;
   for npoint = 1:length(ind_n)
        tic
        rs = 1;
        clear DF MSE SPEED SHIFT B STATS
    
        waves = sensor_waves{npoint}; 
        Ndir = size(waves,1);

        %profile on
        wavesspeed = reshape(cell2mat(waves(:,s)),[Nsens, wavelength, Ndir]); 
        for d = 1:Ndir % for each direction
            tmp = squeeze(wavesspeed(:,:,d)); % matrix (num sensors) by (num time points)
            wavesspeedlin(d,:) = tmp(:); % (num dir) by (num sensors)x(num time points)
        end
        R = size(Data,2)-wavelength; % time points minus the wave length
        range = 1:wavelength;
        for r = 1:R % sliding time window
            tmp = Data(:,range); % time interval with wave
            DataLin = tmp(:);
            lmbd = lambda(wavesspeedlin',DataLin,'CV', 10, 'Alpha', 0.5);
            %[Bestcor{rs}, Bestshift{rs}] = max(corr(DataLin, wavesspeedlin'));
            [B{rs},STATS{rs}] = lasso(wavesspeedlin',DataLin,'CV', 10, 'Alpha', 0.5);%, 'Lambda', lmbd(40)); 
            MSE(rs,:) = sqrt(STATS{rs}.MSE)/norm(DataLin);
            DF(rs,:) = STATS{rs}.DF;
            SPEED(rs) = PARAMS.wave.speeds(s);
            SHIFT(rs,1) = r;
            range = range+1;
            rs = rs+1;
        end;

        [MinMSE MinMSEShift] = min(MSE);
        MinDF = DF(MinMSEShift);

        bestshift(l,npoint) = MinMSEShift;
        bestcoef{npoint, l} = B{1,MinMSEShift};
        bestMSE(npoint,l) = MinMSE;
        
        npoint
        l
    end
    
    l = l+1; 
end

load cortex_infl
surf_hr_inf = cortex_infl;
surf_hr = cortex;
surf_lr = cortex;
vert_inds_lr = [];
vertices = cortex.Vertices;
minmsecort = repmat(min(1-bestMSE),1,size(vertices,1));
minmsecort(ind_n) = 1-bestMSE;
data_lr = minmsecort';

h = plot_brain_cmap2(surf_hr_inf, surf_lr, vert_inds_lr, data_lr)
h = plot_brain_cmap2(surf_hr, surf_lr, vert_inds_lr, data_lr)


% picture
theta = 0:0.01:pi/4;
[val, ind] = sort(abs(bestcoef), 'descend');
bestcoef2 = bestcoef(ind);

figure
for i = 1:5
    if ind(i) ~=5
        rho = sin(2*theta).*cos(2*theta)*2*abs(bestcoef2(i));
        theta2 = theta+pi/8+pi/2*(ind(i)-3);
        if bestcoef2(i)>=0
            col = 'r';
        else col = 'b';
        end
        polar(theta2, rho, col)
        hold on
    else th = linspace(0,2*pi,50);
        r = abs(bestcoef2(i));
        if bestcoef2(i)>=0
            col = 'r';
        else col = 'b';
        end
  
    polar(th,r+zeros(size(th)), col)
    hold on
    end
end


% POSITIVE LASSO
rs = 1;
clear DF MSE SPEED SHIFT B

for s = 1:length(PARAMS.wave.speeds) % for each speed value
    wavesspeed = reshape(cell2mat(waves(:,s)),[Nsens, wavelength, Ndir]); 
    for d = 1:Ndir % for each direction
        tmp = squeeze(wavesspeed(:,:,d)); % matrix (num sensors) by (num time points)
        wavesspeedlin(d,:) = tmp(:); % (num dir) by (num sensors)x(num time points)
    end
    R = size(Data,2)-wavelength; % time points minus the wave length
    range = 1:wavelength;
    for r = 1:R % sliding time window
        tmp = Data(:,range); % time interval with wave
        DataLin = tmp(:);
        [B{rs}, status{rs}] = l1_ls_nonneg(wavesspeedlin([1,3],:)',DataLin,lambda(rs),0.1,1);
        MSE(rs) = sqrt(mean((B{rs}'*wavesspeedlin([1,3],:) - DataLin').^2))/norm(DataLin);
        DF(rs) = sum(B{rs}>10*(-11));
        SPEED(rs) = PARAMS.wave.speeds(s);
        SHIFT(rs,1) = r;
        range = range+1;
        rs = rs+1;
    end;
end

clear MinMSE MinMSEShift MinDF
range  = 1:R;
for r = 1:size(MSE,2)/R
    [MinMSE(r) MinMSEShift(r)] = min(MSE(range));
     MinDF(r) = DF(range(MinMSEShift(r)));
    range= range+R;
end;

figure
subplot(2,1,1)
plot(PARAMS.wave.speeds, MinMSE)
subplot(2,1,2)
plot(PARAMS.wave.speeds, MinDF)

[valsort indsort] = sort(MinMSE);
[valdf inddf] = min(MinDF(indsort(1:3)));
indsorttop = indsort(1:3);
bestspeed = PARAMS.wave.speeds(indsorttop(inddf));
bestshift = MinMSEShift(indsorttop(inddf));
bestcoef = B{1,(indsorttop(inddf)-1)*R+bestshift};
bestMSE = MinMSE(indsorttop(inddf));

% picture
theta = 0:0.01:pi/4;
[val, ind] = sort(abs(bestcoef), 'descend');
bestcoef2 = bestcoef(ind);

figure
for i = 1:2
    rho = sin(2*theta).*cos(2*theta)*2*abs(bestcoef2(i));
    theta2 = theta+pi/8+pi/1*(ind(i)-2);
    if bestcoef2(i)>=0
        col = 'r';
    else col = 'b';
    end
    polar(theta2, rho, col)
    hold on
end

